## inferring D from rejection-sampled paths
library(tidyverse)
library(Rcpp)
source("./funcs_crypt_packing.r")
source("./joining_slide_to_block_to_patient.r")

slide_age_ref = block_to_slide_stag %>% mutate(mark="STAG2") %>% 
   bind_rows(block_to_slide_kdm %>% mutate(mark="KDM6A")) %>% left_join(block_to_patient) %>% 
   rename(Scan_ID = Slide_ID)

## load data
patchdat_10_kdm = load_and_tidy_data("KDM6A", 10)
patchdat_10_stag = load_and_tidy_data("STAG2", 10)

## load large WT patch data as baseline
wt_patches_kdm6a = readxl::read_xlsx("./databases/patches_30_measurements_KDM6A_290520.xlsx") %>% 
   magrittr::set_colnames(c("Name", "Class", "ROI", "centroid_x", "centroid_y",
                            "area", "perimeter", "max_length", "Scan_ID"))
wt_patches_stag2 = readxl::read_xlsx("./databases/pa30_measurements_STAG2_290520.xlsx") %>% 
   magrittr::set_colnames(c("Name", "Class", "ROI", "centroid_x", "centroid_y",
                            "area", "perimeter", "max_length", "Scan_ID"))

## prime data
data_primed_kdm_10 = prime_data(wt_patches_kdm6a, patchdat_10_kdm)
data_primed_stag_10 = prime_data(wt_patches_stag2, patchdat_10_stag)

## define length and area scale
domain_size = data_primed_kdm_10$percrypt_space %>% bind_rows(data_primed_stag_10$percrypt_space) %>%
   summarise(crypt_space_scale = median(space_per_crypt)) %>% pull()
length_scale = 2 * sqrt(domain_size/pi) # domain diameter
area_scale = length_scale*length_scale

stan_ins_ref = generate_stan_ins4(data_primed_stag_10, data_primed_kdm_10, length_scale, area_scale, slide_age_ref, n_grid = 25)
stan_ins_ref$slide_index_tib = stan_ins_ref$slide_index_tib %>% left_join(slide_age_ref) %>% rename(Slide_ID = Scan_ID)


#######################################
RUN = FALSE
if (RUN==TRUE) {
   ## load paths 
   event_time_list = readRDS("./output/paths_sampled_for_all_neighbourhoods.rds")
   
   ## for random samples of the paths from each list element (num neighbourhoods long),
   ## create stan inputs and infer the parameters D & gamma_a
   cmp_model = rstan::stan_model(file = "./stan_models/whitespace_diffusion_model_dripfeed_sampledpaths.stan")
   globiters = 1000
   for (i in 1:globiters) {
      ## sample paths
      evinds = sample(1:nrow(event_time_list[[1]]), length(event_time_list), replace = T)
      event_times = matrix(NA, length(event_time_list), ncol(event_time_list[[1]]))
      for (j in 1:length(evinds)) {
         event_times[j,] = event_time_list[[j]][evinds[j],]
      }
      ## stan inputs
      stan_ins_10s = generate_stan_ins9(data_primed_stag_10, data_primed_kdm_10, length_scale,
                                        area_scale, slide_age_ref, event_times, n_grid = 25)
      ## run sampler
      nchains = 2
      options(mc.cores = nchains)
      iters = 1000
      samp_fit_10s = rstan::sampling(object = cmp_model, data = stan_ins_10s$stan_ins$input, cores=nchains,
                                     init = rep(list(stan_ins_10s$stan_ins$init), nchains), chains=nchains,
                                     iter=iters, thin = 2)
      ## save reduced output for iteration
      fits = rstan::extract(samp_fit_10s)
      fewfits = list()
      fewfits$D = fits$D
      fewfits$gamma_a = fits$gamma_a
      fewfits$a1_gpop = fits$a1_gpop
      fewfits$a2_gpop = fits$a2_gpop
      fewfits$sd_Gacc = fits$sd_Gacc
      fewfits %>% saveRDS(paste0("./stan_fits/diffusion_inference_iter_", paste0(i, ".rds")))
      event_times %>% saveRDS(paste0("./stan_fits/sampled_paths_inference_iter_", paste0(i, ".rds")))
   }
}

# load and combine inference iterations to create master inference
maxiter = length(dir("./stan_fits/", pattern = "diffusion_inference_iter"))
combfits = list()
combfits$D = tibble()
combfits$gamma_a = tibble()
combfits$gamma_a_pop = tibble()
combfits$event_times = tibble()
for (i in (1:maxiter)) {
   fewfits = readRDS(paste0("./stan_fits/diffusion_inference_iter_", paste0(i, ".rds")))
   event_times = readRDS(paste0("./stan_fits/sampled_paths_inference_iter_", paste0(i, ".rds")))
   combfits$D = combfits$D %>% bind_rows(tibble(D = fewfits$D) %>% mutate(iter = 1:n(), run = i))
   combfits$gamma_a_pop = combfits$gamma_a_pop %>% 
      bind_rows(tibble(gamma_a_pop = rbeta(length(fewfits$a1_gpop), fewfits$a1_gpop, fewfits$a2_gpop)) %>% mutate(iter = 1:n(), run = i))
   combfits$gamma_a = combfits$gamma_a %>% 
      bind_rows(fewfits$gamma_a %>% as_tibble() %>% mutate(iter = 1:n()) %>% magrittr::set_colnames(c(1:ncol(fewfits$gamma_a), "iter")) %>% 
                   gather(nbrhood_index, gamma_a, -iter) %>% mutate(nbrhood_index = as.numeric(nbrhood_index), run = i))
   combfits$event_times = combfits$event_times %>% 
      bind_rows(event_times %>% as_tibble() %>% mutate(nbrhood_index = 1:n()) %>% magrittr::set_colnames(c(1:ncol(event_times), "nbrhood_index")) %>% 
                   gather(event_num, event_time, -nbrhood_index) %>% mutate(event_num = as.numeric(event_num), run = i))
}
combfits$D = combfits$D %>% mutate(iter_glob = 1:n())
combfits$gamma_a_pop = combfits$gamma_a_pop %>% mutate(iter_glob = 1:n())
combfits$gamma_a = combfits$gamma_a %>% group_by(nbrhood_index) %>% mutate(iter_glob = 1:n())

combfits$D %>% ggplot() + geom_histogram(aes(D), bins=60)# + xlim(0,10)
combfits$D %>% filter(run<20) %>% ggplot() + geom_line(aes(iter,D)) + facet_wrap(~run, scales="free")
combfits$gamma_a_pop %>% ggplot() + geom_histogram(aes(gamma_a_pop), bins=60)
combfits$gamma_a %>% filter(nbrhood_index<10) %>% ggplot() + geom_histogram(aes(gamma_a), bins=60) + facet_wrap(~nbrhood_index)

combfits$D %>% ggplot() + geom_boxplot(aes(x=1, y=D), width=0.5) + xlim(0,2)
combfits$gamma_a %>% ggplot() + geom_boxplot(aes(nbrhood_index, gamma_a, group=nbrhood_index))

combfits$D %>% ggplot() + geom_density(aes(D)) + xlim(0,10)

## summarise
summary_fits = list()
summary_fits$D = combfits$D %>% summarise(mean_mcmc = mean(D), med_mcmc = median(D), 
                         l95_mcmc = quantile(D, 0.025), u95_mcmc = quantile(D, 0.975))
summary_fits$gamma_a_pop = combfits$gamma_a_pop %>% summarise(mean_mcmc = mean(gamma_a_pop), med_mcmc = median(gamma_a_pop), 
                                          l95_mcmc = quantile(gamma_a_pop, 0.025), u95_mcmc = quantile(gamma_a_pop, 0.975))
summary_fits$gamma_a = combfits$gamma_a %>% group_by(nbrhood_index) %>% 
   summarise(mean_mcmc = mean(gamma_a), med_mcmc = median(gamma_a),
             l95_mcmc = quantile(gamma_a, 0.025), u95_mcmc = quantile(gamma_a, 0.975))
summary_fits$event_times = combfits$event_times %>% group_by(nbrhood_index, event_num) %>% 
   summarise(mean_mcmc = mean(event_time), med_mcmc = median(event_time),
             l95_mcmc = quantile(event_time, 0.025), u95_mcmc = quantile(event_time, 0.975))

## plot summaries
summary_fits$event_times %>% left_join(stan_ins_ref$slide_index_tib) %>% 
   ggplot() + geom_line(aes(event_num, med_mcmc, col=factor(nbrhood_index))) + 
   geom_ribbon(aes(x=event_num, ymin=l95_mcmc, ymax=u95_mcmc, fill=factor(nbrhood_index)), alpha=0.2) +
   geom_hline(aes(yintercept = Age, col=factor(nbrhood_index)), linetype = "dashed") +
   facet_wrap(~nbrhood_index)

summary_fits$event_times %>% left_join(stan_ins_ref$slide_index_tib) %>% 
   ggplot() + geom_line(aes(event_num, med_mcmc, col=factor(mark))) + 
   geom_ribbon(aes(x=event_num, ymin=l95_mcmc, ymax=u95_mcmc, fill=factor(mark)), alpha=0.2) +
   geom_hline(aes(yintercept = Age, col=factor(mark)), linetype = "dashed") +
   facet_wrap(~nbrhood_index)

## get r-dependence predictions for each mutant patch
# new drip feed
library(foreach)
library(doMC)
doMC::registerDoMC(cores = 10)
patient_patch_ref = stan_ins_ref$slide_index_tib %>% arrange(nbrhood_index) %>% left_join(slide_age_ref) %>% 
   left_join(summary_fits$event_times %>% filter(event_num==1)) %>% mutate(pa_mean = Age - mean_mcmc, pa_median = Age - med_mcmc,
                                                                           pa_u95 = Age - l95_mcmc, pa_l95 = Age - u95_mcmc) %>% 
   select(-(event_num:u95_mcmc))
location_by_patch = prepare_r_dep_data_for_plotting(data_primed_kdm_10, data_primed_stag_10, stan_ins_ref)
location_by_patch = location_by_patch %>% mutate(newtype = ifelse(Type=="mut", "Mutant", "Adjacent"))
patch_max_R = location_by_patch %>% group_by(nbrhood_index) %>% summarise(max_R = max(R_mean))
useiters = sample(1:nrow(combfits$D), size = 5000, replace = FALSE)
ref_tib = stan_ins_ref$joinscaledat %>% left_join(slide_age_ref) %>% select(nbrhood_index, a_m, a_wt, Age) %>% distinct()
RUN = FALSE
if (RUN==TRUE) {
   all_nbrhood_runs = tibble()
   for (j in (stan_ins_ref$slide_index_tib %>% pull(nbrhood_index))) {
      thesepaths = tibble(run=combfits$D$run[useiters], iter=useiters) %>% left_join(combfits$event_times %>% filter(nbrhood_index==j))
      thesegamma = tibble(iter_glob=useiters) %>% left_join(combfits$gamma_a %>% filter(nbrhood_index==j))
      thisref = ref_tib %>% filter(nbrhood_index==j)
      maxr = max(10, patch_max_R %>% filter(nbrhood_index==j) %>% pull(max_R))
      rvec = pracma::linspace(0.1, maxr, as.integer(10 * maxr))
      sampleres = foreach (i = useiters)%dopar%{
         get_timeseries_from_paths_real_D(i, thesepaths, thisref$Age, combfits$D$D[i], thisref$a_m, thisref$a_wt, rvec) %>% mutate(iter_glob=i)
      }
      dft_dc = do.call(rbind, sampleres) %>% left_join(thisref %>% rename(t=Age)) %>% left_join(thesegamma)
      dft_dc = dft_dc %>% mutate(rho_r = 1/gamma_a + rho_pert, gamma_r = 1/rho_r) %>% group_by(t, r, nbrhood_index) %>% 
         summarise(med_mcmc = median(gamma_r), l95_mcmc = quantile(gamma_r, 0.025), 
                   u95_mcmc = quantile(gamma_r, 0.975)) %>% ungroup() %>% left_join(thisref %>% rename(t=Age))
      all_nbrhood_runs = all_nbrhood_runs %>% bind_rows(dft_dc)
   }
   all_nbrhood_runs %>% saveRDS("./output/all_nbrhood_runs.rds")   
} else {
   all_nbrhood_runs = readRDS("./output/all_nbrhood_runs.rds")   
}

pl_list_10s = compare_r_dependence_to_data(data_primed_kdm_10, data_primed_stag_10, 
                                           stan_ins_ref, all_nbrhood_runs, WT30_r = NA)
pl_list_10s[[1]]
pl_list_10s[[2]]
pl_list_10s[[3]]
pl_list_10s[[4]]   
pl_list_10s[[5]]
pl_list_10s[[6]]
pl_list_10s[[7]]
pl_list_10s[[8]]


## checking most likely paths given fixed median diffusion coefficient
cmp_model = rstan::stan_model(file = "./stan_models/whitespace_diffusion_model_dripfeed_sampledpaths.stan")
samp_fit_10s = readRDS("./stan_fits/diffusion_model_fit_10ps_hierarchical_normalage.rds")
fits_old = rstan::extract(samp_fit_10s)
rstan::expose_stan_functions(cmp_model)
gamma_a = summary_fits$gamma_a$med_mcmc
D = summary_fits$D$med_mcmc
G_sd = median(fits_old$sd_Gacc)
event_time_list = readRDS("./output/paths_sampled_for_all_neighbourhoods.rds")
likelihood_store = matrix(NA, nrow = nrow(event_time_list[[1]]), ncol = length(data_primed_stag_10$data_primed$Gamma_dat) + length(data_primed_kdm_10$data_primed$Gamma_dat))
for (i in 1:nrow(event_time_list[[1]])) {
   event_times = matrix(NA, length(event_time_list), ncol(event_time_list[[1]]))
   event_times_l = list()
   for (j in 1:length(event_time_list)) {
      event_times[j,] = event_time_list[[j]][i,]
      event_times_l[[j]] = event_time_list[[j]][i,]
   }
   pins = generate_stan_ins9(data_primed_stag_10, data_primed_kdm_10, length_scale,
                             area_scale, slide_age_ref, event_times, n_grid = 25)
   Gam = Gamma_all(gamma_a, D, pins$stan_ins$input$am, pins$stan_ins$input$awt, pins$stan_ins$input$N, 
                   event_times_l, pins$stan_ins$input$pat_ages,
                   pins$stan_ins$input$Rin, pins$stan_ins$input$Rout, pins$stan_ins$input$theta,
                   pins$stan_ins$input$nbrhood_index, pins$stan_ins$input$n_grid, pins$stan_ins$input$num_events)
   likelihood_store[i,] = dnorm(pins$stan_ins$input$Gamma_dat, mean = Gam, sd = G_sd, log = TRUE)
}

# sum likelihoods for neighbourhoods
lik_store_nbrhoods = matrix(NA, nrow(likelihood_store), length(gamma_a))
for (i in 1:length(gamma_a)) {
   lik_store_nbrhoods[,i] = rowSums(likelihood_store[,which(pins$stan_ins$input$nbrhood_index==i)])
}
loglik_ordered_paths_index = lik_store_nbrhoods
for (i in 1:ncol(lik_store_nbrhoods)) loglik_ordered_paths_index[,i] = order(lik_store_nbrhoods[,i], decreasing = T)

all_nbrhood_runs_bestpaths = tibble()
all_nbrhood_pathinfo = tibble()
numbest = 25
for (j in (stan_ins_ref$slide_index_tib %>% pull(nbrhood_index))) {
   thisref = ref_tib %>% filter(nbrhood_index==j)
   thisgamma = gamma_a[j]
   maxr = max(10, patch_max_R %>% filter(nbrhood_index==j) %>% pull(max_R))
   rvec = pracma::linspace(0.1, maxr, as.integer(10 * maxr))
   sampleres = tibble()
   for (i in 1:numbest) {
      path = event_time_list[[j]][loglik_ordered_paths_index[i,j],]
      sampleres = sampleres %>% bind_rows(get_timeseries_from_single_path_real_D(path, thisref$Age, D, thisref$a_m, thisref$a_wt, rvec) %>% 
                                             mutate(index = i, patchage = thisref$Age - path[1], timesincelastfission = thisref$Age - tail(path, 1)))
   }
   dft_dc = sampleres %>% left_join(thisref %>% rename(t=Age))
   all_nbrhood_pathinfo = all_nbrhood_pathinfo %>% bind_rows(dft_dc %>% select(nbrhood_index, t, index, patchage, timesincelastfission) %>% distinct())
   dft_dc = dft_dc %>% mutate(rho_r = 1/thisgamma + rho_pert, gamma_r = 1/rho_r) %>% group_by(t, r, nbrhood_index) %>% 
      summarise(med_mcmc = median(gamma_r), l95_mcmc = quantile(gamma_r, 0.025), 
                u95_mcmc = quantile(gamma_r, 0.975)) %>% ungroup() %>% left_join(thisref %>% rename(t=Age))
   all_nbrhood_runs_bestpaths = all_nbrhood_runs_bestpaths %>% bind_rows(dft_dc)
}

pl_list_10s = compare_r_dependence_to_data(data_primed_kdm_10, data_primed_stag_10, 
                                           stan_ins_ref, all_nbrhood_runs_bestpaths, WT30_r = NA)
pl_list_10s[[1]]
pl_list_10s[[2]]
pl_list_10s[[3]]
pl_list_10s[[4]]   
pl_list_10s[[5]]
pl_list_10s[[6]]
pl_list_10s[[7]]
pl_list_10s[[8]]

all_nbrhood_pathinfo = all_nbrhood_pathinfo %>% left_join(location_by_patch %>% select(Scan_ID, mark, nbrhood_index) %>% distinct()) %>% 
   rename(Age = t)
all_nbrhood_pathinfo %>% ggplot() + geom_boxplot(aes(x=nbrhood_index, y=patchage/Age, group=nbrhood_index, fill=mark))
all_nbrhood_pathinfo %>% ggplot() + geom_boxplot(aes(x=nbrhood_index, y=timesincelastfission, group=nbrhood_index, fill=mark))
patch_time_info = all_nbrhood_pathinfo %>% group_by(nbrhood_index, mark, Age) %>% 
   summarise(med_pa = median(patchage), med_tslf = median(timesincelastfission))

## plotting with Cora's specifications
nbrs_i = location_by_patch %>% pull(nbrhood_index) %>% unique() %>% sort()
for (i in 1:length(nbrs_i)) {
   nbr_i = nbrs_i[i]
   pat_info = location_by_patch %>% filter(nbrhood_index==nbr_i) %>% select(Scan_ID, mark, nbrhood_index) %>% 
      left_join(patient_patch_ref) %>% distinct()
   thisma = location_by_patch %>% filter(nbrhood_index==nbr_i) %>% pull(mark) %>% unique()
   maxr = max(10, location_by_patch %>% filter(nbrhood_index==nbr_i) %>% pull(R_mean) %>% max())
   pl = location_by_patch %>% filter(nbrhood_index==nbr_i) %>%  ggplot() + 
      geom_ribbon(data=all_nbrhood_runs %>% filter(nbrhood_index==nbr_i, r<maxr), 
                  aes(x=r, ymin=l95_mcmc, ymax=u95_mcmc, fill= "Theory" ), alpha=0.2) +
      geom_line(data=all_nbrhood_runs %>% filter(nbrhood_index==nbr_i, r<maxr), 
                aes(x=r, y=med_mcmc, col = "Theory")) + 
      geom_line(data=all_nbrhood_runs_bestpaths %>% filter(nbrhood_index==nbr_i, r<maxr), 
                aes(x=r, y=med_mcmc, col = "Theory"), linetype="dashed") + 
      geom_point(aes(x=R_mean, y=gamma, col=newtype)) + 
      ylab("Stromal fraction") + xlab("r") + theme_bw() + 
      theme(plot.background = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), panel.border = element_blank(),
            axis.line = element_line(color = 'black'),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            axis.text.x= element_text(size = 8),
            axis.text.y= element_text(size = 8))
   if (thisma=="KDM6A") {
      pl = pl + scale_color_manual(aesthetics = c("colour", "fill"), values = c("aquamarine2", "aquamarine4", "black")) +
         guides(fill=FALSE, color=FALSE) + labs(color="")
   } else if (thisma=="STAG2") {
      pl = pl + scale_color_manual(aesthetics = c("colour", "fill"), values = c("tan2", "tan4", "black")) +
         guides(fill=FALSE, color=FALSE) + labs(color="")
   }
   pinfp = patch_time_info %>% slice(i)
   pat_info_string = pat_info %>% 
      mutate(info_string = paste0(DW.Patient.number, 
                                  paste0(paste0("_age", Age), 
                                         paste0("_patchage", 
                                                paste0(round(pinfp$med_pa)),
                                                paste0("_sincelastfission", round(pinfp$med_tslf, digits = 1)))))) %>% pull(info_string)
   ggsave(paste0("./plots/r_dep_gamma_3x2_nbrhood_", 
                 paste(paste(as.character(nbr_i), thisma, sep = "_"), 
                       paste0(pat_info_string, ".pdf"), sep = "_")),
          plot = pl, device = "pdf", width = 6, height = 4, units = "cm")
   ggsave(paste0("./plots/r_dep_gamma_3x1-6_nbrhood_", 
                 paste(paste(as.character(nbr_i), thisma, sep = "_"), 
                       paste0(pat_info_string, ".pdf"), sep = "_")),
          plot = pl, device = "pdf", width = 6, height = 3.2, units = "cm")
}


patchage_dists_stag = readRDS(paste0("./output/patch_age_priorbetafits_", paste0("STAG2", ".rds")))
patchage_dists_kdm = readRDS(paste0("./output/patch_age_priorbetafits_", paste0("KDM6A", ".rds")))
patchage_dists = patchage_dists_stag %>% mutate(mark = "STAG2") %>% 
   bind_rows(patchage_dists_kdm %>% mutate(mark="KDM6A")) %>% rename(Age=age)

patch_time_info = patch_time_info %>% mutate(pa_scaled = med_pa / Age) %>% left_join(patchage_dists) %>% ungroup() %>% 
   left_join(stan_ins_10s$joinscaledat %>% select(Scan_ID, mark, neighbourhood, nbrhood_index, Age) %>% distinct())
bgrid = pracma::linspace(0,1,1000)
bsamps = 100000
allpatch_age_fits = tibble()
for (i in nbrs_i) {
   thispatch = patch_time_info %>% filter(nbrhood_index==i)
   bf = tibble(betafit = rbeta(bsamps, thispatch$alpha, thispatch$beta), nbrhood_index=i,
               neighbourhood=thispatch$neighbourhood, mark=thispatch$mark, 
               Age=thispatch$Age, mostlikely=thispatch$pa_scaled,
               subtitle = paste0("Patch ", paste(thispatch$neighbourhood, paste0(thispatch$Age, " yr"), sep = ", ")))
   allpatch_age_fits = allpatch_age_fits %>% bind_rows(bf)
}
allpatch_age_fits = allpatch_age_fits %>% mutate(subtitle_f = factor(subtitle, levels = (allpatch_age_fits %>% pull(subtitle) %>% unique())))

allpatch_age_fits %>% filter(mark=="STAG2") %>% mutate(mark = "Population distribution") %>% 
   ggplot(aes()) + 
   geom_density(aes(betafit*Age, y=..density.., fill=mark), alpha=0.5) + xlim(0,100) +
   geom_vline(aes(xintercept=mostlikely*Age), linetype="dashed") + facet_wrap(~subtitle_f, ncol=4) +
   theme_bw() + xlab("Patch age") +
   theme(plot.background = element_blank(), panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), panel.border = element_blank(),
         axis.line = element_line(color = 'black'),
         axis.title.x = element_text(size = 8), axis.text.x= element_text(size = 8),
         legend.position = "bottom") + labs(fill="") +
   scale_color_manual(aesthetics = c("colour", "fill"), values = c("tan2", "tan4", "black"))
ggsave("./plots/patch_age_distributions_with_mostlikely_STAG2.pdf",
       device = "pdf", width = 15, height = 21, units = "cm")
allpatch_age_fits %>% filter(mark=="KDM6A") %>% mutate(mark = "Population distribution") %>% 
   ggplot(aes()) + 
   geom_density(aes(betafit*Age, y=..density.., fill=mark), alpha=0.5) + xlim(0,100) +
   geom_vline(aes(xintercept=mostlikely*Age), linetype="dashed") + facet_wrap(~subtitle_f, ncol=4) +
   theme_bw() + xlab("Patch age") +
   theme(plot.background = element_blank(), panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), panel.border = element_blank(),
         axis.line = element_line(color = 'black'),
         axis.title.x = element_text(size = 8), axis.text.x= element_text(size = 8),
         legend.position = "bottom") + labs(fill="") +
   scale_color_manual(aesthetics = c("colour", "fill"), values = c("aquamarine2", "aquamarine4", "black"))
ggsave("./plots/patch_age_distributions_with_mostlikely_KDM6A.pdf",
       device = "pdf", width = 15, height = 18, units = "cm")


###################################################
### numerical experiments 
mark_data = read_csv("../drift_All_Marks/exp_coeff.csv")

# for checking inputs
samp_fit_10s = readRDS("./stan_fits/diffusion_model_fit_10ps_hierarchical_normalage.rds")
fits_10s = rstan::extract(samp_fit_10s)
fit_summaries_10s = extract_param_summaries(fits_10s)

###################### FIG 4 G, H
### Look at range using rho_max, gamma_min?  How many crypts are disturbed by a patch?
av_patch_inputs = stan_ins_ref$joinscaledat %>% select(Scan_ID, mark, neighbourhood, nbrhood_index, epsilon) %>% 
   distinct() %>% group_by(mark) %>% summarise(med_eps = median(epsilon))
gamma_a = summary_fits$gamma_a_pop$med_mcmc #median(fits_10s$gamma_a_pop)
## just define range as the "new homeostasis" assuming diffusion hits some wall
RUN = FALSE
if (RUN==TRUE) {
   for (thresh in c(1.01, 1.05)) {
      # using average mass input
      happy_packing_range = av_patch_inputs %>% mutate(R_H = sqrt(med_eps/(pi/gamma_a * (thresh-1))))
      # or calculate "errorbars"
      indiv_ins = stan_ins_ref$joinscaledat %>% select(Scan_ID, mark, neighbourhood, nbrhood_index, epsilon) %>% 
         distinct() %>% arrange(nbrhood_index)
      outtib = tibble()
      numiters = 10000
      useiters = sample(1:(combfits$gamma_a %>% pull(iter_glob) %>% max()), numiters, replace = FALSE)
      for (ii in 1:numiters) {
         i = useiters[ii]
         outtib = outtib %>% bind_rows(indiv_ins %>% left_join(combfits$gamma_a %>% filter(iter_glob==i) %>% select(-iter)) %>% 
                                          mutate(R_H = sqrt(epsilon/(pi/gamma_a * (thresh-1))), iter=ii))
      }
      outtib %>% group_by(mark) %>% summarise(l95_mcmc = quantile(R_H, 0.025), 
                                              med_mcmc = median(R_H), 
                                              u95_mcmc = quantile(R_H, 0.975)) %>% 
         ggplot() + geom_point(aes(x=mark, y=med_mcmc, col=mark)) +
         geom_errorbar(aes(x=mark, ymin=l95_mcmc, ymax=u95_mcmc, col=mark), width=0.4)
      
      outtib = outtib %>% mutate(N_crypt_footprint = calc_happy_packing_footprint(R_H, stan_ins_ref$joinscaledat %>% 
                                                                                     select(mark, nbrhood_index, s_wt) %>% 
                                                                                     distinct() %>% pull(s_wt) %>% median()))
      outtib %>% ggplot() + geom_boxplot(aes(x=mark, y=R_H, fill=mark))
      outtib %>% ggplot() + geom_boxplot(aes(x=mark, y=N_crypt_footprint, fill=mark))
      outtib %>% saveRDS(paste0("./output/happy_packing_range_10patch_threshold_", 
                                paste0(as.character((thresh - 1)*100), ".rds")))
   }
   rm(outtib); gc()
   range_tib_1pc = readRDS(paste0("./output/happy_packing_range_10patch_threshold_1.rds"))
   range_tib_5pc = readRDS(paste0("./output/happy_packing_range_10patch_threshold_5.rds"))
} else {
   range_tib_1pc = readRDS(paste0("./output/happy_packing_range_10patch_threshold_1.rds"))
   range_tib_5pc = readRDS(paste0("./output/happy_packing_range_10patch_threshold_5.rds"))
}


###################### FIG 4 J
## finding threshold of fission rate for given diffusion rate
marklist = c("MAOA", "STAG2", "KDM6A", "KRAS")
mark_fisrates = mark_data %>% filter(mark%in%marklist) %>% select(mark, rho)
wt_rate = mark_fisrates %>% filter(mark=="MAOA") %>% pull(rho)
mark_fisrates = mark_fisrates %>% mutate(scaled_rho = rho / wt_rate)
mark_fisrates = mark_fisrates %>% bind_rows(tibble(mark = "old_KRAS", rho = 10*wt_rate, scaled_rho = 10))

wtmult_list = 1:20
fis_rate_list = wt_rate * wtmult_list
fis_mult_map = tibble(rho_fi = fis_rate_list, wt_mult = wtmult_list)
rvec = pracma::linspace(0,10,101)
tvec = pracma::linspace(1,50,50)
D = summary_fits$D$med_mcmc #median(fits_10s$D)
gamma_a = summary_fits$gamma_a_pop$med_mcmc #median(fits_10s$gamma_a_pop)
space_tib = stan_ins_ref$joinscaledat %>% select(mark, nbrhood_index, a_m, a_wt, s_wt) %>% distinct()
a_wt = space_tib %>% pull(a_wt) %>% median()
s_wt = space_tib %>% pull(s_wt) %>% median()
a_m = a_wt
iters = 5000
summary_tib = tibble()
patchsize_tib = tibble()
RUN = FALSE
if (RUN==TRUE) {
   for (rho_fi in fis_rate_list) {
      res = iterate_diffusion_time_series(tvec, D, a_m, a_wt, rvec, rho_fi, iters=iters)
      res = res %>% mutate(gamma_r = 1/(1/gamma_a + rho_pert), rho_r = 1/gamma_a + rho_pert)
      sumres_gamma = res %>% group_by(t,r) %>% summarise(med_mcmc = median(gamma_r), l95_mcmc = quantile(gamma_r, 0.025), 
                                                         u95_mcmc = quantile(gamma_r, 0.975)) %>% ungroup()
      sumres_patchsize = res %>% select(t, patchsize, iter) %>% distinct() %>% group_by(t) %>% 
         summarise(med_ps = median(patchsize), l95_ps = quantile(patchsize, 0.025), u95_ps = quantile(patchsize, 0.975)) %>% 
         ungroup()
      summary_tib = summary_tib %>% bind_rows(sumres_gamma %>% mutate(rho_fi = rho_fi))
      patchsize_tib = patchsize_tib %>% bind_rows(sumres_patchsize %>% mutate(rho_fi = rho_fi))
   }
   summary_tib %>% left_join(fis_mult_map) %>% saveRDS("./output/scan_rhofi_for_polyp_threshold_gamma.rds")
   patchsize_tib %>% left_join(fis_mult_map) %>% saveRDS("./output/scan_rhofi_for_polyp_threshold_patchsize.rds")
} else {
   summary_tib = readRDS("./output/scan_rhofi_for_polyp_threshold_gamma.rds")
   patchsize_tib = readRDS("./output/scan_rhofi_for_polyp_threshold_patchsize.rds")
}

## define max/min values from optimal hexagonal packing of equal circles
gamma_min = 1 - pi/sqrt(12)
rho_max = 1/gamma_min

summary_tib %>% 
   ggplot() + geom_line(aes(x=r, y=med_mcmc, col=t, group=t)) +
   geom_ribbon(aes(x=r, ymin=l95_mcmc, ymax=u95_mcmc), alpha=0.2) + 
   geom_hline(aes(yintercept = gamma_min), linetype="dotted") + facet_wrap(~wt_mult)
patchsize_tib %>%
   ggplot() + geom_point(aes(x=t, y=med_ps)) + 
   geom_errorbar(aes(x=t, ymin=l95_ps, ymax=u95_ps), width=0.25) + facet_wrap(~wt_mult) +
   scale_y_log10()

summary_tib %>% 
   ggplot() + geom_line(aes(x=r, y=med_mcmc, col=t, group=t)) +
   geom_ribbon(data = (summary_tib %>% filter(t==50)), aes(x=r, ymin=l95_mcmc, ymax=u95_mcmc), alpha=0.3) + 
   geom_hline(aes(yintercept = gamma_min), linetype="dotted") + facet_wrap(~wt_mult)


######################### FIG 4 K
## If n KRAS mutations at age X, how many will become a polyp in T years?
# rho_fi = mark_fisrates %>% filter(mark=="KRAS") %>% pull(rho) # use new KRAS fission rate of ~17xWT
rho_fi = 17 * wt_rate ## rounding up KRAS mult from 16.5 to nearest integer
iters = 20000
roll_area_smoothed = c(pi*mean(rvec[1:2])^2)
for (i in 2:(length(rvec)-1)) roll_area_smoothed = c(roll_area_smoothed, pi*mean(rvec[(i-1):(i+1)])^2)
roll_area_smoothed = c(roll_area_smoothed, pi*mean(tail(rvec, 2))^2)
rollarea_tib = tibble(r = rvec, roll_area = roll_area_smoothed)
suptimes_list = list()
suptimes_list[[1]] = c(0, 10)
suptimes_list[[2]] = c(10, 20)
suptimes_list[[3]] = c(20, 30)
RUN = TRUE
if (RUN==TRUE) {
   for (st in 1:length(suptimes_list)) {
      suptimes = suptimes_list[[st]]
      suptime_string = paste(as.character(suptimes[1]), as.character(suptimes[2]), sep="-")
      for (suprate in c(0, wt_rate)) {
         res_kras = iterate_diffusion_time_series(tvec, D, a_m, a_wt, rvec, rho_fi, iters=iters, 
                                                  suppressed_rate = suprate, suppression_time = suptimes)
         res_kras = res_kras %>% mutate(gamma_r = 1/(1/gamma_a + rho_pert), rho_r = 1/gamma_a + rho_pert) %>% 
            group_by(t, iter) %>% mutate(Rho_integral = 2*pi*pracma::cumtrapz(r, (1-gamma_r)*r)) %>% 
            group_by(t, iter) %>% mutate(pnt = 1:n(), c_gam = cumsum(gamma_r), roll_av_gam = c_gam / pnt)
         patch_inners = res_kras %>% mutate(ps_am = patchsize * a_m) %>% filter(Rho_integral <= ps_am)
         patch_inners = patch_inners %>% mutate(poly_zone = roll_av_gam < gamma_min) %>% 
            ungroup() %>% left_join(rollarea_tib) %>% 
            group_by(t, iter) %>% mutate(max_area = pi * max(r)^2)
         zone_indicators = patch_inners %>% summarise(inzone = ifelse(sum(poly_zone)>0, 1, 0)) %>% ungroup()
         no_in_zone = zone_indicators %>% filter(inzone == 0)
         in_zone = zone_indicators %>% filter(inzone == 1) %>% left_join(patch_inners %>% ungroup()) %>% 
            group_by(t, iter) %>% filter(poly_zone) %>% summarise(frac_in_zone = max(roll_area) / max_area[1])
         kras_polyp_summedfracs = in_zone %>% group_by(t) %>% summarise(sum_chance = sum(frac_in_zone)/iters)
         kras_polyp_summedfracs = tibble(t = tvec) %>% left_join(kras_polyp_summedfracs) %>% 
            mutate(sum_chance = ifelse(is.na(sum_chance), 0, sum_chance)) %>% 
            mutate(suppressed_rate = suprate, suppression_times = suptime_string)
         kras_polyp_summedfracs %>% saveRDS(paste0("./output/kras17xWT_over_time_forfracpolyp_suppressedrate_", 
                                                   paste0(as.character(round(suprate, 3)), 
                                                          paste0("_suppressedtime", 
                                                                 paste0(suptime_string, ".rds")))))
         rm(res_kras); gc()
      }
   }
}
if (RUN==TRUE) {
   res_kras = iterate_diffusion_time_series(tvec, D, a_m, a_wt, rvec, rho_fi, iters=iters)
   res_kras = res_kras %>% mutate(gamma_r = 1/(1/gamma_a + rho_pert), rho_r = 1/gamma_a + rho_pert) %>% 
      group_by(t, iter) %>% mutate(Rho_integral = 2*pi*pracma::cumtrapz(r, (1-gamma_r)*r)) %>% 
      group_by(t, iter) %>% mutate(pnt = 1:n(), c_gam = cumsum(gamma_r), roll_av_gam = c_gam / pnt)
   patch_inners = res_kras %>% mutate(ps_am = patchsize * a_m) %>% filter(Rho_integral <= ps_am)
   patch_inners = patch_inners %>% mutate(poly_zone = roll_av_gam < gamma_min) %>% 
      ungroup() %>% left_join(rollarea_tib) %>% 
      group_by(t, iter) %>% mutate(max_area = pi * max(r)^2)
   zone_indicators = patch_inners %>% summarise(inzone = ifelse(sum(poly_zone)>0, 1, 0)) %>% ungroup()
   no_in_zone = zone_indicators %>% filter(inzone == 0)
   in_zone = zone_indicators %>% filter(inzone == 1) %>% left_join(patch_inners %>% ungroup()) %>% 
      group_by(t, iter) %>% filter(poly_zone) %>% summarise(frac_in_zone = max(roll_area) / max_area[1])
   kras_polyp_summedfracs = in_zone %>% group_by(t) %>% summarise(sum_chance = sum(frac_in_zone)/iters)
   kras_polyp_summedfracs = tibble(t = tvec) %>% left_join(kras_polyp_summedfracs) %>% 
      mutate(sum_chance = ifelse(is.na(sum_chance), 0, sum_chance)) %>% 
      mutate(suppressed_rate = NA, suppression_times = NA)
   kras_polyp_summedfracs %>% saveRDS("./output/kras17xWT_over_time_forfracpolyp_unsuppressed.rds")
   rm(res_kras); gc()
}

## could instead do as the fraction of the area containing patchsize * a_m that 
## has a rolling average below the threshold?
joinedtib = tibble()
for (st in 1:length(suptimes_list)) {
   suptimes = suptimes_list[[st]]
   suptime_string = paste(as.character(suptimes[1]), as.character(suptimes[2]), sep="-")
   for (suprate in c(0, wt_rate)) {
      kras_polyp_summedfracs = readRDS(paste0("./output/kras17xWT_over_time_forfracpolyp_suppressedrate_", 
                                              paste0(as.character(round(suprate, 3)), 
                                                     paste0("_suppressedtime", 
                                                            paste0(suptime_string, ".rds")))))
      joinedtib = joinedtib %>% bind_rows(kras_polyp_summedfracs)
   }
}
mut_tib = readRDS("./output/kras17xWT_over_time_forfracpolyp_unsuppressed.rds") %>% select(t, sum_chance)

joinedtib %>% 
   ggplot() + geom_line(aes(x=t, y=sum_chance, col=suppression_times)) + 
   facet_wrap(~suppressed_rate) + 
   geom_line(data = mut_tib, aes(x=t, y=sum_chance, col="KRAS"))


####################################
## fit the null model and then compare likelihoods of each neighbourhood under it

## create stan inputs and infer the parameters D & gamma_a
cmp_null_model = rstan::stan_model(file = "./stan_models/whitespace_null_model.stan")
## stan inputs
stanins_null = generate_stan_ins_null(data_primed_stag_10, data_primed_kdm_10, length_scale,
                                      area_scale, slide_age_ref, n_grid = 25)

## run sampler
nchains = 2
options(mc.cores = nchains)
iters = 4000
samp_fit_null = rstan::sampling(object = cmp_null_model, data = stanins_null$stan_ins$input, cores=nchains,
                                init = rep(list(stanins_null$stan_ins$init), nchains), chains=nchains,
                                iter=iters, thin = 2)
rstan::traceplot(samp_fit_null)
samp_fit_null %>% saveRDS("./stan_fits/null_model_fits.rds")

fits_null = rstan::extract(samp_fit_null)
gamma_a_null = fits_null$gamma_a %>% magrittr::set_colnames(1:ncol(fits_null$gamma_a)) %>% as_tibble() %>% mutate(iter = 1:n()) %>% 
   gather(nbrhood_index, val, -iter) %>% mutate(nbrhood_index = as.numeric(nbrhood_index)) %>% group_by(nbrhood_index) %>% 
   summarise(med_mcmc = median(val), l95_mcmc = quantile(val, 0.025), u95_mcmc = quantile(val, 0.975)) %>% pull(med_mcmc)
rstan::expose_stan_functions(cmp_null_model)
G_sd = median(fits_old$sd_Gacc) # use same sd as for diffusion model
pins = stanins_null
Gam = Gamma_all(gamma_a_null, pins$stan_ins$input$N,
                pins$stan_ins$input$Rin, pins$stan_ins$input$Rout, pins$stan_ins$input$theta,
                pins$stan_ins$input$nbrhood_index, pins$stan_ins$input$n_grid)
likelihood_store_null = dnorm(pins$stan_ins$input$Gamma_dat, mean = Gam, sd = G_sd, log = TRUE)

# sum likelihoods for neighbourhoods
lik_store_nbrhoods_null = rep(NA, length(gamma_a_null))
lik_store_diffusion = rep(NA, length(gamma_a_null))
for (i in 1:length(gamma_a)) {
   lik_store_nbrhoods_null[i] = sum(likelihood_store_null[which(pins$stan_ins$input$nbrhood_index==i)])
   lik_store_diffusion[i] = max(lik_store_nbrhoods[,i])
}

diff_loglik_paths = lik_store_nbrhoods %>%  magrittr::set_colnames(1:ncol(lik_store_nbrhoods)) %>% as_tibble() %>% 
   mutate(path = 1:n()) %>% gather(nbrhood_index, loglik, -path) %>% mutate(nbrhood_index = as.numeric(nbrhood_index)) %>% 
   left_join(tibble(null_loglik = lik_store_nbrhoods_null) %>% mutate(nbrhood_index = 1:n()))

diff_loglik_paths %>% ggplot() + geom_histogram(aes(loglik)) + facet_wrap(~nbrhood_index, scales="free") +
   geom_vline(aes(xintercept=null_loglik), linetype="dashed")

diff_loglik_paths %>% group_by(path) %>% summarise(loglik = sum(loglik), null_loglik = sum(null_loglik)) %>% 
   ggplot() + geom_histogram(aes(loglik, y=..density.., fill="Diffusion model"), alpha=0.8) + 
   geom_vline(aes(xintercept=null_loglik, linetype="Null model")) +
   theme_bw() + labs(linetype="", fill="") + xlab("Log likelihood") +
   theme(axis.title.x = element_text(size = 8), axis.text.x= element_text(size = 8),
         axis.title.y = element_text(size = 8), axis.text.y= element_text(size = 8),
         legend.position = "bottom", legend.text = element_text(size = 8))
ggsave("./plots/null_model_loglik.pdf",
       device = "pdf", width = 15, height = 14, units = "cm")
