library(tidyverse)
library(Rcpp)
source("./funcs_crypt_packing.r")
source("./joining_slide_to_block_to_patient.r")

slide_age_ref = block_to_slide_stag %>% mutate(mark="STAG2") %>% 
   bind_rows(block_to_slide_kdm %>% mutate(mark="KDM6A")) %>% left_join(block_to_patient) %>% 
   rename(Scan_ID = Slide_ID)

## load data
patchdat_10_kdm = load_and_tidy_data("KDM6A", 10)
# patchdat_5_kdm = load_and_tidy_data("KDM6A", 5)
patchdat_10_stag = load_and_tidy_data("STAG2", 10)
# patchdat_5_stag = load_and_tidy_data("STAG2", 5)

## load large WT patch data as baseline
wt_patches_kdm6a = readxl::read_xlsx("./databases/patches_30_measurements_KDM6A_290520.xlsx") %>% 
   magrittr::set_colnames(c("Name", "Class", "ROI", "centroid_x", "centroid_y",
                            "area", "perimeter", "max_length", "Scan_ID"))
wt_patches_stag2 = readxl::read_xlsx("./databases/pa30_measurements_STAG2_290520.xlsx") %>% 
   magrittr::set_colnames(c("Name", "Class", "ROI", "centroid_x", "centroid_y",
                            "area", "perimeter", "max_length", "Scan_ID"))

## prime data
# data_primed_kdm_5 = prime_data(wt_patches_kdm6a, patchdat_5_kdm)
# data_primed_stag_5 = prime_data(wt_patches_stag2, patchdat_5_stag)
data_primed_kdm_10 = prime_data(wt_patches_kdm6a, patchdat_10_kdm)
data_primed_stag_10 = prime_data(wt_patches_stag2, patchdat_10_stag)

## define length and area scale
domain_size = data_primed_kdm_10$percrypt_space %>% bind_rows(data_primed_stag_10$percrypt_space) %>%
   summarise(crypt_space_scale = median(space_per_crypt)) %>% pull()
length_scale = 2 * sqrt(domain_size/pi) # domain diameter
area_scale = length_scale*length_scale


#######################################
## stan inputs
stan_ins_10s = generate_stan_ins8(data_primed_stag_10, data_primed_kdm_10, length_scale, 
                                  area_scale, slide_age_ref, n_grid = 25)
## compile and run
RUN = TRUE
if (RUN==TRUE) {
   
   cmp_model = rstan::stan_model(file = "./stan_models/whitespace_diffusion_model.stan")
   nchains = 4
   options(mc.cores = nchains)
   iters = 20000

   samp_fit_10s = rstan::sampling(object = cmp_model, data = stan_ins_10s$stan_ins$input, cores=nchains,
                              init = rep(list(stan_ins_10s$stan_ins$init), nchains), chains=nchains,
                              iter=iters, thin = 5)
   # rstan::traceplot(samp_fit_10s)
} else {
   samp_fit_10s = readRDS("./stan_fits/diffusion_model_fit_10ps_hierarchical_normalage.rds")
   stan_ins_10s = readRDS("./stan_fits/diffusion_model_stanins_10ps_hierarchical_normalage.rds")   
}
rstan::traceplot(samp_fit_10s, pars = c("D"))
rstan::traceplot(samp_fit_10s, pars = c("gamma_a"))
rstan::traceplot(samp_fit_10s, pars = c("a2_gpop"))
rstan::traceplot(samp_fit_10s, pars = c("gamma_a_pop"))
rstan::traceplot(samp_fit_10s, pars = c("event_times"))
pairs(samp_fit_10s, pars = c("gamma_a_pop", "D"))

fits_10s = rstan::extract(samp_fit_10s)
fit_summaries_10s = extract_param_summaries(fits_10s)
pairs(fits_10s$D %>% cbind(fits_10s$gamma_a[,1]) %>% cbind(fits_10s$event_times[,1,]))

fits_10s$event_times[,3,] %>% as_tibble() %>% mutate(iter=1:n()) %>%
   gather(event, time, -iter) %>% ggplot() + geom_line(aes(iter, time, col=event))

tibble(x=fits_10s$D) %>% ggplot() + geom_histogram(aes(x))


###################################################
RUN = FALSE
if (RUN==TRUE) {
   for (ii in 1:10) {
      ## shuffling WT and mutant sources to compare model
      patchdat_10_kdm_shuff = patchdat_10_kdm
      patchdat_10_kdm_shuff$data = patchdat_10_kdm$data %>% select(Mark, Scan_ID, Patch, Type, neighbourhood) %>% distinct() %>% 
         shuffle_multiple_groups("neighbourhood", c("Type")) %>% 
         right_join(patchdat_10_kdm$data %>% select(-Type))
      patchdat_10_stag_shuff = patchdat_10_stag
      patchdat_10_stag_shuff$data = patchdat_10_stag$data %>% select(Mark, Scan_ID, Patch, Type, neighbourhood) %>% distinct() %>% 
         shuffle_multiple_groups("neighbourhood", c("Type")) %>% 
         right_join(patchdat_10_stag$data %>% select(-Type))
      
      data_primed_kdm_10_shuff = prime_data(wt_patches_kdm6a, patchdat_10_kdm_shuff)
      data_primed_stag_10_shuff = prime_data(wt_patches_stag2, patchdat_10_stag_shuff)
      
      stan_ins_10s_shuff = generate_stan_ins2(data_primed_stag_10_shuff, data_primed_kdm_10_shuff, length_scale, 
                                              area_scale, slide_age_ref)
      ## compile and run
      cmp_model = rstan::stan_model(file = "./stan_models/whitespace_diffusion_model.stan")
      nchains = 4
      options(mc.cores = nchains)
      iters = 25000
      samp_fit_shuff = rstan::sampling(object = cmp_model, data = stan_ins_10s_shuff$stan_ins$input, cores=nchains,
                                       init = rep(list(stan_ins_10s_shuff$stan_ins$init), nchains), chains=nchains,
                                       iter=iters, thin = 10)
      samp_fit_shuff %>% saveRDS(paste0("./stan_fits/diffusion_model_fit_10ps_hierarchical_normalage_shuffledsource_", 
                                 paste0(as.character(ii), ".rds")))
   }
}

RUN = FALSE
if (RUN==TRUE) {
   shuff_tib = tibble()
   for (ii in 1:10) {
      samp_fit_shuff = readRDS(paste0("./stan_fits/diffusion_model_fit_10ps_hierarchical_normalage_shuffledsource_", 
                                      paste0(as.character(ii), ".rds")))
      fits_shuff = rstan::extract(samp_fit_shuff)
      shuff_tib = shuff_tib %>% bind_rows(tibble(loglik_shuffle = fits_shuff$lp__, iter = ii))
   }
   raw_tib = tibble(loglik_mut = fits_10s$lp__)
   shuff_tib %>% ggplot() + 
      geom_histogram(aes(x=loglik_shuffle, y=..density.., fill="shuffled"), alpha=0.5) +
      geom_histogram(data = raw_tib, aes(x=loglik_mut, y=..density.., fill="raw"), alpha=0.5)
   shuff_tib %>% saveRDS("./output/sanity_check_shuffled_source_likelihood_draws.rds")
   raw_tib %>% saveRDS("./output/sanity_check_raw_source_likelihood_draws.rds")
   
   ## compare models using LOO
   library(loo)
   log_lik_1 = extract_log_lik(samp_fit_10s, merge_chains = FALSE)
   r_eff_1 = relative_eff(exp(log_lik_1))
   loo_1 = loo(log_lik_1, r_eff = r_eff_1, cores = 4)
   print(loo_1)
   shuff_loo_list = list()
   for (ii in 1:10) {
      samp_fit_shuff = readRDS(paste0("./stan_fits/diffusion_model_fit_10ps_hierarchical_normalage_shuffledsource_", 
                                      paste0(as.character(ii), ".rds")))
      log_lik_2 = extract_log_lik(samp_fit_shuff, merge_chains = FALSE)
      r_eff_2 = relative_eff(exp(log_lik_2))
      loo_2 = loo(log_lik_2, r_eff = r_eff_2, cores = 6)
      print(loo_2)
      shuff_loo_list[[ii]] = loo_2
   }
   
   comp = loo_compare(loo_1, shuff_loo_list[[1]], shuff_loo_list[[2]], shuff_loo_list[[3]],
                      shuff_loo_list[[4]], shuff_loo_list[[5]], shuff_loo_list[[6]], shuff_loo_list[[7]],
                      shuff_loo_list[[8]], shuff_loo_list[[9]], shuff_loo_list[[10]])
   print(comp) # so diffusion model is better
   
   comp %>% as_tibble() %>% mutate(model_type = rownames(comp)) %>% 
      mutate(model_type = ifelse(model_type == "model1", "raw data", model_type)) %>% 
      mutate(model_type = str_replace(model_type, "model", "shuffled ")) %>% 
      saveRDS("./output/LOO_compare_for_10_shuffle_iterations.rds")
   
   comp = readRDS("./output/LOO_compare_for_10_shuffle_iterations.rds")
   comp %>% select(model_type, elpd_diff)   
}

### Explore inference
#####################################################
## compare inferred patch-ages to birth-death processes
ageinfplots = compare_patchage_inference_to_prior(fits_10s, stan_ins_10s, timeatN_stag_10s, timeatN_kdm_10s)
ageinfplots$plot_kdm + xlim(0,100) + ggtitle("Inferred ages, KDM6A")
ageinfplots$plot_stag + xlim(0,100) + ggtitle("Inferred ages, STAG2")
tibble(D = fits_10s$D) %>% ggplot() + geom_histogram(aes(D), bins=100) + xlim(0,0.18) + 
   ggtitle("Diffusion coefficient, in crypt domains per year")
fits_10s$D %>% saveRDS("./output/diffusion_coefficient_draws.rds")

patient_patch_ref = stan_ins_10s$slide_index_tib %>% arrange(nbrhood_index) %>% 
   left_join(slide_age_ref) %>% left_join(fit_summaries_10s$t_p %>% select(param_num, med_mcmc) %>% 
                                             magrittr::set_colnames(c("nbrhood_index", "inferred_patch_age")))
patient_patch_ref %>% write_csv("./output/patient_age_reference.csv")
   
## get r-dependence predictions for each mutant patch
# new drip feed
library(foreach)
library(doMC)
doMC::registerDoMC(cores = 16)
useiters = sample(1:length(fits_10s$D), size = 1000, replace = FALSE)
ref_tib = stan_ins_10s$joinscaledat %>% left_join(slide_age_ref) %>% select(nbrhood_index, a_m, a_wt, Age) %>% distinct()
rvec = pracma::linspace(0.1,20,100)
all_nbrhood_runs = tibble()
for (j in (stan_ins_10s$slide_index_tib %>% pull(nbrhood_index))) {
   thispath = pracma::linspace(stan_ins_10s$stan_ins$input$first_hit_times[j], stan_ins_10s$stan_ins$input$pat_ages[j], stan_ins_10s$stan_ins$input$num_events + 1)
   thispath = thispath[1:stan_ins_10s$stan_ins$input$num_events]
   thisref = ref_tib %>% filter(nbrhood_index==j)
   sampleres = foreach (i = useiters)%dopar%{
      get_timeseries_from_single_path_real_D(i, thispath, thisref %>% pull(Age), fits_10s$D[i], thisref %>% pull(a_m), thisref %>% pull(a_wt), rvec)
   }
   dft_dc = do.call(rbind, sampleres) %>% left_join(thisref %>% rename(t=Age)) %>% group_by(t, r, nbrhood_index) %>% 
      summarise(med_mcmc = median(rho_pert), l95_mcmc = quantile(rho_pert, 0.025), 
                u95_mcmc = quantile(rho_pert, 0.975)) %>% ungroup() %>% left_join(thisref %>% rename(t=Age))
   all_nbrhood_runs = all_nbrhood_runs %>% bind_rows(dft_dc)
}
all_nbrhood_runs = all_nbrhood_runs %>% left_join(fit_summaries_10s$gamma_a %>% select(param_num, med_mcmc) %>% magrittr::set_colnames(c("nbrhood_index", "gamma_a_inf")))
all_nbrhood_runs_gamma = all_nbrhood_runs %>% mutate(med_mcmc = 1/(1/gamma_a_inf + med_mcmc), l95_mcmc = 1/(1/gamma_a_inf + l95_mcmc),
                                                     u95_mcmc = 1/(1/gamma_a_inf + u95_mcmc))
pl_list_10s = compare_r_dependence_to_data(data_primed_kdm_10, data_primed_stag_10, 
                                           stan_ins_10s, all_nbrhood_runs_gamma, WT30_r = NA)
pl_list_10s[[1]]

#############################

## plotting with Cora's specifications
location_by_patch = prepare_r_dep_data_for_plotting(data_primed_kdm_10, data_primed_stag_10, stan_ins_10s)
location_by_patch = location_by_patch %>% mutate(newtype = ifelse(Type=="mut", "Mutant", "Adjacent"))
nbrs_i = location_by_patch %>% pull(nbrhood_index) %>% unique() %>% sort()
for (i in 1:length(nbrs_i)) {
   nbr_i = nbrs_i[i]
   pat_info = location_by_patch %>% filter(nbrhood_index==nbr_i) %>% select(Scan_ID, mark, nbrhood_index) %>% 
      left_join(patient_patch_ref) %>% distinct()
   thisma = location_by_patch %>% filter(nbrhood_index==nbr_i) %>% pull(mark) %>% unique()
   maxr = max(10, location_by_patch %>% filter(nbrhood_index==nbr_i) %>% pull(R_mean) %>% max())
   pl = location_by_patch %>% filter(nbrhood_index==nbr_i) %>%  ggplot() + 
      geom_ribbon(data=r_dep_traj_10s %>% filter(nbrhood_index==nbr_i, r<maxr), 
                  aes(x=r, ymin=l95_mcmc, ymax=u95_mcmc, fill= "Theory" ), alpha=0.2) +
      geom_line(data=r_dep_traj_10s %>% filter(nbrhood_index==nbr_i, r<maxr), 
                aes(x=r, y=med_mcmc, col = "Theory")) + 
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
   pat_info_string = pat_info %>% 
      mutate(info_string = paste0(DW.Patient.number, 
        paste0(paste0("_age", as.character(Age)), 
         paste0("_patchage", as.character(round(inferred_patch_age)))))) %>% pull(info_string)
   ggsave(paste0("./plots/r_dep_gamma_3x2_nbrhood_", 
                 paste(paste(as.character(nbr_i), thisma, sep = "_"), 
                        paste0(pat_info_string, ".pdf"), sep = "_")),
          plot = pl, device = "pdf", width = 6, height = 4, units = "cm")
   ggsave(paste0("./plots/r_dep_gamma_3x1-6_nbrhood_", 
                 paste(paste(as.character(nbr_i), thisma, sep = "_"), 
                       paste0(pat_info_string, ".pdf"), sep = "_")),
          plot = pl, device = "pdf", width = 6, height = 3.2, units = "cm")
}

# Mutant patch in aquamarine4 for KDM6A and tan4 for STAG2, 
# wt patches aquamarine2 for KDM6A WT and tan2 for STAG2 WT?

###################################################

## finding threshold of fission rate for given diffusion rate
mark_data = read_csv("../drift_All_Marks/exp_coeff.csv")
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
D = median(fits_10s$D)
gamma_a = median(fits_10s$gamma_a_pop)
space_tib = stan_ins_10s$joinscaledat %>% select(mark, nbrhood_index, a_m, a_wt, s_wt) %>% distinct()
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
   ggplot() + geom_col(aes(x=t, y=med_ps)) + 
   geom_errorbar(aes(x=t, ymin=l95_ps, ymax=u95_ps), width=0.5) + facet_wrap(~wt_mult) +
   scale_y_log10()

################
## If n KRAS mutations at age X, how many will become a polyp in T years?
mark_data = read_csv("../drift_All_Marks/exp_coeff.csv")
marklist = c("MAOA", "STAG2", "KDM6A", "KRAS")
mark_fisrates = mark_data %>% filter(mark%in%marklist) %>% select(mark, rho)
wt_rate = mark_fisrates %>% filter(mark=="MAOA") %>% pull(rho)
mark_fisrates = mark_fisrates %>% mutate(scaled_rho = rho / wt_rate)
mark_fisrates = mark_fisrates %>% bind_rows(tibble(mark = "old_KRAS", rho = 10*wt_rate, scaled_rho = 10))
# rho_fi = wt_rate * 10 # use "old" published KRAS fission rate of 10xWT
rho_fi = mark_fisrates %>% filter(mark=="KRAS") %>% pull(rho) # use new KRAS fission rate of ~17xWT
gamma_min = 1 - pi/sqrt(12)
rho_max = 1/gamma_min
iters = 20000
D = median(fits_10s$D)
gamma_a = median(fits_10s$gamma_a_pop)
rvec = pracma::linspace(0,10,101)
roll_area_smoothed = c(pi*mean(rvec[1:2])^2)
for (i in 2:(length(rvec)-1)) roll_area_smoothed = c(roll_area_smoothed, pi*mean(rvec[(i-1):(i+1)])^2)
roll_area_smoothed = c(roll_area_smoothed, pi*mean(tail(rvec, 2))^2)
rollarea_tib = tibble(r = rvec, roll_area = roll_area_smoothed)
tvec = pracma::linspace(1,50,50)
space_tib = stan_ins_10s$joinscaledat %>% select(mark, nbrhood_index, a_m, a_wt, s_wt) %>% distinct()
a_wt = space_tib %>% pull(a_wt) %>% median()
s_wt = space_tib %>% pull(s_wt) %>% median()
a_m = a_wt
suptimes_list = list()
suptimes_list[[1]] = c(0, 10)
suptimes_list[[2]] = c(10, 20)
suptimes_list[[3]] = c(20, 30)
rm(fits_10s); rm(samp_fit_10s); gc()
rm(timeatN_stag_10s); rm(timeatN_kdm_10s); gc()
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

         kras_polyp_summedfracs %>% saveRDS(paste0("./output/krasnew_over_time_forfracpolyp_suppressedrate_", 
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
   kras_polyp_summedfracs %>% saveRDS("./output/krasnew_over_time_forfracpolyp_unsuppressed.rds")
   rm(res_kras); gc()
}

## could instead do as the fraction of the area containing patchsize * a_m that 
## has a rolling average below the threshold?
joinedtib = tibble()
for (st in 1:length(suptimes_list)) {
   suptimes = suptimes_list[[st]]
   suptime_string = paste(as.character(suptimes[1]), as.character(suptimes[2]), sep="-")
   for (suprate in c(0, wt_rate)) {
      kras_polyp_summedfracs = readRDS(paste0("./output/krasnew_over_time_forfracpolyp_suppressedrate_", 
                                                paste0(as.character(round(suprate, 3)), 
                                                       paste0("_suppressedtime", 
                                                              paste0(suptime_string, ".rds")))))
      joinedtib = joinedtib %>% bind_rows(kras_polyp_summedfracs)
   }
}
mut_tib = readRDS("./output/krasnew_over_time_forfracpolyp_unsuppressed.rds") %>% select(t, sum_chance)

joinedtib %>% 
   ggplot() + geom_line(aes(x=t, y=sum_chance, col=suppression_times)) + 
   facet_wrap(~suppressed_rate) + 
   geom_line(data = mut_tib, aes(x=t, y=sum_chance, col="KRAS"))

#########################
## Finding the fission rate that would lead to 5% polyp chance
## If n KRAS mutations at age X, how many will become a polyp in T years?
mark_data = read_csv("../drift_All_Marks/exp_coeff.csv")
marklist = c("MAOA", "STAG2", "KDM6A", "KRAS")
mark_fisrates = mark_data %>% filter(mark%in%marklist) %>% select(mark, rho)
wt_rate = mark_fisrates %>% filter(mark=="MAOA") %>% pull(rho)
fission_rate_list = 16:20
gamma_min = 1 - pi/sqrt(12)
rho_max = 1/gamma_min
iters = 10000
D = median(fits_10s$D)
gamma_a = median(fits_10s$gamma_a_pop)
rvec = pracma::linspace(0,10,101)
roll_area_smoothed = c(pi*mean(rvec[1:2])^2)
for (i in 2:(length(rvec)-1)) roll_area_smoothed = c(roll_area_smoothed, pi*mean(rvec[(i-1):(i+1)])^2)
roll_area_smoothed = c(roll_area_smoothed, pi*mean(tail(rvec, 2))^2)
rollarea_tib = tibble(r = rvec, roll_area = roll_area_smoothed)
tvec = pracma::linspace(1,50,50)
space_tib = stan_ins_10s$joinscaledat %>% select(mark, nbrhood_index, a_m, a_wt, s_wt) %>% distinct()
a_wt = space_tib %>% pull(a_wt) %>% median()
s_wt = space_tib %>% pull(s_wt) %>% median()
a_m = a_wt
rm(fits_10s); rm(samp_fit_10s); gc()
rm(timeatN_stag_10s); rm(timeatN_kdm_10s); gc()
RUN = TRUE
if (RUN==TRUE) {
   for (rhomult in fission_rate_list) {
      rho_fi = wt_rate * rhomult
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
      kras_polyp_summedfracs %>% saveRDS(paste0("./output/fisrate_forfracpolyp_", paste0(as.character(rhomult), "xWT.rds")))
      rm(res_kras); gc()  
   }
}

jointib = tibble()
for (rhomult in 10:20) {
   res = readRDS(paste0("./output/fisrate_forfracpolyp_", paste0(as.character(rhomult), "xWT.rds")))
   jointib = jointib %>% bind_rows(res %>% select(t, sum_chance) %>% mutate(wt_mult = rhomult))
}
jointib %>% write_csv("./output/polyp_probs_fissionrates.csv")


######################
### Look at range using rho_max, gamma_min?
av_patch_inputs = stan_ins_10s$joinscaledat %>% select(Scan_ID, mark, neighbourhood, nbrhood_index, epsilon) %>% 
   distinct() %>% group_by(mark) %>% summarise(med_eps = median(epsilon))
gamma_a = median(fits_10s$gamma_a_pop)
## just define range as the "new homeostasis" assuming diffusion hits some wall
RUN = FALSE
if (RUN==TRUE) {
   for (thresh in c(1.01, 1.05)) {
      # using average mass input
      happy_packing_range = av_patch_inputs %>% mutate(R_H = sqrt(med_eps/(pi/gamma_a * (thresh-1))))
      # or calculate "errorbars"
      indiv_ins = stan_ins_10s$joinscaledat %>% select(Scan_ID, mark, neighbourhood, nbrhood_index, epsilon) %>% 
         distinct() %>% arrange(nbrhood_index)
      outtib = tibble()
      for (ii in 1:nrow(fits_10s$gamma_a)) {
         outtib = outtib %>% bind_rows(indiv_ins %>% mutate(gamma_a = fits_10s$gamma_a[ii,]) %>% 
                                          mutate(R_H = sqrt(epsilon/(pi/gamma_a * (thresh-1))), iter=ii))
      }
      outtib %>% group_by(mark) %>% summarise(l95_mcmc = quantile(R_H, 0.025), 
                                              med_mcmc = median(R_H), 
                                              u95_mcmc = quantile(R_H, 0.975)) %>% 
         ggplot() + geom_point(aes(x=mark, y=med_mcmc, col=mark)) +
         geom_errorbar(aes(x=mark, ymin=l95_mcmc, ymax=u95_mcmc, col=mark), width=0.4)
      
      outtib = outtib %>% mutate(N_crypt_footprint = calc_happy_packing_footprint(R_H, stan_ins_10s$joinscaledat %>% 
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


## look at threshold on new homeostatic density for "old patches" in the data
pc_explore = r_dep_traj_10s_rho %>% group_by(mark, Scan_ID, neighbourhood, nbrhood_index) %>% 
   summarise(rho_patch = med_mcmc[1], rho_inf = tail(med_mcmc, 1)) %>% 
   mutate(pc_diff = rho_patch/rho_inf)
pc_explore %>% ggplot() + geom_histogram(aes(x=pc_diff)) + xlim(1,1.1)
pc_explore %>% ungroup() %>% summarise(lq_pc = quantile(pc_diff, 0.25), med_pc = median(pc_diff), uq_pc = quantile(pc_diff, 0.75))
pc_explore %>% ggplot() + geom_boxplot(aes(x=1, y=pc_diff)) + ylim(1, 1.1) + xlim(0, 2)
mult = 1.01 # set at 1%

## pick median patch age and calculate range vs fission rate and range vs mutant crypt size / time
n_mcs = 100
n_fr = 100
space_tib = stan_ins_10s$joinscaledat %>% select(mark, nbrhood_index, a_m, a_wt, s_wt) %>% distinct()
a_wt = space_tib %>% pull(a_wt) %>% median()
s_wt = space_tib %>% pull(s_wt) %>% median()
tvec = pracma::linspace(1,50,n_mcs)
lambda_list = pracma::linspace(wt_rate,10*wt_rate,n_fr)
gamma_a = median(fits_10s$gamma_a_pop)
pc_thresh = 1
patch_age = median(fits_10s$t_p)
out_mat = matrix(0, n_mcs, n_fr)
for (i in 1:n_mcs) {
   out_mat[i,] = calc_happy_packing_range(lambda_list, tvec[i], a_wt, a_wt, gamma_a, pc_thresh)
}
gath_tib = out_mat %>% as_tibble() %>% magrittr::set_colnames(as.character(lambda_list)) %>% 
   mutate(time = tvec) %>% gather(rho, R_H, -time) %>% 
   mutate(R_H = ifelse(is.na(R_H), 0, R_H)) %>% mutate(rho = as.numeric(rho))
gath_tib %>% ggplot(aes(x = time, y = rho, fill = R_H)) +
   labs(x = "age of patch", y = "fission rate", fill = "range") +
   geom_raster()
gath_tib = gath_tib %>% mutate(crypt_footprint = calc_happy_packing_footprint(R_H, s_wt))
gath_tib %>% ggplot(aes(x = time, y = rho, fill = crypt_footprint)) +
   labs(x = "age of patch", y = "fission rate", fill = "n_crypt footprint") +
   geom_raster()

gath_tib %>% ggplot() + geom_contour(aes(x = time, y = rho, z = crypt_footprint, colour = ..level..)) +
   scale_color_viridis_c()

gath_tib %>% ggplot() + 
   aes(x = time, y = rho, z = crypt_footprint, fill = crypt_footprint) + 
   geom_tile() + 
   geom_contour(color = "white", alpha = 0.5) + 
   scale_fill_distiller(palette="Spectral", na.value="white") + 
   theme_bw() + 
   theme(plot.background = element_blank(), panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), panel.border = element_blank(),
         axis.line = element_line(color = 'black'))

gath_tib %>% ggplot() + 
   aes(x = time, y = rho, z = log2(crypt_footprint), fill = log2(crypt_footprint)) + 
   geom_tile() + 
   scale_fill_distiller(palette="Spectral", na.value="white") + 
   theme_bw() + 
   theme(plot.background = element_blank(), panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), panel.border = element_blank(),
         axis.line = element_line(color = 'black'))


### patch ages OBSERVED: save the full patch lifetime
sttimes = c()
entimes = c()
for (i in 1:1e6) {
   ts = generate_fission_event_times(0.01, 100)
   if (length(ts)>=9) {
      sttimes = c(sttimes, ts[9])
      entimes = c(entimes, ts[10])
   }
}

tibble(raw_times = ptimes) %>% mutate(rounded_times = ceiling(raw_times)) %>% 
   group_by(rounded_times) %>% summarise(num_this_time = n()) %>% 
   arrange(rounded_times) %>% mutate(cumulative = cumsum(num_this_time))

sum_p = rep(0, 100)
for (i in 1:length(sttimes)) {
   if (sttimes[i]<=100) {
      if (is.na(entimes[i])) {
         ee = 100
      } else {
         ee = min(100, round(entimes[i]))
      }
      st = round(sttimes[i])
      print(i)
      print(st)
      print(ee)
      sum_p[st:ee] = sum_p[st:ee] + 1   
   }
}
