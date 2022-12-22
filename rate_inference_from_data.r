## Fufi estimations from data
library(readxl)
library(tidyverse)
source("run_funcs.r")

## Estimating the fusion rate given the fission rate
datin = read_xlsx("./databases/Fufi_counts_KDM6A_curated_for_rate_calc.xlsx")
colnames(datin) = c("Image_ID", "monoclonals", "knn", "partials", "not_needed")
datin = datin %>% select(-not_needed)
knn_vec = datin %>% pull(knn)
wts = c()
muts = c()
for (i in 1:length(knn_vec)) {
   elem = knn_vec[i]
   if (elem!="0" && elem!="N/A") {
      splitfufis = strsplit(elem, ", ")[[1]]
      for (j in 1:length(splitfufis)) {
         oneevent = splitfufis[j]
         eventsplit = strsplit(oneevent, " WT ")[[1]]
         wts = c(wts, as.integer(eventsplit[1]))
         muts = c(muts, as.integer(strsplit(eventsplit[2], " mut")[[1]][1]))
      }
   }
}
totN = wts+muts
chi = mean(muts/totN)
# or with perturbation due to unknown neighbours prior to event
#chi = mean(neighbour_perturbation(wts, 0.)) #-0.4

datin = datin %>% mutate(monoclonals = as.integer(monoclonals), partials = as.integer(partials))
datin = datin %>% mutate(monoclonals = ifelse(is.na(monoclonals), 0, monoclonals), partials = ifelse(is.na(partials), 0, partials))
sumsdat = datin %>% summarise(n_mon = sum(monoclonals), n_par = sum(partials))
numfrac = (sumsdat %>% pull(n_mon)) / (sumsdat %>% pull(n_par))
rho_fi = 0.0391560702622755 #/365 # per year or /365 for per day; 0.007/365 is WT value
rho_fu = rho_fus_estimate(rho_fi, numfrac, chi)
rho_fi_WT = 0.007 # /365
rho_fu_WT = rho_fi_WT * rho_fu / rho_fi # if the ratio stays the same from KDM6A!

## NOTE: In Amsterdam Lotte suggested that a large proportion (a third?!) of fufi events are "frozen" -- neither zipping nor unzipping;
# we should take account of this observation in our calculations by looking at how a "freeze rate" betweeen
# 0 and 1 affects the rate and duration calculation (with the frozen fufis randomly chosen from the data and repeated
# multiple times for confidence intervals.)

# calculate error bars and plot
load("./data/fusion_rate_error_estimate_data.Rdata") # output
output = output %>% mutate(fi_rate_real = fi_rate_real*365, fu_rate_real = fu_rate_real*365) # into yearly rate
confwidth = 95
wild_neighbour_perturbation = 0
output = output %>% group_by(run_id) %>% 
   mutate(nP = sum(partials), nM = sum(monoclonals))
rate_error_estimate = output %>% mutate(chi = ifelse(monoclonals==0, 0, neighbour_perturbation(neighbours, wild_neighbour_perturbation)), 
                                        chi_mean = sum(chi)/nM, 
                                        rho_fu_est = rho_fus_estimate(fi_rate_real, nM/nP, chi_mean),
                                        rel_err = rho_fu_est/fu_rate_real - 1,
                                        abs_rel_err = abs(rel_err)) %>% 
   summarise(nP = nP[1], nM = nM[1], abs_rel_err=abs_rel_err[1], 
             rel_err=rel_err[1], rho_fu_est=rho_fu_est[1], chi = chi_mean[1])
average_rate_errs = rate_error_estimate  %>% ungroup() %>% 
   summarise(rate = mean(rho_fu_est),
             rate_95_lwr=quantile(rho_fu_est, probs=((100-confwidth)/2)/100),
             rate_95_upr=quantile(rho_fu_est, probs=(1-((100-confwidth)/2)/100)),
             mean_err=mean(abs_rel_err),
             sd_err= sd(abs_rel_err),
             low_err = quantile(abs_rel_err, probs=((100-confwidth)/2)/100),
             upp_err = quantile(abs_rel_err, probs=(1-((100-confwidth)/2)/100)))
fu_rate_m = average_rate_errs %>% pull(rate)
fu_rate_sd = sd(rate_error_estimate %>% pull(rho_fu_est))
# plotting fusion rate and fission rate with error bars
# fission rate is reported as 0.68 ± 0.15 and 0.72 ± 0.15 (mPAS and MAOA, respectively)
# -> take as 0.7 ± 0.15 (% crypts per year) i.e. 0.007 ± 0.0015 yearly rate
sd_kdm6a = 0.00422329660794263 #/365 #divide by 365 for daily
fistib = tibble(rate = rho_fi, rate_95_lwr = rate-sd_kdm6a, rate_95_upr = rate+sd_kdm6a, mean = 0, sd = 0, low = 0, upp = 0)
fisfus_compare = bind_rows(average_rate_errs %>% mutate(type="fusion") , fistib %>% mutate(type="fission"))
p2 = suppressMessages(fisfus_compare %>% ggplot() + geom_point(aes(x = type, y = rate, col = type), size = 2) + 
                      geom_errorbar(aes(x = type, ymin = rate_95_lwr, ymax = rate_95_upr, col = type), width = 0.1) +
                      ylim(0, 0.05) + ylab("Rate, per year") + xlab("Reaction type (KDM6A)")) + 
                      theme_bw(base_size = 15)
plot(p2)
ggsave(filename = "/home/doran/Work/r_code/crypt_fufi_hexgrid/plots/kdm6a_fusion-fission_yearlyrates.eps")

## Estimating the fufi duration given the rates and the assumption that both events share a similar duration
###################################################################
# (though we are missing crypt counts for one slide)
#datin = read_xlsx("./Fufi_counts_KDM6A_curated_for_duration_calc.xlsx") %>% mutate(Slide_ID = as.integer(Image_ID)) %>% select(-Image_ID)
datin = read_xlsx("./databases/Fufi_counts_KDM6A_curated_for_duration_calc_dropped_missing_slide.xlsx") %>% mutate(Slide_ID = as.integer(Image_ID)) %>% select(-Image_ID)
# datin = read_xlsx("./databases/Fufi_counts_KDM6A_curated_for_duration_calc.xlsx") %>% mutate(Slide_ID = as.integer(Image_ID)) %>% select(-Image_ID)
countsin = read_csv("/home/doran/Work/images/KDM6A_fufi_analysis/Analysed_slides/slide_counts.csv")
# missingcounts = read_csv(file = "./databases/missing_cryptcounts_kdm_COUNTS.csv")
# kdm_slidecounts = kdm_slidecounts %>% bind_rows(missingcounts)
datin = left_join(datin, countsin)
datin = datin %>% mutate_all(list(~replace(., is.na(.), 0)))

# attempting to sort by age; getting block data
################################################
blockdat = read_xlsx("/home/doran/Work/images/Blocks_database.xlsx") %>% select(`Block ID`, Age, SEX)
block_to_slide = read_xlsx("./databases/Image_Locations_KDM6A.xlsx") %>% select(`Block ID`, `Image ID`)
colnames(blockdat) = c("Block_ID", "Age", "Sex")
colnames(block_to_slide) = c("Block_ID", "Image_ID")
blockdat = blockdat %>% mutate(staticID1 = Block_ID) %>% 
   mutate(Block_ID = str_replace(Block_ID, "BL", "")) %>% 
   mutate(Block_ID = str_replace(Block_ID, "PS10:", "")) %>% 
   mutate(Block_ID = str_replace(Block_ID, "PS09:", "")) %>% 
   mutate(Block_ID = str_replace(Block_ID, "HB", "")) %>% 
   separate(Block_ID, c("Block_ID", "slide"), "_") %>% 
   separate(Block_ID, c("Block_ID", "slide"), " ") %>% select(-slide)
block_to_slide = block_to_slide %>% mutate(staticID2 = Block_ID) %>% 
   separate(Block_ID, c("Block_ID", "slide"), "_") %>% select(-slide)

joined_b = block_to_slide %>% left_join(blockdat) %>% select(-staticID1, -staticID2)
joined_b = joined_b %>% mutate(Age = ifelse(is.na(Age), 0, Age)) %>% mutate(Image_ID = as.character(Image_ID))

# kras_data = read_csv("/home/doran/Work/images/all_kras_data.csv")
# joined_b = joined_b %>% left_join(kras_data %>%  rename(Block_ID = Block) %>% mutate(Block_ID = as.character(Block_ID)))
# joined_b %>% write_csv("/home/doran/Work/images/all_kras_data+ImageID.csv", col_names = T)
## And getting clone counts for each slide
clonecounts = read_xlsx("./databases/KDM6A_mutant_counts.xlsx") %>% 
   select(-`DW patient ID`, -Section, -`stain quality`, -Fractions, -Comments)
clonecounts = clonecounts %>% mutate(`2` = `2` * 2,
                       `3` = `3` * 3,
                       `4` = `4` * 4,
                       `5` = `5` * 5,
                       `6` = `6` * 6,
                       `7` = `7` * 7,
                       `8` = `8` * 8,
                       `9` = `9` * 9,
                       `10` = `10` * 10,
                       `11` = `11` * 11,
                       `12` = `12` * 12,
                       `13` = `13` * 13,
                       `14` = `14` * 14)
clonecounts = clonecounts %>% group_by(`Block ID`) %>% mutate(`15+ split` = sum(as.integer(str_split(`15+`, ",")[[1]]))) %>% 
   ungroup() %>% select(-`15+`)
clonecounts = clonecounts %>% gather(patchtype, count, (PPC):(`15+ split`)) %>% group_by(`Block ID`) %>% 
   summarise(NClones = sum(count)) %>% mutate(Block_ID = as.character(`Block ID`)) %>% select(-`Block ID`)
joined_b = joined_b %>% left_join(clonecounts)
# Not all IDs can be matched! We cannot find clone counts for all slides!

# summing all fufis, WTs and Muts
sumpl = datin %>% group_by(Slide_ID, NCrypts) %>% summarise(numfufis = sum(`WT FUFI`) + sum(`Mut FUFI`) + sum(`WT_Mut FUFI`) + sum(`Mut FUFI patch`) +
                                                            sum(`Mut FUFI patch border`) + sum(`WT_Mut FUFI patch border`) + sum(`WT Fufi patch border`)) %>% 
   mutate(frac = numfufis/NCrypts) %>% 
   mutate(ylow  = qbeta(p = 0.025, shape1 = numfufis, shape2 = NCrypts-numfufis), 
          ymid  = qbeta(p = 0.5  , shape1 = numfufis, shape2 = NCrypts-numfufis), 
          yhigh = qbeta(p = 0.975, shape1 = numfufis, shape2 = NCrypts-numfufis))
avs_mean_med = sumpl %>% ungroup() %>% summarise(meanval = mean(frac), medianval = median(frac))
sumpl %>% ggplot() + geom_point(aes(x=as.character(Slide_ID), y=ymid)) + 
   geom_errorbar(aes(x=as.character(Slide_ID), y=ymid, ymin = ylow, ymax = yhigh)) + 
   geom_hline(yintercept = avs_mean_med %>% pull(meanval))

summarycounts = datin %>% summarise(numfufis = sum(`WT FUFI`) + sum(`Mut FUFI`) + sum(`WT_Mut FUFI`) + sum(`Mut FUFI patch`) +
                                       sum(`Mut FUFI patch border`) + sum(`WT_Mut FUFI patch border`) + sum(`WT Fufi patch border`), 
                                    totalcrypts = sum(NCrypts)) %>% mutate(fraction = numfufis/totalcrypts)

## Separating WT, Mut and WT_Mut fufis
sumpl_sep = datin %>% group_by(Slide_ID, NCrypts) %>% summarise(WT = sum(`WT FUFI`) + sum(`WT Fufi patch border`), 
                                                                MUT = sum(`Mut FUFI`) + sum(`Mut FUFI patch`) + sum(`Mut FUFI patch border`),
                                                                WT_MUT = sum(`WT_Mut FUFI`) + sum(`WT_Mut FUFI patch border`),
                                                                MUTall = MUT + WT_MUT)
sumpl_sep = sumpl_sep %>% ungroup() %>%  mutate(Image_ID = as.character(Slide_ID)) %>% select(-Slide_ID)
sumpl_sep = sumpl_sep %>% left_join(joined_b) %>% mutate(NWTCrypts = ifelse(is.na(NClones), NCrypts, NCrypts-NClones))
wt_mut_fufi_sep = sumpl_sep %>% filter(!is.na(NClones)) %>% summarise(wt_fufis = sum(WT) + 0.5*sum(WT_MUT), mut_fufis = sum(MUT) + 0.5*sum(WT_MUT), 
                                          tot_mut_crypts = sum(NClones), tot_wt_crypts = sum(NWTCrypts)) %>% 
   mutate(fraction_fufi = (wt_fufis + mut_fufis) / (tot_wt_crypts + tot_mut_crypts))
frac_wt = (wt_mut_fufi_sep %>% pull(tot_wt_crypts)) / ( (wt_mut_fufi_sep %>% pull(tot_wt_crypts)) + (wt_mut_fufi_sep %>% pull(tot_mut_crypts)) )
# calculate duration using correct rates for WT and Mut fufi events
duration_denominator = (frac_wt*(rho_fi_WT + rho_fu_WT)) + (1-frac_wt)*(rho_fi + rho_fu)
est_duration = -log(1 - wt_mut_fufi_sep %>% pull(fraction_fufi)) / duration_denominator

## Using STAN model for better visualisation of error bars for fufi fractions near zero
######################################################################################
library(rstan)
library(tidybayes)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## fit pooling all fufis
#########################
N = sumpl %>% pull(NCrypts)
n = sumpl %>% pull(numfufis)
fit1 = stan(file = "./hierarchical_model_for_fufi_fraction.stan", data = list(num_slides=length(n), N = N, n = n))
print(fit1, digits=8)
traceplot(fit1)
#fit_ss <- extract(fit1, permuted = TRUE)
draws = as.tibble(as.matrix(fit1)[,1:(length(n)+2)])
population_draws = draws %>% select(mean_p, var_p)
slide_draws = draws %>% select(-mean_p, -var_p)

## plotting per slide fufi fractions and error bars
#####################################################
# All pooled fufis
colnames(slide_draws) = sumpl %>% pull(Slide_ID)
g_slide_draws = slide_draws %>% gather(Slide_ID, iteration_val)
g_slide_draws %>% group_by(Slide_ID) %>% summarise(median = quantile(iteration_val, 0.5),
                                                 lower = quantile(iteration_val, 0.025),
                                                 upper = quantile(iteration_val, 0.975)) %>%
   ggplot() + geom_point(aes(x=as.character(Slide_ID), y=median)) +
   geom_errorbar(aes(x=as.character(Slide_ID), ymin = lower, ymax = upper)) +
   geom_hline(yintercept = mean(population_draws %>% pull(mean_p)), linetype="dotted") +
   ylab("Fraction of crypts fufi-ing") + xlab("Slide") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
# finding outliers
allfufidraws = g_slide_draws %>% group_by(Slide_ID) %>% summarise(median = quantile(iteration_val, 0.5),
                                                   lower = quantile(iteration_val, 0.025),
                                                   upper = quantile(iteration_val, 0.975))

# plotting by age, pooling all fufis
####################################
age_slide_draws = g_slide_draws %>% rename(Image_ID = Slide_ID) %>% left_join(joined_b)
age_slide_draws %>% group_by(Image_ID, Age, Sex) %>% summarise(median = quantile(iteration_val, 0.5),
                                                           lower = quantile(iteration_val, 0.025),
                                                           upper = quantile(iteration_val, 0.975)) %>% 
   ggplot() + geom_point(aes(x=Age, y=median)) + 
   geom_errorbar(aes(x=Age, ymin = lower, ymax = upper)) +
   geom_hline(yintercept = mean(population_draws %>% pull(mean_p)), linetype="dotted") +
   ylab("Fraction of crypts fufi-ing")

##########################################################################################
## to propagate these MCMC confidence intervals forward to the duration calculation, simply calculate 
# duration for every iteration value in the population mean draws
population_draws = population_draws %>% mutate(fufi_duration = -log(1 - mean_p) / duration_denominator)
gpopdraws = population_draws %>% gather(param, val) %>% group_by(param) %>% summarise(middle = quantile(val, 0.5),
                                                                                       lower = quantile(val, 0.025),
                                                                                       upper = quantile(val, 0.975)) 
gpopdraws %>% slice(1) %>% ggplot() + geom_point(aes(x=param, y=middle)) +
   geom_errorbar(aes(x=param, ymin = lower, ymax = upper), width = 0.2)

population_draws %>% gather(param, val) %>% filter(param=="fufi_duration") %>% ggplot() + 
   geom_violin(aes(x=param, y=val)) + geom_point(dat=gpopdraws %>% slice(1), aes(x=param, y=middle)) +
   geom_errorbar(data=gpopdraws %>% slice(1), aes(x=param, ymin = lower, ymax = upper), width = 0.06) +
   xlab("") + ylab("Duration, days")

## Metadata for Cora
agesexcounts_dur = datin %>% select(Slide_ID:NCrypts) %>% left_join(joined_b) %>% select(-NClones)
agesexcounts_dur %>% select(Block_ID:Sex) %>% unique()
# total crypts
agesexcounts_dur %>% summarise(totcrypts = sum(NCrypts))
agesexcounts_dur %>% select(Block_ID:Sex) %>% unique() %>% arrange(desc(Age))
agesexcounts_dur %>% select(Block_ID:Sex) %>% unique() %>% arrange(Age)
agesexcounts_dur %>% filter(Age==0 | is.na(Age))
# output all
kdm6a_agesex %>% select(Block_ID, Slide_ID, Age, Sex) %>% unique() %>% 
   left_join(kdm_slidecounts) %>% 
   write_csv(path = "./fusion_data_used_agesexcounts_kdm6a.csv")

### Extra analysis: splitting WT and Mut fufis (too noisy, really, as not enough Mut data)
##########################################################################################
## fit for only WTs
#########################
# N = sumpl_sep %>% pull(NWTCrypts)
# n = sumpl_sep %>% pull(WT)
# fitWT = stan(file = "./hierarchical_model_for_fufi_fraction.stan", data = list(num_slides=length(n), N = N, n = n))
# traceplot(fitWT)
# drawsWT = as.tibble(as.matrix(fitWT)[,1:(length(n)+2)])
# population_drawsWT = drawsWT %>% select(mean_p, var_p)
# slide_drawsWT = drawsWT %>% select(-mean_p, -var_p)
# 
# ## fit for all MUTs, whether MUT_MUT or WT_MUT
# #########################
# N = sumpl_sep %>% filter(!is.na(NClones), NClones>0) %>% pull(NClones)
# n = sumpl_sep %>% filter(!is.na(NClones), NClones>0) %>% pull(MUTall)
# fitMUTall = stan(file = "./hierarchical_model_for_fufi_fraction.stan", data = list(num_slides=length(n), N = N, n = n))
# traceplot(fitMUTall)
# drawsMUTall = as.tibble(as.matrix(fitMUTall)[,1:(length(n)+2)])
# population_drawsMUTall = drawsMUTall %>% select(mean_p, var_p)
# slide_drawsMUTall = drawsMUTall %>% select(-mean_p, -var_p)
# 
## plotting per slide fufi fractions and error bars
#####################################################
# # WT fufis
# colnames(slide_drawsWT) = sumpl_sep %>% pull(Image_ID)
# g_slide_drawsWT = slide_drawsWT %>% gather(Image_ID, iteration_val) 
# g_slide_drawsWT %>% group_by(Image_ID) %>% summarise(median = quantile(iteration_val, 0.5),
#                                                    lower = quantile(iteration_val, 0.025),
#                                                    upper = quantile(iteration_val, 0.975)) %>% 
#    ggplot() + geom_point(aes(x=as.character(Image_ID), y=median)) + 
#    geom_errorbar(aes(x=as.character(Image_ID), ymin = lower, ymax = upper)) +
#    ylab("Fraction of crypts fufi-ing") + xlab("Slide") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
# 
# # all MUT fufis
# colnames(slide_drawsMUTall) = sumpl_sep %>% filter(!is.na(NClones), NClones>0) %>% pull(Image_ID)
# g_slide_drawsMUTall = slide_drawsMUTall %>% gather(Image_ID, iteration_val)
# g_slide_drawsMUTall %>% group_by(Image_ID) %>% summarise(median = quantile(iteration_val, 0.5),
#                                                          lower = quantile(iteration_val, 0.025),
#                                                          upper = quantile(iteration_val, 0.975)) %>% 
#    ggplot() + geom_point(aes(x=as.character(Image_ID), y=median)) + 
#    geom_errorbar(aes(x=as.character(Image_ID), ymin = lower, ymax = upper)) +
#    ylab("Fraction of crypts fufi-ing") + xlab("Slide") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

####################################

# plotting by age, WT fufis
# age_slide_drawsWT = g_slide_drawsWT %>% left_join(joined_b)
# age_slide_drawsWT %>% group_by(Image_ID, Age, Sex) %>% summarise(median = quantile(iteration_val, 0.5),
#                                                                lower = quantile(iteration_val, 0.025),
#                                                                upper = quantile(iteration_val, 0.975)) %>% 
#    ggplot() + geom_point(aes(x=Age, y=median)) + 
#    geom_errorbar(aes(x=Age, ymin = lower, ymax = upper)) +
#    geom_hline(yintercept = mean(population_draws %>% pull(mean_p)), linetype="dotted") +
#    ylab("Fraction of crypts fufi-ing")
# wtdata_summarised = age_slide_drawsWT %>% group_by(Image_ID, Age, Sex) %>%
#    summarise(median = quantile(iteration_val, 0.5),
#              lower = quantile(iteration_val, 0.025),
#              upper = quantile(iteration_val, 0.975))
# # plotting by age, all MUT fufis
# age_slide_drawsMUTall = g_slide_drawsMUTall %>% left_join(joined_b)
# age_slide_drawsMUTall %>% group_by(Image_ID, Age, Sex) %>% summarise(median = quantile(iteration_val, 0.5),
#                                                                lower = quantile(iteration_val, 0.025),
#                                                                upper = quantile(iteration_val, 0.975)) %>% 
#    ggplot() + geom_jitter(aes(x=Age, y=median)) + 
#    #geom_errorbar(aes(x=Age, ymin = lower, ymax = upper)) +
#    geom_hline(yintercept = mean(population_draws %>% pull(mean_p)), linetype="dotted") +
#    ylab("Fraction of mutant crypts fufi-ing")
# mutdata_summarised = age_slide_drawsMUTall %>% group_by(Image_ID, Age, Sex) %>%
#    summarise(median = quantile(iteration_val, 0.5),
#              lower = quantile(iteration_val, 0.025),
#              upper = quantile(iteration_val, 0.975))

## plotting together
# wtdata_summarised %>% ggplot() + geom_point(aes(x=Age, y=median)) + 
#    geom_errorbar(aes(x=Age, ymin = lower, ymax = upper)) +
#    geom_hline(yintercept = mean(population_draws %>% pull(mean_p)), linetype="dotted") +
#    ylab("Fraction of crypts fufi-ing") +
#    geom_point(data=mutdata_summarised, aes(x=Age, y=median), color = "red") +
#    geom_errorbar(data=mutdata_summarised, aes(x=Age, ymin = lower, ymax = upper), color = "red")
# 
# # Looking at ratio
# get_right_slides = unique(g_slide_drawsMUTall %>% pull(Image_ID))
# mutwtjoined = g_slide_drawsMUTall %>% rename(MUTiter = iteration_val) %>% 
#    bind_cols(g_slide_drawsWT %>% rename(WTiter = iteration_val) %>% filter(Image_ID %in% get_right_slides)) %>% 
#    select(-Image_ID1)
# mutwtjoined = mutwtjoined %>% mutate(ratio = WTiter/MUTiter) %>% group_by(Image_ID) %>% 
#    summarise(median = quantile(ratio, 0.5),
#              lower = quantile(ratio, 0.025),
#              upper = quantile(ratio, 0.975))  %>% left_join(joined_b)
# mutwtjoined %>% ggplot() + geom_point(aes(x=Age, y=median)) + 
#    geom_errorbar(aes(x=Age, ymin = lower, ymax = upper)) + scale_y_log10()
# mutwtjoined %>% ggplot() + geom_point(aes(x=Age, y=median)) + scale_y_log10()
# mutwtjoined %>% ggplot() + geom_jitter(aes(x=Age, y=median))

