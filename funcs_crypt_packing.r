library(tidyverse)

get_timeseries_from_path = function(i, thesepaths, thisage, effective_diff_factor, D, a_m, a_wt, rvec) {
   c_fis_time = (thesepaths %>% filter(iter==i) %>% pull(event_times))[-1] # get rid of first hit event time
   diffusion_time_series_given_event_times(c_fis_time, thisage, effective_diff_factor*D, a_m, a_wt, rvec) %>% 
      rename(dripfeed = rho_pert) %>% mutate(iter=i)
}

create_diffusion_perturbation_with_first_hit_time = function(D, t, a_m, a_wt, rvec, cumulative_fission_times) {
   rho = (a_m - a_wt) / (4*pi*D*(t-cumulative_fission_times[1])) * exp(-rvec^2 / (4*D*(t-cumulative_fission_times[1]))) # initial hit solution
   ps = 1
   for (ft in cumulative_fission_times[-1]) {
      if (t<ft) break
      tau_i = t - ft
      rho = rho + a_m / (4*pi*D*tau_i) * exp(-rvec^2 / (4*D*tau_i)) # subsequent fissions
      ps = ps + 1
   }
   return(list("rho" = rho, "patchsize" = ps))
}

diffusion_time_series_given_event_times_and_first_hit = function(c_fis_time, tvec, D, a_m, a_wt, rvec) {
   timeseries_tib = tibble()
   for (tt in tvec) {
      rho_t = create_diffusion_perturbation_with_first_hit_time(D, tt, a_m, a_wt, rvec, c_fis_time)
      timeseries_tib = timeseries_tib %>% 
         bind_rows(tibble(t = tt, rho_pert = rho_t$rho, r = rvec, patchsize = rho_t$patchsize))
   }
   return(timeseries_tib)
}

get_timeseries_from_paths_real_D = function(i, thesepaths, thisage, D, a_m, a_wt, rvec) {
   c_fis_time = thesepaths %>% filter(iter==i) %>% pull(event_time)
   diffusion_time_series_given_event_times_and_first_hit(c_fis_time, thisage, D, a_m, a_wt, rvec)
}

get_timeseries_from_single_path_real_D = function(path, thisage, D, a_m, a_wt, rvec) {
   diffusion_time_series_given_event_times_and_first_hit(path, thisage, D, a_m, a_wt, rvec)
}

load_and_tidy_data = function(mark, patchsize) {
   if (mark=="KDM6A") {
      if (patchsize==5) {
         patch_pack_dat = read_csv("./databases/packing_problem_KDM6A_pa5_070520.csv")
         new_patch_meas = readxl::read_xlsx("./databases/packing_problem_KDM6A_pa5_280520.xlsx")
         additional_patches = NULL
      }
      if (patchsize==10) {
         ## get rid of dodgy mislabelled KDM6A patch
         patch_pack_dat = read_csv("./databases/packing_problem_KDM6A_020520_nr2.csv") %>% filter(Scan_ID!=627731)
         new_patch_meas = readxl::read_xlsx("./databases/packing_problem_KDM6A_pa10_280520.xlsx") %>% filter(Scan_ID!=627731)
         additional_patches = readxl::read_xlsx("./databases/pa10_radii_060620.xlsx") %>% filter(Mark==mark)
      }
   }
   if (mark=="STAG2") {
      if (patchsize==5) {
         patch_pack_dat = readxl::read_xlsx("./databases/STAG2_pa5_120520.xlsx")
         new_patch_meas = readxl::read_xlsx("./databases/STAG2_pa5_270520.xlsx")
         additional_patches = NULL
      }
      if (patchsize==10) {
         patch_pack_dat = readxl::read_xlsx("./databases/STAG2_pa10_110520.xlsx")
         new_patch_meas = readxl::read_xlsx("./databases/STAG2_pa10_270520.xlsx")
         additional_patches = readxl::read_xlsx("./databases/pa10_radii_060620.xlsx") %>% filter(Mark==mark)

      }
   }
   patch_pack_dat = patch_pack_dat %>% 
      bind_rows(dplyr::setdiff(new_patch_meas %>% select(Scan_ID, Patch), 
                               patch_pack_dat %>% select(Scan_ID, Patch)) %>% 
                   left_join(new_patch_meas)) %>% select(-Total_area) %>% 
      left_join(new_patch_meas %>% select(Scan_ID, Mark, Patch, Type, Total_area))
   if (!is.null(additional_patches)) {
      ## append the new stuff from the following file to patch_pack_dat and retain the Total_area column when joining to new_patch_meas?
      patch_pack_dat = patch_pack_dat %>% 
         bind_rows(dplyr::setdiff(additional_patches %>% select(Scan_ID, Type, crypt1_size, crypt1_xy), 
                                  patch_pack_dat %>% select(Scan_ID, Type, crypt1_size, crypt1_xy)) %>% 
                      left_join(additional_patches)) %>% mutate(Patch = 1:n())
   }
   
   ## getting crypt x, y positions
   crypt_xy_names = paste(paste0(rep("crypt", patchsize), as.character(1:patchsize)), "xy", sep = "_")
   crypt_xy_dat = patch_pack_dat %>% select(Mark, Scan_ID, Patch, Type, crypt_xy_names) %>% 
      gather(crypt, xy, -Mark, -Scan_ID, -Patch, -Type) %>% 
      separate(xy, into = c("x", "y"), sep = ", ") %>% 
      mutate(x=as.numeric(x), y=as.numeric(y)) %>% 
      mutate(crypt = str_remove(crypt, "crypt"), crypt = as.integer(str_remove(crypt, "_xy")))
   
   ## getting crypt sizes
   crypt_size_names = paste(paste0(rep("crypt", patchsize), as.character(1:patchsize)), "size", sep = "_")
   crypt_size_dat = patch_pack_dat %>% select(Mark, Scan_ID, Patch, Type, crypt_size_names) %>% 
      gather(crypt, size, -Mark, -Scan_ID, -Patch, -Type) %>% 
      mutate(crypt = str_remove(crypt, "crypt"), crypt = as.integer(str_remove(crypt, "_size")))
   
   ## join and append total patch size
   patch_pack_dat = crypt_xy_dat %>% arrange(Scan_ID, crypt, Patch, Type) %>% 
      left_join(crypt_size_dat %>% arrange(Scan_ID, crypt, Patch, Type)) %>% 
      left_join(patch_pack_dat %>% select(Mark, Scan_ID, Patch, Type, Total_area))
   
   ## find typos in data reporting
   patch_pack_dat %>% group_by(Scan_ID, Patch) %>% 
      mutate(diff_x = abs(x-median(x)), diff_y = abs(y-median(y))) %>% arrange(desc(diff_x))
   patch_pack_dat %>% group_by(Scan_ID, Patch) %>% 
      mutate(diff_x = abs(x-median(x)), diff_y = abs(y-median(y))) %>% arrange(desc(diff_y))
   
   ## join to new patch area measurements
   # new_patch_meas %>% group_by(Scan_ID, Patch) %>% summarise(n=n()) %>% filter(n>1)
   # patch_pack_dat %>% rename(Total_area_old = Total_area) %>%
   #    left_join(new_patch_meas %>% select(Mark, Scan_ID, Patch, Type, Total_area)) %>%
   #    mutate(Delta_area = Total_area - Total_area_old) %>% View()
   # patch_pack_dat = patch_pack_dat %>% select(-Total_area) %>% 
   #    left_join(new_patch_meas %>% select(Mark, Scan_ID, Patch, Type, Total_area))
   
   ## find the neighbourhoods of each mutant patch
   patch_pack_dat = patch_pack_dat %>% arrange(Scan_ID, Patch, Type, crypt) %>% 
      mutate(ismut = as.integer(Type=="mut")) %>% 
      mutate(bit.change=cumsum(c(1, diff(ismut) != 0))) %>% 
      mutate(bit.change = ifelse(ismut==1, bit.change, bit.change-1)) %>% 
      mutate(neighbourhood = as.integer(as.factor(bit.change))) %>% 
      select(-ismut, -bit.change)
   
   scans = patch_pack_dat %>% pull(Scan_ID) %>% unique()
   nbrhood = 1
   outtib = tibble()
   for (s in scans) {
      assigners = patch_pack_dat %>% filter(Scan_ID==s) %>% filter(Type!="mut") %>% 
         group_by(Patch) %>% summarise(mean_x = mean(x), mean_y = mean(y))
      targets = patch_pack_dat %>% filter(Scan_ID==s) %>% filter(Type=="mut") %>% 
         group_by(Patch) %>% summarise(mean_x = mean(x), mean_y = mean(y))
      nbrhoods = RANN::nn2(targets %>% select(mean_x, mean_y), assigners %>% select(mean_x, mean_y), k = 1)$nn.idx[,1]
      mut_ps = targets %>% pull(Patch)
      match_ps = assigners %>% pull(Patch)
      for (pp in 1:length(mut_ps)) {
         whichmatch = which(nbrhoods==pp)
         outtib = outtib %>% bind_rows(
            patch_pack_dat %>% filter(Scan_ID==s) %>% filter(Type=="mut") %>% filter(Patch==mut_ps[pp]) %>% 
            bind_rows(patch_pack_dat %>% filter(Scan_ID==s) %>% filter(Type!="mut") %>% 
                         filter(Patch%in%match_ps[whichmatch])) %>% 
            mutate(neighbourhood = nbrhood))
         nbrhood = nbrhood + 1
      }
   }
   outtib = outtib %>% arrange(Mark, Scan_ID, neighbourhood, Patch, Type)
   return(list("data" = outtib, "mark" = mark, "patchsize" = patchsize))
}

get_per_crypt_whitespace = function(patch_container) {
   # -- knn on xy within patches,
   # -- approximate as circles,
   # -- find average "white space" left over by distance - radii
   patch_pack_dat = patch_container$data
   scan_ids = patch_pack_dat %>% pull(Scan_ID) %>% unique()
   outtib = tibble()
   for (sc in scan_ids) {
      subdat = patch_pack_dat %>% filter(Scan_ID==sc)
      patch_nums = subdat %>% pull(Patch) %>% unique()
      for (pn in patch_nums) {
         thisdat = subdat %>% filter(Patch==pn)
         knn_obj = thisdat %>% select(x,y) %>% RANN::nn2(., query = ., k = 3)
         radii = thisdat %>% mutate(radius = sqrt(size / pi)) %>% select(crypt, radius)
         knn_obj[["radii_sums"]] = knn_obj$nn.idx
         knn_obj[["white_space"]] = knn_obj$nn.idx
         for (i in 1:nrow(knn_obj$nn.idx)) {
            for (j in 2:ncol(knn_obj$nn.idx)) {
               knn_obj$radii_sums[i,j] = (radii %>% filter(crypt==i) %>% pull(radius)) + 
                  (radii %>% filter(crypt==knn_obj$nn.idx[i,j]) %>% pull(radius))
               knn_obj$white_space[i,j] = max(0, knn_obj$nn.dists[i,j] - knn_obj$radii_sums[i,j])
            }
         }
         # ws = knn_obj$white_space[,2] %>% as.vector()
         ws = knn_obj$white_space[,-1] %>% t() %>% as_tibble(.name_repair = "unique") %>% 
            summarise_all(~min(.)) %>% as.matrix() %>% as.vector()
         outtib = outtib %>% 
            bind_rows(thisdat %>% select(Mark, Scan_ID, Patch, 
                                         Type, Total_area, neighbourhood) %>% 
                         distinct() %>% 
                         mutate(med_ws = median(ws), mean_ws = mean(ws),
                                l95_ws = max(0, quantile(ws, 0.025)), u95_ws = quantile(ws, 0.975),
                                l80_ws = max(0, quantile(ws, 0.1)), u80_ws = quantile(ws, 0.9),
                                l50_ws = max(0, quantile(ws, 0.25)), u50_ws = quantile(ws, 0.75)))
      }
   }
   ## get scale per slide (might not be good if multiple "neighbourhoods" per slide)
   slide_scale = outtib %>% filter(Type=="WT") %>% group_by(Scan_ID, neighbourhood) %>% 
      summarise(scale = median(med_ws))
   outtib = outtib %>% left_join(slide_scale) %>% 
      mutate(scaled_med = med_ws/scale, scaled_mean = mean_ws/scale, 
             scaled_l95 = l95_ws/scale, scaled_u95 = u95_ws/scale, 
             scaled_l80 = l80_ws/scale, scaled_u80 = u80_ws/scale, 
             scaled_l50 = l50_ws/scale, scaled_u50 = u50_ws/scale)
   return(outtib %>% mutate(patchsize = patch_container$patchsize))
}

get_total_whitespace_fraction = function(patch_container) {
   tot_ws_fracs = patch_container$data %>% group_by(Mark, Scan_ID, Type, neighbourhood, Patch) %>% 
      summarise(white_space_frac = 1 - sum(size)/Total_area[1]) %>% ungroup()
   ps_scale = tot_ws_fracs %>% filter(Type=="WT") %>% group_by(Mark, Scan_ID, neighbourhood) %>% 
      summarise(wsf_scale = median(white_space_frac))
   tot_ws_fracs = tot_ws_fracs %>% left_join(ps_scale) %>% 
      mutate(scaled_wsf = white_space_frac / wsf_scale) %>% 
      mutate(patchsize = patch_container$patchsize)
   return(tot_ws_fracs)
}

subsector_params_from_patchdat = function(patchdat) {
   mut_centroids = patchdat$data %>% group_by(Scan_ID, neighbourhood, Type, Patch) %>% 
      summarise(mut_centroid_x = mean(x), mut_centroid_y = mean(y)) %>% ungroup() %>% 
      filter(Type=="mut") %>% select(-Type, -Patch)
   wt_subsectors = patchdat$data %>% group_by(Scan_ID, neighbourhood, Type, Patch) %>% 
      summarise(centroid_x = mean(x), centroid_y = mean(y), area = Total_area[1], n_crypts = n()) %>% ungroup() %>% 
      mutate(Radius = sqrt(area/pi)) %>% left_join(mut_centroids) %>% 
      group_by(Scan_ID, neighbourhood, Type, Patch) %>% 
      mutate(distance_from_mut = sqrt((centroid_x - mut_centroid_x)^2 + (centroid_y - mut_centroid_y)^2),
             theta_subtended = ifelse(distance_from_mut==0, 2*pi, 2 * atan(Radius / distance_from_mut)),
             Rin = ifelse(distance_from_mut==0, 0, max(0, distance_from_mut - 0.5*pi*Radius^2 / (distance_from_mut * theta_subtended))), 
             Rout = sqrt(Rin^2 + 2 * area / theta_subtended),
             subsector_area = 0.5 * theta_subtended * (Rout^2 - Rin^2)) %>% ungroup() %>% 
      select(Scan_ID, neighbourhood, Type, Patch, n_crypts, theta_subtended, Rin, Rout)
   wt_subsectors %>% arrange(Scan_ID, neighbourhood, Type) %>% ungroup()
}

space_params_from_patchdat = function(patchdat) {
   patchdat$data %>% group_by(Scan_ID, Patch, Type, neighbourhood, Total_area) %>% 
      summarise(total_crypt_area = sum(size), n_crypts = n()) %>% 
      mutate(gamma = Total_area - total_crypt_area,
             space_per_crypt = Total_area / n_crypts,
             average_crypt_area = total_crypt_area/n_crypts,
             gamma_per_crypt = space_per_crypt - average_crypt_area) %>% 
      select(Scan_ID, neighbourhood, Type, Patch, n_crypts, Total_area,
             space_per_crypt, average_crypt_area, gamma_per_crypt) %>% ungroup()
}

mutant_patch_area_increase = function(spacedat, percrypt_space) { 
   spacedat %>% left_join(percrypt_space %>% rename(initial_area = space_per_crypt) %>% 
                             select(Scan_ID, initial_area)) %>% 
      mutate(area_increase = Total_area - initial_area, 
             area_increase_frac = Total_area / initial_area) %>% 
      select(Scan_ID, neighbourhood, Type, Patch, n_crypts, 
             initial_area, Total_area, area_increase) %>% ungroup() %>% 
      filter(Type == "mut")
}

u_R = function(r, gamma_a, Dtilde, epsilon) {
   return(epsilon * gamma_a / (pi * Dtilde) * exp(-r^2 / Dtilde))
}

Gamma_perturbation = function(gamma_a, Dtilde, epsilon, Rin, Rout, theta, order="full") {
   ## compute integral around circle radius R of the perturbation of whitespace density
   if (order=="eps") {
      ## just take the order epsilon term, epsilon << 1
      return(-0.5*theta/pi * epsilon * gamma_a^2 * (exp(-Rin^2 / Dtilde) - exp(-Rout^2 / Dtilde)))
   } else if (order=="full") {
      ## take the full perturbation
      xin = u_R(Rin, gamma_a, Dtilde, epsilon)
      xout = u_R(Rout, gamma_a, Dtilde, epsilon)
      return(-0.5*theta*gamma_a*Dtilde * log((1 + xin)/(1  + xout)))
   } else {
      print("order must be either 'full' or 'exp'.")
      return(NA)
   }
}

Gamma_total = function(gamma_a, Dtilde, epsilon, R1, R2, theta, order="full") {
   ## compute integral around circle radius R of the total whitespace density
   Rin = min(R1, R2)
   Rout = max(R1, R2)
   if (Rin==Rout) return(0)
   return(0.5*theta*(Rout^2 - Rin^2) * gamma_a + 
             Gamma_perturbation(gamma_a, Dtilde, epsilon, Rin, Rout, theta, order))
}

gamma_r = function(gamma_a, r, Dtilde, epsilon, order="full") {
   ## compute the white space density as a negative space density diffuses out from r=0
   ## with diffusive scale Dtilde = 4 D t
   ## and negative space density perturbation epsilon = rho_in == Am_incr
   ## Dtilde == <R^2>, mean square displacement of diffusion
   if (order=="eps") {
      ## assumes epsilon << 1
      return(gamma_a + epsilon * gamma_tilde_eps(gamma_a, r, Dtilde, epsilon))
   } else if (order=="full") {
      ## no assumption on size of epsilone
      return(gamma_a / (1 + epsilon*gamma_a/(pi*Dtilde) *  exp(-r^2 / Dtilde)))
   } else {
      print("order must be either 'full' or 'exp'.")
      return(NA)
   }
}

gamma_tilde_eps = function(gamma_a, r, Dtilde, epsilon) {
   ## Dtilde == <R^2>, mean square displacement of diffusion
   ## the first-order term in epsilon, if epsilon << 1
   return(-gamma_a^2 / (pi * Dtilde) * exp(-r^2 / Dtilde))
}

rho_r = function(gamma_a, r, Dtilde, epsilon) {
   ## compute the negative space density 
   ## as a perturbation of size epsilon diffuses out from r=0
   ## with diffusive timescale Dtilde = 4 D t
   ## rho = 1 / gamma
   ## that is, the amount of area per unit white space
   return(1/gamma_a + epsilon * rho_tilde(r, Dtilde))
}

rho_tilde = function(r, Dtilde) {
   ## the r dependence of the negative space density perturbation
   return(1/(pi*Dtilde) *  exp(-r^2 / Dtilde))   
}

outer_circle_radius = function(N_wt, s_wt) {
   return(sqrt(N_wt * s_wt / pi))
}

Dtilde_from_eps_order = function(gamma_a, epsilon, N_wt, s_wt, Am_incr) {
   ## equate the loss of white space due to negative density perturbation
   ## in a circle radius R to the mutant patch area increase, Am_increase,
   ## in the epsilon<<1 limit we can rearrange for Dtilde.
   ## N_wt: number of WT crypts in range
   ## s_wt: the WT crypt space (crypt area + white space)
   ## (possibly have )
   return(-outer_circle_radius(N_wt, s_wt)^2 / log(1 - Am_incr/(epsilon*gamma_a^2)) )
}

epsilon_from_mut_crypt_area = function(a_wt, a_m, N_m) {
   ## approximate the "point source" of extra area we must look to find a 
   ## unit of white space due to the mutant patch of size N_m
   ## a_wt: WT crypt area (not crypt space, just the area)
   ## a_m: mut crypt area
   ## N_m: mut patch size
   return((N_m - 1)*a_m + (a_m - a_wt))
}

Rinfinity = function(gamma_a, Dtilde, epsilon, delta) {
   ## The "range of effect" with tolerance delta
   ## such that gamma(R_inf)/gamma_a ~ delta
   # return(sqrt(Dtilde * log( epsilon*gamma_a/(pi*Dtilde) * (1-delta)/delta)))
   return(sqrt(Dtilde * log(epsilon*gamma_a/(delta*pi*Dtilde))))
}

grad_gammar = function(gamma_a, r, Dtilde, epsilon) {
   expterm = exp(-r^2 / Dtilde)
   preterm = epsilon * gamma_a / (pi * Dtilde)
   return((2 * gamma_a * r / Dtilde  * preterm * expterm) / (1 + preterm * expterm)^2)
}

range_finder = function(gamma_a, D, t, R, lambda, a_m, a_wt, X) {
   A = gamma_a/(4*pi*D*t) * ((exp(lambda*t) - 1)*a_m + (a_m - a_wt))
   B = exp(-R^2 / (4*D*t)) * (exp((1-X)*R^2 / (4*D*t)) - 1)
   C = (1 - exp(-X*R^2 / (4*D*t)))
   return(A*B - C)
}

prime_data = function(wt_patches, patchdat, account_for_edge_fluff=FALSE, mf=2/3) {
   # get ambient summary stats
   patch_outlines = wt_patches %>% filter(ROI=="Polygon")
   patch_crypts = wt_patches %>% filter(ROI=="Ellipse")
   patch_crypt_areas = patch_crypts %>% group_by(Scan_ID) %>% 
      summarise(total_crypt_area = sum(area), n_crypts = n())
   patch_whitespaces = patch_outlines %>% left_join(patch_crypt_areas) %>% 
      mutate(Gamma = area - total_crypt_area) %>% 
      select(Scan_ID, ROI, area, total_crypt_area, Gamma, n_crypts)
   if (account_for_edge_fluff==TRUE) {
      patch_whitespaces = patch_whitespaces %>% group_by(Scan_ID) %>% 
         mutate(real_gamma_pc = Gamma / (n_crypts - mf * hexagonal_patch_edge_length(n_crypts)),
                Gamma = real_gamma_pc * n_crypts,
                area = total_crypt_area + Gamma) %>% select(-real_gamma_pc) %>% ungroup()
   }
   percrypt_space = patch_whitespaces %>% mutate(space_per_crypt = area / n_crypts,
                                                 average_crypt_area = total_crypt_area/n_crypts,
                                                 gamma_per_crypt = space_per_crypt - average_crypt_area) %>% 
      select(Scan_ID, n_crypts, space_per_crypt, average_crypt_area, gamma_per_crypt)
   ambient_values = patch_whitespaces %>% 
      mutate(gamma_a = Gamma / area, s_wt = area / n_crypts, 
             a_wt = total_crypt_area / n_crypts) %>% select(Scan_ID, gamma_a, s_wt, a_wt)
   
   # get slide and patch stats
   subsector_dat = subsector_params_from_patchdat(patchdat)
   space_dat = space_params_from_patchdat(patchdat)
   Amincr_dat = mutant_patch_area_increase(space_dat, percrypt_space) %>% 
      select(Scan_ID, neighbourhood, n_crypts, area_increase) %>% rename(Am_incr = area_increase)
   mutarea_dat = space_dat %>% filter(Type=="mut") %>% 
      select(Scan_ID, neighbourhood, n_crypts, average_crypt_area) %>% rename(a_m = average_crypt_area)
   total_whitespace = space_dat %>% mutate(Gamma_dat = gamma_per_crypt * n_crypts) %>% 
      select(Scan_ID, neighbourhood, Type, Patch, n_crypts, Total_area, Gamma_dat)
   
   # join
   data_primed = subsector_dat %>% left_join(ambient_values) %>% 
      left_join(Amincr_dat) %>% left_join(mutarea_dat) %>% left_join(total_whitespace)
   # calculate mass input epsilon
   data_primed = data_primed %>% mutate(epsilon = epsilon_from_mut_crypt_area(a_wt, a_m, n_crypts))
   # calculate the null model of ambient white space everywhere
   data_primed = data_primed %>% group_by(Scan_ID, neighbourhood, Patch) %>% 
      mutate(Gamma_null = Total_area*gamma_a) %>% ungroup()
   return(list("data_primed" = data_primed, "space_data" = space_dat,
               "percrypt_space" = percrypt_space, "patch_whitespaces" = patch_whitespaces))
}

predict_patch_packing = function(data_primed, Rrmsq) {
   data_primed %>% group_by(Scan_ID, neighbourhood, Patch) %>% 
      mutate(Gamma_pred = Gamma_total(gamma_a, Rrmsq^2, epsilon, Rin, Rout, 
                                      theta_subtended, order="full")) %>%
      mutate(abs_err = abs(Gamma_dat - Gamma_pred)) %>% ungroup()
}

find_min_MSE_Rrmsq = function(data_primed) {
   Rrmsq_l = pracma::linspace(1, 2000, 5000)
   outtib = tibble()
   outtib2 = tibble()
   for (Rrmsq in Rrmsq_l) {
      ## also predict the mutant patch packing itself! This will prevent the 
      ## minimum from being R=1 before the mass can spread out to the WT patches....
      ## actually if R~size of mutant patch, only check prediction error on mutant patch itself?
      this_err = predict_patch_packing(data_primed %>% ungroup(), Rrmsq) %>% 
         summarise(RMSE = mean(abs_err), RMedSE = median(abs_err)) %>% mutate(Rrms = Rrmsq)
      this_err2 = predict_patch_packing(data_primed %>% ungroup(), Rrmsq) %>% 
         filter(abs_err<quantile(abs_err, 0.9)) %>% 
         summarise(RMSE = mean(abs_err), RMedSE = median(abs_err)) %>% mutate(Rrms = Rrmsq)
      outtib = outtib %>% bind_rows(this_err)
      outtib2 = outtib2 %>% bind_rows(this_err2)
      # do something less sensitive to outliers but smoother than median: mean of bottom 90%?
   }
   outtib %>% ggplot() + geom_line(aes(x=Rrms, y=RMSE))
   outtib %>% ggplot() + geom_line(aes(x=Rrms, y=RMedSE))
   outtib2 %>% ggplot() + geom_line(aes(x=Rrms, y=RMSE))
   outtib2 %>% ggplot() + geom_line(aes(x=Rrms, y=RMedSE))
   
   ## take minimum MSE and predict packing problem, plot comparison with data
   # Rrmsq_min = outtib %>% arrange(RMSE) %>% slice(1) %>% pull(Rrms)
   Rrmsq_min = outtib2 %>% arrange(RMSE) %>% slice(1) %>% pull(Rrms)
   return(Rrmsq_min)
}

calc_diffusion_coefficient = function(Rrmsq_min, mark) {
   ## calculate diffusion coefficient from Rrmsq and patch life times
   paste0("./output/first_passage_time_",
          paste0(mark, "_1-5.rds"))
   fpt_1_5 = readRDS(paste0("./output/first_passage_time_",
                            paste0(mark, "_1-5.rds")))
   fpt_5_6 = readRDS(paste0("./output/first_passage_time_",
                            paste0(mark, "_5-6.rds")))
   fpt_5_10 = readRDS(paste0("./output/first_passage_time_",
                             paste0(mark, "_5-10.rds")))
   patch_lifetime_5 = fpt_1_5$fpt_summary %>% pull(mean_fpt) + 0.5 * fpt_5_6$fpt_summary %>% pull(mean_fpt)
   patch_lifetime_5_lower = fpt_5_10$fpt_summary %>% pull(mean_fpt) + 0.5 * fpt_5_6$fpt_summary %>% pull(mean_fpt)
   
   D = Rrmsq_min^2 / (4 * patch_lifetime_5)
   D_lower = Rrmsq_min^2 / (4 * patch_lifetime_5_lower)
   return(list("D_1_5" = D, "D_5_10" = D_lower))
}

calculate_new_Rrms_from_D = function(Ds, mark) {
   fpt_1_10 = readRDS(paste0("./output/first_passage_time_",
                             paste0(mark, "_1-10.rds")))
   fpt_10_11 = readRDS(paste0("./output/first_passage_time_",
                              paste0(mark, "_10-11.rds")))
   patch_lifetime_10 = fpt_1_10$fpt_summary %>% pull(mean_fpt) + 0.5 * fpt_10_11$fpt_summary %>% pull(mean_fpt)
   Rrmsq_new = sqrt(Ds$D_1_5 * 4 * patch_lifetime_10)
   Rrmsq_new_lower = sqrt(Ds$D_5_10 * 4 * patch_lifetime_10)
   return(list("Rrmsq_1_5" = Rrmsq_new, "Rrmsq_5_10" = Rrmsq_new_lower))
}

generate_stan_ins = function(data_primed_stag, data_primed_kdm, length_scale, 
                             area_scale, t_gam_pars_stag, t_gam_pars_kdm, slide_age_ref) {
   joinscaledat = data_primed_kdm$data_primed %>% mutate(mark = "KDM6A") %>% 
      bind_rows(data_primed_stag$data_primed %>% mutate(mark = "STAG2")) %>% 
      mutate(Rin = Rin / length_scale, Rout = Rout / length_scale,
             s_wt = s_wt / area_scale, a_wt = a_wt / area_scale, Am_incr = Am_incr / area_scale,
             a_m = a_m / area_scale, Total_area = Total_area / area_scale, 
             Gamma_dat = Gamma_dat / area_scale, epsilon = epsilon / area_scale, 
             Gamma_null = Gamma_null / area_scale)
   
   joinscaledat = joinscaledat %>% arrange(mark, Scan_ID, neighbourhood, Type) %>% 
      mutate(ismut = as.integer(Type=="mut")) %>% 
      mutate(bit.change=cumsum(c(1, diff(ismut) != 0))) %>% 
      mutate(bit.change = ifelse(ismut==1, bit.change, bit.change-1)) %>% 
      mutate(nbrhood_index = as.integer(as.factor(bit.change))) %>% 
      select(-ismut, -bit.change)
   
   # join ages
   joinscaledat = joinscaledat %>% left_join(slide_age_ref %>% select(Scan_ID, Age))
   age_M = joinscaledat %>% select(nbrhood_index, Age) %>% distinct() %>% arrange(nbrhood_index) %>% pull(Age)
   
   n_nbrhoods = length(joinscaledat %>% pull(nbrhood_index) %>% unique())
   slide_index_tib = joinscaledat %>% select(mark, Scan_ID, nbrhood_index, gamma_a) %>% distinct() %>% 
      mutate(slide_index = as.integer(as.factor(Scan_ID)),
             mark_index = ifelse(mark=="STAG2", 1, 2))
   wt30_gamma_a_tib = slide_index_tib %>% select(Scan_ID, gamma_a, slide_index) %>% distinct()
   n_slides = nrow(wt30_gamma_a_tib)
   
   tgam_shapes = c(t_gam_pars_stag$alpha %>% pull(map_est), t_gam_pars_kdm$alpha %>% pull(map_est))
   tgam_rates = c(t_gam_pars_stag$beta %>% pull(map_est), t_gam_pars_kdm$beta %>% pull(map_est))
   
   input_list = list(N = nrow(joinscaledat), 
                     Gamma_dat = joinscaledat %>% pull(Gamma_dat),
                     Rin = joinscaledat %>% pull(Rin),
                     Rout = joinscaledat %>% pull(Rout),
                     theta = joinscaledat %>% pull(theta_subtended),
                     M = n_nbrhoods,
                     epsilon = joinscaledat %>% select(epsilon, nbrhood_index) %>% distinct() %>% pull(epsilon),
                     nbrhood_index = joinscaledat %>% pull(nbrhood_index),
                     slide_index = slide_index_tib %>% pull(slide_index),
                     S = n_slides,
                     mark_index = slide_index_tib %>% pull(mark_index),
                        a_gam_t = tgam_shapes, b_gam_t = tgam_rates,
                     pat_ages = age_M)
                     # gamma_a_wt30 = wt30_gamma_a_tib %>% pull(gamma_a),
   
   ## initialisation
   # t_p_init = slide_index_tib %>% group_by(nbrhood_index) %>% 
   #    mutate(tp_samp = rgamma(1, shape = tgam_shapes[mark_index], rate = tgam_rates[mark_index])) %>% 
   #    pull(tp_samp)
   t_p_init = rbeta(n_nbrhoods, 1, 1)
   init_list = list(gamma_a = joinscaledat %>% select(gamma_a, nbrhood_index) %>% 
                       distinct() %>% pull(gamma_a),
                    D = fdrtool::rhalfnorm(1, theta=fdrtool::sd2theta(1)),
                    t_p = t_p_init,
                    sd_Gacc = fdrtool::rhalfnorm(1, fdrtool::sd2theta(1)),
                    a1_gpop = rgamma(1, 10, 1),
                    a2_gpop = rgamma(1, 10, 1))
                    # a_gpat = rgamma(n_slides, 10, 1), 
                    # b_gpat = rgamma(n_slides, 10, 1),
                    # b1_gpop = rgamma(1, 1, 1),
                    # b2_gpop = rgamma(1, 1, 1))
   stan_ins = list(input = input_list, init = init_list)
   return(list("stan_ins" = stan_ins, "joinscaledat" = joinscaledat, "slide_index_tib" = slide_index_tib))
}

generate_stan_ins2 = function(data_primed_stag, data_primed_kdm, length_scale, 
                             area_scale, slide_age_ref) {
   joinscaledat = data_primed_kdm$data_primed %>% mutate(mark = "KDM6A") %>% 
      bind_rows(data_primed_stag$data_primed %>% mutate(mark = "STAG2")) %>% 
      mutate(Rin = Rin / length_scale, Rout = Rout / length_scale,
             s_wt = s_wt / area_scale, a_wt = a_wt / area_scale, Am_incr = Am_incr / area_scale,
             a_m = a_m / area_scale, Total_area = Total_area / area_scale, 
             Gamma_dat = Gamma_dat / area_scale, epsilon = epsilon / area_scale, 
             Gamma_null = Gamma_null / area_scale)
   
   joinscaledat = joinscaledat %>% arrange(mark, Scan_ID, neighbourhood, Type) %>% 
      mutate(ismut = as.integer(Type=="mut")) %>% 
      mutate(bit.change=cumsum(c(1, diff(ismut) != 0))) %>% 
      mutate(bit.change = ifelse(ismut==1, bit.change, bit.change-1)) %>% 
      mutate(nbrhood_index = as.integer(as.factor(bit.change))) %>% 
      select(-ismut, -bit.change)
   
   # join ages
   joinscaledat = joinscaledat %>% left_join(slide_age_ref %>% select(Scan_ID, Age))
   age_M = joinscaledat %>% select(nbrhood_index, Age) %>% distinct() %>% arrange(nbrhood_index) %>% pull(Age)
   
   n_nbrhoods = length(joinscaledat %>% pull(nbrhood_index) %>% unique())
   slide_index_tib = joinscaledat %>% select(mark, Scan_ID, nbrhood_index, gamma_a) %>% distinct() %>% 
      mutate(slide_index = as.integer(as.factor(Scan_ID)),
             mark_index = ifelse(mark=="STAG2", 1, 2))
   wt30_gamma_a_tib = slide_index_tib %>% select(Scan_ID, gamma_a, slide_index) %>% distinct()
   n_slides = nrow(wt30_gamma_a_tib)
   
   # min_ages = min(c(timeatN_stag, timeatN_kdm))
   
   input_list = list(N = nrow(joinscaledat), 
                     Gamma_dat = joinscaledat %>% pull(Gamma_dat),
                     Rin = joinscaledat %>% pull(Rin),
                     Rout = joinscaledat %>% pull(Rout),
                     theta = joinscaledat %>% pull(theta_subtended),
                     M = n_nbrhoods,
                     epsilon = joinscaledat %>% select(epsilon, nbrhood_index) %>% distinct() %>% pull(epsilon),
                     nbrhood_index = joinscaledat %>% pull(nbrhood_index),
                     slide_index = slide_index_tib %>% pull(slide_index),
                     S = n_slides,
                     mark_index = slide_index_tib %>% pull(mark_index),
                     pat_ages = age_M)
   # min_ages = min_ages,
   # gamma_a_wt30 = wt30_gamma_a_tib %>% pull(gamma_a),
   
   ## initialisation
   # t_p_init = slide_index_tib %>% group_by(nbrhood_index) %>% 
   #    mutate(tp_samp = rgamma(1, shape = tgam_shapes[mark_index], rate = tgam_rates[mark_index])) %>% 
   #    pull(tp_samp)
   t_p_init = age_M/2
   init_list = list(gamma_a = joinscaledat %>% select(gamma_a, nbrhood_index) %>% 
                       distinct() %>% pull(gamma_a),
                    D = fdrtool::rhalfnorm(1, theta=fdrtool::sd2theta(2)),
                    t_p = t_p_init,
                    sd_Gacc = fdrtool::rhalfnorm(1, fdrtool::sd2theta(1)),
                    a1_gpop = rgamma(1, 10, 1),
                    a2_gpop = rgamma(1, 10, 1))
   # a_gpat = rgamma(n_slides, 10, 1), 
   # b_gpat = rgamma(n_slides, 10, 1),
   # b1_gpop = rgamma(1, 1, 1),
   # b2_gpop = rgamma(1, 1, 1))
   stan_ins = list(input = input_list, init = init_list)
   return(list("stan_ins" = stan_ins, "joinscaledat" = joinscaledat, "slide_index_tib" = slide_index_tib))
}

generate_stan_ins3 = function(data_primed_stag, data_primed_kdm, length_scale, area_scale, slide_age_ref) {
   joinscaledat = data_primed_kdm$data_primed %>% mutate(mark = "KDM6A") %>% 
      bind_rows(data_primed_stag$data_primed %>% mutate(mark = "STAG2")) %>% 
      mutate(Rin = Rin / length_scale, Rout = Rout / length_scale,
             s_wt = s_wt / area_scale, a_wt = a_wt / area_scale, Am_incr = Am_incr / area_scale,
             a_m = a_m / area_scale, Total_area = Total_area / area_scale, 
             Gamma_dat = Gamma_dat / area_scale, epsilon = epsilon / area_scale, 
             Gamma_null = Gamma_null / area_scale)
   
   joinscaledat = joinscaledat %>% arrange(mark, Scan_ID, neighbourhood, Type) %>% 
      mutate(ismut = as.integer(Type=="mut")) %>% 
      mutate(bit.change=cumsum(c(1, diff(ismut) != 0))) %>% 
      mutate(bit.change = ifelse(ismut==1, bit.change, bit.change-1)) %>% 
      mutate(nbrhood_index = as.integer(as.factor(bit.change))) %>% 
      select(-ismut, -bit.change)
   
   # join ages
   joinscaledat = joinscaledat %>% left_join(slide_age_ref %>% select(Scan_ID, Age))
   age_M = joinscaledat %>% select(nbrhood_index, Age) %>% distinct() %>% arrange(nbrhood_index) %>% pull(Age)
   
   n_nbrhoods = length(joinscaledat %>% pull(nbrhood_index) %>% unique())
   slide_index_tib = joinscaledat %>% select(mark, Scan_ID, nbrhood_index, gamma_a) %>% distinct() %>% 
      mutate(slide_index = as.integer(as.factor(Scan_ID)),
             mark_index = ifelse(mark=="STAG2", 1, 2))
   wt30_gamma_a_tib = slide_index_tib %>% select(Scan_ID, gamma_a, slide_index) %>% distinct()
   n_slides = nrow(wt30_gamma_a_tib)
   
   # get inputs for beta priors on [0,1] patch age multipliers
   tp_raw_prior_tib = tibble()
   for (thismark in (slide_index_tib %>% select(mark, mark_index) %>% distinct() %>% arrange(mark_index) %>% pull(mark))) {
      fit_tib = readRDS(paste0("./output/patch_age_priorbetafits_", paste0(thismark, ".rds")))
      sum_fit_tib = fit_tib %>% summarise(alpha=median(alpha), beta=median(beta))
      tp_raw_prior_tib = tp_raw_prior_tib %>% bind_rows(sum_fit_tib)
   }
   
   input_list = list(N = nrow(joinscaledat), 
                     Gamma_dat = joinscaledat %>% pull(Gamma_dat),
                     Rin = joinscaledat %>% pull(Rin),
                     Rout = joinscaledat %>% pull(Rout),
                     theta = joinscaledat %>% pull(theta_subtended),
                     M = n_nbrhoods,
                     epsilon = joinscaledat %>% select(epsilon, nbrhood_index) %>% distinct() %>% pull(epsilon),
                     nbrhood_index = joinscaledat %>% pull(nbrhood_index),
                     slide_index = slide_index_tib %>% pull(slide_index),
                     S = n_slides,
                     mark_index = slide_index_tib %>% pull(mark_index),
                     pat_ages = age_M,
                     alpha_patchage = tp_raw_prior_tib %>% pull(alpha),
                     beta_patchage = tp_raw_prior_tib %>% pull(beta))
   
   ## initialisation
   t_p_raw_init = rep(0.5, n_nbrhoods)
   init_list = list(gamma_a = joinscaledat %>% select(gamma_a, nbrhood_index) %>% 
                       distinct() %>% pull(gamma_a),
                    D = fdrtool::rhalfnorm(1, theta=fdrtool::sd2theta(2)),
                    t_p_raw = t_p_raw_init,
                    sd_Gacc = fdrtool::rhalfnorm(1, fdrtool::sd2theta(1)),
                    a1_gpop = rgamma(1, 10, 1),
                    a2_gpop = rgamma(1, 10, 1))
   stan_ins = list(input = input_list, init = init_list)
   return(list("stan_ins" = stan_ins, "joinscaledat" = joinscaledat, "slide_index_tib" = slide_index_tib))
}

generate_stan_ins_simdata = function(full_data, maxpatchsize, hitmark, n_grid=50, Nstem=7, Ncrypts=1e7) {
   ## stan inputs
   joinscaledat = full_data
   joinscaledat = joinscaledat %>% mutate(nbrhood_index = neighbourhood)
   age_M = joinscaledat %>% select(nbrhood_index, age) %>% distinct() %>% arrange(nbrhood_index) %>% pull(age)
   n_nbrhoods = length(joinscaledat %>% pull(nbrhood_index) %>% unique())
   mark_index = rep(1, n_nbrhoods)
   
   # inputs for event time generation and numerical integration
   num_events = maxpatchsize
   mark_data = read_csv("../drift_All_Marks/exp_coeff.csv")
   mark_data = mark_data %>% filter(mark==hitmark)
   fisrates = c(mark_data %>% pull(rho), 0)
   alf = mark_data %>% pull(alpha)
   deltaC = mark_data %>% pull(mean_slope)
   alf_Nstem_deltaC_Ncrypts = c(alf * deltaC * Nstem * Ncrypts, 0)
   
   input_list = list(N = nrow(joinscaledat), 
                     Gamma_dat = joinscaledat %>% pull(Gamma_dat),
                     Rin = joinscaledat %>% pull(Rin),
                     Rout = joinscaledat %>% pull(Rout),
                     theta = joinscaledat %>% pull(theta_subtended),
                     M = n_nbrhoods,
                     am = joinscaledat %>% select(a_m, nbrhood_index) %>% distinct() %>% pull(a_m),
                     awt = joinscaledat %>% select(a_wt, nbrhood_index) %>% distinct() %>% pull(a_wt),
                     nbrhood_index = joinscaledat %>% pull(nbrhood_index),
                     mark_index = mark_index,
                     pat_ages = age_M,
                     fisrate = fisrates,
                     alf_Nstem_deltaC_Ncrypts = alf_Nstem_deltaC_Ncrypts,
                     n_grid = n_grid,
                     num_events = num_events)
   
   ## initialisation
   # simulate initial good paths for each mutant patch
   # event_times = matrix(NA, n_nbrhoods, num_events)
   event_times = matrix(NA, n_nbrhoods, num_events+1)
   for (i in 1:n_nbrhoods) {
      done = FALSE
      while (done==FALSE) {
         # hit = 999
         # while (hit>(age_M[i])) hit = next_fission_time(alf_Nstem_deltaC_Ncrypts[mark_index[i]], 1, runif(1))
         hit = age_M[i] - runif(1, 15, min(100, age_M[i]))
         posspath = generate_fission_event_times_given_hit(fisrates[mark_index[i]], age_M[i], hittime = hit)
         if (length(posspath)==num_events) {
            done = TRUE
            # truepath = rep(0, num_events)
            posspath = c(posspath, age_M[i] + runif(1, 1, 15))
            truepath = rep(0, num_events+1)
            truepath[1] = hit
            # for (k in 2:num_events) truepath[k] = posspath[k] - posspath[k-1]
            for (k in 2:(num_events+1)) truepath[k] = posspath[k] - posspath[k-1]
            event_times[i,] = truepath
         }
      }
   }
   init_list = list(gamma_a = joinscaledat %>% select(gamma_a, nbrhood_index) %>% 
                       distinct() %>% pull(gamma_a),
                    D = fdrtool::rhalfnorm(1, theta=fdrtool::sd2theta(2)),
                    sd_Gacc = fdrtool::rhalfnorm(1, fdrtool::sd2theta(1)),
                    a1_gpop = rgamma(1, 10, 1),
                    a2_gpop = rgamma(1, 10, 1),
                    event_times = event_times)
   stan_ins = list(input = input_list, init = init_list)
   return(list("stan_ins" = stan_ins, "joinscaledat" = joinscaledat))
}

generate_stan_ins_simdata2 = function(full_data, maxpatchsize, hitmark, n_grid=50, Nstem=7, Ncrypts=1e7) {
   ## stan inputs
   joinscaledat = full_data
   joinscaledat = joinscaledat %>% mutate(nbrhood_index = neighbourhood)
   age_M = joinscaledat %>% select(nbrhood_index, age) %>% distinct() %>% arrange(nbrhood_index) %>% pull(age)
   n_nbrhoods = length(joinscaledat %>% pull(nbrhood_index) %>% unique())
   mark_index = rep(1, n_nbrhoods)
   
   # inputs for event time generation and numerical integration
   num_events = maxpatchsize
   mark_data = read_csv("../drift_All_Marks/exp_coeff.csv")
   mark_data = mark_data %>% filter(mark==hitmark)
   fisrates = c(mark_data %>% pull(rho), 0)
   alf = mark_data %>% pull(alpha)
   deltaC = mark_data %>% pull(mean_slope)
   alf_Nstem_deltaC_Ncrypts = c(alf * deltaC * Nstem * Ncrypts, 0)
   
   input_list = list(N = nrow(joinscaledat), 
                     Gamma_dat = joinscaledat %>% pull(Gamma_dat),
                     Rin = joinscaledat %>% pull(Rin),
                     Rout = joinscaledat %>% pull(Rout),
                     theta = joinscaledat %>% pull(theta_subtended),
                     M = n_nbrhoods,
                     am = joinscaledat %>% select(a_m, nbrhood_index) %>% distinct() %>% pull(a_m),
                     awt = joinscaledat %>% select(a_wt, nbrhood_index) %>% distinct() %>% pull(a_wt),
                     nbrhood_index = joinscaledat %>% pull(nbrhood_index),
                     mark_index = mark_index,
                     pat_ages = age_M,
                     alf_Nstem_deltaC_Ncrypts = alf_Nstem_deltaC_Ncrypts,
                     n_grid = n_grid,
                     num_events = num_events)
   
   ## initialisation
   # simulate initial good paths for each mutant patch
   event_times = rep(0, n_nbrhoods)
   for (i in 1:n_nbrhoods) {
      hit = 999
      while (hit>age_M[i]) hit = next_fission_time(alf_Nstem_deltaC_Ncrypts[mark_index[i]], 1, runif(1))
      event_times[i] = hit
   }
   init_list = list(gamma_a = joinscaledat %>% select(gamma_a, nbrhood_index) %>% 
                       distinct() %>% pull(gamma_a),
                    D = rgamma(1, shape = 1, rate = 0.25),
                    sd_Gacc = fdrtool::rhalfnorm(1, fdrtool::sd2theta(1)),
                    a1_gpop = rgamma(1, 10, 1),
                    a2_gpop = rgamma(1, 10, 1),
                    event_times = event_times)
   stan_ins = list(input = input_list, init = init_list)
   return(list("stan_ins" = stan_ins, "joinscaledat" = joinscaledat))
}

generate_stan_ins_simdata3 = function(full_data, maxpatchsize, hitmark, first_hit_times, n_grid=50) {
   ## stan inputs
   joinscaledat = full_data
   joinscaledat = joinscaledat %>% mutate(nbrhood_index = neighbourhood)
   age_M = joinscaledat %>% select(nbrhood_index, age) %>% distinct() %>% arrange(nbrhood_index) %>% pull(age)
   n_nbrhoods = length(joinscaledat %>% pull(nbrhood_index) %>% unique())
   mark_index = rep(1, n_nbrhoods)
   
   # inputs for event time generation and numerical integration
   num_events = maxpatchsize
   input_list = list(N = nrow(joinscaledat), 
                     Gamma_dat = joinscaledat %>% pull(Gamma_dat),
                     Rin = joinscaledat %>% pull(Rin),
                     Rout = joinscaledat %>% pull(Rout),
                     theta = joinscaledat %>% pull(theta_subtended),
                     M = n_nbrhoods,
                     am = joinscaledat %>% select(a_m, nbrhood_index) %>% distinct() %>% pull(a_m),
                     awt = joinscaledat %>% select(a_wt, nbrhood_index) %>% distinct() %>% pull(a_wt),
                     nbrhood_index = joinscaledat %>% pull(nbrhood_index),
                     mark_index = mark_index,
                     pat_ages = age_M,
                     n_grid = n_grid,
                     num_events = num_events,
                     first_hit_times = first_hit_times)
   
   ## initialisation
   init_list = list(gamma_a = joinscaledat %>% select(gamma_a, nbrhood_index) %>% 
                       distinct() %>% pull(gamma_a),
                    D = rgamma(1, shape = 1, rate = 0.25),
                    sd_Gacc = fdrtool::rhalfnorm(1, fdrtool::sd2theta(1)),
                    a1_gpop = rgamma(1, 10, 1),
                    a2_gpop = rgamma(1, 10, 1))
   stan_ins = list(input = input_list, init = init_list)
   return(list("stan_ins" = stan_ins, "joinscaledat" = joinscaledat))
}

generate_stan_ins4 = function(data_primed_stag, data_primed_kdm, length_scale, area_scale, slide_age_ref, n_grid=50, Nstem=7, Ncrypts=1e7) {
   joinscaledat = data_primed_kdm$data_primed %>% mutate(mark = "KDM6A") %>% 
      bind_rows(data_primed_stag$data_primed %>% mutate(mark = "STAG2")) %>% 
      mutate(Rin = Rin / length_scale, Rout = Rout / length_scale,
             s_wt = s_wt / area_scale, a_wt = a_wt / area_scale, Am_incr = Am_incr / area_scale,
             a_m = a_m / area_scale, Total_area = Total_area / area_scale, 
             Gamma_dat = Gamma_dat / area_scale, epsilon = epsilon / area_scale, 
             Gamma_null = Gamma_null / area_scale)
   
   joinscaledat = joinscaledat %>% arrange(mark, Scan_ID, neighbourhood, Type) %>% 
      mutate(ismut = as.integer(Type=="mut")) %>% 
      mutate(bit.change=cumsum(c(1, diff(ismut) != 0))) %>% 
      mutate(bit.change = ifelse(ismut==1, bit.change, bit.change-1)) %>% 
      mutate(nbrhood_index = as.integer(as.factor(bit.change))) %>% 
      select(-ismut, -bit.change)
   
   # join ages
   joinscaledat = joinscaledat %>% left_join(slide_age_ref %>% select(Scan_ID, Age))
   age_M = joinscaledat %>% select(nbrhood_index, Age) %>% distinct() %>% arrange(nbrhood_index) %>% pull(Age)
   
   n_nbrhoods = length(joinscaledat %>% pull(nbrhood_index) %>% unique())
   slide_index_tib = joinscaledat %>% select(mark, Scan_ID, nbrhood_index, gamma_a) %>% distinct() %>% 
      mutate(slide_index = as.integer(as.factor(Scan_ID)),
             mark_index = ifelse(mark=="STAG2", 1, 2))
   wt30_gamma_a_tib = slide_index_tib %>% select(Scan_ID, gamma_a, slide_index) %>% distinct()
   n_slides = nrow(wt30_gamma_a_tib)
   mark_index = slide_index_tib %>% pull(mark_index)

   # inputs for event time generation and numerical integration
   num_events = joinscaledat %>% pull(n_crypts) %>% unique()
   mark_data = read_csv("../drift_All_Marks/exp_coeff.csv")
   mark_data = slide_index_tib %>% select(mark, mark_index) %>% distinct() %>% arrange(mark_index) %>% left_join(mark_data)
   fisrates = mark_data %>% pull(rho)
   alf = mark_data %>% pull(alpha)
   deltaC = mark_data %>% pull(mean_slope)
   alf_Nstem_deltaC_Ncrypts = alf * deltaC * Nstem * Ncrypts
   
   input_list = list(N = nrow(joinscaledat), 
                     Gamma_dat = joinscaledat %>% pull(Gamma_dat),
                     Rin = joinscaledat %>% pull(Rin),
                     Rout = joinscaledat %>% pull(Rout),
                     theta = joinscaledat %>% pull(theta_subtended),
                     M = n_nbrhoods,
                     am = joinscaledat %>% select(a_m, nbrhood_index) %>% distinct() %>% pull(a_m),
                     awt = joinscaledat %>% select(a_wt, nbrhood_index) %>% distinct() %>% pull(a_wt),
                     nbrhood_index = joinscaledat %>% pull(nbrhood_index),
                     mark_index = slide_index_tib %>% pull(mark_index),
                     pat_ages = age_M,
                     fisrate = fisrates,
                     alf_Nstem_deltaC_Ncrypts = alf_Nstem_deltaC_Ncrypts,
                     n_grid = n_grid,
                     num_events = num_events)
   
   ## initialisation
   # simulate initial good paths for each mutant patch
   event_times = matrix(NA, n_nbrhoods, num_events)
   for (i in 1:n_nbrhoods) {
      done = FALSE
      while (done==FALSE) {
         hit = 999
         while (hit>age_M[i]) hit = next_fission_time(alf_Nstem_deltaC_Ncrypts[mark_index[i]], 1, runif(1))
         posspath = generate_fission_event_times_given_hit(fisrates[mark_index[i]], age_M[i], hittime = hit)
         if (length(posspath)==num_events) {
            done = TRUE
            truepath = rep(0, num_events)
            truepath[1] = hit
            for (k in 2:num_events) truepath[k] = posspath[k] - posspath[k-1]
            event_times[i,] = truepath
         }
      }
   }
   init_list = list(gamma_a = joinscaledat %>% select(gamma_a, nbrhood_index) %>% 
                       distinct() %>% pull(gamma_a),
                    D = fdrtool::rhalfnorm(1, theta=fdrtool::sd2theta(2)),
                    sd_Gacc = fdrtool::rhalfnorm(1, fdrtool::sd2theta(1)),
                    a1_gpop = rgamma(1, 10, 1),
                    a2_gpop = rgamma(1, 10, 1),
                    event_times = event_times)
   stan_ins = list(input = input_list, init = init_list)
   return(list("stan_ins" = stan_ins, "joinscaledat" = joinscaledat, "slide_index_tib" = slide_index_tib))
}

generate_stan_ins5 = function(data_primed_stag, data_primed_kdm, length_scale, area_scale, slide_age_ref, n_grid=50, Nstem=7, Ncrypts=1e7) {
   joinscaledat = data_primed_kdm$data_primed %>% mutate(mark = "KDM6A") %>% 
      bind_rows(data_primed_stag$data_primed %>% mutate(mark = "STAG2")) %>% 
      mutate(Rin = Rin / length_scale, Rout = Rout / length_scale,
             s_wt = s_wt / area_scale, a_wt = a_wt / area_scale, Am_incr = Am_incr / area_scale,
             a_m = a_m / area_scale, Total_area = Total_area / area_scale, 
             Gamma_dat = Gamma_dat / area_scale, epsilon = epsilon / area_scale, 
             Gamma_null = Gamma_null / area_scale)
   
   joinscaledat = joinscaledat %>% arrange(mark, Scan_ID, neighbourhood, Type) %>% 
      mutate(ismut = as.integer(Type=="mut")) %>% 
      mutate(bit.change=cumsum(c(1, diff(ismut) != 0))) %>% 
      mutate(bit.change = ifelse(ismut==1, bit.change, bit.change-1)) %>% 
      mutate(nbrhood_index = as.integer(as.factor(bit.change))) %>% 
      select(-ismut, -bit.change)
   
   # join ages
   joinscaledat = joinscaledat %>% left_join(slide_age_ref %>% select(Scan_ID, Age))
   age_M = joinscaledat %>% select(nbrhood_index, Age) %>% distinct() %>% arrange(nbrhood_index) %>% pull(Age)
   
   n_nbrhoods = length(joinscaledat %>% pull(nbrhood_index) %>% unique())
   slide_index_tib = joinscaledat %>% select(mark, Scan_ID, nbrhood_index, gamma_a) %>% distinct() %>% 
      mutate(slide_index = as.integer(as.factor(Scan_ID)),
             mark_index = ifelse(mark=="STAG2", 1, 2))
   wt30_gamma_a_tib = slide_index_tib %>% select(Scan_ID, gamma_a, slide_index) %>% distinct()
   n_slides = nrow(wt30_gamma_a_tib)
   mark_index = slide_index_tib %>% pull(mark_index)
   
   # inputs for event time generation and numerical integration
   num_events = joinscaledat %>% pull(n_crypts) %>% unique()
   mark_data = read_csv("../drift_All_Marks/exp_coeff.csv")
   mark_data = slide_index_tib %>% select(mark, mark_index) %>% distinct() %>% arrange(mark_index) %>% left_join(mark_data)
   fisrates = mark_data %>% pull(rho)
   alf = mark_data %>% pull(alpha)
   deltaC = mark_data %>% pull(mean_slope)
   alf_Nstem_deltaC_Ncrypts = alf * deltaC * Nstem * Ncrypts
   
   input_list = list(N = nrow(joinscaledat), 
                     Gamma_dat = joinscaledat %>% pull(Gamma_dat),
                     Rin = joinscaledat %>% pull(Rin),
                     Rout = joinscaledat %>% pull(Rout),
                     theta = joinscaledat %>% pull(theta_subtended),
                     M = n_nbrhoods,
                     am = joinscaledat %>% select(a_m, nbrhood_index) %>% distinct() %>% pull(a_m),
                     awt = joinscaledat %>% select(a_wt, nbrhood_index) %>% distinct() %>% pull(a_wt),
                     nbrhood_index = joinscaledat %>% pull(nbrhood_index),
                     mark_index = slide_index_tib %>% pull(mark_index),
                     pat_ages = age_M,
                     alf_Nstem_deltaC_Ncrypts = alf_Nstem_deltaC_Ncrypts,
                     n_grid = n_grid,
                     num_events = num_events)
   
   ## initialisation
   # simulate initial good paths for each mutant patch
   event_times = rep(0, n_nbrhoods)
   for (i in 1:n_nbrhoods) {
      hit = 999
      while (hit>age_M[i]) hit = next_fission_time(alf_Nstem_deltaC_Ncrypts[mark_index[i]], 1, runif(1))
      event_times[i] = hit
   }
   init_list = list(gamma_a = joinscaledat %>% select(gamma_a, nbrhood_index) %>% 
                       distinct() %>% pull(gamma_a),
                    D = fdrtool::rhalfnorm(1, theta=fdrtool::sd2theta(2)),
                    sd_Gacc = fdrtool::rhalfnorm(1, fdrtool::sd2theta(1)),
                    a1_gpop = rgamma(1, 10, 1),
                    a2_gpop = rgamma(1, 10, 1),
                    event_times = event_times)
   stan_ins = list(input = input_list, init = init_list)
   return(list("stan_ins" = stan_ins, "joinscaledat" = joinscaledat, "slide_index_tib" = slide_index_tib))
}

generate_stan_ins6 = function(data_primed_stag, data_primed_kdm, length_scale, area_scale, slide_age_ref, n_grid=50, betatightness=50) {
   joinscaledat = data_primed_kdm$data_primed %>% mutate(mark = "KDM6A") %>% 
      bind_rows(data_primed_stag$data_primed %>% mutate(mark = "STAG2")) %>% 
      mutate(Rin = Rin / length_scale, Rout = Rout / length_scale,
             s_wt = s_wt / area_scale, a_wt = a_wt / area_scale, Am_incr = Am_incr / area_scale,
             a_m = a_m / area_scale, Total_area = Total_area / area_scale, 
             Gamma_dat = Gamma_dat / area_scale, epsilon = epsilon / area_scale, 
             Gamma_null = Gamma_null / area_scale)
   
   joinscaledat = joinscaledat %>% arrange(mark, Scan_ID, neighbourhood, Type) %>% 
      mutate(ismut = as.integer(Type=="mut")) %>% 
      mutate(bit.change=cumsum(c(1, diff(ismut) != 0))) %>% 
      mutate(bit.change = ifelse(ismut==1, bit.change, bit.change-1)) %>% 
      mutate(nbrhood_index = as.integer(as.factor(bit.change))) %>% 
      select(-ismut, -bit.change)
   
   # join ages
   joinscaledat = joinscaledat %>% left_join(slide_age_ref %>% select(Scan_ID, Age))
   age_M = joinscaledat %>% select(nbrhood_index, Age) %>% distinct() %>% arrange(nbrhood_index) %>% pull(Age)
   
   n_nbrhoods = length(joinscaledat %>% pull(nbrhood_index) %>% unique())
   slide_index_tib = joinscaledat %>% select(mark, Scan_ID, nbrhood_index, gamma_a) %>% distinct() %>% 
      mutate(slide_index = as.integer(as.factor(Scan_ID)),
             mark_index = ifelse(mark=="STAG2", 1, 2))
   wt30_gamma_a_tib = slide_index_tib %>% select(Scan_ID, gamma_a, slide_index) %>% distinct()
   n_slides = nrow(wt30_gamma_a_tib)
   mark_index = slide_index_tib %>% pull(mark_index)
   
   # inputs for event time generation and numerical integration
   num_events = (joinscaledat %>% pull(n_crypts) %>% unique()) - 1 # not counting first hit
   mark_data = read_csv("../drift_All_Marks/exp_coeff.csv")
   mark_data = slide_index_tib %>% select(mark, mark_index) %>% distinct() %>% arrange(mark_index) %>% left_join(mark_data)
   fisrates = mark_data %>% pull(rho)

   tp_raw_prior_tib = tibble()
   for (thismark in (slide_index_tib %>% select(mark, mark_index) %>% distinct() %>% arrange(mark_index) %>% pull(mark))) {
      fit_tib = readRDS(paste0("./output/patch_age_priorbetafits_", paste0(thismark, ".rds")))
      sum_fit_tib = fit_tib %>% summarise(alpha=median(alpha)*betatightness, beta=median(beta)*betatightness)
      tp_raw_prior_tib = tp_raw_prior_tib %>% bind_rows(sum_fit_tib)
   }
   
   input_list = list(N = nrow(joinscaledat), 
                     Gamma_dat = joinscaledat %>% pull(Gamma_dat),
                     Rin = joinscaledat %>% pull(Rin),
                     Rout = joinscaledat %>% pull(Rout),
                     theta = joinscaledat %>% pull(theta_subtended),
                     M = n_nbrhoods,
                     am = joinscaledat %>% select(a_m, nbrhood_index) %>% distinct() %>% pull(a_m),
                     awt = joinscaledat %>% select(a_wt, nbrhood_index) %>% distinct() %>% pull(a_wt),
                     nbrhood_index = joinscaledat %>% pull(nbrhood_index),
                     mark_index = slide_index_tib %>% pull(mark_index),
                     pat_ages = age_M,
                     fisrate = fisrates,
                     n_grid = n_grid,
                     num_events = num_events,
                     a_hit = tp_raw_prior_tib$alpha,
                     b_hit = tp_raw_prior_tib$beta)
   
   ## initialisation
   # simulate initial good paths for each mutant patch
   event_times = matrix(NA, n_nbrhoods, num_events)
   tp_init = rep(NA, n_nbrhoods)
   for (i in 1:n_nbrhoods) {
      done = FALSE
      while (done==FALSE) {
         tp_init[i] = rbeta(1, tp_raw_prior_tib$alpha[mark_index[i]], tp_raw_prior_tib$beta[mark_index[i]])
         hit = age_M[i] * (1 - tp_init[i])
         posspath = generate_fission_event_times_given_hit(fisrates[mark_index[i]], age_M[i], hittime = hit)
         if (length(posspath)==num_events) {
            done = TRUE
            truepath = rep(0, num_events)
            truepath[1] = hit
            for (k in 2:num_events) truepath[k] = posspath[k] - posspath[k-1]
            event_times[i,] = truepath
         }
      }
   }
   init_list = list(gamma_a = joinscaledat %>% select(gamma_a, nbrhood_index) %>% 
                       distinct() %>% pull(gamma_a),
                    D = fdrtool::rhalfnorm(1, theta=fdrtool::sd2theta(2)),
                    sd_Gacc = fdrtool::rhalfnorm(1, fdrtool::sd2theta(1)),
                    a1_gpop = rgamma(1, 10, 1),
                    a2_gpop = rgamma(1, 10, 1),
                    event_times = event_times,
                    t_p_raw = tp_init)
   stan_ins = list(input = input_list, init = init_list)
   return(list("stan_ins" = stan_ins, "joinscaledat" = joinscaledat, "slide_index_tib" = slide_index_tib))
}

generate_stan_ins7 = function(data_primed_stag, data_primed_kdm, length_scale, area_scale, slide_age_ref, n_grid=50) {
   joinscaledat = data_primed_kdm$data_primed %>% mutate(mark = "KDM6A") %>% 
      bind_rows(data_primed_stag$data_primed %>% mutate(mark = "STAG2")) %>% 
      mutate(Rin = Rin / length_scale, Rout = Rout / length_scale,
             s_wt = s_wt / area_scale, a_wt = a_wt / area_scale, Am_incr = Am_incr / area_scale,
             a_m = a_m / area_scale, Total_area = Total_area / area_scale, 
             Gamma_dat = Gamma_dat / area_scale, epsilon = epsilon / area_scale, 
             Gamma_null = Gamma_null / area_scale)
   
   joinscaledat = joinscaledat %>% arrange(mark, Scan_ID, neighbourhood, Type) %>% 
      mutate(ismut = as.integer(Type=="mut")) %>% 
      mutate(bit.change=cumsum(c(1, diff(ismut) != 0))) %>% 
      mutate(bit.change = ifelse(ismut==1, bit.change, bit.change-1)) %>% 
      mutate(nbrhood_index = as.integer(as.factor(bit.change))) %>% 
      select(-ismut, -bit.change)
   
   # join ages
   joinscaledat = joinscaledat %>% left_join(slide_age_ref %>% select(Scan_ID, Age))
   age_M = joinscaledat %>% select(nbrhood_index, Age) %>% distinct() %>% arrange(nbrhood_index) %>% pull(Age)
   
   n_nbrhoods = length(joinscaledat %>% pull(nbrhood_index) %>% unique())
   slide_index_tib = joinscaledat %>% select(mark, Scan_ID, nbrhood_index, gamma_a) %>% distinct() %>% 
      mutate(slide_index = as.integer(as.factor(Scan_ID)),
             mark_index = ifelse(mark=="STAG2", 1, 2))
   wt30_gamma_a_tib = slide_index_tib %>% select(Scan_ID, gamma_a, slide_index) %>% distinct()
   n_slides = nrow(wt30_gamma_a_tib)
   mark_index = slide_index_tib %>% pull(mark_index)
   
   # inputs for event time generation and numerical integration
   num_events = (joinscaledat %>% pull(n_crypts) %>% unique()) - 1 # not counting first hit
   mark_data = read_csv("../drift_All_Marks/exp_coeff.csv")
   mark_data = slide_index_tib %>% select(mark, mark_index) %>% distinct() %>% arrange(mark_index) %>% left_join(mark_data)
   fisrates = mark_data %>% pull(rho)
   
   tp_raw_prior_tib = tibble()
   for (thismark in (slide_index_tib %>% select(mark, mark_index) %>% distinct() %>% arrange(mark_index) %>% pull(mark))) {
      fit_tib = readRDS(paste0("./output/patch_age_priorbetafits_", paste0(thismark, ".rds")))
      sum_fit_tib = fit_tib %>% mutate(mean_pa =  1/(1 + beta/alpha)) %>% summarise(mean_tpa = mean(mean_pa))
      tp_raw_prior_tib = tp_raw_prior_tib %>% bind_rows(sum_fit_tib)
   }
   
   hit_times = rep(NA, n_nbrhoods)
   for (i in 1:n_nbrhoods) {
      hit_times[i] = age_M[i] * (1 - tp_raw_prior_tib$mean_tpa[mark_index[i]])
   }
   
   input_list = list(N = nrow(joinscaledat), 
                     Gamma_dat = joinscaledat %>% pull(Gamma_dat),
                     Rin = joinscaledat %>% pull(Rin),
                     Rout = joinscaledat %>% pull(Rout),
                     theta = joinscaledat %>% pull(theta_subtended),
                     M = n_nbrhoods,
                     am = joinscaledat %>% select(a_m, nbrhood_index) %>% distinct() %>% pull(a_m),
                     awt = joinscaledat %>% select(a_wt, nbrhood_index) %>% distinct() %>% pull(a_wt),
                     nbrhood_index = joinscaledat %>% pull(nbrhood_index),
                     mark_index = slide_index_tib %>% pull(mark_index),
                     pat_ages = age_M,
                     fisrate = fisrates,
                     n_grid = n_grid,
                     num_events = num_events,
                     first_hit_times = hit_times)
   
   ## initialisation
   # simulate initial good paths for each mutant patch
   init_list = list(gamma_a = joinscaledat %>% select(gamma_a, nbrhood_index) %>% 
                       distinct() %>% pull(gamma_a),
                    D = fdrtool::rhalfnorm(1, theta=fdrtool::sd2theta(2)),
                    sd_Gacc = fdrtool::rhalfnorm(1, fdrtool::sd2theta(1)),
                    a1_gpop = rgamma(1, 10, 1),
                    a2_gpop = rgamma(1, 10, 1))
   stan_ins = list(input = input_list, init = init_list)
   return(list("stan_ins" = stan_ins, "joinscaledat" = joinscaledat, "slide_index_tib" = slide_index_tib))
}

generate_stan_ins8 = function(data_primed_stag, data_primed_kdm, length_scale, area_scale, slide_age_ref, n_grid=50) {
   joinscaledat = data_primed_kdm$data_primed %>% mutate(mark = "KDM6A") %>% 
      bind_rows(data_primed_stag$data_primed %>% mutate(mark = "STAG2")) %>% 
      mutate(Rin = Rin / length_scale, Rout = Rout / length_scale,
             s_wt = s_wt / area_scale, a_wt = a_wt / area_scale, Am_incr = Am_incr / area_scale,
             a_m = a_m / area_scale, Total_area = Total_area / area_scale, 
             Gamma_dat = Gamma_dat / area_scale, epsilon = epsilon / area_scale, 
             Gamma_null = Gamma_null / area_scale)
   
   joinscaledat = joinscaledat %>% arrange(mark, Scan_ID, neighbourhood, Type) %>% 
      mutate(ismut = as.integer(Type=="mut")) %>% 
      mutate(bit.change=cumsum(c(1, diff(ismut) != 0))) %>% 
      mutate(bit.change = ifelse(ismut==1, bit.change, bit.change-1)) %>% 
      mutate(nbrhood_index = as.integer(as.factor(bit.change))) %>% 
      select(-ismut, -bit.change)
   
   # join ages
   joinscaledat = joinscaledat %>% left_join(slide_age_ref %>% select(Scan_ID, Age))
   age_M = joinscaledat %>% select(nbrhood_index, Age) %>% distinct() %>% arrange(nbrhood_index) %>% pull(Age)
   
   n_nbrhoods = length(joinscaledat %>% pull(nbrhood_index) %>% unique())
   slide_index_tib = joinscaledat %>% select(mark, Scan_ID, nbrhood_index, gamma_a) %>% distinct() %>% 
      mutate(slide_index = as.integer(as.factor(Scan_ID)),
             mark_index = ifelse(mark=="STAG2", 1, 2))
   wt30_gamma_a_tib = slide_index_tib %>% select(Scan_ID, gamma_a, slide_index) %>% distinct()
   n_slides = nrow(wt30_gamma_a_tib)
   mark_index = slide_index_tib %>% pull(mark_index)
   
   # inputs for event time generation and numerical integration
   num_events = (joinscaledat %>% pull(n_crypts) %>% unique())
   mark_data = read_csv("../drift_All_Marks/exp_coeff.csv")
   mark_data = slide_index_tib %>% select(mark, mark_index) %>% distinct() %>% arrange(mark_index) %>% left_join(mark_data)
   fisrates = mark_data %>% pull(rho)
   
   tp_raw_prior_tib = tibble()
   for (thismark in (slide_index_tib %>% select(mark, mark_index) %>% distinct() %>% arrange(mark_index) %>% pull(mark))) {
      fit_tib = readRDS(paste0("./output/patch_age_priorbetafits_", paste0(thismark, ".rds")))
      sum_fit_tib = fit_tib %>% mutate(mean_pa =  1/(1 + beta/alpha)) %>% summarise(mean_tpa = mean(mean_pa))
      tp_raw_prior_tib = tp_raw_prior_tib %>% bind_rows(sum_fit_tib)
   }
   
   hit_times = rep(NA, n_nbrhoods)
   for (i in 1:n_nbrhoods) {
      hit_times[i] = age_M[i] * (1 - tp_raw_prior_tib$mean_tpa[mark_index[i]])
   }
   
   input_list = list(N = nrow(joinscaledat), 
                     Gamma_dat = joinscaledat %>% pull(Gamma_dat),
                     Rin = joinscaledat %>% pull(Rin),
                     Rout = joinscaledat %>% pull(Rout),
                     theta = joinscaledat %>% pull(theta_subtended),
                     M = n_nbrhoods,
                     am = joinscaledat %>% select(a_m, nbrhood_index) %>% distinct() %>% pull(a_m),
                     awt = joinscaledat %>% select(a_wt, nbrhood_index) %>% distinct() %>% pull(a_wt),
                     nbrhood_index = joinscaledat %>% pull(nbrhood_index),
                     mark_index = slide_index_tib %>% pull(mark_index),
                     pat_ages = age_M,
                     fisrate = fisrates,
                     n_grid = n_grid,
                     num_events = num_events,
                     first_hit_times = hit_times)
   
   ## initialisation
   # simulate initial good paths for each mutant patch
   event_times = matrix(NA, n_nbrhoods, num_events)
   for (i in 1:n_nbrhoods) {
      done = FALSE
      while (done==FALSE) {
         hit = hit_times[i]
         posspath = generate_fission_event_times_given_hit(fisrates[mark_index[i]], age_M[i], hittime = hit)
         if (length(posspath)==num_events) {
            done = TRUE
            truepath = rep(0, num_events)
            truepath[1] = hit
            for (k in 2:num_events) truepath[k] = posspath[k] - posspath[k-1]
            event_times[i,] = truepath
         }
      }
   }
   init_list = list(gamma_a = joinscaledat %>% select(gamma_a, nbrhood_index) %>% 
                       distinct() %>% pull(gamma_a),
                    D = fdrtool::rhalfnorm(1, theta=fdrtool::sd2theta(2)),
                    sd_Gacc = fdrtool::rhalfnorm(1, fdrtool::sd2theta(1)),
                    a1_gpop = rgamma(1, 10, 1),
                    a2_gpop = rgamma(1, 10, 1),
                    event_times = event_times)
   stan_ins = list(input = input_list, init = init_list)
   return(list("stan_ins" = stan_ins, "joinscaledat" = joinscaledat, "slide_index_tib" = slide_index_tib))
}

generate_stan_ins9 = function(data_primed_stag, data_primed_kdm, length_scale, 
                              area_scale, slide_age_ref, event_times, n_grid=50) {
   joinscaledat = data_primed_kdm$data_primed %>% mutate(mark = "KDM6A") %>% 
      bind_rows(data_primed_stag$data_primed %>% mutate(mark = "STAG2")) %>% 
      mutate(Rin = Rin / length_scale, Rout = Rout / length_scale,
             s_wt = s_wt / area_scale, a_wt = a_wt / area_scale, Am_incr = Am_incr / area_scale,
             a_m = a_m / area_scale, Total_area = Total_area / area_scale, 
             Gamma_dat = Gamma_dat / area_scale, epsilon = epsilon / area_scale, 
             Gamma_null = Gamma_null / area_scale)
   
   joinscaledat = joinscaledat %>% arrange(mark, Scan_ID, neighbourhood, Type) %>% 
      mutate(ismut = as.integer(Type=="mut")) %>% 
      mutate(bit.change=cumsum(c(1, diff(ismut) != 0))) %>% 
      mutate(bit.change = ifelse(ismut==1, bit.change, bit.change-1)) %>% 
      mutate(nbrhood_index = as.integer(as.factor(bit.change))) %>% 
      select(-ismut, -bit.change)
   
   # join ages
   joinscaledat = joinscaledat %>% left_join(slide_age_ref %>% select(Scan_ID, Age))
   age_M = joinscaledat %>% select(nbrhood_index, Age) %>% distinct() %>% arrange(nbrhood_index) %>% pull(Age)
   
   n_nbrhoods = length(joinscaledat %>% pull(nbrhood_index) %>% unique())
   slide_index_tib = joinscaledat %>% select(mark, Scan_ID, nbrhood_index, gamma_a) %>% distinct() %>% 
      mutate(slide_index = as.integer(as.factor(Scan_ID)),
             mark_index = ifelse(mark=="STAG2", 1, 2))
   wt30_gamma_a_tib = slide_index_tib %>% select(Scan_ID, gamma_a, slide_index) %>% distinct()
   n_slides = nrow(wt30_gamma_a_tib)
   mark_index = slide_index_tib %>% pull(mark_index)
   
   # inputs for event time generation and numerical integration
   num_events = joinscaledat %>% pull(n_crypts) %>% unique()
   
   input_list = list(N = nrow(joinscaledat), 
                     Gamma_dat = joinscaledat %>% pull(Gamma_dat),
                     Rin = joinscaledat %>% pull(Rin),
                     Rout = joinscaledat %>% pull(Rout),
                     theta = joinscaledat %>% pull(theta_subtended),
                     M = n_nbrhoods,
                     am = joinscaledat %>% select(a_m, nbrhood_index) %>% distinct() %>% pull(a_m),
                     awt = joinscaledat %>% select(a_wt, nbrhood_index) %>% distinct() %>% pull(a_wt),
                     nbrhood_index = joinscaledat %>% pull(nbrhood_index),
                     mark_index = slide_index_tib %>% pull(mark_index),
                     pat_ages = age_M,
                     n_grid = n_grid,
                     num_events = num_events,
                     event_times = event_times)
   
   ## initialisation
   init_list = list(gamma_a = joinscaledat %>% select(gamma_a, nbrhood_index) %>% 
                       distinct() %>% pull(gamma_a),
                    D = rgamma(1, shape = 1, rate = 0.25),
                    sd_Gacc = fdrtool::rhalfnorm(1, fdrtool::sd2theta(1)),
                    a1_gpop = rgamma(1, 10, 1),
                    a2_gpop = rgamma(1, 10, 1))
   stan_ins = list(input = input_list, init = init_list)
   return(list("stan_ins" = stan_ins, "joinscaledat" = joinscaledat, "slide_index_tib" = slide_index_tib))
}

generate_stan_ins_null = function(data_primed_stag, data_primed_kdm, length_scale, 
                                  area_scale, slide_age_ref, n_grid=50) {
   joinscaledat = data_primed_kdm$data_primed %>% mutate(mark = "KDM6A") %>% 
      bind_rows(data_primed_stag$data_primed %>% mutate(mark = "STAG2")) %>% 
      mutate(Rin = Rin / length_scale, Rout = Rout / length_scale,
             s_wt = s_wt / area_scale, a_wt = a_wt / area_scale, Am_incr = Am_incr / area_scale,
             a_m = a_m / area_scale, Total_area = Total_area / area_scale, 
             Gamma_dat = Gamma_dat / area_scale, epsilon = epsilon / area_scale, 
             Gamma_null = Gamma_null / area_scale)
   
   joinscaledat = joinscaledat %>% arrange(mark, Scan_ID, neighbourhood, Type) %>% 
      mutate(ismut = as.integer(Type=="mut")) %>% 
      mutate(bit.change=cumsum(c(1, diff(ismut) != 0))) %>% 
      mutate(bit.change = ifelse(ismut==1, bit.change, bit.change-1)) %>% 
      mutate(nbrhood_index = as.integer(as.factor(bit.change))) %>% 
      select(-ismut, -bit.change)
   
   # join ages
   joinscaledat = joinscaledat %>% left_join(slide_age_ref %>% select(Scan_ID, Age))
   
   n_nbrhoods = length(joinscaledat %>% pull(nbrhood_index) %>% unique())
   slide_index_tib = joinscaledat %>% select(mark, Scan_ID, nbrhood_index, gamma_a) %>% distinct() %>% 
      mutate(slide_index = as.integer(as.factor(Scan_ID)),
             mark_index = ifelse(mark=="STAG2", 1, 2))
   wt30_gamma_a_tib = slide_index_tib %>% select(Scan_ID, gamma_a, slide_index) %>% distinct()
   n_slides = nrow(wt30_gamma_a_tib)
   mark_index = slide_index_tib %>% pull(mark_index)
   
   input_list = list(N = nrow(joinscaledat), 
                     Gamma_dat = joinscaledat %>% pull(Gamma_dat),
                     Rin = joinscaledat %>% pull(Rin),
                     Rout = joinscaledat %>% pull(Rout),
                     theta = joinscaledat %>% pull(theta_subtended),
                     M = n_nbrhoods,
                     nbrhood_index = joinscaledat %>% pull(nbrhood_index),
                     mark_index = slide_index_tib %>% pull(mark_index),
                     n_grid = n_grid)
   
   ## initialisation
   init_list = list(gamma_a = joinscaledat %>% select(gamma_a, nbrhood_index) %>% 
                       distinct() %>% pull(gamma_a),
                    sd_Gacc = fdrtool::rhalfnorm(1, fdrtool::sd2theta(1)),
                    a1_gpop = rgamma(1, 10, 1),
                    a2_gpop = rgamma(1, 10, 1))
   stan_ins = list(input = input_list, init = init_list)
   return(list("stan_ins" = stan_ins, "joinscaledat" = joinscaledat, "slide_index_tib" = slide_index_tib))
}

prepare_r_dep_data_for_plotting = function(data_primed_kdm, data_primed_stag, stan_ins) {
   density_by_patch = data_primed_kdm$space_data %>% mutate(mark="KDM6A") %>% 
      bind_rows(data_primed_stag$space_data %>% mutate(mark="STAG2")) %>% 
      mutate(Rho = (average_crypt_area * n_crypts) / Total_area,
             gamma = (Total_area - average_crypt_area * n_crypts) / Total_area) %>% 
      select(mark, Scan_ID, neighbourhood, Type, Patch, Rho, gamma)
   just_mut_density = density_by_patch %>% filter(Type=="mut") %>% 
      mutate(Rho_mut = Rho, gamma_mut = gamma) %>% 
      select(mark, Scan_ID, neighbourhood, Rho_mut, gamma_mut)
   density_by_patch = density_by_patch %>% left_join(just_mut_density) %>% 
      mutate(Rho_sc = Rho / Rho_mut, gamma_sc = gamma / gamma_mut)
   location_by_patch = stan_ins$joinscaledat %>% left_join(density_by_patch) %>% 
      mutate(R_mean = ifelse(Type!="mut", 0.5*(Rin+Rout), 0))
   return(location_by_patch)
}

compare_r_dependence_to_data = function(data_primed_kdm, data_primed_stag, stan_ins, 
                                        r_dep_traj, WT30_r=10, useinds=NULL) {
   location_by_patch = prepare_r_dep_data_for_plotting(data_primed_kdm, data_primed_stag, stan_ins)
   # add in predicted trajectories and attach the WT 30 patches at some large r value
   slide_to_nbri_map = location_by_patch %>% select(mark, Scan_ID, nbrhood_index) %>% distinct()
   rinf_patches = data_primed_kdm$patch_whitespaces %>% mutate(mark = "KDM6A") %>% 
      bind_rows(data_primed_stag$patch_whitespaces %>% mutate(mark = "STAG2")) %>% 
      mutate(gamma = Gamma / area) %>% select(Scan_ID, gamma) %>% mutate(R_mean = WT30_r) %>% 
      right_join(slide_to_nbri_map)
   nbrs_i = location_by_patch %>% pull(nbrhood_index) %>% unique() %>% sort()
   pl_list = list()
   if (is.null(useinds)) {
      for (i in 1:ceiling(length(nbrs_i)/6)) {
         nbr_i = nbrs_i[((i-1)*6+1):min(i*6, length(nbrs_i))]
         pl = location_by_patch %>% filter(nbrhood_index%in%nbr_i) %>%  ggplot() + 
            geom_point(aes(x=R_mean, y=gamma, col=Type)) + 
            geom_errorbarh(aes(xmin=Rin, xmax=Rout, y=gamma, col=Type), height=0.005, alpha=0.6) + 
            geom_line(data=r_dep_traj %>% filter(nbrhood_index%in%nbr_i), aes(x=r, y=med_mcmc)) + 
            geom_ribbon(data=r_dep_traj %>% filter(nbrhood_index%in%nbr_i), aes(x=r, ymin=l95_mcmc, ymax=u95_mcmc), alpha=0.2) + 
            facet_wrap(~nbrhood_index, scales="free_y")
         if (!is.na(WT30_r)) pl = pl + geom_point(dat=rinf_patches %>% filter(nbrhood_index%in%nbr_i), aes(x=R_mean, y=gamma, col="WT30"))
         pl_list[[i]] = pl
      }   
   } else {
      pl = location_by_patch %>% filter(nbrhood_index%in%useinds) %>%  ggplot() + 
         geom_point(aes(x=R_mean, y=gamma, col=Type)) + 
         geom_errorbarh(aes(xmin=Rin, xmax=Rout, y=gamma, col=Type), height=0.005, alpha=0.6) +
         geom_line(data=r_dep_traj %>% filter(nbrhood_index%in%useinds), aes(x=r, y=med_mcmc)) + 
         geom_ribbon(data=r_dep_traj %>% filter(nbrhood_index%in%useinds), aes(x=r, ymin=l95_mcmc, ymax=u95_mcmc), alpha=0.2) + 
         facet_wrap(~nbrhood_index, scales="free_y")
      if (!is.na(WT30_r)) pl = pl + geom_point(dat=rinf_patches %>% filter(nbrhood_index%in%useinds), aes(x=R_mean, y=gamma, col="WT30"))
      pl_list[[1]] = pl
   }
   return(pl_list)
}

gather_patch_age_fits = function(fits, stan_ins) {
   stag_inds = which((stan_ins$slide_index_tib %>% pull(mark_index))==1)
   kdm_inds = which((stan_ins$slide_index_tib %>% pull(mark_index))==2)
   
   stag_tp_tib = fits$t_p[,stag_inds] %>% as_tibble() %>% mutate(iter = 1:n()) %>% 
      gather(patch, age, -iter) %>% 
      mutate(patch = as.integer(str_remove(patch, "V"))) %>% mutate(type="inferred")
   kdm_tp_tib = fits$t_p[,kdm_inds] %>% as_tibble() %>% mutate(iter = 1:n()) %>% 
      gather(patch, age, -iter) %>% 
      mutate(patch = as.integer(str_remove(patch, "V"))) %>% mutate(type="inferred")
   return(list("stag_tib" = stag_tp_tib, "kdm_tib" = kdm_tp_tib))
}

compare_patchage_inference_to_prior = function(fits, stan_ins, timeatN_stag, timeatN_kdm) {
   stag_bd = tibble(bd_ages = timeatN_stag) %>% sample_n(4000, replace = FALSE) %>% 
      mutate(type="fission/fusion\nbirth-death")
   kdm_bd = tibble(bd_ages = timeatN_kdm) %>% sample_n(4000, replace = FALSE) %>% 
      mutate(type="fission/fusion\nbirth-death")

   mark_tibs = gather_patch_age_fits(fits, stan_ins)
   
   pl_s = ggplot() + 
      geom_histogram(data=mark_tibs$stag_tib, aes(x=age, y=..density.., fill=type), alpha=0.5) + 
      geom_histogram(data=stag_bd, aes(x=bd_ages, y=..density.., fill=type), alpha=0.5) +
      facet_wrap(~patch)
   pl_k = ggplot() + 
      geom_histogram(data=mark_tibs$kdm_tib, aes(x=age, y=..density.., fill=type), alpha=0.5) + 
      geom_histogram(data=kdm_bd, aes(x=bd_ages, y=..density.., fill=type), alpha=0.5) +
      facet_wrap(~patch)
   return(list("plot_stag" = pl_s, "plot_kdm" = pl_k))
}

generate_average_gamma_r_evolution = function(fit_summaries, stan_ins, RUN=FALSE, rs = c(0,15), ns = 200) {
   N = stan_ins$joinscaledat %>% pull(n_crypts) %>% unique()
   if (RUN==TRUE) {
      calc_tib = stan_ins$joinscaledat %>% filter(Type=="mut") %>%
         select(mark, Scan_ID, neighbourhood, nbrhood_index, n_crypts, epsilon) %>% 
         left_join(fit_summaries$gamma_a %>% rename(nbrhood_index = param_num, 
                                                    gamma_a = med_mcmc) %>% 
                      select(nbrhood_index, gamma_a)) %>% 
         left_join(fit_summaries$Rdiff %>% rename(nbrhood_index = param_num, 
                                                  Rdiff = med_mcmc) %>% 
                      select(nbrhood_index, Rdiff))
      outtib = tibble()
      rspace = pracma::linspace(rs[1], rs[2], ns)
      for (rr in rspace) {
         this_res = calc_tib %>% group_by(mark, Scan_ID, nbrhood_index) %>%
            mutate(gamma_r = gamma_r(gamma_a, rr, Rdiff^2, epsilon, order="full"),
                   r = rr) %>% ungroup() %>% select(nbrhood_index, r, gamma_r)
         outtib = outtib %>% bind_rows(this_res)
      }
      r_dep_traj = calc_tib %>% select(mark, Scan_ID, neighbourhood, nbrhood_index, n_crypts) %>%
         left_join(outtib)
      r_dep_traj %>% saveRDS(paste0("./output/diffusion_average_gammar_", paste0(as.character(N), "patches.rds")))
   } else {
      r_dep_traj = readRDS(paste0("./output/diffusion_average_gammar_", paste0(as.character(N), "patches.rds")))   
   }
   return(r_dep_traj)
}

generate_gamma_r_evolution = function(fits, stan_ins, RUN=FALSE, rs = c(0,10), ns = 20) {
   N = stan_ins$joinscaledat %>% pull(n_crypts) %>% unique()
   if (RUN==TRUE) {
      calc_tib = stan_ins$joinscaledat %>% filter(Type=="mut") %>%
         select(mark, Scan_ID, neighbourhood, nbrhood_index, n_crypts, epsilon)
      outtib = tibble()
      rspace = pracma::linspace(rs[1], rs[2], ns)
      iters = 500
      iter_inds = sort(sample(1:dim(fits$gamma_a)[1], size = iters, replace = FALSE))
      for (i in 1:iters) {
         k = iter_inds[i]
         for (rr in rspace) {
            this_res = calc_tib %>% group_by(mark, Scan_ID, nbrhood_index) %>%
               mutate(gamma_r = gamma_r(fits$gamma_a[k,nbrhood_index], rr,
                                        fits$Rdiff[k,nbrhood_index]^2, epsilon, order="full"),
                      r = rr, iter = i) %>% ungroup() %>% select(nbrhood_index, r, gamma_r, iter)
            outtib = outtib %>% bind_rows(this_res)
         }
      }
      
      sumout = outtib %>% group_by(nbrhood_index, r) %>%
         summarise(med_mcmc = median(gamma_r),
                   l95_mcmc = quantile(gamma_r, 0.025),
                   u95_mcmc = quantile(gamma_r, 0.975)) %>%
         ungroup()
      r_dep_traj = calc_tib %>% select(mark, Scan_ID, neighbourhood, nbrhood_index, n_crypts) %>%
         left_join(sumout)
      r_dep_traj %>% saveRDS(paste0("./output/diffusion_gammar_", paste0(as.character(N), "patches.rds")))
   } else {
      r_dep_traj = readRDS(paste0("./output/diffusion_gammar_", paste0(as.character(N), "patches.rds")))   
   }
   return(r_dep_traj)
}

generate_gradgamma_r_evolution = function(fits, stan_ins, RUN=FALSE, rs = c(0,10), ns = 20) {
   N = stan_ins$joinscaledat %>% pull(n_crypts) %>% unique()
   calc_tib = stan_ins$joinscaledat %>% filter(Type=="mut") %>%
      select(mark, Scan_ID, neighbourhood, nbrhood_index, n_crypts, epsilon)
   if (RUN==TRUE) {
      outtib = tibble()
      rspace = pracma::linspace(rs[1], rs[2], ns)
      iters = 500
      iter_inds = sort(sample(1:dim(fits$gamma_a)[1], size = iters, replace = FALSE))
      for (i in 1:iters) {
         k = iter_inds[i]
         for (rr in rspace) {
            this_res = calc_tib %>% group_by(mark, Scan_ID, nbrhood_index) %>%
               mutate(gradgamma_r = grad_gammar(fits$gamma_a[k,nbrhood_index], rr,
                                        fits$Rdiff[k,nbrhood_index]^2, epsilon),
                      r = rr, iter = i) %>% ungroup() %>% 
               select(mark, Scan_ID, nbrhood_index, r, gradgamma_r, iter)
            outtib = outtib %>% bind_rows(this_res)
         }
      }
      outtib %>% saveRDS(paste0("./output/diffusion_gradgammar_", paste0(as.character(N), "patches.rds")))
   } else {
      outtib = readRDS(paste0("./output/diffusion_gradgammar_", paste0(as.character(N), "patches.rds")))   
   }
   sumout = outtib %>% group_by(nbrhood_index, r) %>%
      summarise(med_mcmc = median(gradgamma_r),
                l95_mcmc = quantile(gradgamma_r, 0.025),
                u95_mcmc = quantile(gradgamma_r, 0.975)) %>%
      ungroup()
   r_dep_traj = calc_tib %>% select(mark, Scan_ID, neighbourhood, nbrhood_index, n_crypts) %>%
      left_join(sumout)
   return(list("raw_out" = outtib, "r_dep_traj" = r_dep_traj))
}

shuffle_joint_columns = function(datin, cols_to_shuff) {
   shuffinds = sample(1:nrow(datin), replace=FALSE)
   datin %>% mutate_at(cols_to_shuff, ~.[shuffinds])
}

shuffle_multiple_groups = function(datin, groupname, cols_to_shuff) {
   the_groups = datin %>% pull(!!as.symbol(groupname)) %>% unique()
   outtib = tibble()
   for (i in 1:length(the_groups)) {
      outtib = outtib %>% bind_rows(datin %>% filter(!!as.symbol(groupname)==the_groups[i]) %>% 
                                       shuffle_joint_columns(cols_to_shuff))
   }
   return(outtib)
}

summarise_mcmc = function(mcmc_data, map_ind=NULL) {
   dims = dim(mcmc_data)
   if (length(dims)==1) thistib = mcmc_data %>% as_tibble() %>% magrittr::set_colnames(c("V1"))
   if (length(dims)==2) thistib = mcmc_data %>% as_tibble() %>% magrittr::set_colnames(paste0(rep("V", dims[2]), as.character(1:dims[2])))
   if (length(dims)>2) thistib = mcmc_data %>% as_tibble()
   if (!is.null(map_ind)) {
      return(thistib %>% mutate(iter = 1:n()) %>% 
                gather(param_num, val, -iter) %>% 
                group_by(param_num) %>% arrange(iter) %>% 
                summarise(mean_mcmc = mean(val, na.rm=T), med_mcmc = median(val, na.rm=T), 
                          u95_mcmc = quantile(val, 0.975, na.rm=T), l95_mcmc = quantile(val, 0.025, na.rm=T),
                          map_est = val[map_ind]) %>% 
                mutate(param_num=ifelse(str_detect(param_num, "V"), as.numeric(str_remove(param_num, "V")), param_num)) %>% 
                arrange(param_num))
   } else {
      return(thistib %>% mutate(iter = 1:n()) %>% 
                gather(param_num, val, -iter) %>% 
                group_by(param_num) %>% arrange(iter) %>% 
                summarise(mean_mcmc = mean(val, na.rm=T), med_mcmc = median(val, na.rm=T), 
                          u95_mcmc = quantile(val, 0.975, na.rm=T), l95_mcmc = quantile(val, 0.025, na.rm=T)) %>% 
                mutate(param_num=ifelse(str_detect(param_num, "V"), as.numeric(str_remove(param_num, "V")), param_num)) %>% 
                arrange(param_num))
   }
}

extract_param_summaries = function(fit, paramnames=NULL) {
   map_est = which.max(fit$lp__)
   if (is.null(paramnames)) {
      paramnames = names(fit)
   }
   paramnames = setdiff(paramnames, c("lp__"))
   outlist = list()
   for (thisname in paramnames) {
      outlist[[thisname]] = fit[[thisname]] %>% summarise_mcmc(map_ind = map_est) %>% mutate(param = thisname)
   }
   return(outlist)
}

plot_beta_prior = function(shapes, lwrbnd, uprbnd, npts) {
   alpha = shapes[1]
   beta = shapes[2]
   grid = pracma::linspace(lwrbnd,uprbnd,npts)
   samps = dbeta(x = grid, shape1 = alpha, shape2 = beta)
   plot(grid, samps)
}

plot_gamma_prior = function(shape, rate, lwrbnd, uprbnd, npts) {
   grid = pracma::linspace(lwrbnd,uprbnd,npts)
   samps = dgamma(x = grid, shape = shape, rate = rate)
   plot(grid, samps)
}

theor_range_plot = function(fits, av_patch_inputs, ages_tib, nr=50, rs=c(0,10), iters=300, get_errbars=FALSE) {
   outtib = tibble()
   rspace = pracma::linspace(rs[1], rs[2], nr)
   if (get_errbars==TRUE) {
      iter_inds = sort(sample(1:dim(fits$gamma_a)[1], size = iters, replace = FALSE))
      for (i in 1:iters) {
         k = iter_inds[i]
         for (rr in rspace) {
            this_res = av_patch_inputs %>% left_join(ages_tib) %>% group_by(mark) %>%
               mutate(gamma_r = gamma_r(fits$gamma_a_pop[k], rr,
                                        4*fits$D[k]*age, med_eps, order="full"),
                      r = rr, iter = i, age = age) %>% ungroup() %>% 
               select(mark, r, gamma_r, iter, age)
            outtib = outtib %>% bind_rows(this_res)
         }
      }
      sumout = outtib %>% group_by(mark, r) %>%
         summarise(med_mcmc = median(gamma_r),
                   l95_mcmc = quantile(gamma_r, 0.025),
                   u95_mcmc = quantile(gamma_r, 0.975)) %>%
         ungroup()
      return(sumout)
   } else {
      for (rr in rspace) {
         this_res = av_patch_inputs %>% left_join(ages_tib) %>% 
            group_by(mark) %>% mutate(Dtilde = 4*median(fits$D)*age) %>% 
            mutate(gamma_r = gamma_r(median(fits$gamma_a_pop), rr,
                                     Dtilde, med_eps, order="full"), 
                   r = rr, age = age) %>% 
            ungroup() %>% select(mark, r, gamma_r, age)
         outtib = outtib %>% bind_rows(this_res)
      }
      return(outtib)
   }
}

single_r_trajectory = function(gamma_a, D, patch_input, time, nr=500, rs=c(0,20)) {
   outtib = tibble()
   rspace = pracma::linspace(rs[1], rs[2], nr)
   Rdiffsq = 4*D*time
   for (rr in rspace) {
      this_res = tibble(t = time, 
                        r = rr,
                        gamma_r = gamma_r(gamma_a, rr, Rdiffsq, patch_input, order="full")) 
      outtib = outtib %>% bind_rows(this_res)
   }
   return(outtib)
}

get_kras_rdep = function(gamma_a, D, epsilon, age) {
   rspace = pracma::linspace(0,10,200)
   gamma_r_kras = c()
   for (rr in rspace) {
      Dtilde = 4*D*age
      gamma_r_kras = c(gamma_r_kras, gamma_r(gamma_a, rr, Dtilde, epsilon, order="full"))
   }
   return(tibble(r = rspace, gamma_r = gamma_r_kras))
}

reticulate::conda_list()
reticulate::use_condaenv('spring', required=TRUE)
reticulate::source_python("./integrate_polyinterp.py")

next_fission_time = function(rate, numcrypts, U01) {
   return(-log(1 - U01) / (rate*numcrypts))  
}

create_diffusion_perturbation = function(D, t, a_m, a_wt, rvec, cumulative_fission_times) {
   rho = (a_m - a_wt) / (4*pi*D*t) * exp(-rvec^2 / (4*D*t)) # initial hit solution
   ps = 1
   for (ft in cumulative_fission_times) {
      if (t<ft) break
      tau_i = t - ft
      rho = rho + a_m / (4*pi*D*tau_i) * exp(-rvec^2 / (4*D*tau_i)) # subsequent fissions
      ps = ps + 1
   }
   return(list("rho" = rho, "patchsize" = ps))
}

generate_fission_event_times_given_hit = function(fission_rate, maxt, hittime=0) {
   fis_time = hittime
   c_fis_time = c(fis_time)
   n = 0
   while (fis_time<maxt) {
      n = n + 1
      fis_time = tail(c_fis_time, 1) + next_fission_time(fission_rate, n, runif(1))
      if (fis_time<maxt) c_fis_time = c(c_fis_time, fis_time)
   }
   return(c_fis_time)
}

generate_fission_event_times = function(fission_rate, maxt, suppressed_rate=NULL, suppression_time=NULL) {
   fis_time = 0
   c_fis_time = c(fis_time)
   n = 0
   if (is.null(suppressed_rate)) {
      while (fis_time<maxt) {
         n = n + 1
         fis_time = tail(c_fis_time, 1) + next_fission_time(fission_rate, n, runif(1))
         if (fis_time<maxt) c_fis_time = c(c_fis_time, fis_time)
      }
      c_fis_time = c_fis_time[-1] # remove leading zero
      return(c_fis_time)
   } else {
      while (TRUE) {
         n = n + 1
         new_event_time = tail(c_fis_time, 1) + next_fission_time(fission_rate, n, runif(1))
         if (new_event_time>suppression_time[1]) {
            n = n - 1
            break
         } else {
            c_fis_time = c(c_fis_time, new_event_time)
         }
      }
      c_fis_time = c_fis_time[-1]
      c_fis_time_supp = c(suppression_time[1])
      while (TRUE) {
         n = n + 1
         new_event_time = tail(c_fis_time_supp, 1) + next_fission_time(suppressed_rate, n, runif(1))
         if (new_event_time>suppression_time[2]) {
            n = n - 1
            break
         } else {
            c_fis_time_supp = c(c_fis_time_supp, new_event_time)
         }
      }
      c_fis_time_supp = c_fis_time_supp[-1]
      c_fis_time2 = c(suppression_time[2])
      while (TRUE) {
         n = n + 1
         new_event_time = tail(c_fis_time2, 1) + next_fission_time(fission_rate, n, runif(1))
         if (new_event_time>maxt) {
            n = n - 1
            break
         } else {
            c_fis_time2 = c(c_fis_time2, new_event_time)
         }
      }
      c_fis_time2 = c_fis_time2[-1]
      return(c(c_fis_time, c_fis_time_supp, c_fis_time2))
   }
}

generate_diffusion_time_series = function(tvec, D, a_m, a_wt, rvec, fission_rate,
                                          suppressed_rate=NULL, suppression_time=NULL) {
   c_fis_time = generate_fission_event_times(fission_rate, tail(tvec,1), suppressed_rate=suppressed_rate, suppression_time=suppression_time)
   timeseries_tib = diffusion_time_series_given_event_times(c_fis_time, tvec, D, a_m, a_wt, rvec)
   return(timeseries_tib)
}

diffusion_time_series_given_event_times = function(c_fis_time, tvec, D, a_m, a_wt, rvec) {
   timeseries_tib = tibble()
   for (tt in tvec) {
      rho_t = create_diffusion_perturbation(D, tt, a_m, a_wt, rvec, c_fis_time)
      timeseries_tib = timeseries_tib %>% 
         bind_rows(tibble(t = tt, rho_pert = rho_t$rho, r = rvec, patchsize = rho_t$patchsize))
   }
   return(timeseries_tib)
}

iterate_diffusion_time_series = function(tvec, D, a_m, a_wt, rvec, fission_rate, iters=100,
                                         suppressed_rate=NULL, suppression_time=NULL) {
   library(doParallel)
   registerDoParallel()
   res = foreach(i=1:iters) %dopar% {
      generate_diffusion_time_series(tvec, D, a_m, a_wt, rvec, fission_rate, 
                                     suppressed_rate=suppressed_rate, suppression_time=suppression_time) %>% 
         mutate(iter = i)
   }
   return(do.call("rbind", res))
}

integrate_Gamma_parallel = function(res) {
   library(doParallel)
   registerDoParallel(cores=8)
   nruns = res %>% pull(iter) %>% unique()
   res_list = foreach(i=1:length(nruns)) %dopar% {
      res %>% filter(iter==nruns[i]) %>% group_by(t) %>% mutate(Gamma_integral = 2*pi*integrate(gamma_r*r, r))
   }
   return(do.call("rbind", res_list))
}

integrate_Rho_parallel = function(res, method="simple") {
   library(doParallel)
   registerDoParallel(cores=8)
   nruns = res %>% pull(iter) %>% unique()
   if (method=="simple") {
      res_list = foreach(i=1:length(nruns)) %dopar% {
         res %>% filter(iter==nruns[i]) %>% group_by(t) %>% mutate(Rho_integral = 2*pi*pracma::cumtrapz(r, (1-gamma_r)*r))
      }
   } else {
      res_list = foreach(i=1:length(nruns)) %dopar% {
         res %>% filter(iter==nruns[i]) %>% group_by(t) %>% mutate(Rho_integral = 2*pi*integrate((1-gamma_r)*r, r))
      }
   }
   return(do.call("rbind", res_list))
}

transform_to_realspace = function(timepoint, rho_thresh) {
   # can we use this to find "range" by setting max density to "just above" ambient?   
   if (nrow(timepoint %>% filter(rho_r>rho_thresh))==0) {
      return(timepoint)
   } else {
      old_r = timepoint %>% filter(rho_r>rho_thresh) %>% pull(r) %>% max()
      high_bit = timepoint %>% filter(rho_r>rho_thresh)
      if (nrow(high_bit)<2) {
         ## approximate integral around point value
         new_r = timepoint %>% arrange(r) %>% slice(1:2) %>%
            mutate(Rho_integral = 2*pi*(r[2]-r[1])*mean(r)*mean(rho_r) ) %>% 
            summarise(r1 = sqrt(max(Rho_integral)/(pi*rho_thresh))) %>% pull()
      } else if (nrow(high_bit)<10) {
         new_r = high_bit %>% mutate(Rho_integral = 2*pi*(pracma::cumtrapz(r, r*rho_r) %>% as.vector())) %>% 
            summarise(r1 = sqrt(max(Rho_integral)/(pi*rho_thresh))) %>% pull()
      } else {
         new_r = timepoint %>% filter(rho_r>rho_thresh) %>% mutate(Rho_integral = 2*pi*integrate(r*rho_r, r)) %>% 
            summarise(r1 = sqrt(max(Rho_integral)/(pi*rho_thresh))) %>% pull()   
      }
      shifted = timepoint %>% filter(r>old_r) %>% mutate(r = r + (new_r - old_r))
      flattened = timepoint %>% filter(r<=new_r) %>% mutate(rho_r = rho_thresh)
      new_timepoint = flattened %>% bind_rows(shifted)
      return(new_timepoint)
   }
   timepoint %>% ggplot() + geom_line(aes(x=r, y=rho_r, col = "raw")) +
      geom_line(data = new_timepoint, aes(x=r, y=rho_r, col = "transformed"))
}

calc_happy_packing_range = function(lambda, t, a_m, a_wt, gamma_a, pchresh) {
   # lambda == fission rate
   # thresh == percentage above ambient crypt density at "new happy" level (e.g. 1%, 5%)
   mass_in = (exp(lambda*t) - 1)*a_m + (a_m - a_wt)
   R_H = sqrt(mass_in / (pi/gamma_a * pchresh/100))
   return(R_H)
}

calc_happy_packing_footprint = function(R_H, s_wt) {
   ## calculate the number of crypts normally found
   ## in the happy packing range
   return((R_H^2 * pi) / s_wt)
}
