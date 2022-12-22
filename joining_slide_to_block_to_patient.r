library(readxl)
library(tidyverse)

## define a vector of bad patients, e.g. radiation therapy
# bad_pats = c("CRA_00370")
bad_pats = c("CRA_00370", "CRA_00299", "CRA_00349")

trim_block_ID = function(tib) {
   tib %>% mutate(staticID1 = Block_ID) %>% 
      mutate(Block_ID = str_replace(Block_ID, "BL", "")) %>% 
      mutate(Block_ID = str_replace(Block_ID, "NAUT", "")) %>% 
      mutate(Block_ID = str_replace(Block_ID, "PS[:digit:][:digit:][:punct:]", "")) %>% 
      mutate(Block_ID = str_replace(Block_ID, "PS[:digit:][:punct:]", "")) %>% 
      mutate(Block_ID = str_replace(Block_ID, "HB", "")) %>% 
      separate(Block_ID, c("Block_ID", "section1", "section2"), "[:punct:](?=[:digit:])") %>% 
      mutate(Block_ID = ifelse(is.na(section2), Block_ID, paste(Block_ID, section1, sep = "_"))) %>% 
      mutate(section = as.numeric(ifelse(is.na(section2), section1, section2))) %>% 
      select(-staticID1, -section1, -section2)
}

## global databases
block_to_patient = read_xls("./databases/Block_database.xls") %>% select(`Block ID`, Age, SEX, `DW Patient number`)
colnames(block_to_patient) = c("Block_ID", "Age", "Sex", "DW.Patient.number")
block_to_patient = trim_block_ID(block_to_patient) %>% select(-section)
block_to_patient = block_to_patient %>% filter(!is.na(Block_ID), !is.na(DW.Patient.number)) %>% 
   distinct()

block_to_slide_kdm = read_xlsx("./databases/Image_locations_all_marks.xlsx", sheet = "KDM6A") %>% 
   select(`Block ID`, `Image ID`)
block_to_slide_mpas = read_xlsx("./databases/Image_locations_all_marks.xlsx", sheet = "mPAS") %>% 
   select(`Block ID`, `Image ID`)
block_to_slide_stag = read_xlsx("./databases/Image_locations_all_marks.xlsx", sheet = "STAG2") %>% 
   select(`Block ID`, `Image ID`)
colnames(block_to_slide_kdm) = c("Block_ID", "Slide_ID")
colnames(block_to_slide_mpas) = c("Block_ID", "Slide_ID")
colnames(block_to_slide_stag) = c("Block_ID", "Slide_ID")
block_to_slide_kdm = trim_block_ID(block_to_slide_kdm %>% distinct())
block_to_slide_mpas = trim_block_ID(block_to_slide_mpas %>% distinct())
block_to_slide_stag = trim_block_ID(block_to_slide_stag %>% distinct())
block_to_slide_kdm = block_to_slide_kdm %>% filter(!is.na(Block_ID)) %>% distinct()
block_to_slide_mpas = block_to_slide_mpas %>% filter(!is.na(Block_ID)) %>% distinct()
block_to_slide_stag = block_to_slide_stag %>% filter(!is.na(Block_ID)) %>% distinct()

## KDM6A data
###################
dat_raw_kdm6a = read_xlsx("./databases/Fufi_counts_KDM6A.xlsx", sheet = 2)
colnames(dat_raw_kdm6a) = c("Slide_ID", "monoclonals", "knn1", "partials", "knn3",
                             "monoclonals_patch", "knn2", "partials_patch", "knn4")
dat_raw_kdm6a = dat_raw_kdm6a %>% mutate(monoclonals = as.integer(monoclonals), 
                                           partials = as.integer(partials),
                                           monoclonals_patch = as.integer(monoclonals_patch), 
                                           partials_patch = as.integer(partials_patch)) %>% 
   mutate(monoclonals = ifelse(is.na(monoclonals), 0, monoclonals), 
          partials = ifelse(is.na(partials), 0, partials),
          monoclonals_patch = ifelse(is.na(monoclonals_patch), 0, monoclonals_patch),
          partials_patch = ifelse(is.na(partials_patch), 0, partials_patch))

# attach WT fufis
wtfufis_kdm6a = read_xlsx("./databases/Fufi_counts_KDM6A.xlsx", sheet = 1) %>% select(Image_ID, `WT FUFI`) %>% 
   rename(WT_Fufi = `WT FUFI`, Slide_ID = Image_ID) %>% mutate(WT_Fufi = as.numeric(WT_Fufi))
dat_raw_kdm6a = dat_raw_kdm6a %>% left_join(wtfufis_kdm6a)

# drop any missing neighbours
dat_raw_kdm6a = dat_raw_kdm6a %>% mutate(flag = (monoclonals_patch>0 & is.na(knn2))) %>% filter(flag==FALSE) %>%
   mutate(flag = (partials_patch>0 & is.na(knn4))) %>% filter(flag==FALSE) %>% select(-flag)

# fufi data joined across patient/block/slide
data_used_join_kdm = dat_raw_kdm6a %>% 
   left_join(block_to_slide_kdm) %>% 
   left_join(block_to_patient)

# missing slide -> Block/patient conversion?
data_used_join_kdm %>% filter(is.na(Block_ID) | is.na(DW.Patient.number))

## mPAS data
###################
dat_raw_mpas = read_csv("./databases/Fufi_counts_mPAS.csv", col_types = cols(.default = "d", WT.Neighbours.mut.mut = "c",
                                                                              Mut.Neighbours.mut.mut = "c"))
colnames(dat_raw_mpas) = c("Slide_ID", "WT_Fufi", "monoclonals_patch", "wts", "muts", "partials_patch", "wts2", "muts2", 
                            "monoclonals", "partials", "partialclonefufi", "not_needed")

dat_raw_mpas = dat_raw_mpas %>% select(-not_needed)
dat_raw_mpas = dat_raw_mpas %>% mutate(monoclonals = as.integer(monoclonals), 
                                         partials = as.integer(partials),
                                         monoclonals_patch = as.integer(monoclonals_patch), 
                                         partials_patch = as.integer(partials_patch),
                                         partialclonefufi = as.integer(partialclonefufi)) %>% 
   mutate(monoclonals = ifelse(is.na(monoclonals), 0, monoclonals), 
          partials = ifelse(is.na(partials), 0, partials),
          monoclonals_patch = ifelse(is.na(monoclonals_patch), 0, monoclonals_patch),
          partials_patch = ifelse(is.na(partials_patch), 0, partials_patch),
          partialclonefufi = ifelse(is.na(partialclonefufi), 0, partialclonefufi))

# remove repeats with no section
block_to_slide_mpas = block_to_slide_mpas %>% group_by(Slide_ID) %>% 
   mutate(num = n()) %>% ungroup() %>% filter(!(num>1 & is.na(section))) %>% select(-num)
# remove repeats with missing 1 from start of Block_ID
block_to_slide_mpas = block_to_slide_mpas %>% group_by(Slide_ID) %>% 
   mutate(num = n()) %>% ungroup() %>% filter(!(num>1 & Block_ID=="23490")) %>% select(-num)
# remove repeat with misslabelled block 97051/97050 (or wrong repeated slide id)
block_to_slide_mpas = block_to_slide_mpas %>% group_by(Slide_ID) %>% 
   mutate(num = n()) %>% ungroup() %>% filter(!(num>1 & Block_ID=="97051")) %>% select(-num)

# drop any missing neighbours
dat_raw_mpas = dat_raw_mpas %>% mutate(flag = (monoclonals_patch>0 & is.na(wts))) %>% filter(flag==FALSE) %>%
   mutate(flag = (partials_patch>0 & is.na(wts2))) %>% filter(flag==FALSE) %>% select(-flag)

# fufi data joined across patient/block/slide
data_used_join_mpas = dat_raw_mpas %>% 
   left_join(block_to_slide_mpas) %>% 
   left_join(block_to_patient)

# missing slide -> Block/patient conversion?
# data_used_join_mpas %>% filter(is.na(Block_ID) | is.na(DW.Patient.number)) %>% 
#    write_csv("./databases/missing_mpas_slideblock_conversion.csv")

# solving above missing issues
missingdat_mpas = read_csv("./databases/missing_mpas_slideblock_conversion_sorted.csv") %>% 
   mutate(Block_ID=as.character(Block_ID), 
          wts=as.character(wts), muts=as.character(muts))
data_used_join_mpas = data_used_join_mpas %>% 
   filter(!(is.na(Block_ID) | is.na(DW.Patient.number))) %>% bind_rows(missingdat_mpas)

data_used_join_mpas %>% filter(is.na(Block_ID) | is.na(DW.Patient.number))

## make neighbour data follow KDM6A / STAG2 format
data_used_join_mpas = data_used_join_mpas %>% group_by(Slide_ID) %>% 
   mutate(knn1 = paste(paste(str_split(wts, pattern = ",")[[1]], 
                             rep("WT", length(str_split(wts, pattern = ",")[[1]])), sep = " "), 
                        paste(str_split(muts, pattern = ",")[[1]], sep = " "), 
                       rep("mut", length(str_split(muts, pattern = ",")[[1]])),sep = " ", collapse = ", ")) %>% 
   select(Slide_ID, WT_Fufi, monoclonals_patch, knn1, everything()) %>% 
   mutate(knn2 = paste(paste(str_split(wts2, pattern = ",")[[1]], 
                             rep("WT", length(str_split(wts2, pattern = ",")[[1]])), sep = " "), 
                       paste(str_split(muts2, pattern = ",")[[1]], sep = " "), 
                       rep("mut", length(str_split(muts2, pattern = ",")[[1]])),sep = " ", collapse = ", ")) %>% 
   select(Slide_ID:muts, partials_patch, knn2, everything())

## STAG2 data
###################
# dat_raw_stag2 = read_csv("./databases/STAG2_fufi_data_for_Doran_021019.csv")
dat_raw_stag2 = read_csv("./databases/STAG2_fufi_data_for_Doran_110520.csv")
stagdontuse = read_csv("./databases/stag2_bad_slides_dontuse.csv")
dat_raw_stag2 = dat_raw_stag2 %>% filter(!(Slide_ID %in% (stagdontuse %>% pull(Slide_ID))))
# if filtering out slides with "higher levels" of fufis
# dat_raw_stag2 = dat_raw_stag2 %>% filter(Slide_ID %in% (joined_marks_filt %>% filter(mark=="STAG2") %>% pull(Slide_ID)))
dat_raw_stag2 = dat_raw_stag2 %>% rename(monoclonals = mut_fufi,
                                           monoclonals_patch = mut_fufi.patch.border,
                                           partials = mut_WT.fufi,
                                           partials_patch = mut_WT_fufi.patch.border,
                                           WT_Fufi = Fufi_QC)
dat_raw_stag2 = dat_raw_stag2 %>% mutate(monoclonals = as.integer(monoclonals), 
                                           partials = as.integer(partials),
                                           monoclonals_patch = as.integer(monoclonals_patch), 
                                           partials_patch = as.integer(partials_patch)) %>% 
   mutate(monoclonals = ifelse(is.na(monoclonals), 0, monoclonals), 
          partials = ifelse(is.na(partials), 0, partials),
          monoclonals_patch = ifelse(is.na(monoclonals_patch), 0, monoclonals_patch),
          partials_patch = ifelse(is.na(partials_patch), 0, partials_patch)) %>% 
   select(-NFufis, -NMutantCrypts, -NClones, -NPatches)

# drop any missing neighbours
dat_raw_stag2 = dat_raw_stag2 %>% mutate(flag = (monoclonals_patch>0 & is.na(neighbours))) %>% filter(flag==FALSE) %>%
   mutate(flag = (partials_patch>0 & is.na(neighbours.1))) %>% filter(flag==FALSE) %>% select(-flag)

# fufi data joined across patient/block/slide
data_used_join_stag = dat_raw_stag2 %>% 
   left_join(block_to_slide_stag) %>% 
   left_join(block_to_patient)

# missing slide -> Block/patient conversion?
data_used_join_stag %>% filter(is.na(Block_ID) | is.na(DW.Patient.number))

## parse knns
data_used_join_stag = data_used_join_stag %>% mutate(neighbours = str_replace_all(neighbours, ", ", ","), neighbours.1 = str_replace_all(neighbours.1, ", ", ",")) %>% 
   mutate(neighbours = str_replace_all(neighbours, ",", ", "), neighbours.1 = str_replace_all(neighbours.1, ",", ", ")) %>% 
   mutate(neighbours = str_replace_all(neighbours, " , ", ", "), neighbours.1 = str_replace_all(neighbours.1, " , ", ", "))

###################################
### join clone counts and patch edge lengths to patients
take_max_shared_count = function(cnt1, cnt2) {
   temp = cnt1 %>% filter((Slide_ID %in% (cnt2 %>% pull(Slide_ID)))) %>% 
      left_join(cnt2, by = "Slide_ID") %>% group_by(Slide_ID) %>% 
      mutate(indicator = ifelse(NMutantCrypts.x>NMutantCrypts.y, 1, 2))
   temp %>% mutate(patch_edge_length = ifelse(indicator==1, patch_edge_length.x, patch_edge_length.y), 
                   singlet_clones = ifelse(indicator==1, singlet_clones.x, singlet_clones.y), 
                   patch_clones = ifelse(indicator==1, patch_clones.x, patch_clones.y), 
                   NMutantCrypts = ifelse(indicator==1, NMutantCrypts.x, NMutantCrypts.y)) %>% 
      select(-patch_edge_length.x, -patch_edge_length.y, -NMutantCrypts.x, -NMutantCrypts.y, 
             -singlet_clones.x, -singlet_clones.y, -patch_clones.x, -patch_clones.y) %>% ungroup()
}

condense_and_join = function(data_used_join, cnt1, cnt2) {
   # sort out any discrepancies
   sub_info1 = take_max_shared_count(cnt1, cnt2)
   sub_info2 = cnt1 %>% filter(!(Slide_ID %in% (cnt2 %>% pull(Slide_ID))))
   sub_info3 = cnt2 %>% filter(!(Slide_ID %in% c(sub_info1 %>% pull(Slide_ID), sub_info2 %>% pull(Slide_ID))))
   
   outtib = tibble()
   if (nrow(sub_info1)>0) {
      outtib = outtib %>% bind_rows(data_used_join %>% filter(Slide_ID %in% (sub_info1 %>% pull(Slide_ID))) %>% 
         left_join(sub_info1 %>% select(Slide_ID, NMutantCrypts, patch_edge_length, singlet_clones, patch_clones)))
   }
   if (nrow(sub_info2)>0) {
      outtib = outtib %>% bind_rows(data_used_join %>% filter(Slide_ID %in% (sub_info2 %>% pull(Slide_ID))) %>% 
         left_join(sub_info2 %>% select(Slide_ID, NMutantCrypts, patch_edge_length, singlet_clones, patch_clones)))
   }
   if (nrow(sub_info3)>0) {
      outtib = outtib %>% bind_rows(data_used_join %>% filter(Slide_ID %in% (sub_info3 %>% pull(Slide_ID))) %>% 
         left_join(sub_info3 %>% select(Slide_ID, NMutantCrypts, patch_edge_length, singlet_clones, patch_clones)))
   }
   return(outtib)
}

## kdm6a
patchedgelengths_kdm = read_tsv("./databases/patch_edge_lengths_kdm6a.tsv", 
                                col_types = cols(.default = "d", patch_clones = "c")) %>% 
   distinct() %>% mutate(patch_clones = ifelse(is.na(patch_clones), "0", patch_clones)) # DeCryptICS scores
patchedgelengths2_kdm = read_tsv("./databases/patch_edge_lengths2_KDM6A.tsv", 
                                 col_types = cols(.default = "d", patch_clones = "c")) %>% 
   mutate(Block_ID = as.character(Block_ID)) %>% distinct() %>% 
   mutate(patch_clones = str_remove(patch_clones, ",0"))
patchedgelengths2_kdm = patchedgelengths2_kdm %>% left_join(block_to_slide_kdm)
# sort out any discrepancies, priorities manual scoring
data_used_join_kdm = condense_and_join(data_used_join_kdm, patchedgelengths2_kdm, 
                                       patchedgelengths_kdm %>% filter(!(Slide_ID%in%(patchedgelengths2_kdm %>% pull(Slide_ID))))) %>% distinct()

## mpas
patchedgelengths_mpas = read_tsv("./databases/patch_edge_lengths_mpas.tsv",
                                 col_types = cols(.default = "d", patch_clones = "c", 
                                                  DW.Patient.number="c", Block_ID="c")) %>% 
   distinct() %>% mutate(patch_clones = ifelse(is.na(patch_clones), "0", patch_clones))
# patchedgelengths2_mpas = read_tsv("./databases/patch_edge_lengths2_mpas.tsv") %>% distinct() # just patient-labelled
patchedgelengths3_mpas = read_tsv("./databases/patch_edge_lengths_mpas3.tsv", 
                                  col_types = cols(.default = "d", patch_clones = "c")) %>% 
   distinct() %>% mutate(patch_clones = ifelse(is.na(patch_clones), "0", patch_clones))
patchedgelengths_mpas = trim_block_ID(patchedgelengths_mpas) %>% select(-section)
# sort out any discrepancies
data_used_join_mpas = condense_and_join(data_used_join_mpas, patchedgelengths_mpas, patchedgelengths3_mpas) %>% distinct()

## stag2
patchedgelengths_stag = read_tsv("./databases/patch_edge_lengths_stag2.tsv", 
                                 col_types = cols(.default = "d", patch_clones = "c")) %>% 
   mutate(Block_ID = as.character(Block_ID)) %>% distinct() %>% 
   mutate(patch_clones = str_remove(patch_clones, ",0")) %>% 
   left_join(block_to_slide_stag)
patchedgelengths2_stag = read_tsv("./databases/patch_edge_lengths2_stag2.tsv", 
                                  col_types = cols(.default = "d", patch_clones = "c")) %>% 
   mutate(Block_ID = as.character(Block_ID)) %>% distinct() %>% 
   mutate(patch_clones = str_remove(patch_clones, ",0"))
patchedgelengths2_stag = patchedgelengths2_stag %>% filter(!is.nan(NMutantCrypts) & !is.na(NMutantCrypts) & !is.nan(Slide_ID) & !is.na(Slide_ID))
patchedgelengths_stag = patchedgelengths_stag %>% filter(!is.nan(NMutantCrypts) & !is.na(NMutantCrypts) & !is.nan(Slide_ID) & !is.na(Slide_ID))
# sort out any discrepancies
data_used_join_stag = condense_and_join(data_used_join_stag, patchedgelengths_stag, patchedgelengths2_stag) %>% distinct()

##################
## comparing data sources
data_used_join_mpas %>% select(monoclonals, partials, monoclonals_patch, partials_patch) %>% summarise_all(~sum(.))
data_used_join_mpas %>% filter(DW.Patient.number!="CRA_00370") %>% 
   select(monoclonals, partials, monoclonals_patch, partials_patch) %>% summarise_all(~sum(.))
data_used_join_kdm %>% select(monoclonals, partials, monoclonals_patch, partials_patch) %>% summarise_all(~sum(.))
data_used_join_stag %>% select(monoclonals, partials, monoclonals_patch, partials_patch) %>% summarise_all(~sum(.))
data_used_join_stag %>% filter(DW.Patient.number!="CRA_00370") %>% 
   select(monoclonals, partials, monoclonals_patch, partials_patch) %>% summarise_all(~sum(.))

## looking at patch sizes and partial fufis
data_used_join_mpas %>% filter(partials>0 | partials_patch>0) %>% select(Slide_ID, partials_patch, partials, 
                                                                         singlet_clones, patch_clones, patch_edge_length)
data_used_join_kdm %>% filter(partials>0 | partials_patch>0) %>% select(Slide_ID, partials_patch, partials, 
                                                                        singlet_clones, patch_clones, patch_edge_length)
data_used_join_stag %>% filter(partials>0 | partials_patch>0) %>% select(Slide_ID, partials_patch, partials, 
                                                                         singlet_clones, patch_clones, patch_edge_length)
data_used_join_stag %>% filter(partials>0 | partials_patch>0) %>% select(Slide_ID, partials_patch, partials, 
                                                                         singlet_clones, patch_clones, patch_edge_length)

data_used_join_stag %>% summarise(sum_partials_patch = sum(partials_patch), sum_partials = sum(partials),
                                  sum_singlet_clones = sum(singlet_clones), sum_patch_edge = sum(patch_edge_length-singlet_clones)) %>% 
   mutate(singlet_partial_fraction = sum_partials / sum_singlet_clones,
          patch_partial_fraction = sum_partials_patch / sum_patch_edge)

data_used_join_kdm %>% summarise(sum_partials_patch = sum(partials_patch), sum_partials = sum(partials),
                                 sum_singlet_clones = sum(singlet_clones), sum_patch_edge = sum(patch_edge_length-singlet_clones)) %>% 
   mutate(singlet_partial_fraction = sum_partials / sum_singlet_clones,
          patch_partial_fraction = sum_partials_patch / sum_patch_edge)

data_used_join_mpas %>% summarise(sum_partials_patch = sum(partials_patch), sum_partials = sum(partials),
                                  sum_singlet_clones = sum(singlet_clones), sum_patch_edge = sum(patch_edge_length-singlet_clones)) %>% 
   mutate(singlet_partial_fraction = sum_partials / sum_singlet_clones,
          patch_partial_fraction = sum_partials_patch / sum_patch_edge)


## outputting data used in raw calculations
data_used_join_kdm %>% filter(!(DW.Patient.number%in%bad_pats)) %>% 
   select(DW.Patient.number, Sex, Age, Block_ID, section, Slide_ID) %>% mutate(mark = "KDM6A") %>% 
   bind_rows(data_used_join_mpas %>% filter(!(DW.Patient.number%in%bad_pats)) %>% 
                select(DW.Patient.number, Sex, Age, Block_ID, section, Slide_ID) %>% mutate(mark = "mPAS")) %>% 
   bind_rows(data_used_join_stag %>% filter(!(DW.Patient.number%in%bad_pats)) %>% 
                select(DW.Patient.number, Sex, Age, Block_ID, section, Slide_ID) %>% mutate(mark = "STAG2")) %>%
   write_csv("./output/compiled_patient_list_fusionrate.csv")

# if (all((data_used_join_stag %>% 
#          pull(Slide_ID) %>% 
#          sort()) == (dat_raw_stag2 %>% 
#                      pull(Slide_ID) %>% 
#                      sort()))) { ## check we have the same data set
#    data_used_join_stag %>% 
#       select(DW.Patient.number, Age, Sex, Block_ID, section, Slide_ID) %>% 
#       mutate(mark = "STAG2") %>% 
#       write_csv("./output/fusion_raw_patient_data_STAG2.csv")
#    data_used_join_stag %>% 
#       select(DW.Patient.number, Age, Sex) %>% distinct() %>% 
#       mutate(mark = "STAG2") %>% 
#       write_csv("./output/fusion_raw_patient_data_pooledslides_STAG2.csv")
# }
# if (all((data_used_join_kdm %>% 
#          pull(Slide_ID) %>% 
#          sort()) == (dat_raw_kdm6a %>% 
#                      pull(Slide_ID) %>% 
#                      sort()))) { ## check we have the same data set
#    data_used_join_kdm %>% 
#       select(DW.Patient.number, Age, Sex, Block_ID, section, Slide_ID) %>% 
#       mutate(mark = "KDM6A") %>% 
#       write_csv("./output/fusion_raw_patient_data_KDM6A.csv")
#    data_used_join_kdm %>% 
#       select(DW.Patient.number, Age, Sex) %>% distinct() %>% 
#       mutate(mark = "KDM6A") %>% 
#       write_csv("./output/fusion_raw_patient_data_pooledslides_KDM6A.csv")
# }
# if (all((data_used_join_mpas %>% 
#          pull(Slide_ID) %>% 
#          sort()) == (dat_raw_mpas %>% 
#                      pull(Slide_ID) %>% 
#                      sort()))) { ## check we have the same data set
#    data_used_join_mpas %>%
#       select(DW.Patient.number, Age, Sex, Block_ID, section, Slide_ID) %>% 
#       mutate(mark = "mPAS") %>% 
#       write_csv("./output/fusion_raw_patient_data_mPAS.csv")
#    data_used_join_mpas %>% 
#       select(DW.Patient.number, Age, Sex) %>% distinct() %>% 
#       mutate(mark = "mPAS") %>% 
#       write_csv("./output/fusion_raw_patient_data_pooledslides_mPAS.csv")
# }

## getting image locations
# block_to_slide_mpas = read_xlsx("~/Work/images/Image_locations_all_marks.xlsx", sheet = "mPAS") %>%
#    select(`Block ID`, `Image ID`, `File Location`)
# colnames(block_to_slide_mpas) = c("Block_ID", "Slide_ID", "File_Location")
# block_to_slide_mpas = trim_block_ID(block_to_slide_mpas %>% distinct())
# block_to_slide_mpas = block_to_slide_mpas %>% filter(!is.na(Block_ID)) %>% distinct()
# # datin_mpas %>% select(Slide_ID, Block_ID, section, DW.Patient.number, Age, Sex) %>% left_join(block_to_slide_mpas) %>% write_csv("./output/mPAS_slide_locations.csv")
# # datin_mpas %>% select(Slide_ID, Block_ID, section, DW.Patient.number, Age, Sex) %>% left_join(block_to_slide_mpas) %>% View()
# datin_mpas %>% mutate(Mut_Fufi_frac = Mut_Fufi / (NMutantCrypts+1e-14), 
#                       WT_Fufi_frac = WT_Fufi / (NCrypts - NMutantCrypts)) %>% 
#    filter(Mut_Fufi_frac>0.01) %>% 
#    select(Slide_ID, Block_ID, section, DW.Patient.number, Age, Sex, Mut_Fufi_frac) %>% 
#    left_join(block_to_slide_mpas) %>% arrange(desc(Mut_Fufi_frac)) %>% 
#    write_csv("./output/mPAS_slide_locations_highfufi.csv")
