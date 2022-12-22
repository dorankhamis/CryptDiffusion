import glob
import pandas as pd
import xlrd

def hexagonal_patch_edge_length(patchsize):
   ''' Hexagonal (perfect) patch sizes go up as 1 + 1*6 + 2*6 + 3*6 + ...
   where the last part gives the patch edge length
   we use the partial sum of the sum of integers:
   patchsize = 1 + 6 * sum_{n=1}^{n=N}(n) ,
   where sum_{n=1}^{n=N}(n) = N(N+1)/2,
   so we solve the quadratic for N and then multiply
   by 6 for the patch edge length.
   If less than 7, all crypts are on the edge anyway.
   '''
   if patchsize<7:
      return patchsize
   else:
      partialsum_N = -0.5 + np.sqrt( (patchsize-1)/3. + 0.25 )
      return int(np.ceil(partialsum_N * 6))

#####################
## KDM6A
basepath = "/home/doran/Work/images/KDM6A_fufi_analysis/Analysed_slides/"
slide_folders = glob.glob(basepath + "Analysed_*")
slide_counts = np.loadtxt(basepath + "slide_counts.csv", skiprows=1, delimiter = ",")      
outfolder = "/home/doran/Work/r_code/crypt_fufi_hexgrid/databases/"
with open(outfolder + "patch_edge_lengths_kdm6a.tsv", 'w') as fo:
   fo.write("Slide_ID\tpatch_edge_length\tNMutantCrypts\n")   
   for i in range(slide_counts.shape[0]):
      patchedgelength = 0
      N_mutant_crypts = 0
      slide_id = str(int(slide_counts[i,0]))
      slidepath = basepath + 'Analysed_' + slide_id + '/'
      patchsizes = np.loadtxt(slidepath + 'patch_sizes.txt', ndmin = 1)
      # find single mutant crypts
      singlemutants = slide_counts[i,3] - sum(patchsizes)
      patchedgelength += singlemutants
      N_mutant_crypts += singlemutants
      # account for patches
      for j in range(len(patchsizes)):
         patchedgelength += hexagonal_patch_edge_length(patchsizes[j])
         N_mutant_crypts += patchsizes[j]
      fo.write("%s\t%d\t%d\n" % (slide_id, patchedgelength, N_mutant_crypts) )

basefolder2 = "/home/doran/immwork/r_code/crypt_fufi_hexgrid/databases/"
slide_counts = pd.read_excel(basefolder2 + "KDM6A_scoring_summary.xlsx")
with open(basefolder2 + "patch_edge_lengths2_KDM6A.tsv", 'w') as fo:
   fo.write("Block_ID\tsection\tpatch_edge_length\tNMutantCrypts\n")
   for i in range(slide_counts.shape[0]):
      patchedgelength = 0
      N_mutant_crypts = 0
      if not np.isnan(slide_counts['Block ID'][i]):
         block_id = str(int(slide_counts['Block ID'][i]))
         sec_id = int(slide_counts['section'][i])
         # find single mutant crypts
         singlemutants = slide_counts['Total WPC'][i]
         patchedgelength += singlemutants
         N_mutant_crypts += singlemutants
         # account for patches
         for j in range(2,15):
            numthissize = slide_counts[j][i]
            if not np.isnan(numthissize):
               numthissize = int(numthissize)
               for k in range(numthissize):
                  patchedgelength += hexagonal_patch_edge_length(j)
                  N_mutant_crypts += j
         # combined scorings above 14
         numsthissize = slide_counts['15+'][i]
         if type(numsthissize)==str:
            numsthissize = [int(s) for s in numsthissize.split(',')]
            for k in range(len(numsthissize)):
               patchedgelength += hexagonal_patch_edge_length(numsthissize[k])
               N_mutant_crypts += numsthissize[k]
         else:
            if not np.isnan(numsthissize):
               onepatchthissize = int(numsthissize)
               patchedgelength += hexagonal_patch_edge_length(onepatchthissize)
               N_mutant_crypts += onepatchthissize
         fo.write("%s\t%d\t%1.0f\t%1.0f\n" % (block_id, sec_id, patchedgelength, N_mutant_crypts) )
      
#####################
## STAG2
basefolder = "/home/doran/Work/r_code/crypt_fufi_hexgrid/databases/"
slide_counts = pd.read_excel(basefolder + "Stag2_scoring_summary.xlsx")
#Index([      'section', 'DW patient ID',      'Block ID',           'Box',
#       'stain quality',           'PPC',     'Fractions',     'Total WPC',
#                     2,               3,               4,               5,
#                     6,               7,               8,               9,
#                    10,              11,              12,              13,
#                    14,           '15+',      'Comments'],
#      dtype='object')
      
with open(basefolder + "patch_edge_lengths_stag2.tsv", 'w') as fo:
   fo.write("Block_ID\tsection\tpatch_edge_length\tNMutantCrypts\n")
   for i in range(slide_counts.shape[0]):
      patchedgelength = 0
      N_mutant_crypts = 0
      if not np.isnan(slide_counts['Block ID'][i]):
         block_id = str(int(slide_counts['Block ID'][i]))
         sec_id = int(slide_counts['section'][i])
         # find single mutant crypts
         singlemutants = slide_counts['Total WPC'][i]
         patchedgelength += singlemutants
         N_mutant_crypts += singlemutants
         # account for patches
         for j in range(2,15):
            numthissize = slide_counts[j][i]
            if not np.isnan(numthissize):
               numthissize = int(numthissize)
               for k in range(numthissize):
                  patchedgelength += hexagonal_patch_edge_length(j)
                  N_mutant_crypts += j
         # combined scorings above 14
         numsthissize = slide_counts['15+'][i]
         if type(numsthissize)==str:
            numsthissize = [int(s) for s in numsthissize.split(',')]
            for k in range(len(numsthissize)):
               patchedgelength += hexagonal_patch_edge_length(numsthissize[k])
               N_mutant_crypts += numsthissize[k]
         else:
            if not np.isnan(numsthissize):
               onepatchthissize = int(numsthissize)
               patchedgelength += hexagonal_patch_edge_length(onepatchthissize)
               N_mutant_crypts += onepatchthissize
         print("i = %d" % i)
         print(N_mutant_crypts)
         fo.write("%s\t%d\t%1.0f\t%1.0f\n" % (block_id, sec_id, patchedgelength, N_mutant_crypts))

basefolder = "/home/doran/Work/r_code/crypt_fufi_hexgrid/databases/"
slide_counts = pd.read_csv(basefolder + "missing_stag2_clonecounts_sorted.csv")
with open(basefolder + "patch_edge_lengths2_stag2.tsv", 'w') as fo:
   fo.write("Slide_ID\tBlock_ID\tsection\tpatch_edge_length\tNMutantCrypts\n")
   for i in range(slide_counts.shape[0]):
      patchedgelength = 0
      N_mutant_crypts = 0
      if not np.isnan(slide_counts['Block_ID'][i]):
         block_id = str(int(slide_counts['Block_ID'][i]))
         slide_id = int(slide_counts['Slide_ID'][i])
         sec_id = str(int(slide_counts['section'][i]))
         # find single mutant crypts
         singlemutants = slide_counts['Total WPC'][i]
         patchedgelength += singlemutants
         N_mutant_crypts += singlemutants
         # account for patches
         for j in range(2,15):
            numthissize = slide_counts[str(j)][i]
            if not np.isnan(numthissize):
               numthissize = int(numthissize)
               for k in range(numthissize):
                  patchedgelength += hexagonal_patch_edge_length(j)
                  N_mutant_crypts += j
         # combined scorings above 14
         numsthissize = slide_counts['15+'][i]
         if type(numsthissize)==str:
            numsthissize = [int(s) for s in numsthissize.split(',')]
            for k in range(len(numsthissize)):
               patchedgelength += hexagonal_patch_edge_length(numsthissize[k])
               N_mutant_crypts += numsthissize[k]
         else:
            if not np.isnan(numsthissize):
               onepatchthissize = int(numsthissize)
               patchedgelength += hexagonal_patch_edge_length(onepatchthissize)
               N_mutant_crypts += onepatchthissize     
         fo.write("%d\t%s\t%s\t%1.0f\t%1.0f\n" % (slide_id, block_id, sec_id, patchedgelength, N_mutant_crypts) )

#####################
## mPAS
outfolder = "/home/doran/Work/r_code/crypt_fufi_hexgrid/databases/"
#basefolder = "/home/doran/Work/r_code/drift_All_Marks/Clone_data/mPAS/"
#slide_counts = pd.read_csv(basefolder + "mPAS_pat_table.csv")
# better data:
slide_counts = pd.read_csv(outfolder + "all_mpas_slides.csv")

#Index(['DW.Patient.number', 'Age', 'SEX', 'mutant_hom', 'num_crypts',
#       'anyClone', 'single_cell', 'partial', 'full', 'fullx2', 'fullx3',
#       'fullx4', 'fullx5', 'fullx6', 'fullx7', 'fullx8', 'fullx9', 'fullx10',
#       'fullxmore10', 'monoclonal', 'mPAS_het_Anna', 'freq_mono', 'freq_part',
#       'freq_single'],
#      dtype='object')

patchlist = ['full', 'fullx2', 'fullx3', 'fullx4', 'fullx5', 
             'fullx6', 'fullx7', 'fullx8', 'fullx9', 'fullx10']

with open(outfolder + "patch_edge_lengths_mpas.tsv", 'w') as fo:
   fo.write("DW.Patient.number\tSlide_ID\tBlock_ID\tpatch_edge_length\tNMutantCrypts\n")
   for i in range(slide_counts.shape[0]):
      patchedgelength = 0
      N_mutant_crypts = 0
      if type(slide_counts['DW.Patient.number'][i])==str:
         pat_id = slide_counts['DW.Patient.number'][i]
         slide_id = int(slide_counts['Slide'][i])
         block_id = str(slide_counts['Block_ID'][i])
         # account for patches
         for j in range(len(patchlist)):
            col_id = patchlist[j]
            numthissize = slide_counts[col_id][i]
            if not np.isnan(numthissize):
               numthissize = int(numthissize)
               for k in range(numthissize):
                  patchedgelength += hexagonal_patch_edge_length(j+1)
                  N_mutant_crypts += j+1
         ## combined scorings above 10         
         numsthissize = slide_counts['fullxmore10'][i]
#         # if string of actual size
#         if type(numsthissize)==str:
#            numsthissize = [int(s) for s in numsthissize.split(',')]
#            for k in range(len(numsthissize)):
#               patchedgelength += hexagonal_patch_edge_length(numsthissize[k])
#         else:
#            if not np.isnan(numsthissize):
#               onepatchthissize = int(numsthissize)
#               patchedgelength += hexagonal_patch_edge_length(onepatchthissize)
         # if just a count
         if not np.isnan(numthissize):
            numthissize = int(numthissize)
            for k in range(numthissize):
               patchedgelength += hexagonal_patch_edge_length(11+k)
               N_mutant_crypts += 11+k
         print(slide_id)
         print(N_mutant_crypts)
         fo.write("%s\t%d\t%s\t%1.0f\t%1.0f\n" % (pat_id, slide_id, block_id, patchedgelength, N_mutant_crypts) )

# more
outfolder = "/home/doran/Work/r_code/crypt_fufi_hexgrid/databases/"
basefolder = "/home/doran/Work/r_code/drift_All_Marks/Clone_data/mPAS/"
slide_counts = pd.read_csv(basefolder + "mPAS_pat_table.csv")

#Index(['DW.Patient.number', 'Age', 'SEX', 'mutant_hom', 'num_crypts',
#       'anyClone', 'single_cell', 'partial', 'full', 'fullx2', 'fullx3',
#       'fullx4', 'fullx5', 'fullx6', 'fullx7', 'fullx8', 'fullx9', 'fullx10',
#       'fullxmore10', 'monoclonal', 'mPAS_het_Anna', 'freq_mono', 'freq_part',
#       'freq_single'],
#      dtype='object')

patchlist = ['full', 'fullx2', 'fullx3', 'fullx4', 'fullx5', 
             'fullx6', 'fullx7', 'fullx8', 'fullx9', 'fullx10']

with open(outfolder + "patch_edge_lengths2_mpas.tsv", 'w') as fo:
   fo.write("DW.Patient.number\tpatch_edge_length\tNMutantCrypts\n")
   for i in range(slide_counts.shape[0]):
      patchedgelength = 0
      N_mutant_crypts = 0
      if type(slide_counts['DW.Patient.number'][i])==str:
         pat_id = slide_counts['DW.Patient.number'][i]
         # account for patches
         for j in range(len(patchlist)):
            col_id = patchlist[j]
            numthissize = slide_counts[col_id][i]
            if not np.isnan(numthissize):
               numthissize = int(numthissize)
               for k in range(numthissize):
                  patchedgelength += hexagonal_patch_edge_length(j+1)
                  N_mutant_crypts += j+1
         ## combined scorings above 10         
         numsthissize = slide_counts['fullxmore10'][i]
         # if string of actual size
         if type(numsthissize)==str:
            numsthissize = [int(s) for s in numsthissize.split(',')]
            for k in range(len(numsthissize)):
               patchedgelength += hexagonal_patch_edge_length(numsthissize[k])
               N_mutant_crypts += numsthissize[k]
         else:
            if not np.isnan(numsthissize):
               onepatchthissize = int(numsthissize)
               patchedgelength += hexagonal_patch_edge_length(onepatchthissize)
               N_mutant_crypts += onepatchthissize
         print(pat_id)
         print(N_mutant_crypts)
         fo.write("%s\t%1.0f\t%1.0f\n" % (pat_id, patchedgelength, N_mutant_crypts))         
         
# even more
basefolder = "/home/doran/Work/r_code/crypt_fufi_hexgrid/databases/"
slide_counts = pd.read_csv(basefolder + "slide_counts_mPAS_combined_041119.csv")
#Index(['Slide_ID', 'NCrypts', 'NFufis', 'NMutantCrypts', 'NClones', 'NPatches',
#       'Partials', 'Fractions', 'Pa2', 'Pa3', 'Pa4', 'Pa5', 'Pa6', 'Pa7',
#       'Pa8', 'Pa9', 'Pa10', 'Pa11', 'comments'],
#      dtype='object')
      
with open(basefolder + "patch_edge_lengths_mpas3.tsv", 'w') as fo:
   fo.write("Slide_ID\tpatch_edge_length\tNMutantCrypts\n")
   for i in range(slide_counts.shape[0]):
      patchedgelength = 0
      N_mutant_crypts = 0
      if not np.isnan(slide_counts['Slide_ID'][i]):
         sld_id = str(int(slide_counts['Slide_ID'][i]))
         # account for patches
         totalmuts = 0
         for j in range(2,12):
            numthissize = slide_counts['Pa'+str(j)][i]
            if not np.isnan(numthissize):
               numthissize = int(numthissize)
               for k in range(numthissize):
                  totalmuts += j
                  patchedgelength += hexagonal_patch_edge_length(j)
                  N_mutant_crypts += j
         # find single mutant crypts
         singlemutants = slide_counts['NMutantCrypts'][i] - totalmuts
         patchedgelength += singlemutants
         N_mutant_crypts += singlemutants
         # sort out partials
         if not np.isnan(slide_counts['Partials'][i]):
            patchedgelength += slide_counts['Partials'][i]
            N_mutant_crypts += slide_counts['Partials'][i]
         print(sld_id)
         print(N_mutant_crypts)
         fo.write("%s\t%1.0f\t%1.0f\n" % (sld_id, patchedgelength, N_mutant_crypts))
         
         
