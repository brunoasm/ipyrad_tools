#This script contains a function to write a HapMap file for migrate-n from an adegenet genind object
#Written by B. de Medeiros, starting on 19-Nov-2018
library(adegenet)
library(stringr)
library(dplyr)

genpop2hapmap = function(genpop_obj, project_title = 'MIGRATE'){
  
  snp1_counts = genpop_obj %>%
    genind2genpop %>%
    as.data.frame() %>%
    .[,seq(1,dim(.)[2]-1,by = 2)]
  
  snp2_counts = genpop_obj %>%
    genind2genpop %>%
    as.data.frame() %>%
    .[,seq(2,dim(.)[2],by = 2)]
  
  sample_sizes = snp1_counts + snp2_counts
  
  in_all_pops = sample_sizes %>%
    apply(.,2,function(x){sum(x!=0) >= length(x)})
  
  snp1_counts_filtered = snp1_counts[,in_all_pops]
  snp2_counts_filtered = snp2_counts[,in_all_pops]
  sample_sizes_filtered = sample_sizes[,in_all_pops]
  
  #these translations work because I made genpop files myself from sequences
  #they are not universal
  nuc_translate = c('001' = 'A',
                    '004' = 'T',
                    '003' = 'G',
                    '002' = 'C')
  
  snp1_nucs = colnames(snp1_counts_filtered) %>%
    str_remove('^.+\\.') %>%
    nuc_translate[.]
  
  snp2_nucs = colnames(snp2_counts_filtered) %>%
    str_remove('^.+\\.') %>%
    nuc_translate[.]
  
  output_string = '#HapMap input for migrate-n made from adegenet file'
  output_string = c(output_string, '#Scripts written by B. de Medeiros, starting on Nov 2018')
  
  output_string = c(output_string,str_c('H',
                                        nrow(snp1_counts_filtered), 
                                        ncol(snp1_counts_filtered),
                                        project_title,
                                        sep = ' '))
  
  for (i in 1:nrow(snp1_counts_filtered)){
    n_samples = max(sample_sizes_filtered[i,])
    output_string = c(output_string, str_c(n_samples,
                                           rownames(sample_sizes_filtered)[1],
                                           sep = ' '))
    for (j in 1:ncol(snp1_counts_filtered)){
      output_string = c(output_string, str_c(j,
                                             snp1_nucs[j],
                                             snp1_counts_filtered[i,j],
                                             snp2_nucs[j],
                                             snp2_counts_filtered[i,j],
                                             sample_sizes_filtered[i,j],
                                             sep = '\t'))
    }
  }
  
  return(output_string)
}
