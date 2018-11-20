Making an input file for migrate-n (HapMap format)
================

Here we will create an input file for migrate-n in the HapMap format.

See manual for more information:
<https://peterbeerli.com/programs/migrate/distribution_4.x/migratedoc4.x.pdf>

We will start from a genpop file, filter alleles to only biallelic, then
to those under hardy-weinberg equilibrium

We will then use a custom function to write the gening object as an
input file for migrate-n. We will only include in this file sites that
were sequenced for at least one individual per population.

First, let’s load the function that writes a HapMap string:

``` r
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
                                           rownames(sample_sizes_filtered)[i],
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
```

This function can also be loaded from a script:

``` r
source('genpop2hapmap.R')
```

    ## Loading required package: ade4

    ## 
    ##    /// adegenet 2.1.1 is loaded ////////////
    ## 
    ##    > overview: '?adegenet'
    ##    > tutorials/doc/questions: 'adegenetWeb()' 
    ##    > bug reports/feature requests: adegenetIssues()

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

Now, let’s load required libraries:

``` r
library(adegenet)
library(pegas)
```

    ## Loading required package: ape

    ## 
    ## Attaching package: 'pegas'

    ## The following object is masked from 'package:ape':
    ## 
    ##     mst

    ## The following object is masked from 'package:ade4':
    ## 
    ##     amova

``` r
library(dplyr)
library(stringr)
```

Now we read the input file:

``` r
snps = read.genepop('test_data.gen',ncode=3)
```

    ## 
    ##  Converting data from a Genepop .gen file to a genind object... 
    ## 
    ## 
    ## File description:  test data created from an ipyrad unlinked SNPS file 
    ## 
    ## ...done.

We will filter out non-biallelic SNPs and sites not under HW
equilibrium:

``` r
biallelic_snps = which(snps@loc.n.all == 2) %>% 
  names %>%
  match(snps@loc.fac,.,nomatch = NA) %>%
  is.na(.) %>%
  `!` %>%
  snps[,.]

hwe_res = hw.test(biallelic_snps)

snps_filtered = which(hwe_res[,"Pr.exact"] > 0.05) %>%
  names %>%
  match(biallelic_snps@loc.fac,.,nomatch = NA) %>%
  is.na(.) %>%
  `!` %>%
  biallelic_snps[,.]
```

Finally, let’s write a string in the HapMap format using our custom
function. This function automatically filters out sites not sequenced
for all populations:

``` r
out_string = genpop2hapmap(genpop_obj = snps_filtered,
                           project_title = 'TEST')
```

    ## 
    ##  Converting data from a genind to a genpop object... 
    ## 
    ## ...done.
    ## 
    ## 
    ##  Converting data from a genind to a genpop object... 
    ## 
    ## ...done.

Let’s now write out the file:

``` r
write(out_string,file = 'test.hapmap')
```
