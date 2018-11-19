# ipyrad_tools

A collection of scripts that I use to change ipyrad output into other formats:

## adegenet2migrate

Scripts to create an input file for migrate-n (https://peterbeerli.com/migrate-html5/index.html) from a genepop file.

It uses adegenet (http://adegenet.r-forge.r-project.org) to load the genepop file into R, filter non-biallelic SNPs, those not under Hardy-Weiberg equilibrium, and those not sequenced for at least one individual per population.

The genepop file was created using unlinked_snps_to_genepop.py, which will be uploaded to this repository.
