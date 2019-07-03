#!/bin/bash

# 32 vCPUs 208GB memory machine for Verneri's phenotype file.
# start with 4 vCPUs w 32gb mem each... maybe need more?
# Boot disk - change default disk size from 10GB to 200GB. (Debian GNU/Linux - though this doesn't matter too much).
#start off with not much memory while getting things set up, then scale up for the actual phenotyping.
# load up a terminal window, under 'Remote access' clicked ssh and then 'open in browser window'

sudo apt-get update
HOME=/home/eatkinso
INSTALL_DIR=/home/eatkinso/R/

# Install a bunch of packages to get R up and running on the cluster.
PKGS="bzip2 build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libpcre3-dev gfortran openjdk-8-jdk"
for p in $PKGS; do
    sudo apt-get install -y $p
done

# Install R 3.4.
wget https://cran.r-project.org/src/base/R-3/R-3.4.1.tar.gz
tar -xvzf R-3.4.1.tar.gz
rm R-3.4.1.tar.gz

cd R-3.4.1
./configure --with-readline=no --with-x=no --prefix=$INSTALL_DIR
make
make install

# Add R to the path. 
#export PATH=$INSTALL_DIR/bin:$PATH
export PATH="/home/eatkinso/R-3.4.1/bin:$PATH" 

#set up an R mirror
#loaded up R, then selected 65, a USA one in TN. Put into my general R profile for future use
vi ~/.Rprofile:

local({r <- getOption("repos")
   r["CRAN"] <- "https://mirrors.nics.utk.edu/cran" 
   options(repos=r)
})

# Install a bunch of packages:
#Rscript -e 'install.packages(c("data.table", repos="http://cran.rstudio.com")'
Rscript -e 'install.packages(c("data.table"))'
Rscript -e 'install.packages(c("bit64"))'

# Install git.
sudo apt-get install git
cd $HOME

# clone PHESANT library
git clone https://github.com/astheeggeggs/PHESANT.git

# Move the required reengineering_phenofile.r over.
# It just renames some columns and has some rules about which of the visits we include
##### Don't have access... maybe ask Sam about this
#gsutil cp gs://ukbb_association/reengineering_phenofile_neale_lab.r $HOME
#Duncan moved to a location I do have access to hopefully:
gsutil cp gs://phenotype_31063/reengineering_phenofile_neale_lab.r $HOME

##still no access. Try setting up locally, with glcoud
alias gsutil=/Users/elizabeth/Desktop/google-cloud-sdk/bin/gsutil 
alias gcloud=/Users/elizabeth/Desktop/google-cloud-sdk/bin/gcloud
###setting up Hail
export HAIL_HOME=/Users/elizabeth/Desktop/hail
export PATH=$PATH:$HAIL_HOME/bin/
alias gsutil=/Users/elizabeth/Desktop/google-cloud-sdk/bin/gsutil 
PATH=$PATH:/Users/elizabeth/Desktop/google-cloud-sdk/bin
alias gcloud=/Users/elizabeth/Desktop/google-cloud-sdk/bin/gcloud
alias l='ls -lhtr'
#gcloud init

#change the new configuration name: ukbb-diverse-pops

#see if can copy over outside the VM which might not have access
#gsutil cp gs://phenotype_31063/reengineering_phenofile_neale_lab.r $HOME
gsutil cp gs://phenotype_31063/reengineering_phenofile_neale_lab.r gs://ukb-diverse-pops/Phenotypes/reengineering_phenofile_neale_lab.r


# Move the raw phenotype file.
#gsutil cp gs://phenotype_31063/ukb11214.csv $HOME
gsutil cp gs://phenotype_31063/ukb11214.csv gs://ukb-diverse-pops/Phenotypes/ukb11214.csv

#copy these all into the VM also from my Google Cloud Phenotypes directory
gsutil cp gs://ukb-diverse-pops/Phenotypes/reengineering_phenofile_neale_lab_ega.r $HOME
gsutil cp gs://ukb-diverse-pops/Phenotypes/reengineering_phenofile_neale_lab_ega2.r $HOME

gsutil cp gs://ukb-diverse-pops/Phenotypes/ukb11214.csv $HOME
#needed to crank the boot disk size up to 400GB to get this to work. Need 25 for this file, another 25 for the output fil, and then R will need maybe 10fold that to be able to run
#also made a tmp directory, which it wanted to use for the transfer 

# Run the reegineering_phenofile.r on the raw phenotype data.
export PATH="/home/eatkinso/R-3.4.1/bin:$PATH" 

#errored out running reengineering script. need to install.packages('bit64') to let it handle the integer64 style entries properly...
#Warning message: In require_bit64() : Some columns are type 'integer64' but package bit64 is not installed. Those columns will print as strange looking floating point data. There is no need to reload the data. Simply install.packages('bit64') to obtain the integer64 print method and print the data again.
#Error in is.data.frame(x) : object 'bd' not found
#Calls: gsub -> colnames -> is.data.frame
#Execution halted
Rscript -e 'install.packages(c("bit64"))'

Rscript reengineering_phenofile_neale_lab_ega.r 
Rscript reengineering_phenofile_neale_lab_ega2.r 

# This will write a .tsv file that we then execute a PHESANT script on.
#copy over the new formatted pheno file to my google cloud bucket
gsutil cp neale_lab_parsed.tsv gs://ukb-diverse-pops/Phenotypes/neale_lab_parsed.tsv


# Next, we wish to restrict to the subset of samples from our QC analysis.
#gsutil cp gs://ukb31063-mega-gwas/qc/ukb31063.keep_samples.txt $HOME
#I am not starting off restricting to anyone, include all samples. 
#there's also another raw pheno file in tsv format 
#gs://phenotype_31063/ukb31063.raw_phenotypes.tsv.bgz	Raw phenotype file provided by UK Biobank with application #31063 sample IDs; converted from csv to tsv and quotation marks surrounding fields removed (502,616 samples)
#just do the same pipeline as what Duncan did to make sure that it works on the raw csv from ukbb.

#Rscript restrict_to_QC_samples.r
#I will include everyone in this first pass, so don't need to restrict to QC samples yet...
#also needed to install the package optparseoptparse
Rscript -e 'install.packages(c("optparse"))'
export PATH="/home/eatkinso/R-3.4.1/bin:$PATH" 
screen #run in a screen session in case I lose internet at any point
#launched part 2 at 11:40. Part 1 took maybe 1-1.5 hours to run.

cd PHESANT/WAS
# Execute the PHESANT script on the parsed phenotype file.
NUMPARTS=4

#now just start at part 2 since part 1 finished before
for i in `seq 2 $NUMPARTS`;
do
    Rscript phenomeScan.r \
        --phenofile="../../neale_lab_parsed.tsv" \
        --variablelistfile="../variable-info/outcome_info_final_round2.tsv" \
        --datacodingfile="../variable-info/data-coding-ordinal-info.txt" \
        --userId="userId" \
        --resDir="../../" \
        --out="uk_round2_allSamples_phenos_phesant" \
        --partIdx="$i" \
        --numParts="${NUMPARTS}"
done


#at the end of part 4, threw one error message:
#Warning message:
In fread(opt$phenofile, header = TRUE, sep = "\t") :
  Found and resolved improper quoting out-of-sample. First healed line 29855: <<1298545 1       1       0       1       1957   3
7       44      96      105     189     144     3       2008-11-15      11007           -0.067                                 1
0                                                                                                                              5
5       58      9       387500  435500  2132043 0       0       1       0       -4.60767        3       4               2      3
        3       6               0       5               125     489             2       3       3       6       1       0      1
        3       6       7       2       5       0       3       3       3       9       0       4       3       9       8      2
        5       0       1       1       0       0       1       0       1       0       0       1       1               469    4
53                      391             437                     391     344     1       2       12      4       3       4      1
4       38      5       6       1       1       1               5       10      0               5       45      3       2      1
        2       5       3                       3       1       -10     -10     -10     1       2       4       2       1      0
        2       7       3       -1      1       1       1       1       0       3       0       0       -1      3       1      -
10      0       1       1       3       2       2       1       1       4       2       0       6       3       7       4      1
        2       2       3       1       1       0       2       1       1       0       5       0       0       >>. If the field
s are not quoted (e.g. field separator does not appear within any field), try quote="" to avoid this warning.

#when done, copy over all the files produced to my google cloud bucket
gsutil cp uk_round2_allSamples_phenos_phesant.1.log gs://ukb-diverse-pops/Phenotypes/uk_round2_allSamples_phenos_phesant.1.log
gsutil cp uk_round2_allSamples_phenos_phesant.1.tsv gs://ukb-diverse-pops/Phenotypes/uk_round2_allSamples_phenos_phesant.1.tsv
gsutil cp uk_round2_allSamples_phenos_phesant.2.tsv gs://ukb-diverse-pops/Phenotypes/uk_round2_allSamples_phenos_phesant.2.tsv
gsutil cp uk_round2_allSamples_phenos_phesant.3.tsv gs://ukb-diverse-pops/Phenotypes/uk_round2_allSamples_phenos_phesant.3.tsv
gsutil cp uk_round2_allSamples_phenos_phesant.4.tsv gs://ukb-diverse-pops/Phenotypes/uk_round2_allSamples_phenos_phesant.4.tsv

#gsutil cp uk_round2_allSamples_phenos_phesant.* gs://ukb-diverse-pops/Phenotypes/. 
#that put them all into a folder called ./; just specify the filenames for each for now
#gsutil cp neale_lab_parsed.tsv gs://ukb-diverse-pops/Phenotypes/neale_lab_parsed_allSamples.tsv
 