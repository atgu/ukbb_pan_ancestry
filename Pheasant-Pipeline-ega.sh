#!/bin/bash

######  Running PHESANT on a Google Could VM
# started with 4 vCPUs w 32gb mem each
# Boot disk - change default disk size from 10GB to 200GB. (Debian GNU/Linux).
#start off with not much memory while getting things set up, then scale up for the actual phenotyping.
# load up a terminal window, under 'Remote access' click ssh and then 'open in browser window'

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
Rscript -e 'install.packages(c("optparse"))'
export PATH="/home/eatkinso/R-3.4.1/bin:$PATH" 

# Install git.
sudo apt-get install git
cd $HOME

# clone PHESANT library
git clone https://github.com/astheeggeggs/PHESANT.git

# Move the required reengineering_phenofile.r over and the phenotypes
gsutil cp gs://phenotype_31063/ukb11214.csv gs://ukb-diverse-pops/Phenotypes/ukb11214.csv

#copy these all into the VM also from my Google Cloud Phenotypes directory
gsutil cp gs://ukb-diverse-pops/Phenotypes/reengineering_phenofile_neale_lab_ega2.r $HOME
gsutil cp gs://ukb-diverse-pops/Phenotypes/ukb11214.csv $HOME
#needed to crank the boot disk size up to 400GB to get this to work. Need 25 for this file, another 25 for the output fil, and then R will need maybe 10fold that to be able to run
#also made a tmp directory to use for the transfer 

# Run the reegineering_phenofile.r on the raw phenotype data.
Rscript reengineering_phenofile_neale_lab_ega.r 

# This will write a .tsv file that we then execute a PHESANT script on.
#copy over the new formatted pheno file to my google cloud bucket
gsutil cp neale_lab_parsed.tsv gs://ukb-diverse-pops/Phenotypes/neale_lab_parsed.tsv

# Next, we wish to restrict to the subset of samples from our QC analysis.
#gsutil cp gs://ukb31063-mega-gwas/qc/ukb31063.keep_samples.txt $HOME
#Rscript restrict_to_QC_samples.r

screen   #run in a screen session in case I lose internet at any point
cd PHESANT/WAS
# Execute the PHESANT script on the parsed phenotype file.
NUMPARTS=4

for i in `seq 1 $NUMPARTS`;
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


#when done, copy over all the files produced to my google cloud bucket
gsutil cp uk_round2_allSamples_phenos_phesant.* gs://ukb-diverse-pops/Phenotypes/. 
#gsutil cp neale_lab_parsed.tsv gs://ukb-diverse-pops/Phenotypes/neale_lab_parsed_allSamples.tsv
 
