# Workflow

This working flow runs on `epyc` node on `clust1-sub`

# Call the function
check_install_packages(packages)
```

# UNITE workflow


## index bam files (when needed)
```bash
for i in *.bam; do sbatch --time=0-1:00:00 --mem=4G --partition=epyc -J index --wrap="source ~/.bashrc; conda activate projectx; samtools index ${i}"; done
```

## calculate bam stats (e.g., number of reads, and coverage.)
```bash
for i in *.bam
do
sbatch --time=0-10 --mem=128G  -J "count" -p "epyc" --wrap="source ~/.bashrc; conda activate R4_2; Rscript --vanilla /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/summariseBam.R  ${i}  .  FALSE"
```


## downsample the bam file to target_depth (modify the variables below when needed) 
```bash
downsampled_dir="../downsampled_bam"
#the unit of target_depth is "x"
target_depth=0.1
mkdir -p ${downsampled_dir}
for i in *.bam.summary.csv
do
sbatch --time=0-20 --mem=32G  -J "downsamp" -p "epyc" --wrap="source ~/.bashrc; conda activate R4_1; Rscript --vanilla /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/downsample_gatk.R ${i}  ${target_depth}  ${downsampled_dir}  /home/nrlab/tools/anaconda3/envs/projectx/bin/gatk"
done
```

## run UNITE
```bash
working_dir=""
/home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/run_ulyses.sh ${working_dir}
```

## UNITE plots

In last step, for each bam file, there will be "*ulyses_image.rds" files generated.
The code below will generate corresponding plots for each rds file generated.
```bash
for i in *.ulyses_image.rds
do 
sbatch --time=0-00:15:00 \
  --mem=1G --partition=epyc \
  -J "plot" \
  --wrap="source ~/.bashrc; \
    conda activate R4_2; \
    Rscript /home/nrlab/wang04/ulyses/explore_packaged/functions/ulyses_image_to_features.R \
      --image_file $i \
      --specimen plasma \
      --isize_filter_min 100 \
      --isize_filter_max 220 "
done
```

## Quality Control the UNITE results
This step checks if the result files for each samples are successfully generated.
Enter the directory containing the UNITE result files, then run the code below:
```bash
sbatch --time=1-0 --mem=8G -J QC --partition=epyc --wrap="bash /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/QC_ulyses_results.sh"
```

#################################################################################
#archived codes below
#################################################################################


## Check and build R env

```R
check_install_packages <- function(packages) {
  missing_packages <- setdiff(packages, installed.packages()[, "Package"])
  
  if (length(missing_packages) > 0) {
    message("Installing missing packages:")
    message(missing_packages)
  
    # Install missing packages
    install.packages(missing_packages, dependencies = TRUE)
  
    message("*******************************************")
    message("Installed missing packages.")
  } else {
    message("All required packages are already installed.")
  }
}

# List of packages to check and install if missing
packages <- c(
  "parallel", "future", "furrr", "progressr", "parallelly",
  "tidyverse", "Biobase", "QDNAseq", "Homo.sapiens", "plyranges",
  "GenomicRanges", "GenomicAlignments", "matrixStats", "AnnotationHub",
  "BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Hsapiens.UCSC.hg38",
  "Rsamtools", "rslurm", "abind", "tictoc", "patchwork", "EBImage",
  "gridExtra", "grid", "reshape2", "purrr", "fastDummies", "ggalluvial",
  "MultiAssayExperiment", "nord", "plotly", "htmlwidgets", "ggpubr"
)



## Demisfying ref genomes

https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use

https://gatk.broadinstitute.org/hc/en-us/articles/360035890711-GRCh37-hg19-b37-humanG1Kv37-Human-Reference-Discrepancies

http://assets.thermofisher.com/TFS-Assets/LSG/Vector-Information/human-genome-38-faq.pdf

## finaleDB tsv to GRanges.rds

```bash


#------------------------------------------------------------------------------
# hg19
#------------------------------------------------------------------------------



for i in $(pwd)/*hg19.frag.tsv.bgz; do if [ ! -f ${i}.GRanges.rds ]; then echo ${i}; fi; done  |\
 while read -r filename || [ -n "$filename" ]; do sbatch --time=0-5 --mem=64G  -J "bgz-rds" -p "general" --wrap="source ~/.bashrc; conda activate R4_1; Rscript --vanilla /home/nrlab/wang04/ulyses/finaleDB/read_finaledb_frag_tsv_to_granges.r  ${filename}"; done

for i in $(pwd)/*hg19.frag.tsv.bgz; do if [ ! -f ${i}.GRanges.rds ]; then echo ${i}; fi; done  |\
 while read -r filename || [ -n "$filename" ]; do sbatch --time=0-3 --mem=128G  -J "bgz-rds" -p "general" --wrap="source ~/.bashrc; conda activate R4_1; Rscript --vanilla /home/nrlab/wang04/ulyses/finaleDB/read_finaledb_frag_tsv_to_granges.r  ${filename}"; done


# following samples :

/scratchb/nrlab/wang04/finaledb/hg19/EE86274.hg19.frag.tsv.bgz
/scratchb/nrlab/wang04/finaledb/hg19/EE86859.hg19.frag.tsv.bgz
/scratchb/nrlab/wang04/finaledb/hg19/EE86965.hg19.frag.tsv.bgz
/scratchb/nrlab/wang04/finaledb/hg19/EE87484.hg19.frag.tsv.bgz
/scratchb/nrlab/wang04/finaledb/hg19/EE87560.hg19.frag.tsv.bgz
/scratchb/nrlab/wang04/finaledb/hg19/EE87749.hg19.frag.tsv.bgz

# * ? some samples don't have chrMT. to check why.

#------------------------------------------------------------------------------
# hg38
#------------------------------------------------------------------------------


for i in $(pwd)/*hg38.frag.tsv.bgz; do if [ ! -f ${i}.GRanges.rds ]; then echo ${i}; fi; done |\
 while read -r filename || [ -n "$filename" ]; do sbatch --time=0-5 --mem=64G  -J "bgz-rds" -p "general" --wrap="source ~/.bashrc; conda activate R4_1; Rscript --vanilla /home/nrlab/wang04/ulyses/finaleDB/read_finaledb_frag_tsv_to_granges.r  ${filename}"; done

# following samples failed:

/scratchb/nrlab/wang04/finaledb/hg38/EE86965.hg38.frag.tsv.bgz
/scratchb/nrlab/wang04/finaledb/hg38/EE87484.hg38.frag.tsv.bgz
/scratchb/nrlab/wang04/finaledb/hg38/EE87560.hg38.frag.tsv.bgz
/scratchb/nrlab/wang04/finaledb/hg38/EE87749.hg38.frag.tsv.bgz

```


## Downsample finaleDB rds files

```bash

sbatch --time=0-1  --job-name='rds_ds_igniter' --partition="general" --mem=1G --wrap="source ~/.bashrc; conda activate R4_1; Rscript /home/nrlab/wang04/ulyses/finaleDB/finaleDB_rds_downsample.r  /scratchb/nrlab/wang04/finaledb/hg19  1 /scratchb/nrlab/wang04/finaledb/hg19_downsampled_rds general"



sbatch --time=0-1  --job-name='rds_ds_igniter' --partition="general" --mem=1G --wrap="source ~/.bashrc; conda activate R4_1; Rscript /home/nrlab/wang04/ulyses/finaleDB/finaleDB_rds_downsample.r  /scratchb/nrlab/wang04/finaledb/hg38  1 /scratchb/nrlab/wang04/finaledb/hg38_downsampled_rds general"


```

## Summarize bam stats

```shell


for i in *.bam

do

sbatch --time=0-4 --mem=4G  -J "count" -p "epyc" --wrap="source ~/.bashrc; conda activate R4_1; Rscript --vanilla /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/summariseBam.R  ${i}  .  FALSE"

done

```

## Based on the stats, downsample the bam to targeted depth, e.g. 0.1x

### stm

```shell

for i in *.csv

do

sbatch --time=0-20 --mem=32G  -J "downsamp" -p "epyc" --wrap="source ~/.bashrc; conda activate R4_1; Rscript --vanilla /home/nrlab/wang04/ulyses/scripts/downsample_gatk.R ${i}  0.1  /scratchc/nrlab/wang04/ulyses/data/internal/stm/downsampled_bam  /home/nrlab/tools/anaconda3/envs/projectx/bin/gatk"

done

```

### liquorice

```shell
for i in *.csv

do

sbatch --time=0-20 --mem=128G  -J "downsamp" -p "epyc" --wrap="source ~/.bashrc; conda activate R4_1; Rscript --vanilla /home/nrlab/wang04/ulyses/scripts/downsample_gatk.R ${i}  0.1  /scratchc/nrlab/wang04/ulyses/data/external/liquorice/downsampled_bam  /home/nrlab/tools/anaconda3/envs/projectx/bin/gatk"

done



```

### delfi

```shell

for i in *.csv

do

sbatch --time=0-20 --mem=32G  -J "downsamp" -p "general" --wrap="source ~/.bashrc; conda activate R4_1; Rscript --vanilla /home/nrlab/wang04/ulyses/scripts/downsample_gatk.R ${i}  0.1  /scratchb/nrlab/wang04/4_delfi/downsampled_bam  /home/nrlab/tools/anaconda3/envs/projectx/bin/gatk"

done

```

## Run Ulyses

### liquorice 0.1x

```shell

# the bam file need to be fullpath filename

WRAPPER_FILE="/home/nrlab/wang04/ulyses/explore_packaged/functions/wrapper_func2.r"


for i in $(pwd)/*DL{001..010}-DL*.bam
do 
sbatch --time=3-0 --job-name='uys' --partition="epyc" --mem=64G --wrap="source ~/.bashrc; conda activate R4_1; Rscript ${WRAPPER_FILE}  ${i} ${i%.bam}.bai 1000000  epyc  /home/nrlab/wang04/ulyses/explore_packaged"
done


for i in $(pwd)/*DL{011..030}-DL*.bam
do 
sbatch --time=3-0 --job-name='uys' --partition="epyc" --mem=64G --wrap="source ~/.bashrc; conda activate R4_1; Rscript ${WRAPPER_FILE}  ${i} ${i%.bam}.bai 1000000  epyc  /home/nrlab/wang04/ulyses/explore_packaged"
done


for i in $(pwd)/*DL{031..060}-DL*.bam
do 
sbatch --time=3-0 --job-name='uys' --partition="epyc" --mem=64G --wrap="source ~/.bashrc; conda activate R4_1; Rscript ${WRAPPER_FILE}  ${i} ${i%.bam}.bai 1000000  epyc  /home/nrlab/wang04/ulyses/explore_packaged"
done


for i in $(pwd)/*DL{061..100}-DL*.bam
do 
sbatch --time=3-0 --job-name='uys' --partition="epyc" --mem=64G --wrap="source ~/.bashrc; conda activate R4_1; Rscript ${WRAPPER_FILE}  ${i} ${i%.bam}.bai 1000000  epyc  /home/nrlab/wang04/ulyses/explore_packaged"
done


for i in $(pwd)/*DL{101..140}-DL*.bam
do 
sbatch --time=3-0 --job-name='uys' --partition="epyc" --mem=64G --wrap="source ~/.bashrc; conda activate R4_1; Rscript ${WRAPPER_FILE}  ${i} ${i%.bam}.bai 1000000  epyc  /home/nrlab/wang04/ulyses/explore_packaged"
done



for i in $(pwd)/*DL{141..180}-DL*.bam
do 
sbatch --time=3-0 --job-name='uys' --partition="epyc" --mem=64G --wrap="source ~/.bashrc; conda activate R4_1; Rscript /home/nrlab/wang04/ulyses/explore_packaged/functions/wrapper_func2.r  ${i} ${i%.bam}.bai 1000000  epyc  /home/nrlab/wang04/ulyses/explore_packaged"
done



for i in $(pwd)/*DL{181..220}-DL*.bam
do 
sbatch --time=3-0 --job-name='uys' --partition="epyc" --mem=64G --wrap="source ~/.bashrc; conda activate R4_1; Rscript /home/nrlab/wang04/ulyses/explore_packaged/functions/wrapper_func2.r  ${i} ${i%.bam}.bai 1000000  epyc  /home/nrlab/wang04/ulyses/explore_packaged"
done

for i in $(pwd)/*DL{221..263}-DL*.bam
do 
sbatch --time=3-0 --job-name='uys' --partition="epyc" --mem=64G --wrap="source ~/.bashrc; conda activate R4_1; Rscript /home/nrlab/wang04/ulyses/explore_packaged/functions/wrapper_func2.r  ${i} ${i%.bam}.bai 1000000  epyc  /home/nrlab/wang04/ulyses/explore_packaged"
done




```

### delfi 0.1x

```shell

sbatch --time=0-00:30:00  --job-name='u-w-igniter' --partition="epyc" --mem=1G --wrap="source ~/.bashrc; conda activate R4_1; Rscript /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/sbatch_ulyses.R  /scratchc/nrlab/wang04/ulyses/data/external/delfi/downsampled_bams"





```

### stm 0.1x

```shell


sbatch --time=0-00:30:00  --job-name='u-w-igniter' --partition="epyc" --mem=1G --wrap="source ~/.bashrc; conda activate R4_1; Rscript /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/sbatch_ulyses.R  /scratchc/nrlab/wang04/ulyses/data/internal/stm/downsampled_bam"



```

### finaleDB 1M

```shell

sbatch --time=0-00:30:00  --job-name='u-w-igniter' --partition="epyc" --mem=1G --wrap="source ~/.bashrc; conda activate R4_1; Rscript /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/sbatch_ulyses_rds_as_input.R  /scratchc/nrlab/wang04/ulyses/data/external/finaledb/hg19_downsampled_rds"






```

## Ulyses plot

### finaleDB

```bash

sintr --time=3-0 --mem=8G --partition="epyc"

conda activate R4_2

cd /scratchc/nrlab/wang04/ulyses/data/external/finaledb/hg19_downsampled_rds

for i in *.ulyses_image.rds; do Rscript /home/nrlab/wang04/ulyses/explore_packaged/functions/ulyses_image_to_features.R   $i; done

```

### liquorice

```bash
cd /scratchc/nrlab/wang04/ulyses/data/external/liquorice/downsampled_bam

for i in *.ulyses_image.rds; do Rscript /home/nrlab/wang04/ulyses/explore_packaged/functions/ulyses_image_to_features.R   $i; done

```

### delfi

```bash
cd /scratchc/nrlab/wang04/ulyses/data/external/delfi/downsampled_bams
for i in *.ulyses_image.rds; do Rscript /home/nrlab/wang04/ulyses/explore_packaged/functions/ulyses_image_to_features.R   $i; done


```

### stm

```bash
cd /scratchc/nrlab/wang04/ulyses/data/internal/stm/downsampled_bams
for i in *.ulyses_image.rds; do Rscript /home/nrlab/wang04/ulyses/explore_packaged/functions/ulyses_image_to_features.R   $i; done


```

## Run ichorCNA TF%

```R

# convert finaleDB RDS to bam file 

# for hg19

for i in /scratchc/nrlab/wang04/ulyses/data/external/finaledb/hg19_downsampled_rds/*.hg19.frag.tsv.bgz.GRanges.rds.1M.rds
do
sbatch --time=0-0:30:00 --mem=1G --partition="epyc" -J "bam" --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/galp_to_bam.R ${i} /scratchc/nrlab/wang04/ulyses/data/external/finaledb/hg19_downsampled_rds FALSE"
done

# for hg38

for i in /scratchc/nrlab/wang04/ulyses/data/external/finaledb/hg38_downsampled_rds/*.hg38.frag.tsv.bgz.GRanges.rds.1M.rds
do
sbatch --time=0-0:30:00 --mem=1G --partition="epyc" -J "bam" --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/galp_to_bam.R ${i} /scratchc/nrlab/wang04/ulyses/data/external/finaledb/hg38_downsampled_rds FALSE"
done


# run ichorCNA on all 0.1x bams
# be careful finaleDB is single end while internal and ega data is paired end bams

run `/home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/sbatch_ichorCNA.R` mannualy.

```

## Build MAE

```R

cd /scratchc/nrlab/wang04/ulyses/data

sintr --time=20-0 --mem=300G -J 'all_mae' --partition="epyc"

conda activate R4_2

Rscript /home/nrlab/wang04/ulyses/multiassayexperiment/ulyses_mae.R




```
