
## download, trim and align using TAP pipeline

```bash
spack load openjdk@17

 /home/bioinformatics/software/pipelines/kickstart/current/bin/kickstart   -f "Library Type" -f "Index Type" -l=SLX-11034 -l=SLX-11868

```

## calculate bam stats

```
cd /mnt/scratchb/nrlab/wang04/csf/processed

for i in *.bam

do

sbatch --time=0-10 --mem=128G  -J "count" -p "general" --wrap="source ~/.bashrc; conda activate R4_1; Rscript --vanilla /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/summariseBam.R  ${i}  .  FALSE"

done

```

## downsample to 0.1x (1Million fragments)

```bash

cd /mnt/scratchb/nrlab/wang04/csf/processed
mkdir -p /mnt/scratchb/nrlab/wang04/csf/downsampled_bam
downsampled_dir=/mnt/scratchb/nrlab/wang04/csf/downsampled_bam

for i in *.csv

do

sbatch --time=0-20 --mem=32G  -J "downsamp" -p "general" --wrap="source ~/.bashrc; conda activate R4_1; Rscript --vanilla /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/downsample_gatk.R ${i}  0.1  ${downsampled_dir}  /home/nrlab/tools/anaconda3/envs/projectx/bin/gatk"

done



```


## run ulyses: generate ulyses_image.rds files

```bash
downsampled_dir=/scratchc/nrlab/wang04/ulyses/data/internal/csf/downsampled_bam

sbatch --time=0-00:30:00  --job-name='u-w-igniter' --partition="epyc" --mem=1G --wrap="source ~/.bashrc; conda activate R4_1; Rscript /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/sbatch_ulyses.R  ${downsampled_dir}"
```

## ulyses plot

```bash
cd /scratchc/nrlab/wang04/ulyses/data/internal/csf/downsampled_bam

conda activate R4_2

for i in *.ulyses_image.rds
do 
Rscript /home/nrlab/wang04/ulyses/explore_packaged/functions/ulyses_image_to_features.R   ${i}
done

```







