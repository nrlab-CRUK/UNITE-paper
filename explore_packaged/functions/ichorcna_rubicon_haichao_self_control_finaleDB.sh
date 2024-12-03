#!/bin/bash
#activate the conda ichorcna environment
source /home/${USER}/.bashrc
conda activate ichorcna

BAM=$1
wigfile=$2
outdir=$3
bam_basename=$4

#run readCounter
readCounter --build ${BAM}
readCounter  --window 1000000 \
-c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22  \
${BAM}  > ${wigfile}

#run ichorCNA
Rscript /home/nrlab/tools/anaconda3/envs/ichorcna/ichorCNA/scripts/runIchorCNA.R --id ${bam_basename}.ichor \
  --WIG ${wigfile} --ploidy "c(2,3)" --normal "c(0.95, 0.99, 0.995, 0.999)" --maxCN 5 \
  --gcWig /home/nrlab/tools/anaconda3/envs/ichorcna/ichorCNA/inst/extdata/gc_hg19_1000kb.wig \
  --mapWig /home/nrlab/tools/anaconda3/envs/ichorcna/ichorCNA/inst/extdata/map_hg19_1000kb.wig \
  --centromere /home/nrlab/tools/anaconda3/envs/ichorcna/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
  --includeHOMD False --chrs "c(1:22)" --chrTrain "c(1:22)" \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
  --scStates "c(1, 3)" --txnE 0.9999 --txnStrength 10000 --outDir ${outdir} 

