

python3 "/Users/wang04/.nextflow/assets/nrlab-CRUK/NRLAB_TAP/python/FastqSplit.py" \
    --source="/mnt/scratchb/PGDX5883P1_WGS.sorted_processed.namesort.r1.fastq.gz" \
    --prefix="DL-00001.DL421-DL959.H5WNTBBXX.8.r_1" \
    --reads=10000000





DL-00001.DL421-DL959.H5WNTBBXX.8.r_1-C000000.fq.gz
DL-00001.DL421-DL959.H5WNTBBXX.8.r_2-C000000.fq.gz

DL-00001.DL458-DL996.H5WNTBBXX.8.r_2-C000000.fq.gz
DL-00001.DL458-DL996.H5WNTBBXX.8.r_1-C000000.fq.gz




---
so it seems these two files are broken:

[wang04@clust1-headnode bam2fq]$ zcat PGDX6865P1_WGS.sorted_processed.namesort.r1.fastq.gz | wc -l

gzip: PGDX6865P1_WGS.sorted_processed.namesort.r1.fastq.gz: unexpected end of file
138643325





PGDX5883P1_WGS.sorted_processed.namesort.r1.fastq.gz




----


/scratchb/nrlab/wang04/external_dataset_download/delfi/delfi_ega_data/bam_debug






for i in *.bam

do

sbatch --time=0-4 --mem=4G  -J "count" -p "epyc" --wrap="source ~/.bashrc; conda activate R4_1; Rscript --vanilla /scratchc/nrlab/wang04/ulyses/scripts/summariseBam.R ${i}  .  FALSE"

done



sbatch --time=5-0 --mem=16G  -J "" -p "epyc" --wrap="source ~/.bashrc; conda activate R4_1; Rscript --vanilla /scratchc/nrlab/wang04/ulyses/scripts/finaledb_download.R "


