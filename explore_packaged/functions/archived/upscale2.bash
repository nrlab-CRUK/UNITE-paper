
WRAPPER_FILE="/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_packaged/functions/wrapper_func2.r"

for i in *_filtered.bam
do 
sbatch --time=1-0 --job-name='c_ulyses' --mem=32G --wrap="source ~/.bashrc; conda activate R4_1; /home/nrlab/tools/anaconda3/envs/R4_1/bin/Rscript ${WRAPPER_FILE}  $(pwd)/${i} 1000000"
done



# try list(time = '10:00:00',  mem = '24G')

 SLX-13900.D708_D504.tagtrim.bwamem2_mrkdup_filtered.bam

i="SLX-13900_D708_D504.sorted_mrkdup_filtered.bam"
sbatch --time=1-0 --job-name='c_ulyses' --mem=32G --wrap="source ~/.bashrc; conda activate R4_1; /home/nrlab/tools/anaconda3/envs/R4_1/bin/Rscript ${WRAPPER_FILE}  $(pwd)/${i} 1000000"


sbatch --nodes=1 --output=/dev/null --time=00:01:00 --dependency=singleton --job-name="test" --wait --wrap="hostname"



#   slurm_options = list(time = '1-0',  mem = '80G'), ELAPSED 4:48:00 
 SLX-13900.D708_D504.tagtrim.bwamem2_mrkdup_filtered.bam

i="SLX-13900_D708_D504.sorted_mrkdup_filtered.bam"
sbatch --time=1-0 --job-name='c_ulyses' --mem=32G --wrap="source ~/.bashrc; conda activate R4_1; /home/nrlab/tools/anaconda3/envs/R4_1/bin/Rscript ${WRAPPER_FILE}  $(pwd)/${i} 1000000"




### try again slurm_options = list(time = '1-0',  mem = '80G') 
WRAPPER_FILE="/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_packaged/functions/wrapper_func2.r"

for i in *_filtered.bam
do 
sbatch --time=1-0 --job-name='d_ulyses' --mem=16G --wrap="source ~/.bashrc; conda activate R4_1; /home/nrlab/tools/anaconda3/envs/R4_1/bin/Rscript ${WRAPPER_FILE}  $(pwd)/${i} 1000000"
done



 #run a urine file 
WRAPPER_FILE="/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_packaged/functions/wrapper_func2.r"

for i in SLX-10991_iPCRtagT026_SLX-10340_iPCRtagT043.merge.bam.markduplicates.bam.downsamp_0.5x.mrkdup.bam
do 
sbatch --time=1-0 --job-name='d_ulyses' --mem=16G --wrap="source ~/.bashrc; conda activate R4_1; /home/nrlab/tools/anaconda3/envs/R4_1/bin/Rscript ${WRAPPER_FILE}  $(pwd)/${i} 1000000"
done


# check status 

cat /scratcha/nrlab/wang04/urine/delfi/down_sample/bams/0.5x/slurm-20858162.out 


# m0 m1 m2 in tagseq artefact

WRAPPER_FILE="/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_packaged/functions/wrapper_func2.r"

for i in *.bam
do 
sbatch --time=1-0 --job-name='d_ulyses' --mem=16G --wrap="source ~/.bashrc; conda activate R4_1; /home/nrlab/tools/anaconda3/envs/R4_1/bin/Rscript ${WRAPPER_FILE}  $(pwd)/${i} 1000000"
done


# run ulyses for the liquorice pilot samples

WRAPPER_FILE="/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_packaged/functions/wrapper_func2.r"

for i in *.t40.mean167.sd100_bwa1.bam.mrkdup.bam
do 
sbatch --time=2-0 --job-name='d_ulyses' --mem=32G --wrap="source ~/.bashrc; conda activate R4_1; /home/nrlab/tools/anaconda3/envs/R4_1/bin/Rscript ${WRAPPER_FILE}  $(pwd)/${i} 5000000"
done



# re-run ulyses control sample due to mem allocation errors
WRAPPER_FILE="/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_packaged/functions/wrapper_func2.r"

for i in liquorice.Ctrl_8.t40.mean167.sd100_bwa1.bam.mrkdup.bam
do 
sbatch --time=2-0 --job-name='d_ulyses' --mem=64G --wrap="source ~/.bashrc; conda activate R4_1; /home/nrlab/tools/anaconda3/envs/R4_1/bin/Rscript ${WRAPPER_FILE}  $(pwd)/${i} 5000000"
done


# run high TF sample,expect noisy CN profiles

WRAPPER_FILE="/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_packaged/functions/wrapper_func2.r"

for i in liquorice.EwS_7_1.t40.mean167.sd100_bwa1.bam.mrkdup.bam
do 
sbatch --time=2-0 --job-name='d_ulyses' --mem=64G --wrap="source ~/.bashrc; conda activate R4_1; /home/nrlab/tools/anaconda3/envs/R4_1/bin/Rscript ${WRAPPER_FILE}  $(pwd)/${i} 5000000"
done
