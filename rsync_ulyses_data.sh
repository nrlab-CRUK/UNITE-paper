#!/bin/bash
# login to clust1-headnode.cri.camres.org
ssh wang04@clust1-sub.cri.camres.org



# stm
stm_from=/scratchc/nrlab/wang04/ulyses/data/stm_paper_data/downsampled_bam/
stm_to=/mnt/nas-data/nrlab/group_folders/wang04/ulyses_datasets/stm_paper_data/downsampled_bam/
screen -S rsync_stm -d -m
screen -S rsync_stm -X stuff "rsync -av --progress --stats ${stm_from} ${stm_to}^M"
screen -S rsync_stm -X detach


# liquorice
liquorice_from=/scratchc/nrlab/wang04/ulyses/data/liquorice/downsampled_bam/
liquorice_to=/mnt/nas-data/nrlab/group_folders/wang04/ulyses_datasets/liquorice/downsampled_bam/
screen -S rsync_liquorice -d -m
screen -S rsync_liquorice -X stuff "rsync -av --progress --stats ${liquorice_from} ${liquorice_to}^M"
screen -S rsync_liquorice -X detach

# finaleDB
finaleDB_from=/scratchc/nrlab/wang04/ulyses/data/finaledb/hg19_downsampled_rds/
finaleDB_to=/mnt/nas-data/nrlab/group_folders/wang04/ulyses_datasets/finaledb/hg19_downsampled_rds/
screen -S rsync_finaledb -d -m
screen -S rsync_finaledb -X stuff "rsync -av --progress --stats ${finaleDB_from} ${finaleDB_to}^M"
screen -S rsync_finaledb -X detach

# delfi
delfi_from=/scratchc/nrlab/wang04/ulyses/data/delfi/downsampled_bams/
delfi_to=/mnt/nas-data/nrlab/group_folders/wang04/ulyses_datasets/delfi/downsampled_bams/
screen -S rsync_delfi -d -m
screen -S rsync_delfi -X stuff "rsync -av --progress --stats ${delfi_from} ${delfi_to}^M"
screen -S rsync_delfi -X detach


# urine
urine_from=/scratchc/nrlab/wang04/ulyses/data/urine/downsampled_bam/
urine_to=/mnt/nas-data/nrlab1/group_folders/wang04/ulyses_datasets/urine/downsampled_bam/
screen -S rsync_urine -d -m
screen -S rsync_urine -X stuff "rsync -av --progress --stats ${urine_from} ${urine_to}^M"
screen -S rsync_urine -X detach

# dbs
dbs_from=/scratchc/nrlab/wang04/ulyses/data/dbs/ulyses_analysis/dbs/downsampled_bam/
dbs_to=/mnt/nas-data/nrlab1/group_folders/wang04/ulyses_datasets/dbs/ulyses_analysis/dbs/downsampled_bam/
screen -S rsync_dbs -d -m
screen -S rsync_dbs -X stuff "rsync -av --progress --stats ${dbs_from} ${dbs_to}^M"
screen -S rsync_dbs -X detach


# csf
csf_from=/scratchc/nrlab/wang04/ulyses/data/csf/ulyses_analysis/csf/downsampled_bam/
csf_to=/mnt/nas-data/nrlab1/group_folders/wang04/ulyses_datasets/csf/ulyses_analysis/csf/downsampled_bam/
screen -S rsync_csf -d -m
screen -S rsync_csf -X stuff "rsync -av --progress --stats ${csf_from} ${csf_to}^M"
screen -S rsync_csf -X detach


################################################################################
# second batch
################################################################################

# pbcp 
pbcp_from=/scratchc/nrlab/wang04/ulyses_second_batch/data/pbcp/downsampled_bam/
pbcp_to=/mnt/nas-data/nrlab1/group_folders/wang04/ulyses_datasets/plasma/pbcp/downsampled_bam/
screen -S rsync_pbcp -d -m
screen -S rsync_pbcp -X stuff "rsync -av --progress --stats ${pbcp_from} ${pbcp_to}^M"
screen -S rsync_pbcp -X detach


# danlandau

dan_from=/scratchc/nrlab/wang04/ulyses_second_batch/data/dan_landau/downsampled_bam/
dan_to=/mnt/nas-data/nrlab1/group_folders/wang04/ulyses_datasets/plasma/dan_landau/downsampled_bam/
screen -S rsync_dan -d -m
screen -S rsync_dan -X stuff "rsync -av --progress --stats ${dan_from} ${dan_to}^M"
screen -S rsync_dan -X detach

# lucid

lucid_from=/scratchc/nrlab/wang04/ulyses_second_batch/data/lucid/downsampled_bam/
lucid_to=/mnt/nas-data/nrlab/group_folders/wang04/ulyses_datasets/lucid/downsampled_bam/
screen -S rsync_lucid -d -m
screen -S rsync_lucid -X stuff "rsync -av --progress --stats ${lucid_from} ${lucid_to}^M"
screen -S rsync_lucid -X detach

# oee

oee_from=/scratchc/nrlab/wang04/ulyses_second_batch/data/oee/downsampled_bam/
oee_to=/mnt/nas-data/nrlab/group_folders/wang04/ulyses_datasets/oee/downsampled_bam/
screen -S rsync_oee -d -m
screen -S rsync_oee -X stuff "rsync -av --progress --stats ${oee_from} ${oee_to}^M"
screen -S rsync_oee -X detach