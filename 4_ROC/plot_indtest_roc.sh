#!/bin/bash
source ~/.bashrc; conda activate R4_3; Rscript /home/nrlab/wang04/ulyses/4_ROC/indtest_roc.R --tf_min 0 --tf_max 0.03

source ~/.bashrc;conda activate R4_3 ; Rscript /home/nrlab/wang04/ulyses/4_ROC/indtest_roc.R --tf_min 0.03 --tf_max 0.1

source ~/.bashrc;conda activate R4_3 ; Rscript /home/nrlab/wang04/ulyses/4_ROC/indtest_roc.R --tf_min 0.1 --tf_max 1

source ~/.bashrc;conda activate R4_3 ; Rscript /home/nrlab/wang04/ulyses/4_ROC/indtest_roc.R --tf_min 0 --tf_max 1

