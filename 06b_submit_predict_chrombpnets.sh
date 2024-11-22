#!/bin/bash

InputIDs=("AF04" "AF06" "AF08" "AF10" "AF12" "AF14" "AF16" "AF18" "AF20" "AF22" "AF24" "AF26" "AF28" "AF30" "AF34" "AF36" "AF38" "EU03" "EU05" "EU07" "EU09" "EU13" "EU15" "EU17" "EU19" "EU21" "EU25" "EU27" "EU29" "EU33" "EU37" "EU39" "EU41" "EU43" "EU47")

# sbatch file
SBATCH_FILE="/project/lbarreiro/USERS/daniel/deeplearn_pred/github/06a_predict_chrombpnet.sbatch"

# log file to keep track of submitted jobs
LOG_FILE="/project/lbarreiro/USERS/daniel/deeplearn_pred/scripts/submitted_jobs.log"

# function to count currently running jobs
count_running_jobs() {
    squeue -u $USER | grep " gpu " | wc -l
}

# initialize the log file if it doesn't exist
if [ ! -f $LOG_FILE ]; then
    touch $LOG_FILE
fi

# read the TSV file and process each indv
for ID in ${InputIDs[@]}; do
    # check if the indv has already been processed
    if grep -q "^$ID$" $LOG_FILE; then
        echo "Job for $ID has already been submitted. Skipping..."
        continue
    fi

    # check the number of running jobs and wait if it's 6 or more
    while [ $(count_running_jobs) -ge 6 ]; do
        echo "Waiting for jobs to finish."
        sleep 600  # wait for 10 minutes before checking again
    done

    # submit the job using batch
    sbatch $SBATCH_FILE $ID

    echo "Submitted job for: $ID"

    # log the submitted gene
    echo $ID >> $LOG_FILE
done
