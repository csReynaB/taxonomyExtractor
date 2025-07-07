#!/bin/bash
#SBATCH --job-name=LINEAGES        # Job name
#SBATCH --output=/lisc/scratch/ccr/add_Taxonomy_oldTables/log_%A.out
#SBATCH --error=/lisc/scratch/ccr/add_Taxonomy_oldTables/log_%A.err
#SBATCH --time=30:00:00                # Maximum runtime (hh:mm:ss)
#SBATCH --ntasks=1                    # Number of tasks (1 task for a single script)
#SBATCH --nodes 1
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=1G                     # Memory allocation per node

module load conda
conda activate /lisc/scratch/ccr/conda_envs/taxonomics_env

# Arguments
metadata="$1"
out="$2"
organism_col="$3"
rep_id_col="$4"
prot_id_col="$5"
taxid_col="$6"
priority="$7"

# command
cmd="python3 extract_lineages.py --input \"$metadata\" --output \"$out\""

# Append arguments if given
[ -n "$organism_col" ] && cmd+=" --organism_col \"$organism_col\""
[ -n "$rep_id_col" ] && cmd+=" --rep_id_col \"$rep_id_col\""
[ -n "$prot_id_col" ] && cmd+=" --prot_id_col \"$prot_id_col\""
[ -n "$taxid_col" ] && cmd+=" --taxid_col \"$taxid_col\""
[ -n "$priority" ] && cmd+=" --method_priority \"$priority\""

# Echo for logging
echo "Running command:"
echo $cmd

# Run the command
eval $cmd
