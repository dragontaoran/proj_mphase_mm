#!/bin/bash

dir="res"

mkdir -p ${dir}

echo "#!/bin/bash" > tmp.slurm
echo "" >> tmp.slurm
echo "#SBATCH --mail-user=r.tao@vanderbilt.edu  # email address" >> tmp.slurm
echo "#SBATCH --mail-type=NONE  # Alerts sent when job begins, ends, or aborts" >> tmp.slurm
echo "#SBATCH --nodes=1   # Number of nodes required" >> tmp.slurm
echo "#SBATCH --ntasks=1   # Number of nodes required" >> tmp.slurm
echo "#SBATCH --mem=8G  # Total memory (RAM) required, per node" >> tmp.slurm
echo "#SBATCH --time=04-00:00:00  # Wall Clock time (dd-hh:mm:ss) [max of 14 days]" >> tmp.slurm
echo "#SBATCH --array=1-200" >> tmp.slurm
echo "#SBATCH --output=${dir}/%A_%a.slog  # output and error messages go to this file" >> tmp.slurm
echo "#SBATCH --job-name=aODS # job name" >> tmp.slurm
echo "#SBATCH --constraint=\"sandybridge|haswell|skylake\"" >> tmp.slurm
echo "" >> tmp.slurm
echo "echo \"SLURM_JOBID: \" \$SLURM_JOBID" >> tmp.slurm
echo "echo \"SLURM_ARRAY_TASK_ID: \" \$SLURM_ARRAY_TASK_ID" >> tmp.slurm
echo "echo \"SLURM_ARRAY_JOB_ID: \" \$SLURM_ARRAY_JOB_ID" >> tmp.slurm
echo "" >> tmp.slurm
echo "Rscript analysis.R \$SLURM_ARRAY_TASK_ID ${dir} ACML" >> tmp.slurm

# cat tmp.slurm
sbatch tmp.slurm
rm -rf tmp.slurm

