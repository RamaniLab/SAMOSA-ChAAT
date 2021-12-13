#!/bin/env bash
#
#$ -pe smp 1
#$ -l mem_free=3G
#$ -l h_rt=24:00:00

# >>> conda initialize >>>
# !! Contents within this block are managed by conda init !!
__conda_setup="$(/wynton/home/goodarzi/cpmcnally/bin/miniconda3/bin/conda shell.bash hook 2> /dev/null)"
if [ $? -eq 0 ]; then
        eval "$__conda_setup"
else
        if [ -f "/wynton/home/goodarzi/cpmcnally/bin/miniconda3/etc/profile.d/conda.sh" ]; then
                . "/wynton/home/goodarzi/cpmcnally/bin/miniconda3/etc/profile.d/conda.sh"
        else
                export PATH="/wynton/home/goodarzi/cpmcnally/bin/miniconda3/bin:$PATH"
        fi
fi
unset __conda_setup
# <<< conda initialize <<<

cd /wynton/group/goodarzilab/ramanilab/results/pacbio
conda activate procSAMOSA

python /wynton/home/goodarzi/cpmcnally/code/scripts/runAccessibilityHMM.py $1 $2
