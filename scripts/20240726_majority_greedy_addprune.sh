#!/bin/bash
#$ -S /bin/bash
#$ -l intel
#$ -l s_vmem=12G
#$ -l mem_req=12G
#$ -l d_rt=192:00:00
#$ -l s_rt=192:00:00
#$ -t 1-100:1
signal=low
#$ -e ./addprune_signal_low/
#$ -o ./addprune_signal_low/
#$ -cwd
#$ -N majority_greedy_addprune

function jobfinish_kubo(){
curl -X POST -H "Content-type: application/json" --data "{\"text\":\"<@UNHH
M0F6U> from nig ID: ${1} - ${2}\"}" https://hooks.slack.com/services/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
}

source ~/.bashrc
source ../ENV/bin/activate
export PYTHONPATH=/home/atsu1217/ConsensusProj/src
num=$SGE_TASK_ID
python3 20240726_transfer_greedy_addprune.py $num --signal $signal
deactivate

jobfinish_kubo majority_low ${signal}_$num
