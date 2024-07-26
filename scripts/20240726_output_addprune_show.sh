#ASTRAL_signal_high
echo "Creating_DIST,Creating_TD,iterate_count,add-prune,all" > majority_greedy_addprune_signal_high.csv
for i in {1..100}; do cat addprune_signal_high/majority_greedy_addprune.o*.$i|grep '^TIME'| tr '　' ' '|awk '{print $2}'|xargs|tr ' ' ',' ;done >> majority_greedy_addprune_signal_high.csv
awk -F, 'NR>1 {for (i = 1; i <= NF; i++) {sum[i] += $i;} count++;} END {for (i = 1; i <= NF; i++) {printf "%.2f", sum[i] / count; if (i < NF) {printf ",";}} print "";}' majority_greedy_addprune_signal_high.csv
#ASTRAL_signal_low
echo "Creating_DIST,Creating_TD,iterate_count,add-prune,all" > majority_greedy_addprune_signal_low.csv
for i in {1..100}; do cat addprune_signal_low/majority_greedy_addprune.o*.$i|grep '^TIME'| tr '　' ' '|awk '{print $2}'|xargs|tr ' ' ',' ;done >> majority_greedy_addprune_signal_low.csv
awk -F, 'NR>1 {for (i = 1; i <= NF; i++) {sum[i] += $i;} count++;} END {for (i = 1; i <= NF; i++) {printf "%.2f", sum[i] / count; if (i < NF) {printf ",";}} print "";}' majority_greedy_addprune_signal_low.csv
