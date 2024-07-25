#ASTRAL_signal_high
echo "iterate,first_calc,pruning,all" > ASTRAL_greedy_pruning_signal_high.csv
for i in {1..100}; do cat signal_high/ASTRAL_greedy_pruning.o26284632.$i|grep '^TIME'| tr '　' ' '|awk '{print $2}'|xargs|tr ' ' ',' ;done >> ASTRAL_greedy_pruning_signal_high.csv
awk -F, 'NR>1 {for (i = 1; i <= NF; i++) {sum[i] += $i;} count++;} END {for (i = 1; i <= NF; i++) {printf "%.2f", sum[i] / count; if (i < NF) {printf ",";}} print "";}' ASTRAL_greedy_pruning_signal_high.csv
#ASTRAL_signal_low
echo "iterate,first_calc,pruning,all" > ASTRAL_greedy_pruning_signal_low.csv
for i in {1..100}; do cat signal_low/ASTRAL_greedy_pruning.o*.$i|grep '^TIME'| tr '　' ' '|awk '{print $2}'|xargs|tr ' ' ',' ;done >> ASTRAL_greedy_pruning_signal_low.csv
awk -F, 'NR>1 {for (i = 1; i <= NF; i++) {sum[i] += $i;} count++;} END {for (i = 1; i <= NF; i++) {printf "%.2f", sum[i] / count; if (i < NF) {printf ",";}} print "";} ASTRAL_greedy_pruning_signal_low.csv
