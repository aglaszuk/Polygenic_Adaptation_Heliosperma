#!/bin/bash

# Extract top10 model estimates from any directory for both the 2-origins and the 1-origin models

for dir in *or* 
 do
  cd $dir/HELIO
  cat run*/best_* >> best_likelihoods.txt
  cut -f 21 best_likelihoods.txt | grep -v "MaxEstLhood" | sort | head -n 1
  cut -f 21 best_likelihoods.txt | grep -v "MaxEstLhood" | sort | awk '{sum+=$1} END {print "Average: " sum/ NR}'
  cut -f 21 best_likelihoods.txt | grep -v "MaxEstLhood" | sort | head -n 10 | sed 's/-//g' > top10_MaxEstL.txt
  head -n 1 best_likelihoods.txt > top10_results.txt
   while read line
    do
     grep $line best_likelihoods.txt >> top10_results.txt
    done < top10_MaxEstL.txt 
  rm top10_MaxEstL.txt
  cd ../..
done 
