#!/bin/bash
for p in 4 20
	do
    for n in 500 2000 5000
	do
	for alpha_id in {1..9}
		do
	    for gamma_id in {1..5}
		do 
		for seed in {1..100} 
			do
  Rscript ../simulations/pred_pac.R $p $n $alpha_id $gamma_id $seed
  Rscript ../simulations/pred_pac_est.R $p $n $alpha_id $gamma_id $seed
done
done
done
done
done