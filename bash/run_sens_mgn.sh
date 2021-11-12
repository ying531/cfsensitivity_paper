#!/bin/bash
for gamma_id in {1..6}
	do
    for diff_id in {1.5}
	do
	for ite in {1..2}
	    do
	    for seed in {1..100} 
		do
  Rscript ../simulations/sens_mgn.R $gamma_id $diff_id $ite $seed
  Rscript ../simulations/sens_mgn_est.R $gamma_id $diff_id $ite $seed
done
done
done
done