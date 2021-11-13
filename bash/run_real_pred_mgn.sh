#!/bin/bash
for alpha_id in {1..9}
	do
	for gamma_id in {1..5}
		do 
		for seed in {1..100} 
			do
  Rscript ../realdata/syn_pred_mgn.R $alpha_id $gamma_id $seed
done
done
done