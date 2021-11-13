#!/bin/bash
for pos in {1..2}
	do
	for alpha_id in {1..9}
		do 
		for seed in {1..10} 
			do
  Rscript ../realdata/sens_pac.R $pos $alpha_id $seed
done
done
done