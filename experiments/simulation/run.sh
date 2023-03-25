#!/bin/bash

for N in 100  400 800 1200
do 
	for ID in 1 2 3 4 5
	do
		Rscript main.r $N  $ID 'df1_gen' &
	done
done

wait 
echo -e "df1 done"


for N in 100  400 800 1200
do 
	for ID in 1 2 3 4 5
	do
		Rscript main.r $N  $ID 'df2_gen' &
	done
done

wait 
echo -e "df2 done"

for N in 100  400 800  1200
do 
	for ID in 1 2 3 4 5
	do
		Rscript main.r $N  $ID 'df3_gen' &
	done
done

wait 
echo -e "df3 done"


