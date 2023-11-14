#!/bin/bash

therm=10000
B_spacing=50
B_max=5000

mkdir -p results_MOM_tcorr_scalar

Nt=50
Ns=30

#Remember to adjust for your lattice dimension
kappa_max=$((Ns/2))
echo "Analysis started correctly!"

for kappa in $(seq 0 1 ${kappa_max}); do
                
		echo "Parte analisi kappa = ${kappa}"
		./script_MOM_analysis.sh ${therm} ${kappa} ${B_spacing} ${B_max}
                
done
echo "Ora scrivo i risultati"
./print_MOM_scalar_results.sh
