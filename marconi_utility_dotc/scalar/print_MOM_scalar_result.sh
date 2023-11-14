#!/bin/bash

# block size
block=200

Nt=50
Ns=30

kappa_max=$((Ns/2))

# print G(tT)/T^5 vs tT for several values of n_cool and kappa

for kappa in $( seq 0 1 ${kappa_max} );do
                num_line=$( wc -l < results_MOM_tcorr_scalar/kappa_${kappa}/mean_tcorr_kappa_${kappa}.dat)

                out_file="results_MOM_tcorr_scalar/kappa_${kappa}/tcorr_MOM_scalar_kappa_${kappa}_Re.dat"
                echo '#t       Re{G(t)}   err_Re{G(t)}' > ${out_file}
		tcorr=$( awk -v B=${block} -v Nt=${Nt} -v num_line=${num_line} '(NR<=num_line/2 && $1==B && $2<Nt){printf("%0.15lg\t%0.15lg\t%0.15lg\n", $2/Nt, $3, $4)}' results_MOM_tcorr_scalar/kappa_${kappa}/mean_tcorr_kappa_${kappa}.dat )
                echo "${tcorr}" >> ${out_file}
                echo '#-----------------------------------#' >> ${out_file}
                out_file="results_MOM_tcorr_scalar/kappa_${kappa}/tcorr_MOM_scalar_kappa_${kappa}_Im.dat"
                echo '#t Im{G(t)} err_Im{G(t)}' > ${out_file}
		tcorr=$( awk -v B=${block} -v Nt=${Nt} -v num_line=${num_line} '(NR>num_line/2 && $1==B && $2<Nt){printf("%.15lg %.15lg %.15lg\n", $2/Nt,  $3, $4)}' results_MOM_tcorr_scalar/kappa_${kappa}/mean_tcorr_kappa_${kappa}.dat)
                echo "${tcorr}" >> ${out_file}
done

echo "Analisi completata con successo!"
