#!/bin/bash

if [ "$#" -ne 4 ]; then
        echo "ERROR! Correct usage: $0 thermalization B_spacing B_max"
        exit 1
fi

therm=$1
kappa=$2
B_spacing=$3
B_max=$4

save_folder="results_MOM_tcorr_scalar/kappa_${kappa}"
mkdir -p ${save_folder}

data_file="data_MOM_tcorr_kappa_${kappa}.dat"
out_file="mean_tcorr_kappa_${kappa}.dat"
rm -f ${data_file}

L=30
T=50


infile="../prova2D.txt"

awk -v t=${therm} -v T=${T} -v kappa=${kappa} '(NR> 1 && $1>t  && $3==kappa){for(i=4;i<=2*T+3;i++){printf "%s ", $i}; printf "\n"}' ${infile} >> ${data_file}

echo "Sto creando: ${out_file}"
./analysis_scalar_MOM_tcorr.x ${data_file} ${out_file} ${B_max} ${B_spacing} ${L} ${T}
echo "Analisi per ${kappa} eseguita correttamente"

num_line=$(wc -l < ${out_file})


#Nota al posto del 9 ci va T-1
tmax=$((T-1))
for t in $( seq 0 1 ${tmax} ); do

        out_file_t="scalar_MOM_Retcorr_${t}_vs_block_size.dat"
        awk -v t=${t} -v num_line=${num_line} '(NR<=num_line/2 && $2==t){print($1, $3, $4, $5, $6)}' ${out_file} > ${out_file_t}
        mv ${out_file_t} ${save_folder}
        out_file_t="scalar_MOM_Imtcorr_${t}_vs_block_size.dat"
        awk -v t=${t} -v num_line=${num_line} '(NR>num_line/2 && $2==t){print($1, $3, $4, $5, $6)}' ${out_file} > ${out_file_t}
        mv ${out_file_t} ${save_folder}
done

mv ${out_file} ${save_folder}
mv ${data_file} ${save_folder}
