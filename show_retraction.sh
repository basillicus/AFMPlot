#!/bin/bash

function get_distance_crash {
# $1: retraction number
# $2: number of atoms
file="../CONTCAR.$1"
alat=$(head -n 5  $file | tail -n 1  | awk '{print $3}')
min_z="10"
max_z="0"

# Adjusting the lines where the atomic coordinates are
let N=$2+9

for j in $(seq 10 $N); do
#for j in $(seq 10 11); do
    #z=$(sed -n "${j}p" < $file | awk '{print $3}')
    z=$(head -n ${j} $file | tail -n1 | awk '{print $3}')
    FLAG=$(echo "${z} < ${min_z}" | bc -l)
    if [[ $FLAG -eq 1 ]]; then 
        min_z=$z
    fi
    if (( $(echo "${z} > ${max_z}" | bc -l) )); then 
        max_z=$z
    fi
done
distance=$(echo "($max_z - $min_z) * $alat" | bc -l)
echo ${distance}
}

function get_distance_retract {
# $1: retraction number
# $2: number of atoms
file="CONTCAR.$1"
alat=$(head -n 5  $file | tail -n 1  | awk '{print $3}')
min_z="10"
max_z="0"

# Adjusting the lines where the atomic coordinates are
let N=$2+9

for j in $(seq 10 $N); do
    #z=$(sed -n "${j}p" < $file | awk '{print $3}')
    z=$(head -n ${j} $file | tail -n1 | awk '{print $3}')
    FLAG=$(echo "${z} < ${min_z}" | bc -l)
    if [[ $FLAG -eq 1 ]]; then 
        min_z=$z
    fi
    if (( $(echo "${z} > ${max_z}" | bc -l) )); then 
        max_z=$z
    fi
done
distance=$(echo "($max_z - $min_z) * $alat" | bc -l)
echo ${distance}
}

CWD=$(pwd)
PATH_APPROACHING=${CWD%/*}
FOLDER_RETRACTION=${CWD##*/}
NUM_RET=${FOLDER_RETRACTION/ret_}

echo $CWD
echo $FOLDER_RETRACTION
echo $NUM_RET


# Clean temporary files in case they exist
if [ -f temp.ener ]; then
    rm temp.ener
fi
if [ -f tmp.xyz ]; then
    rm tmp.xyz
fi
# test function

NATOMS=$(head -n1 evolution.xyz)

# Coger las geometrias del crashing (evolution .xyz) hasta (NUM_RET - 1) == temp_approach.xyz
echo NATOMS $NATOMS
let lala=$NATOMS+2

NLINES=$(echo "${lala} * ${NUM_RET}" | bc)
head -n $NLINES ../evolution.xyz > tmp.xyz
# Anadir las geoemtrias del retracting temp_approach.xyz
cat ./evolution.xyz >> tmp.xyz
mv tmp.xyz crash_reatraction_${NUM_RET}.xyz


# # Coger las energias  del crashing y retracting
# for i in $(seq 1 ${NUM_RET}); do 
#     ener=$(grep 'ee  e' ../OUTCAR.${i} | tail -n 1 | awk '{print $5}')
#     dist=$(get_distance_crash ${i} ${NATOMS})
#     #echo dist: $dist
#     echo $dist   $ener >> temp.ener
# done
# 
# echo "" >> temp.ener
# echo "" >> temp.ener
# for i in $(ls -v ./OUTCAR.*); do 
#     num=${i##*OUTCAR.}
#     ener=$(grep 'ee  e' ${i} | tail -n 1 | awk '{print $5}')
#     dist=$(get_distance_retract ${num} ${NATOMS})
#     echo  $dist  $ener  >> temp.ener
# done
# 
# mv temp.ener ener_crash_retracting_${NUM_RET}.dat
# 
# # Hacer magia, vamos, plottear
