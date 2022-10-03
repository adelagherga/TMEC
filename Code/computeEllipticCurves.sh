#!/bin/bash
# computeEllipticCurves.sh
#
# This function initiates the Thue and Thue--Mahler solver in Magma in parallel
# and amalgamates the elliptic curves in seperate files, by conductor.
#
# Authors
#    Adela Gherga <adelagherga@gmail.com>
# Created
#    24 August 2022

while getopts ":l:" opt; do
    case $opt in
	l)
	    # List of conductors.
	    list+=("$OPTARG")
	    name="${list[*]}"
	    name="["${data_string//${IFS:0:1}/,}"]"
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG." >&2
	    exit 1 ;;
	:)
	    echo "Option -$OPTARG requires an argument." >&2
	    exit 1 ;;
    esac
done
if [ -z "${list}" ]; then
    if [ $# -eq 0 ]; then
	echo "Argument required." >&2
	exit 1
    fi
    if [ $# -eq 1 ]; then
	list=($1)
	name="[$1]"
    elif [ $# -eq '2' ]; then
	list=($(seq $1 $2))
	name="[$1..$2]"
    else
	echo "Invalid argument."
	exit 1
    fi
fi

# Generate a directory Data/${name} for all output.
# If such a directory already exists in Data/, generate a directory
# Data/${name}i, where Data/${name}j exists for all j < i.
if [ ! -d "Data/${name}" ]; then
    mkdir "Data/${name}"
else
    iter=1
    name="${name}${iter}"
    while [ -d "Data/${name}" ]; do
	iter=$(( $iter + 1 ))
	name="${name}${iter}"
    done
    mkdir "Data/${name}"
fi

# Generate necessary subdirectories
mkdir Data/${name}/EllipticCurves
mkdir Data/${name}/TMOutfiles
mkdir Data/${name}/TMLogfiles

# Generate files for each conductor and populate each file with the
# corresponding elliptic curves.
for N in "${list[@]}"; do
    touch "Data/${name}/EllipticCurves/${N}.csv"
done

# Generate all required Thue--Mahler forms in parallel, applying all necessary
# local tests in the process.
echo "Generating all required cubic forms for conductors in $name."
(for N in "${list[@]}"; do echo "$N"; done) | parallel -j20 magma -b N:={} name:=${name} Code/findForms.m 2>&1

# Amalgamate all Thue--Mahler forms into a single document.
for N in "${list[@]}"; do
    F1="Data/${name}/${N}Forms.csv"
    F2="Data/${name}/${N}tmp.txt"
    [ -f "$F1" ] && cat "$F1" >> "Data/${name}/${name}TMForms.csv"
    rm -f "$F1"
    rm -f "$F2"
done

# Remove redundant Thue--Mahler equations.
# chmod +x Code/gatherFormRedundancy.py # DO WE NEED THIS?
echo "Removing redundant cubic forms."
python Code/gatherFormRedundancy.py "Data/${name}/${name}TMForms.csv" "Data/${name}/${name}SortedTMForms.csv"

# Generate optimal Thue--Mahler forms and all S-unit equations.
cat Data/${name}/${name}TMForms.csv | parallel -j20 magma -b set:={} name:=${name} Code/optimalForm.m 2>&1

# Amalgamate all S-unit equations into a single document.
while IFS= read -r line; do
    F="Data/${name}/${line}.csv"
    F2="Data/${name}/${line}tmp.txt"
    [ -f "$F" ] && cat "$F" >> "Data/${name}/${name}SUnitTMForms.csv"
    rm -f "$F"
    rm -f "$F2"
done < "Data/${name}/${name}TMForms.csv"
mv Data/${name}/${name}SUnitTMForms.csv Data/${name}/${name}TMForms.csv

# Run Thue--Mahler code in parallel.
# That is, for each line "set" of Data/TMForms/${name}Forms.csv, run
# magma set:="set" Code/computeEllipticCurvesTM.m &.
# The following code runs these jobs using GNU parallel, running no more than
# 20 (-j20) jobs at once, and storing GNU parallel's progress in the logfile
# Data/${name}TMLog (--joblog Data/${name}TMLog).
cat Data/${name}/${name}TMForms.csv | parallel -j20 --joblog Data/${name}/${name}TMLog magma -b set:={} name:=${name} Code/computeEllipticCurvesTM.m 2>&1

# Amalgamate all logfiles pertaining to the same Thue--Mahler equation.



#for F in "Data/${name}/TMLogfiles"/*; do
 #   filename="${F##*/}"
  #  filename="$(echo ${filename} | cut -d']' -f -3)""]"
#    cat "$F" >> "Data/${name}/TMLogfiles/${filename}Log.csv"
#    rm -f "$F"
#done
#for F in "Data/${name}/TMOutfiles"/*; do
#    filename="${F##*/}"
#    filename="$(echo ${filename} | cut -d']' -f -3)""]"
#    Ns="$(echo ${filename} | cut -d']' -f -1)"
#    Ns="${Ns:1}"
#    readarray -td '' Ns < <(awk '{ gsub(/,+/,"\0"); print; }' <<<"$Ns")
#    cat "$F" >> "Data/${name}/TMOutfiles/${filename}Out.csv"

#    for N in "${Ns[@]}"; do
#	EC="${N}"
#	cat "$F" >> "Data/${name}/TMOutfiles/${filename}Out.csv"

#    done
#    rm -f "$F"
done
