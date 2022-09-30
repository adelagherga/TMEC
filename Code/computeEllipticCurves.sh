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




# Generate files for each conductor and populate each file with the
# corresponding elliptic curves.
for N in "${list[@]}"; do
    touch "Data/${name}/EllipticCurves/${N}.csv"
done

# Generate all required Thue--Mahler forms in parallel, applying all necessary
# local tests in the process.
echo "Generating all required cubic forms for conductors in $name."
(for N in "${list[@]}"; do echo "$N"; done) | parallel -j20 magma N:={} name:=${name} Code/findForms.m 2>&1

# Amalgamate all Thue--Mahler forms into a single document.
for N in "${list[@]}"; do
    F="Data/${name}/${N}Forms.csv"
    [ -f "$F" ] && cat "$F" >> "Data/${name}/${name}TMForms.csv"
    rm -f "$F"
done

# Remove redundant Thue--Mahler equations.
chmod +x Code/gatherFormRedundancy.py
python Code/gatherFormRedundancy.py "Data/${name}/${name}TMForms.csv" "Data/${name}/${name}SortedTMForms.csv"
mv Data/${name}/${name}SortedTMForms.csv Data/${name}/${name}TMForms.csv

# Generate optimal Thue--Mahler forms and all S-unit equations.
cat Data/${name}/${name}TMForms.csv | parallel -j20 magma set:={} name:=${name} Code/optimalForm.m 2>&1

# Amalgamate all S-unit equations into a single document.
while IFS= read -r line; do
    F0="Data/TMForms/$line.csv"
    [ -f "$F" ] && cat "$F" >> "Data/TMForms/${name}SUnitForms.csv"
    rm -f "$F"
done < "Data/TMForms/${name}Forms.csv"

# Run Thue--Mahler code in parallel.
# That is, for each line "set" of Data/TMForms/${name}Forms.csv, run
# magma set:="set" Code/computeEllipticCurvesTM.m &.
# The following code runs these jobs using GNU parallel, running no more than
# 20 (-j20) jobs at once, and storing GNU parallel's progress in the logfile
# Data/${name}TMLog (--joblog Data/${name}TMLog).
cat Data/TMForms/${name}SUnitForms.csv | parallel -j20 --joblog Data/${name}TMLog magma set:={} Code/computeEllipticCurvesTM.m 2>&1

## LEFT OFF HERE

# Amalgamate all elliptic curves from Thue--Mahler output.
while IFS= read -r line; do
    F="Data/TMOutfiles/$line.csv"
    [ -f "$F" ] && cat "$F" >> "Data/TMForms/${name}SortedForms.csv"
    rm -f "$F"
done < "Data/TMForms/${name}Forms.csv"
mv Data/TMForms/${name}SortedForms.csv Data/TMForms/${name}Forms.csv



for F in Data/TMOutfiles/*; do
    N=$(echo $F | grep -o -E '[0-9]+' | head -1 | sed -e 's/^0\+//')
    cat "$F" >> "Data/EllipticCurves/$N.csv"
done
# Amalgamate all elliptic curves from Thue output.
for F in Data/ThueOutfiles/*; do
    N=$(echo $F | grep -o -E '[0-9]+' | head -1 | sed -e 's/^0\+//')
    cat "$F" >> "Data/EllipticCurves/$N.csv"
done
