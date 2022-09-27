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
    elif [ $# -eq '2' ]; then
	list=($(seq $1 $2))
    else
	echo "Invalid argument."
	exit 1
    fi
fi

mkdir Data
mkdir Data/EllipticCurves
# Generate files for each conductor and populate each file with the
# corresponding elliptic curves.
for N in "${list[@]}"; do
    touch "Data/EllipticCurves/${N}.csv"
done

echo "Generating cubic forms"
(for N in "${list[@]}"; do echo "$N"; done) | parallel -j20 magma N:={} Code/findForms.m 2>&1
(for N in "${list[@]}"; do echo "$N"; done) | magma N:=$N Code/findForms.m 2>&1


# Run ThueMahler code in parallel.
# That is, for each line "set" of Data/Forms/TMTestForms.csv, run
# magma set:="set" Code/computeEllipticCurvesTM.m &.
# The following code runs these jobs using GNU parallel, running no more than
# 20 (-j20) jobs at once, and storing GNU parallel's progress in the logfile
# Data/TMTest1Log (--joblog Data/TMTest1Log).
cat Data/Forms/TMTestForms.csv | parallel -j20 --joblog Data/TMTest1Log magma set:={} Code/computeEllipticCurvesTM.m 2>&1

# Generate files for each conductor and populate each file with the
# corresponding elliptic curves.
END=500999
START=500000
for ((N=START;N<=END;N++)); do
    touch Data/EllipticCurves/$N.csv
done

# Amalgamate all elliptic curves from Thue--Mahler output.
for F in Data/TMOutfiles/*; do
    N=$(echo $F | grep -o -E '[0-9]+' | head -1 | sed -e 's/^0\+//')
    cat "$F" >> "Data/EllipticCurves/$N.csv"
done
# Amalgamate all elliptic curves from Thue output.
for F in Data/ThueOutfiles/*; do
    N=$(echo $F | grep -o -E '[0-9]+' | head -1 | sed -e 's/^0\+//')
    cat "$F" >> "Data/EllipticCurves/$N.csv"
done
