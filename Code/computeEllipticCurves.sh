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
    tmpname="${name}${iter}"
    while [ -d "Data/${tmpname}" ]; do
	iter=$(( $iter + 1 ))
	tmpname="${name}${iter}"
    done
    name=$tmpname
    mkdir "Data/${name}"
fi

# Generate necessary subdirectories
Dir="Data/${name}"
ECDir="${Dir}/EllipticCurves"
TMOutDir="${Dir}/TMOutfiles"
TMLogDir="${Dir}/TMLogfiles"
mkdir "${ECDir}"
mkdir "${TMOutDir}"
mkdir "${TMLogDir}"

# Generate files for each conductor and populate each file with the
# corresponding elliptic curves.
for N in "${list[@]}"; do
    touch "${ECDir}/${N}.csv"
done

# Generate all required Thue--Mahler forms in parallel, applying all necessary
# local tests in the process.
printf "Generating all required cubic forms for conductors in"
printf " $(echo ${name} | cut -d']' -f -1)]..."
(for N in "${list[@]}"; do echo "$N"; done) | parallel -j20 magma -b N:={} name:=${name} Code/findForms.m 2>&1

# Amalgamate all Thue--Mahler forms into a single document.
for N in "${list[@]}"; do
    F1="${Dir}/${N}Forms.csv"
    F2="${Dir}/${N}tmp.txt"
    [ -f "$F1" ] && cat "$F1" >> "${Dir}/TMForms.csv"
    rm -f "$F1"
    rm -f "$F2"
done
printf "Done.\n"

# Remove redundant Thue--Mahler equations.
# chmod +x Code/gatherFormRedundancy.py
printf "Removing redundant cubic forms..."
python Code/gatherFormRedundancy.py "${Dir}/TMForms.csv" "${Dir}/tmpTMForms.csv"
printf "Done.\n"

# Generate optimal Thue--Mahler forms and all S-unit equations.
printf "Generating optimal GL2(Z)-equivalent cubic forms..."
cat ${Dir}/TMForms.csv | parallel -j20 magma -b set:={} name:=${name} Code/optimalForm.m 2>&1

# Amalgamate all S-unit equations into a single document.
while IFS= read -r line; do
    F="${Dir}/${line}.csv"
    F2="${Dir}/${line}tmp.txt"
    [ -f "$F" ] && cat "$F" >> "${Dir}/tmpTMForms.csv"
    rm -f "$F"
    rm -f "$F2"
done < "${Dir}/TMForms.csv"
mv "${Dir}/tmpTMForms.csv" "${Dir}/TMForms.csv"
printf "Done.\n"

# Run Thue--Mahler code in parallel.
# That is, for each line "set" of Data/TMForms/${name}Forms.csv, run
# magma set:="set" Code/computeEllipticCurvesTM.m &.
# The following code runs these jobs using GNU parallel, running no more than
# 20 (-j20) jobs at once, and storing GNU parallel's progress in the logfile
# Data/${name}TMLog (--joblog Data/${name}TMLog).
printf "Solving the Thue--Mahler equations..."
cat ${Dir}/TMForms.csv | parallel -j20 --joblog ${Dir}/TMLog magma -b set:={} name:=${name} Code/computeEllipticCurvesTM.m 2>&1
printf "Done.\n"

# Amalgamate all logfiles pertaining to the same Thue--Mahler equation.
printf "Gathering all files..."
while IFS= read -r form_i; do
    form="$(echo ${form_i} | cut -d']' -f -3)""]"
    if [ "${form_i}" != "${form}" ]; then
	IFOut="${TMOutDir}/${form_i}Out.csv"
	OFOut="${TMOutDir}/${form}Out.csv"
	IFLog="${TMLogDir}/${form_i}Log.txt"
	OFLog="${TMLogDir}/${form}Log.txt"
	[ -f "${IFOut}" ] && cat "${IFOut}" >> "${OFOut}"
	cat "${IFLog}" >> "${OFLog}"
	rm -f "${IFOut}"
	rm -f "${IFLog}"
    fi
done < "${Dir}/TMForms.csv"
printf "Done.\n"

for F in "${TMOutDir}"/*; do
    cat "$F" >> "${ECDir}/AllCurves.csv"
done
sort -t, -k1,1n -k2,2n -k3,3n -k4,4n -k5,5n -k6,6n "${ECDir}/AllCurves.csv" > "${ECDir}/tmpAllCurves.csv"
mv "${ECDir}/tmpAllCurves.csv" "${ECDir}/AllCurves.csv"
awk -F, 'BEGIN {OFS=FS} {print $1,$2,$3,$4,$5,$6}' "${ECDir}/AllCurves.csv" > "${ECDir}/tmpAllCurves.csv"
mv "${ECDir}/tmpAllCurves.csv" "${ECDir}/AllCurves.csv"
awk -i inplace '!seen[$0]++' "${ECDir}/AllCurves.csv"

while IFS= read -r Ncurve; do
    N="$(echo ${Ncurve} | cut -d',' -f -1)"
    echo ${Ncurve} >> "${ECDir}/${N}.csv"
done < "${ECDir}/AllCurves.csv"
printf "Done.\n"
printf "Finished computing all elliptic curves of conductor"
printf " $(echo ${name} | cut -d']' -f -1)]."
