#!/bin/bash
# computeEllipticCurves.sh
#
# This function computes all elliptic curves of conductors N1,...,Nn and
# sorts the curves into seperate files by conductor.
#
# Parameters
#     [OPTIONAL] -l

# Authors
#    Adela Gherga <adelagherga@gmail.com>
# Created
#    24 August 2022

usage() {
    echo "usage: "
    echo "  $0 N1 [N2]"
    echo "  $0 [-l N1] N2..."
    echo ""
    echo "  Generate elliptic curves having conductors in the range N1 to N2."
    echo "  If N2 is omitted, generate elliptic curves having conductor N1."
    echo ""
    echo "  -l  Generate elliptic curves having conductors in the list"\
	 "[N1,N2,...]."
    exit 1
}

getNRange() {

    # Parses terminal input to determine all elliptic curves of conductors
    # N1,...,N2, or if N2 is omitted, of conductor N1. This function determines
    # the list of conductors to be resolved, along with an appropriate directory
    # name.
    #
    # Parameters
    #     N1 [N2]
    #         A single conductor, N1, or the range [N1,...,N2].
    #     [-l N1] [N2...]: [OPTIONAL]
    #         A finite, arbitrary list of conductors to generate [N1,N2,...].
    # Returns
    #     list
    #         The list of conductors to be resolved.
    #     name
    #         The name of the directory in which computations will take place,
    #         in the format "[N1,..,N2]".

    local N1
    local N2

    if [ $# -eq 0 ]; then
	echo "Invalid input: argument required." >&2
	usage >&2
    elif [ $# -eq 1 ]; then
	N1=$((10#"$1"))
	list=("${N1}")
	name="[""${N1}""]"
    elif [ $# -eq 2 ]; then
	if [ "$1" -eq "$2" ]; then
	    N1=$((10#"$1"))
	    list=("${N1}")
	    name="[""${N1}""]"
	elif [ "$1" -lt "$2" ]; then
	    N1=$((10#"$1"))
	    N2=$((10#"$2"))
	    list=($(seq "${N1}" "${N2}"))
	    name="[""${N1}""..""${N2}""]"
	else
	    echo "Invalid input: N1 must be less than or equal to N2." >&2
	    usage >&2
	fi
    else
	echo "Invalid input: too many arguments: use option -l." >&2
	usage >&2
    fi
}

getNList() {

    # Parses terminal input and generates the list of conductors to be resolved,
    # along with an appropriate directory name. This function handles a single
    # conductor N, a range of conductors [N1,...,N2], or, with the flag -l, an
    # arbitrary finite list of conductors.
    #
    # Parameters
    #     N1 [N2]
    #         A single conductor, N1, or the range [N1,...,N2].
    #     [-l N1] [N2...]: [OPTIONAL]
    #         A finite, arbitrary list of conductors to generate [N1,N2,...].
    # Returns
    #     list
    #         The list of conductors to be resolved.
    #     name
    #         The name of the directory in which computations will take place,
    #         in the format "[N1,..,N2]" or "[N1,N2,..,Nn]".

    local OPTIND
    local N
    local Nlist

    while getopts ":l:" opt; do
	case $opt in
	    l)
		list+=("$((10#${OPTARG}))")
		;;
	    \?)
		echo "Invalid option: -${OPTARG}." >&2
		usage >&2
		;;
	    :)
		echo "Option -$OPTARG requires an argument." >&2
		usage >&2
		;;
	esac
    done
    shift $(($OPTIND - 1))

    for N in "$@"; do
	if ! [[ "$N" =~ ^[0-9]+$ ]]; then
	    echo "Invalid input: positive integers only." >&2
	    usage >&2
	fi
    done
    if [ -z "${list}" ]; then
	# The input is not in list format.
	getNRange "$@"
    else
	for N in "$@"; do
	    list+=("$((10#$N))")
	done
	printf -v Nlist '%s,' "${list[@]}"
	name="[""${Nlist%,}""]"
    fi
    list="$(printf '%s\n' "${list[@]}")"
}

generateDirectories() {

    # Generates a directory Data/${name} for all output, as well as all
    # necessary subdirectories and files. If such a directory
    # already exists in Data/, generates a directory Data/${name}i, where
    # Data/${name}j exists for all j < i.
    #
    # Parameters
    #     list
    #         The list of conductors to be resolved.
    #     name
    #         The name of the directory in which computations will take place,
    #         in the format "[N1,..,N2]" or "[N1,N2,..,Nn]".
    # Returns
    #     Dir
    #         The directory Data/${name}, storing all output, logs, and errors.
    #     ECDir
    #         The directory Data/${name}/EllipticCurves, storing all resulting
    #         elliptic curves.
    #     TMOutDir
    #         The directory Data/${name}/TMOutfiles, storing all output from the
    #         Thue--Mahler computation. Files will be created in this directory
    #         only when the corresponding Thue--Mahler equation produces relevant
    #         elliptic curves.
    #     TMLogDir
    #         The directory Data/${name}/TMLogfiles, storing all logs from the
    #         Thue--Mahler computation.
    #     Dir/Errors.txt
    #         A file tracking all errors from Magma. REDIRECT TO STDERR?!
    #     ECDir/N.csv
    #         A file storing all elliptic curves of conductor N, generated for
    #         each N in list.

    local iter
    local tmpname
    local N

    # Generate Data directory, if it does not already exist.
    if [ ! -d "Data/" ]; then
	mkdir "Data/"
    fi

    # Generate Data/${name} directory.
    if [ ! -d "Data/${name}" ]; then
	Dir="Data/${name}"
    else
	iter=1
	tmpname="${name}${iter}"
	while [ -d "Data/${tmpname}" ]; do
	    iter=$(( "${iter}" + 1 ))
	    tmpname="${name}${iter}"
	done
	Dir="Data/${tmpname}"
    fi

    # Generate necessary subdirectories and files.
    mkdir "${Dir}"
    ECDir="${Dir}/EllipticCurves"
    TMOutDir="${Dir}/TMOutfiles"
    TMLogDir="${Dir}/TMLogfiles"
    mkdir "${ECDir}"
    mkdir "${TMOutDir}"
    mkdir "${TMLogDir}"
    touch "${Dir}/Errors.txt"

    # Generate files for each conductor.
    for N in "${list[@]}"; do
	touch "${ECDir}/${N}.csv"
    done
}

runParallel() {

    # Runs magma in parallel. That is, for each entry of parInput, arg, runs
    # magma -b arg:="arg" dir:=${Dir} Code/funcFile. The following code runs
    # these jobs using GNU parallel, running no more than
    # 20 (-j20) jobs at once, and storing GNU parallel's progress in the logfile
    # ${Dir}/TMLog (--joblog ${Dir}/TMLog).
    #
    # Parameters
    #     parInput
    #         A list of entries, arg, for the parameter arg in magma, to be parsed
    #         in parallel.
    #     arg
    #         A variable to be defined in magma.
    #     funcFile
    #         The magma startup file.

    echo "$1" | parallel -j20 --joblog ${Dir}/"$2" "$3"

}

gatherRedundancy() {

    local N
    local F

    for N in "${list[@]}"; do
	F="${Dir}/"$1".csv"
	[ -f "${F}" ] && cat "${F}" >> "${Dir}/TMForms.csv"
	rm -f "${F}"
    done

    # Remove redundant Thue--Mahler equations.
    python Code/gatherFormRedundancy.py "${Dir}/TMForms.csv" \
	   "${Dir}/tmpTMForms.csv"

}

amalgamateFormFiles() {

    # Amalgamates all Thue--Mahler forms into a single document.

    # Parameters
    #     $1
    #         An output file LEFT OFF HERE

    local line
    local F1
    local F2

    while IFS= read -r line; do
	F1="${Dir}/""${line}"".csv"
	F2="${Dir}/""${line}""tmp.txt"
	[ -f "${F1}" ] && cat "${F1}" >> "${Dir}/tmpTMForms.csv"
	rm -f "${F1}"
	if grep -q "error" "${F2}"; then
	    echo "${line}:optimalForm.m" >> "${Dir}/Errors.txt"
	fi
	rm -f "${F2}"
    done < "${Dir}/TMForms.csv"
    mv "${Dir}/tmpTMForms.csv" "${Dir}/TMForms.csv"

    if ! [ -s "${Dir}/TMForms.csv" ]; then
	printf "Finished computing all elliptic curves of conductor ${name}.\n"
	exit 0
    fi

}

amalgamateTMFiles() {

    # Amalgamate all logfiles and outfiles pertaining to the same Thue--Mahler
    # equation.

    local form_i
    local form
    local IFOUT
    local OFOUT
    local IFLog
    local OFLog

    while IFS= read -r form_i; do
	form="$(echo ${form_i} | cut -d']' -f -3)""]"
	if [ "${form_i}" != "${form}" ]; then
	    IFOut="${TMOutDir}/${form_i}Out.csv"
	    OFOut="${TMOutDir}/${form}Out.csv"
	    IFLog="${TMLogDir}/${form_i}Log.txt"
	    OFLog="${TMLogDir}/${form}Log.txt"
	    if grep -q "error" "$IFLog"; then
		echo "${form_i}:computeEllipticCurvesTM.m" >> "${Dir}/Errors.txt"
	    fi
	    [ -f "${IFOut}" ] && cat "${IFOut}" >> "${OFOut}"
	    cat "${IFLog}" >> "${OFLog}"
	    rm -f "${IFOut}"
	    rm -f "${IFLog}"
	fi
    done < "${Dir}/TMForms.csv"

}

sortEllipticCurves() {
    # Amagamate all elliptic curves.
    local F
    local Ncurve
    local N

    for F in "${TMOutDir}"/*; do
	cat "$F" >> "${ECDir}/AllCurves.csv"
    done
    sort -t, -k1,1n -k2,2n -k3,3n -k4,4n -k5,5n -k6,6n "${ECDir}/AllCurves.csv" \
	 > "${ECDir}/tmpAllCurves.csv"
    mv "${ECDir}/tmpAllCurves.csv" "${ECDir}/AllCurves.csv"
    # Output only N,aInvariants for each elliptic curve.
    awk -F, 'BEGIN {OFS=FS} {print $1,$2,$3,$4,$5,$6}' "${ECDir}/AllCurves.csv" \
	> "${ECDir}/tmpAllCurves.csv"
    mv "${ECDir}/tmpAllCurves.csv" "${ECDir}/AllCurves.csv"
    # Remove duplicate elliptic curves.
    awk -i inplace '!seen[$0]++' "${ECDir}/AllCurves.csv"

    # Sort elliptic curves by conductor into seperate files.
    while IFS= read -r Ncurve; do
	N="$(echo ${Ncurve} | cut -d',' -f -1)"
	echo ${Ncurve} >> "${ECDir}/${N}.csv"
    done < "${ECDir}/AllCurves.csv"
}

printErrors() {
    # Raise warning if errors were detected.
    if [ -s "${Dir}/Errors.txt" ]; then
	# The file is not-empty.
	printf "Errors detected:\n"
	cat "${Dir}/Errors.txt"
    fi
}

main () {

    # Establishes run order.

    local N
    local line
    getNList "$@"
    generateDirectories

    # Generate all required Thue--Mahler forms in parallel, applying all
    # necessary local tests in the process.
    printf "Generating all required cubic forms for conductors in ${name}..."
    runParallel "${list}" Log "Code/N2TME {} '${Dir}' > /dev/null"
    printf "Done.\n"

    printf "Removing redundant cubic forms..."
    gatherRedundancy
    printf "Done.\n"

    # Generate optimal Thue--Mahler forms and all S-unit equations.
    printf "Generating optimal GL2(Z)-equivalent cubic forms..."
    runParallel "$(cat ${Dir}/TMForms.csv)" TMLog \
	       "magma set:={} dir:=${Dir} Code/optimalForm.m 2>&1"
    printf "Done.\n"
    amalgamateFormFiles

    printf "Solving the Thue--Mahler equations..."
    runParallel "$(cat ${Dir}/TMForms.csv)" TMLog \
		"magma set:={} dir:=${Dir} Code/computeEllipticCurvesTM.m 2>&1"
    printf "Done.\n"
    amalgamateTMFiles

    sortEllipticCurves
    printErrors
    printf "Finished computing all elliptic curves of conductor ${name}.\n"

}

main "$@"
