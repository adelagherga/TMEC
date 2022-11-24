#!/bin/bash
# computeEllipticCurves.sh
#
# These functions compute all elliptic curves of given conductor(s), taking as
# input a single conductor N, a range of conductors [N1,...,N2], or, with the
# flag -l, an arbitrary finite list of conductors, with various computations
# executed in parallel.
#
# Parameters
#     N1 [N2]
#         A single conductor, N1, or the range [N1,...,N2].
#     [-l N1] [N2...]: [OPTIONAL]
#         A finite, arbitrary list of conductors to generate [N1,N2,...].
# Returns
#     Dir
#         A subdirectory or Data storing all output, logs, errors, and elliptic
#         curves.
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

validate () {

    # Verifies terminal input is a positive integer, otherwise prints usage
    # statement and terminates the program.
    #
    # Parameters
    #     N
    #         Terminal input, excluding any function flags.

    re='^[1-9][0-9]*(,[1-9][0-9]*)*$'
    if ! [[ "$1" =~ $re ]]; then
        echo "Invalid input: positive integers only." >&2
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
		validate "${OPTARG}"
		list+=("${OPTARG}")
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
	validate "$N"
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
    #     TMDir
    #         The directory Data/${name}/TM, storing all data related to the
    #         Thue--Mahler computation.
    #     TMOutDir
    #         The directory Data/${name}/TM/Outfiles, storing all output from
    #         the Thue--Mahler computation. Files will be created in this
    #         directory only when the corresponding Thue--Mahler equation
    #         produces relevant elliptic curves.
    #     TMLogDir
    #         The directory Data/${name}/TM/Logfiles, storing all logs from the
    #         Thue--Mahler computation.
    #     XYZ2Dir
    #         The directory Data/${name}/XYZ2, storing all data relateed to the
    #         X+Y=Z^2 computation.
    #     XYZ2OutDir
    #         The directory Data/${name}/XYZ2/Outfiles, storing all output from
    #         the X+Y=Z^2 computation. Files will be created in this
    #         directory only when the corresponding S-unit equation produces
    #         relevant elliptic curves.
    #     XYZ2LogDir
    #         The directory Data/${name}/XYZ2/Logfiles, storing all logs from the
    #         X+Y=Z^2 computation.
    #     Dir/Errors.txt
    #         A file tracking all errors from magma.
    #     ECDir
    #         The directory Data/${name}/EllipticCurves, storing all resulting
    #         elliptic curves.
    #     ECDir/${N}.csv
    #         A file storing all elliptic curves of conductor N, generated for
    #         each N in list, in the format N [a1,a2,a3,a4,a6].

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
    TMDir="${Dir}/TM"
    TMOutDir="${TMDir}/Outfiles"
    TMLogDir="${TMDir}/Logfiles"
    XYZ2Dir="${Dir}/XYZ2"
    XYZ2OutDir="${XYZ2Dir}/Outfiles"
    XYZ2LogDir="${XYZ2Dir}/Logfiles"
    mkdir "${ECDir}"
    mkdir "${TMDir}"
    mkdir "${TMOutDir}"
    mkdir "${TMLogDir}"
    mkdir "${XYZ2Dir}"
    mkdir "${XYZ2OutDir}"
    mkdir "${XYZ2LogDir}"
    touch "${Dir}/Errors.txt"

    # Generate files for each conductor.
    for N in "${list[@]}"; do
	touch "${ECDir}/${N}.csv"
    done
}

runParallel() {

    # Runs GNU parallel across 20 cores. That is, runs
    # echo <inFile> | parallel -j20 --joblog <logFile> <program>
    # taking, as input for <program>, each line of <inFile>, and storing GNU
    # parallel's progress in <logFile>.
    #
    # Parameters
    #     inFile
    #         A file containing line-seperated entries to be run with <program>.
    #     logFile
    #         A file storing GNU parallel's progress.
    #     program
    #         The program to be run in parallel, including any redirected stdout,
    #         stderr.

    echo "$1" | parallel -j20 --joblog ${Dir}/"$2" "$3"

}

verifyNonEmpty() {
# FIX THIS!!! AND VERIFY!
    # Verifies whether the file TM/TMForms.csv is empty. When true, this function
    # populates and sorts all elliptic curve files before terminating the
    # program. This function also generates a single file containing all
    # resulting curves, to be used for final comparisons.
    #
    # Parameters
    #     Dir/TM/TMForms.csv
    #         The file Data/${name}/TM/TMForms.csv containing, in each line, the
    #         Thue--Mahler form to be solved.
    # Returns
    #     Dir/EllipticCurves/AllCurves.csv
    #         The file containing all elliptic curves generated in the program,
    #         sorted and free of duplicates, in the format N [a1,a2,a3,a4,a6].

    local msg

    if [ ! -s "$1" ]; then
	msg="Finished computing all elliptic curves of conductor ${name}"
	msg="${msg} corresponding to the $2 solver.\n"
	printf "${msg}"
	sortCurvesByN
	exit 0
    fi
}

gatherRedundantForms() {
# FIX DESCRIOTION AND TEST
    # Amalgamates all Thue--Mahler form files into a single file, TM/TMForms.csv,
    # and removes redundant forms.
    #
    # Parameters
    #     Dir/TM/${N}Forms.csv
    #         The file containing, in each line, the Thue--Mahler form to be
    #         solved, pertaining to conductor N. One such file exists for every N
    #         in ${list}.
    # Returns
    #     Dir/TM/TMForms.csv
    #         The file Data/${name}/TM/TMForms.csv containing, in each line, the
    #         Thue--Mahler form to be solved.

    local forms
    local N
    local F

    if [[ "$1" == "TM" ]]; then
	forms="${TMDir}/TMForms.csv"
	for N in "${list[@]}"; do
	    F="${TMDir}/${N}Forms.csv"
	    if [ -f "${F}" ]; then
		cat "${F}" >> "${forms}"
		rm "${F}"
	    fi
	done
    elif [[ "$1" == "XYZ2" ]]; then
	forms="${XYZ2Dir}/XYZ2Forms.csv"
    else
	echo "Invalid input: file type must be either TM or XYZ2." >&2
	exit 1
    fi

    verifyNonEmpty "${forms}" "$1"

    printf "Removing redundant cases..."
    #FIX THIS FUNCTION AS WELL - ADD COMMENTS AND CHECK determineOF FUNC
    python Code/gatherRedundancy.py "${forms}" "$1"
    printf "Done.\n"
}

amalgamateFormFiles() {

    # Amalgamates all S-unit equations corresponding to Thue--Mahler forms into
    # a single document, TM/TMForms.csv, cleaning up any additional temporary
    # files in the process and forwarding any magma errors to the file
    # Errors.txt.
    #
    # Parameters
    #     Dir/TM/${line}.csv
    #         The output file from Code/TM/optimalForm.m containing, for each
    #         line, the optimal Thue--Mahler S-unit equations to be solved.
    #     Dir/TM/${line}tmp.txt
    #         The magma logfile from Code/TM/optimalForm.m tracking any magma
    #         errors.
    # Returns
    #     Dir/TM/TMForms.csv
    #         The file Data/${name}/TM/TMForms.csv containing, in each line, the
    #         optimal Thue--Mahler S-unit equations to be solved.

    local line
    local F1
    local F2
    local rerun

    while IFS= read -r line; do
	F1="${TMDir}/${line}.csv"
	F2="${TMDir}/${line}tmp.txt"
	if [ -f "${F1}" ]; then
	    cat "${F1}" >> "${TMDir}/tmpTMForms.csv"
	    rm "${F1}"
	fi
	if grep -q "error" "${F2}"; then
	    rerun="magma -b set:=${line} dir:='${TMDir}' "
	    rerun="${rerun} Code/TM/optimalForm.m 2>&1"
	    echo "${rerun}" >> "${Dir}/Errors.txt"
	fi
	rm "${F2}"
    done < "${TMDir}/TMForms.csv"
    mv "${TMDir}/tmpTMForms.csv" "${TMDir}/TMForms.csv"

    verifyNonEmpty "${TMDir}/TMForms.csv"
}

moveTMCurves() {

    # For each line of inFile, extracts the conductor N and copies the line to
    # the file ${ECDir}/${N}.csv.
    #
    # Parameters
    #     inFile
    #         The output file from Code/TM/runTMEC.m containing, for each line,
    #         the elliptic curve obtained by solving the Thue--Mahler S-unit
    #         equation in the filename.

    local line
    local N

    while IFS= read -r line; do
	N=$(echo "${line}" | cut -d' ' -f -1)
	echo "${line}" >> "${ECDir}/${N}.csv"
    done < "$1"
}

sortCurvesByN() {

    # Sorts all elliptic curve files and generates a single file containing all
    # curves.
    #
    # Returns
    #     Dir/EllipticCurves/AllCurves.csv
    #         The file containing all elliptic curves generated in the program,
    #         sorted and free of duplicates, in the format N [a1,a2,a3,a4,a6].

    local N
    local F

    for N in "${list[@]}"; do
	F="${ECDir}/${N}.csv"
	sort -t, -k2,2n -k3,3n -k4,4n -k5,5n -k6,6n "$F" > "${ECDir}/tmp${N}.csv"
	mv "${ECDir}/tmp${N}.csv" "${ECDir}/${N}.csv"
	# Remove duplicate elliptic curves.
	awk -i inplace '!seen[$0]++' "${ECDir}/${N}.csv"
	cat "${ECDir}/${N}.csv" >> "${ECDir}/AllCurves.csv"
    done
}

sortCurves() {

    # Amalgamates all logfiles pertaining to the same Thue--Mahler equation,
    # forwards any magma errors to the file Errors.txt, and populates and sorts
    # all elliptic curve files. This function also generates a single file
    # containing all resulting curves, to be used for final comparisons.
    #
    # Parameters
    #     Dir/TM/${line}Out.csv
    #         The output file from Code/TM/runTMEC.m containing, for each line,
    #         the elliptic curve obtained by solving the Thue--Mahler S-unit
    #         equation in ${line}.
    #     Dir/TM/${line}Log.txt
    #         The magma logfile from Code/TM/runTMEC.m tracking the code's
    #         progress, solutions, and any magma errors that arise.
    # Returns
    #     Dir/TM/${form}Log.txt
    #         The file containing the amalgamated logfiles Dir/TM/${line}Log.txt
    #         pertaining to the same Thue--Mahler form.
    #     Dir/EllipticCurves/AllCurves.csv
    #         The file containing all elliptic curves generated in the program,
    #         sorted and free of duplicates, in the format N [a1,a2,a3,a4,a6].

    local line
    local form
    local Out
    local Log
    local allLog
    local rerun

    while IFS= read -r line; do
	form="$(echo ${line} | cut -d']' -f -3)""]"
	if [ "${line}" != "${form}" ]; then
	    Out="${TMOutDir}/${line}Out.csv"
	    Log="${TMLogDir}/${line}Log.txt"
	    allLog="${TMLogDir}/${form}Log.txt"
	    if grep -q "error" "${Log}"; then
		rerun="magma -b set:='${line}' dir:='${TMDir}'"
		rerun="${rerun} Code/TM/runTMEC.m 2>&1"
		echo "${rerun}" >> "${Dir}/Errors.txt"
	    fi
	    if [ -f "${Out}" ]; then
		moveTMCurves "${Out}"
	    fi
	    cat "${Log}" >> "${allLog}"
	    rm "${Log}"
	fi
    done < "${TMDir}/TMForms.csv"
    rm -r "${TMOutDir}"

    sortCurvesByN
}


    getNList "$@"
    generateDirectories

    # Generate all elliptic curves of j-invariant 0 in conductors list, in
    # parallel. That is, run
    # Code/CurvesNj0 N > ${Dir}/EllipticCurves/N.csv
    # in parallel, with N an entry of ${conductors}, storing GNU parallel's
    # progress in the file ${Dir}/j0Log.
    conductors=$(printf '%s\n' "${list[@]}")
    program="Code/CurvesNj0 {} > '${ECDir}'/{}.csv"
    printf "Generating all j-invariant 0 curves for conductors in ${name}..."
    runParallel "${conductors}" j0Log "${program}"
    printf "Done.\n"

    # Generate all required Thue--Mahler forms in parallel, applying all
    # necessary local tests in the process. That is, run
    # Code/N2TME N '${TMDir}' > /dev/null
    # in parallel, with N an entry of ${conductors}, storing GNU parallel's
    # progress in the file ${Dir}/formLog.
    program="Code/N2TME {} '${TMDir}' > /dev/null"
    printf "Generating all required cubic forms for conductors in ${name}..."
    runParallel "${conductors}" formLog "${program}"
    printf "Done.\n"

    # Remove redundant Thue--Mahler equations.
    gatherRedundantForms "TM"
