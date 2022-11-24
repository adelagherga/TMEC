#!/bin/bash
# computeEllipticCurves.sh

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
    #         The directory Data/${name}/TM/TMOutfiles, storing all output from
    #         the Thue--Mahler computation. Files will be created in this
    #         directory only when the corresponding Thue--Mahler equation
    #         produces relevant elliptic curves.
    #     TMLogDir
    #         The directory Data/${name}/TMLogfiles, storing all logs from the
    #         Thue--Mahler computation.
    #     XYZ2Dir
    #         The directory Data/${name}/XYZ2, storing all data relateed to the
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
    TMOutDir="${TMDir}/TMOutfiles"
    TMLogDir="${TMDir}/TMLogfiles"
    XYZ2Dir="${Dir}/XYZ2"
    XYZ2OutDir="${XYZ2Dir}/XYZ2Outfiles"
    XYZ2LogDir="${XYZ2Dir}/XYZ2Logfiles"
    mkdir "${ECDir}"
    mkdir "${TMDir}"
    mkdir "${XYZ2Dir}"
    mkdir "${XYZ2OutDir}"
    mkdir "${XYZ2LogDir}"
    mkdir "${TMOutDir}"
    mkdir "${TMLogDir}"
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

    if [ ! -s "$1" ]; then
	printf "Finished computing all elliptic curves of conductor ${name}.\n"
	sortCurvesByN
	exit 0
    fi
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


main() {


    getNList "$@"
    generateDirectories
    conductors=$(printf '%s\n' "${list[@]}")

printf "Generating all required prime lists for conductors in ${name}..."
echo "${conductors}" | factor > ${XYZ2Dir}/XYZ2Forms.csv
printf "Done.\n"

verifyNonEmpty "${XYZ2Dir}/XYZ2Forms.csv"


# Remove redundant prime sets.
python Code/gatherRedundancy.py "${XYZ2Dir}/XYZ2Forms.csv" "${XYZ2Dir}/tmpXYZ2Forms.csv" "XYZ2"

# Generate optimal Thue--Mahler forms and all S-unit equations in parallel.
# That is, run
# magma -b set:=<line> dir:=${TMDir} Code/TM/optimalForm.m 2>&1
# in parallel, for each <line> of ${TMForms}, storing GNU parallel's
# progress in the file ${Dir}/optimalLog.
XYZ2Forms=$(cat ${XYZ2Dir}/XYZ2Forms.csv)
program="magma -b set:={} dir:='${XYZ2Dir}' Code/XYZ2/seperateForm.m 2>&1"
printf "Preparing prime list data for parallelization..."
runParallel "${XYZ2Forms}" seperateLog "${program}"
printf "Done.\n"


    local line
    local F1
    local F2
    local rerun

    while IFS= read -r line; do
	F1="${XYZ2Dir}/${line}.csv"
	F2="${XYZ2Dir}/${line}tmp.txt"
	if [ -f "${F1}" ]; then
	    cat "${F1}" >> "${XYZ2Dir}/tmpXYZ2Forms.csv"
	    rm "${F1}"
	fi
	if grep -q "error" "${F2}"; then
	    rerun="magma -b set:=${line} dir:='${XYZ2Dir}' "
	    rerun="${rerun} Code/XYZ2/seperateForm.m 2>&1"
	    echo "${rerun}" >> "${Dir}/Errors.txt"
	fi
	rm "${F2}"
    done < "${XYZ2Dir}/XYZ2Forms.csv"

    mv "${XYZ2Dir}/tmpXYZ2Forms.csv" "${XYZ2Dir}/XYZ2Forms.csv"



    # Solve all Thue--Mahler S-unit equations in parallel.
    # That is, run
    # magma -b set:=<line> dir:=${TMDir} Code/TM/runTMEC.m 2>&1
    # in parallel, for each <line> of ${TMForms}, storing GNU parallel's
    # progress in the file ${Dir}/optimalLog.
    XYZ2Forms=$(cat ${XYZ2Dir}/XYZ2Forms.csv)
    program="magma -b set:={} dir:='${XYZ2Dir}' Code/XYZ2/runXYZ2EC.m 2>&1"
    printf "Solving the Thue--Mahler equations..."
    runParallel "${XYZ2Forms}" XYZ2Log "${program}"
    printf "Done.\n"



    while IFS= read -r line; do
	form="$(echo ${line} | cut -d']' -f -2)""]" #FYI CHANGED THIS FROM 3 TO 2!!
	if [ "${line}" != "${form}" ]; then
	    Out="${XYZ2OutDir}/${line}Out.csv"
	    Log="${XYZ2LogDir}/${line}Log.txt"
	    allLog="${XYZ2LogDir}/${form}Log.txt"
	    if grep -q "error" "${Log}"; then
		rerun="magma -b set:='${line}' dir:='${XYZ2Dir}'"
		rerun="${rerun} Code/XYZ2/runXYZ2EC.m 2>&1"
		echo "${rerun}" >> "${Dir}/Errors.txt"
	    fi
	    if [ -f "${Out}" ]; then
		moveTMCurves "${Out}"
	    fi
	    cat "${Log}" >> "${allLog}"
	    rm "${Log}"
	fi
    done < "${XYZ2Dir}/XYZ2Forms.csv"
    rm -r "${XYZ2OutDir}"

    printf "Sorting all elliptic curves..."
    sortCurvesByN
    printf "Done.\n"
    printf "Finished computing all elliptic curves of conductor ${name}.\n"
}

main "$@"
