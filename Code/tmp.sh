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

    # Parse terminal input to determine all elliptic curves of conductors
    # N1,...,N2, or if N2 is omitted, of conductor N1. Determine the list of
    # conductors to be resolved, along with an appropriate directory name.
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

validate() {

    # Verify terminal input is a positive integer, otherwise print usage
    # statement and terminate the program.
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

    # Parse terminal input and generate the list of conductors to be resolved,
    # along with an appropriate directory name. Input can be a single conductor
    # N, a range of conductors [N1,...,N2], or, with the flag -l, an arbitrary
    # finite list of conductors.
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

    # Generate a directory Data/${name} for all output, as well as all necessary
    # subdirectories and files. If such a directory already exists in Data/,
    # generate a directory Data/${name}i, where Data/${name}j exists for all
    # j < i.
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

    # Run GNU parallel across 20 cores. That is, run
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

    # Verify whether the file Data/coresponding to <fileType> is empty. When
    # true, terminate the current program.
    #
    # Parameters
    #     Dir/TM/TMForms.csv
    #         The file Data/${name}/TM/TMForms.csv containing, in each line, the
    #         Thue--Mahler form to be solved.
    # Returns
    #     Dir/EllipticCurves/AllCurves.csv
    #         The file containing all elliptic curves generated in the program,
    #         sorted and free of duplicates, in the format N [a1,a2,a3,a4,a6].

    local formFile
    local logFile
    local msg

    # how can we print this message in the respective log?
    formFile="${Dir}/$1/$1Forms.csv"
    logFile="${Dir}/$1Status.txt"

    if [ ! -s "${formFile}" ]; then
	msg="Finished computing all elliptic curves of conductor ${name}"
	msg="${msg} corresponding to the $1 solver.\n"
	printf "${msg}" >> "${logFile}"
	exit 0
    fi
}

gatherRedundantForms() {

    # Remove redundant equations across conductors for the Thue--Mahler solver
    # and the XYZ2 solver. In the Thue--Mahler case, amalgamate all Thue--Mahler
    # form files into a single file, TM/TMForms.csv, before removing redundant
    # forms.
    #
    # Parameters
    #     fileType
    #         The string denoting the input file, either "TM" or "XYZ2",
    #         corresponding to the Thue--Mahler or XYZ2 solver, respectively.
    # Returns
    #     outFile
    #         The output file containing, in each line, the equation to be
    #         resolved. In the Thue--Mahler case, this file is
    #         Data/${name}/TM/TMForms.csv, and in the XYZ2 case, this file is
    #         Data/${name}/XYZ2/XYZ2Forms.csv

    local formFile
    local logFile
    local N
    local F

    formFile="${Dir}/$1/$1Forms.csv"
    logFile="${Dir}/$1Status.txt"
    printf "This is the form and log file. \n"
    echo ${formFile}
    echo ${logFile}

    if [[ "$1" == "TM" ]]; then
	printf "We've got TM forms! \n"
	printf "time to sort TM stuff! \n"
	for N in "${list[@]}"; do
	    F="${Dir}/TM/${N}Forms.csv"
	    if [ -f "${F}" ]; then
		cat "${F}" >> "${formFile}"
		rm "${F}"
	    fi
	done
    elif [[ "$1" == "XYZ2" ]]; then
	printf "We've got XYZ2 forms! \n"
	printf "Do nothing here \n!"
    fi

    printf "Verifying nonempty ..." >> ${logFile}
    verifyNonEmpty "$1"
    printf "Done verifying nonempty.\n" >> ${logFile}

    printf "Removing redundant cases..."
    python Code/gatherRedundancy.py "${forms}" "$1"
    printf "Done.\n"
}

runTM() (

    # in a subshell so that we can run veriyNonEmpty without terminating the whole program
    # LEFT OFF HERE

    local TMLog
    local program

    printf "Inside TM. \n"
    TMLog="${Dir}/TMStatus.txt"
    touch "${TMLog}"

    # Generate all required Thue--Mahler forms in parallel, applying all
    # necessary local tests in the process. That is, run
    # Code/N2TME N '${TMDir}' > /dev/null
    # in parallel, with N an entry of ${conductors}, storing GNU parallel's
    # progress in the file ${Dir}/formJoblog.
    program="Code/N2TME {} '${TMDir}' > /dev/null"
    printf "Generating all required cubic forms for conductors in ${name}..." \
	   >> "${TMLog}"
    runParallel "${conductors}" formJoblog "${program}"
    printf "Done.\n" >> "${TMLog}"

    # Remove redundant Thue--Mahler equations.
    gatherRedundantForms "TM"

    printf "Done TM.\n"
)

runXYZ2() (

    local XYZ2Log

    printf "Inside XYZ2. \n"
    XYZ2Log="${Dir}/XYZ2Status.txt"
    touch "${XYZ2Log}"

    printf "Generating all required prime lists for conductors in ${name}..." \
	   >> "${XYZ2Dir}"
    echo "${conductors}" | factor > ${XYZ2Dir}/XYZ2Forms.csv
    printf "Done.\n" >> "${XYZ2Dir}"

    # Remove redundant prime lists.
    gatherRedundantForms "XYZ2"

    printf "Done XYZ2.\n"
)

main() {

    # Generates all elliptic curves of given conductor(s), taking as input a
    # single conductor N, a range of conductors [N1,...,N2], or, with the
    # flag -l, an arbitrary finite list of conductors.
    #
    # Parameters
    #     N1 [N2]
    #         A single conductor, N1, or the range [N1,...,N2].
    #     [-l N1] [N2...]: [OPTIONAL]
    #         A finite, arbitrary list of conductors to generate [N1,N2,...].
    # Returns
    #     Dir
    #         A subdirectory or Data storing all output, logs, errors, and
    #         elliptic curves.

    local conductors
    local program
    local TMForms

    printf "We're in main, not yet in subroutines...\n"

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

    printf "Entering TM .\n"
    runTM
    printf "Now out of TM and entering XYZ2 .\n"
    runXYZ2
    printf "Out and back and done!\n"
}
