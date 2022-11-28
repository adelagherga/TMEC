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
    #     ECDir
    #         The directory Data/${name}/EllipticCurves, storing all resulting
    #         elliptic curves.
    #     ECDir/${N}.csv
    #         A file storing all elliptic curves of conductor N, generated for
    #         each N in list, in the format N [a1,a2,a3,a4,a6].
    #     JLDir
    #         The directory Data/${name}Joblogs, storing all files detailing GNU
    #         parallel's progress.

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
    ECDir="${Dir}/EllipticCurves"
    JLDir="${Dir}/Joblogs"
    mkdir "${Dir}"
    mkdir "${ECDir}"
    mkdir "${JLDir}"

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

    echo "$1" | parallel -j20 --joblog ${JLDir}/"$2" "$3"
}

verifyNonEmpty() {

    # Verify whether the files TMForms.csv and XYZ2Forms are empty. When true,
    # terminate the corresponding program.
    #
    # Parameters
    #     fileType
    #         The string denoting the input file, either "TM" or "XYZ2",
    #         corresponding to the Thue--Mahler or XYZ2 solver, respectively.

    local formFile
    local logFile
    local msg

    formFile="${Dir}/$1/$1Forms.csv"
    logFile="${Dir}/$1Status.txt"
    if [ ! -s "${formFile}" ]; then
	rm -r "${Dir}/$1/Outfiles"
	msg="Finished computing all elliptic curves of conductor ${name}"
	msg="${msg} corresponding to the $1 solver.\n"
	printf "${msg}" >> "${logFile}"
	exit 0
    fi
}

gatherRedundantCases() {

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

    if [[ "$1" == "TM" ]]; then
	for N in "${list[@]}"; do
	    F="${Dir}/TM/${N}Forms.csv"
	    if [ -f "${F}" ]; then
		cat "${F}" >> "${formFile}"
		rm "${F}"
	    fi
	done
    fi

    verifyNonEmpty "$1"

    printf "Removing redundant cases..." >> "${logFile}"
    python3 Code/gatherRedundancy.py "${formFile}" "$1"
    printf "Done.\n" >> "${logFile}"
}

amalgamateSUnitFiles() {

    # Amalgamate all S-unit equations, corresponding to either the Thue--Mahler
    # solver or the XYZ2 solver, into a single document. Clean up any additional
    # temporary files in the process and forward any magma errors to the files
    # TM/Errors.txt or XYZ2/Errors.txt, respectively.
    #
    # Parameters
    #     fileType
    #         The string denoting the input file, either "TM" or "XYZ2",
    #         corresponding to the Thue--Mahler or XYZ2 solver, respectively.
    # Returns
    #     outFile
    #         The output file containing, in each line, the S-unit equation to be
    #         resolved. In the Thue--Mahler case, these S-unit equations pertain
    #         to the optimal Thue--Mahler form. In this case, this file is
    #         Data/${name}/TM/TMForms.csv, and in the XYZ2 case, this file is
    #         Data/${name}/XYZ2/XYZ2Forms.csv

    local formFile
    local logFile

    local line
    local F1
    local F2
    local rerun

    formFile="${Dir}/$1/$1Forms.csv"
    logFile="${Dir}/$1Status.txt"

    while IFS= read -r line; do
	F1="${Dir}/$1/${line}.csv"
	F2="${Dir}/$1/${line}tmp.txt"
	if [ -f "${F1}" ]; then
	    cat "${F1}" >> "${Dir}/$1/tmp$1Forms.csv"
	    rm "${F1}"
	fi
	if grep -q "error" "${F2}"; then
	    rerun="magma -b set:=${line} dir:='${Dir}/$1' "
	    if [[ "$1" == "TM" ]]; then
		rerun="${rerun} Code/TM/optimalForm.m 2>&1"
	    else
		rerun="${rerun} Code/XYZ2/seperatelForm.m 2>&1"
	    fi
	    echo "${rerun}" >> "${Dir}/$1/Errors.txt"
	fi
	rm "${F2}"
    done < "${formFile}"
    mv "${Dir}/$1/tmp$1Forms.csv" "${formFile}"

    verifyNonEmpty "$1"
}

moveCurves() {

    # For each line of inFile, extract the conductor N and copy the line to the
    # file ${ECDir}/${N}.csv.
    #
    # Parameters
    #     inFile
    #         The output file from Code/TM/runTMEC.m or Code/XYZ2/runXYZ2EC.m
    #         containing, for each line, the elliptic curve obtained by solving
    #         the Thue--Mahler or XYZ2 S-unit equation in the filename.

    local line
    local N

    while IFS= read -r line; do
	N=$(echo "${line}" | cut -d' ' -f -1)
	echo "${line}" >> "${ECDir}/${N}.csv"
    done < "$1"
}

amalgamateSolFiles() {

    # Amalgamate all logfiles pertaining to either the same Thue--Mahler
    # equation or the same XYZ2 equation, and forward any magma errors to the
    # files TM/Errors.txt or XYZ2/Errors.txt, respectively. Finally, populate
    # the elliptic curve files ${ECDir}/${N}.csv.
    #
    # Parameters
    #     fileType
    #         The string denoting the input file, either "TM" or "XYZ2",
    #         corresponding to the Thue--Mahler or XYZ2 solver, respectively.
    # Returns
    #     logFile
    #         The file containing the amalgamated logfiles pertaining to the same
    #         Thue--Mahler or XYZ2 equation. In the Thue--Mahler case, this file
    #         is Data/${name}/TM/${line}Log.txt, and in the XYZ2 case, this file
    #         is Data/${name}/XYZ2/${line}Log.txt.

    local formFile
    local logFile
    local line
    local form
    local out
    local log
    local allLog
    local rerun

    formFile="${Dir}/$1/$1Forms.csv"
    logFile="${Dir}/$1Status.txt"

    while IFS= read -r line; do
	form="$(echo ${line} | sed 's|\(.*\),\[.*|\1|')"
	if [ "${line}" != "${form}" ]; then
	    out="${Dir}/$1/Outfiles/${line}Out.csv"
	    log="${Dir}/$1/Logfiles/${line}Log.txt"
	    allLog="${Dir}/$1/Logfiles/${form}Log.txt"
	    if grep -q "error" "${log}"; then
		rerun="magma -b set:='${line}' dir:='${Dir}/$1'"
		rerun="${rerun} Code/$1/run$1EC.m 2>&1"
		echo "${rerun}" >> "${Dir}/$1/Errors.txt"
	    fi
	    if [ -f "${out}" ]; then
		moveCurves "${out}"
	    fi
	    cat "${log}" >> "${allLog}"
	    rm "${log}"
	fi
    done < "${formFile}"
    rm -r "${Dir}/$1/Outfiles"
}

runTM() (

    # Generate all necessary and optimal Thue--Mahler forms, without redundancy,
    # and solve all respective S-unit equations. This process is launched in a
    # subshell to run in parallel with the XYZ2 solver.
    #
    # Returns
    #     Dir/TM
    #         The directory Data/${name}/TM, storing all solutions and data
    #         related to the Thue--Mahler computation.

    local logFile
    local program
    local TMForms

    logFile="${Dir}/TMStatus.txt"
    mkdir "${Dir}/TM"
    mkdir "${Dir}/TM/Outfiles"
    mkdir "${Dir}/TM/Logfiles"
    touch "${Dir}/TM/Errors.txt"

    # Generate all required Thue--Mahler forms in parallel, applying all
    # necessary local tests in the process. That is, run
    # Code/N2TME N '${Dir}/TM' > /dev/null
    # in parallel, with N an entry of ${conductors}, storing GNU parallel's
    # progress in the file ${JLDir}/formJoblog.
    program="Code/N2TME {} '${Dir}/TM' > /dev/null"
    printf "Generating all required cubic forms for conductors in ${name}..." \
	   >> "${logFile}"
    runParallel "${conductors}" formJoblog "${program}"
    printf "Done.\n" >> "${logFile}"

    # Remove redundant Thue--Mahler equations.
    gatherRedundantCases "TM"

    # Generate optimal Thue--Mahler forms and all S-unit equations in parallel.
    # That is, run
    # magma -b set:=<line> dir:=${Dir}/TM Code/TM/optimalForm.m 2>&1
    # in parallel, for each <line> of ${TMForms}, storing GNU parallel's
    # progress in the file ${JLDir}/optimalJoblog.
    TMForms=$(cat ${Dir}/TM/TMForms.csv)
    program="magma -b set:={} dir:='${Dir}/TM' Code/TM/optimalForm.m 2>&1"
    printf "Generating optimal GL2(Z)-equivalent cubic forms..." \
	   >> "${logFile}"
    runParallel "${TMForms}" optimalJoblog "${program}"
    printf "Done.\n" >> "${logFile}"

    # Clean up directory and amalgamate all Thue--Mahler S-unit equations into
    # the single file ${Dir}/TM/TMForms.csv.
    amalgamateSUnitFiles "TM"

    # Solve all Thue--Mahler S-unit equations in parallel.
    # That is, run
    # magma -b set:=<line> dir:=${Dir}/TM Code/TM/runTMEC.m 2>&1
    # in parallel, for each <line> of ${TMForms}, storing GNU parallel's
    # progress in the file ${JLDir}/TMJoblog.
    TMForms=$(cat ${Dir}/TM/TMForms.csv)
    program="magma -b set:={} dir:='${Dir}/TM' Code/TM/runTMEC.m 2>&1"
    printf "Solving the Thue--Mahler equations..." >> "${logFile}"
    runParallel "${TMForms}" TMJoblog "${program}"
    printf "Done.\n" >> "${logFile}"
    exit 0
)

runXYZ2() (

    # Generate all necessary prime lists for the XYZ2 solver, without redundancy,
    # and solve all respective S-unit equations. This process is launched in a
    # subshell to run in parallel with the Thue--Mahler solver.
    #
    # Returns
    #     Dir/XYZ2
    #         The directory Data/${name}/XYZ2, storing all solutions and data
    #         related to the XYZ2 computation.

    local logFile
    local program
    local XYZ2Forms

    logFile="${Dir}/XYZ2Status.txt"
    mkdir "${Dir}/XYZ2"
    mkdir "${Dir}/XYZ2/Outfiles"
    mkdir "${Dir}/XYZ2/Logfiles"
    touch "${Dir}/XYZ2/Errors.txt"

    # Generate all required prime lists for the XYZ2 solver by determining the
    # prime factorization of all relevant conductors.
    printf "Generating all required prime lists for conductors in ${name}..." \
	   >> "${logFile}"
    echo "${conductors}" | factor > "${Dir}/XYZ2/XYZ2Forms.csv"
    printf "Done.\n" >> "${logFile}"

    # Remove redundant prime lists.
    gatherRedundantCases "XYZ2"

    # Generate all XYZ2 S-unit equations in parallel.
    # That is, run
    # magma -b set:=<line> dir:=${Dir}/XYZ2 Code/XYZ2/seperateForm.m 2>&1
    # in parallel, for each <line> of ${XYZ2Forms}, storing GNU parallel's
    # progress in the file ${JLDir}/seperateJoblog.
    XYZ2Forms=$(cat ${Dir}/XYZ2/XYZ2Forms.csv)
    program="magma -b set:={} dir:='${Dir}/XYZ2' Code/XYZ2/seperateForm.m 2>&1"
    printf "Preparing prime list data for parallelization..." >> "${logFile}"
    runParallel "${XYZ2Forms}" seperateJoblog "${program}"
    printf "Done.\n" >> "${logFile}"

    # Clean up directory and amalgamate all XYZ2 S-unit equations into the single
    # file ${Dir}/XYZ2/XYZ2Forms.csv.
    amalgamateSUnitFiles "XYZ2"

    # Solve all XYZ2 S-unit equations in parallel.
    # That is, run
    # magma -b set:=<line> dir:=${Dir}/XYZ2 Code/TM/runTMEC.m 2>&1
    # in parallel, for each <line> of ${XYZ2Forms}, storing GNU parallel's
    # progress in the file ${JLDir}/XYZ2Joblog.
    XYZ2Forms=$(cat ${Dir}/XYZ2/XYZ2Forms.csv)
    program="magma -b set:={} dir:='${Dir}/XYZ2' Code/XYZ2/runXYZ2EC.m 2>&1"
    printf "Solving the XYZ2 S-unit equations..." >> "${logFile}"
    runParallel "${XYZ2Forms}" XYZ2Joblog "${program}"
    printf "Done.\n" >> "${logFile}"
    exit 0
)

sortCurves() {

    # Sort all elliptic curve files and generate a single file containing all
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
    local missing

    getNList "$@"
    generateDirectories

    # Generate all elliptic curves of j-invariant 0 in conductors list, in
    # parallel. That is, run
    # Code/CurvesNj0 N > ${Dir}/EllipticCurves/N.csv
    # in parallel, with N an entry of ${conductors}, storing GNU parallel's
    # progress in the file ${JLDir}/j0Joblog.
    conductors=$(printf '%s\n' "${list[@]}")
    program="Code/CurvesNj0 {} > '${ECDir}'/{}.csv"
    printf "Generating all j-invariant 0 curves for conductors in ${name}..."
    runParallel "${conductors}" j0Joblog "${program}"
    printf "Done.\n"

    printf "Generating all remaining curves for conductors in ${name}..."
    runTM &
    runXYZ2 &
    wait
    printf "Done.\n"

    # If TM/TMForms.csv is not empty, clean up the ${Dir}/TM directory and
    # amalgamate all logfiles pertaining to the same Thue--Mahler form,
    # populating the elliptic curve files ${ECDir}/${N}.csv in the process.
    if [ -s "${Dir}/TM/TMForms.csv" ]; then
	amalgamateSolFiles "TM"
    fi
    # If XYZ2/XYZ2Forms.csv is not empty, clean up the ${Dir}/XYZ2 directory and
    # amalgamate all logfiles pertaining to the same XYZ2 equation,
    # populating the elliptic curve files ${ECDir}/${N}.csv in the process.
    if [ -s "${Dir}/XYZ2/XYZ2Forms.csv" ]; then
	amalgamateSolFiles "XYZ2"
    fi

    # Sort all elliptic curve files and generate a single file containing all
    # curves.
    printf "Sorting all elliptic curves..."
    sortCurves
    printf "Done.\n"
    printf "Finished computing all elliptic curves of conductor ${name}.\n"

    # If there are no errors in the Thue--Mahler and XYZ2 computation, run a
    # comparison with the LMFDB, where possible, and print any missing curves in
    # ${ECDir}/MissingCurves.csv
    missing="${Dir}/EllipticCurves/MissingCurves.csv"
    if [ ! -s "${Dir}/TM/Errors.txt" ] && [ ! -s "${Dir}/TM/Errors.txt" ]; then
	printf "Running comparison with LMFDB, where possible..."
	python3 Code/compareLMFDB.py "${Dir}"
	printf "Done.\n"
	if [ ! -s "${missing}" ] && [ -e "${missing}" ]; then
	    printf "There are no missing elliptic curves.\n"
	elif [ -s "${missing}" ]; then
	    printf "Error: there are missing elliptic curves.\n"
	fi
    else
	printf "Errors in Thue--Mahler or XYZ2 code; "
	printf "cannot run LMFDB comparison.\n"
    fi
}

main "$@"
