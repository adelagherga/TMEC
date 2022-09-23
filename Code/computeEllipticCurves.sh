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

# Run Thue code in parallel.
# That is, for each line "set" of Data/Forms/ThueTestForms.csv, run
# magma set:="set" Code/computeEllipticCurvesThue.m &.
# The following code runs these jobs using GNU parallel, running no more than
# 20 jobs at once (-j20), and storing GNU parallel's progress in the logfile
# Data/ThueTest1Log (--joblog Data/ThueTest1Log).
cat Data/Forms/ThueTestForms.csv | parallel -j20 --joblog Data/ThueTest1Log magma set:={} Code/computeEllipticCurvesThue.m 2>&1

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
