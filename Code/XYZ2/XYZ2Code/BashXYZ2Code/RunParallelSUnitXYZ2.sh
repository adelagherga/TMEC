#!/bin/bash
# RunParallelSUnitXYZ2.sh

# INPUT:
#    N/A 
#        
# OUTPUT:
#     S_i.out:= the file containing the set of all solutions X + Y = Z^2 with good redeuction outside the set of primes S_i, where 
#        S_i:= [2,p_1,...,p_{n_i}], as listed in the file nSets.out
#        
# COMMENTS:
#    This algorithm is a bash function. It runs through all sets S as printed in nSets.out and executes SUnitXYZ2.m on each set. The output is printed in a file S.out, where S:= [ p_1, p_2, ..., p_n ] 
#    This alrogithm executes SUnitXYZ2.m on each set S in parallel
#    To run in the background, enter the following in a terminal window: 
#        $ nohup $ABSOLUTE_PATH/SUnitCode/XYZ2Code/BashXYZ2Code/RunSUnitXYZ2.sh $ABSOLUTE_PATH &
#                
# EXAMPLE:
#    N/A


ABSOLUTE_PATH=$1        # selects the first item as the input, $ABSOLUTE_PATH; this step is necessary as the shell does not recognize the variable $ABSOLUTE_PATH otherwise
LINE=$2         # selects the second item as the input, $LINE
MAGMA=$3        # selects the third item as the input, $MAGMA; this step is necessary as the shell does not recognize the variable $MAGMA

cd $ABSOLUTE_PATH/SUnitCode     # changes the directory to SUnitCode so that all called-upon Magma functions are loaded properly
line=$(sed -n ${LINE}p $ABSOLUTE_PATH/SUnitOutput/nSets.out)    # selects the indicated line, LINE, from nSets.out (ie. the indicated set S)
line2="$(echo -e "${line}" | tr -d '[:space:]')"        # removes white space between the characters of the set S; ie. writes S:= [ p_1, p_2, ..., p_n ] as [p_1,p_2,...,p_n]
                                                        # this step is necessary as Magma input does not allow for empty spaces
$MAGMA S:=$line2 $ABSOLUTE_PATH/SUnitCode/XYZ2Code/BashToMagmaRunSUnit/SUnitXYZ2FromCommandLineInput.m > $ABSOLUTE_PATH/SUnitOutput/$line2.out      
        # runs the set SUnitXYZ2.m on the set S:= [ p_1, p_2, ..., p_n ] and prints the results in a file S.out
# rm -r $ABSOLUTE_PATH/Requested  # removes unnecessary empty files created by GNU parallel
