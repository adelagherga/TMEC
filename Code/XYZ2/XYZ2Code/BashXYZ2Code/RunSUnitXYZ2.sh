#!/bin/bash
# RunSUnitXYZ2.sh

# INPUT:
#    N/A 
#        
# OUTPUT:
#     S_i.out:= the file containing the set of all solutions X + Y = Z^2 with good redeuction outside the set of primes S_i, where 
#        S_i:= [2,p_1,...,p_{n_i}], as listed in the file nSets.out
#        
# COMMENTS:
#    This algorithm is a bash function. It runs through all sets S as printed in nSets.out and executes SUnitXYZ2.m on each set. The output is printed in a file S.out, where S:= [ p_1, p_2, ..., p_n ] 
#    This alrogithm executes SUnitXYZ2.m on each set S in sequence, ie. on 1 core
#    To run in the background, enter the following in a terminal window: 
#        $ nohup $ABSOLUTE_PATH/SUnitCode/XYZ2Code/BashXYZ2Code/RunSUnitXYZ2.sh $ABSOLUTE_PATH &
#                
# EXAMPLE:
#    N/A


ABSOLUTE_PATH=$1        # selects the first item as the input, $ABSOLUTE_PATH; this step is necessary as the shell does not recognize the variable $ABSOLUTE_PATH otherwise
cd $ABSOLUTE_PATH/SUnitCode     # changes the directory to SUnitCode so that all called-upon Magma functions are loaded properly
while read line         # iterates through all lines (ie. sets S) of nSets.out
    do line2="$(echo -e "${line}" | tr -d '[:space:]')"         # removes white space between the characters of the set S; ie. writes S:= [ p_1, p_2, ..., p_n ] as [p_1,p_2,...,p_n]
                                                                # this step is necessary as Magma input does not allow for empty spaces
    magma S:=$line2 $ABSOLUTE_PATH/SUnitCode/XYZ2Code/BashToMagmaRunSUnit/SUnitXYZ2FromCommandLineInput.m > $ABSOLUTE_PATH/SUnitOutput/$line2.out       # runs the set SUnitXYZ2.m on the set 
                                                                                                                                                        # S:= [ p_1, p_2, ..., p_n ] and prints the results in a file S.out
done < $ABSOLUTE_PATH/SUnitOutput/nSets.out     # iterates through all lines (ie. sets S) of nSets.out
