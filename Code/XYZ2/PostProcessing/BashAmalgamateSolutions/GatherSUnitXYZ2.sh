#!/bin/bash
# GatherSUnitXYZ2.sh

# INPUT:
#    N/A 
#        
# OUTPUT:
#     S_i.out:= the file containing the set of all solutions X + Y = Z^2 with good redeuction outside the set of primes S_i, where 
#        S_i:= [2,p_1,...,p_{n_i}], as listed in the file nSets.out
#        
# COMMENTS:
#    This algorithm is a bash function. It runs through all sets S.out where S:= [ p_1, p_2, ..., p_n ] as printed in nSets.out and output by SUnitXYZ2.m  
#    This alrogithm sorts through the files S.out and amalgamates all corresponding solutions into one file, S=[p_1,...,p_n]
#    If errors exist in any of the files, an error is printed, along with the filename S.out
#    
# EXAMPLE:
#    N/A

mkdir $ABSOLUTE_PATH/SUnitOutput/Errors         # generates a folder to gather all files containing errors
mkdir $ABSOLUTE_PATH/SUnitOutput/SetSolutions   # generates a folder to gather all files with amalgamated solutions correspoinding to S=[p_1,...,p_s]
# mkdir $ABSOLUTE_PATH/SUnitOutput/CompletedJobs     # generates a folder to gather all files as originally output in the folder $ABSOLUTE_PATH/SUnitOutput as output by SUnitXYZ2.m


FILES=$ABSOLUTE_PATH/SUnitOutput/*      # sets the variable $FILES as the set of all files in the folder $ABSOLUTE_PATH/SUnitOutput; this includes $ABSOLUTE_PATH/SUnitOutput/nSets.out

for S in $FILES; do     # runs through all files in $FILES
    if [ "$S" != "$ABSOLUTE_PATH/SUnitOutput/nSets.out" ] && [ "$S" != "$ABSOLUTE_PATH/SUnitOutput/SetSolutions" ] && [ "$S" != "$ABSOLUTE_PATH/SUnitOutput/Errors" ]; then
            # $S runs through all files except $ABSOLUTE_PATH/SUnitOutput/nSets.out, and $ABSOLUTE_PATH/SUnitOutput/SetSolutions, adn $ABSOLUTE_PATH/SUnitOutput/Errors  
        FileS="${S##*/}"        # truncates the name of the file S=[p_1,...,p_s] without the parent folders; ie. removes "$ABSOLUTE_PATH/SUnitOutput"
        echo $FileS
        E=$(grep "error" $S)    # sets the variable $E to the line containing the string "error" in $S, if it exists; ie. looks for errors in each file $S

        if [[ $E = *[!\ ]* ]]; then     # if $E is not an empty string; ie. if $S contains an error, meaning that there is an error in the set associated to $S
            echo "Errors abound! Oh nooooo"  
            echo $FileS  # prints the name of the file containing the error without the parent folders
            mv $S $ABSOLUTE_PATH/SUnitOutput/Errors/$FileS      # copies the file containing an error to the folder "Error"
            # mv $S $ABSOLUTE_PATH/SUnitOutput/CompletedJobs/$FileS        # moves $S to the folder Jobs 
        else    # if there are no errors 
            s=$(grep -n '^\[$' $S | grep -Eo '^[^:]+')  # computes the line of text in $S containing only the character "["; the set of solutions S lie on line $s+1
            f=$(grep -n '^\]$' $S | grep -Eo '^[^:]+')  # computes the line of text in $S containing only the character "]"; the set of solutions S lie on line $f-1
            nS=$(sed -n "$((s-1)) p " $S)       # finds and generates the string "S = [...], D = ..."; this string lies on the line $s-1. 
            nS=${nS%D*}         # removes everything on the line nS after "D = ..."
            nS=${nS%,*}         # removes the last "," after the set S; nS now represents the set S that is being iterated over in $S
            nS="$(echo -e "${nS}" | tr -d '[:space:]')"         # removes white space between the characters of the set S; ie. writes S:= [ p_1, p_2, ..., p_n ] as [p_1,p_2,...,p_n]
                                                                    # this step is necessary as Magma input does not allow for empty spaces 
            sed -n "$((s+1))","$((f-1)) p" $S >> $ABSOLUTE_PATH/SUnitOutput/SetSolutions/$nS    # moves all relevant text (ie. solutions) to the file $nS (ie. S=[p_1,...,p_n])
            # mv $S $ABSOLUTE_PATH/SUnitOutput/CompletedJobs/$FileS  # moves $S to the folder Jobs
        fi
    fi
done


FILES=$ABSOLUTE_PATH/SUnitOutput/SetSolutions/*         # sets the variable $FILES as the set of all files in the folder $ABSOLUTE_PATH/SUnitOutput/SetSolutions
for S in $FILES; do
    FILENAME="${S##*/}"         # truncates the name of the file S=[p_1,...,p_s] without the parent folders; ie. removes "$ABSOLUTE_PATH/SUnitOutput/SetSolutions"
    echo $FILENAME
    sed 's/\,$//'  $ABSOLUTE_PATH/SUnitOutput/SetSolutions/$FILENAME > temp.txt ; mv temp.txt $ABSOLUTE_PATH/SUnitOutput/SetSolutions/$FILENAME   
        # removes all the characters "," at the end of each line of the file S=[FILENAME], where it exists
    sed 's/$/\,/'  $ABSOLUTE_PATH/SUnitOutput/SetSolutions/$FILENAME > temp.txt ; mv temp.txt $ABSOLUTE_PATH/SUnitOutput/SetSolutions/$FILENAME
        # adds the character "," at the end of each line of the file S=[FILENAME]; ensures that each line is separated by "," for processing in Magma
    sort -g -k 2 $ABSOLUTE_PATH/SUnitOutput/SetSolutions/$FILENAME | uniq > temp.txt ; mv temp.txt $ABSOLUTE_PATH/SUnitOutput/SetSolutions/$FILENAME
        # sorts the solutions in the file S=[FILENAME] and removes duplicate solutions
    sed '1s/^/Sol:=[/' $ABSOLUTE_PATH/SUnitOutput/SetSolutions/$FILENAME > temp.txt ;  mv temp.txt $ABSOLUTE_PATH/SUnitOutput/SetSolutions/$FILENAME
        # adds "Sol:=[" at the start of each file $FILENAME; necessary for parsing by MAGMA to generate elliptic curves
    sed '$ s/.$/];/' $ABSOLUTE_PATH/SUnitOutput/SetSolutions/$FILENAME > temp.txt ;  mv temp.txt $ABSOLUTE_PATH/SUnitOutput/SetSolutions/$FILENAME
        # adds "];" at the end of each file $FILENAME; necessary for parsing by MAGMA to generate elliptic curves
done