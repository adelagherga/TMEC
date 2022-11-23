#!/bin/bash
# GenerateSets.sh

# INPUT:
#    m:= the upper bound on the product of the primes of S:= [2,p_1,...,p_{n-1}]
#        This input is entered in a bash terminal window 
#        
# OUTPUT:
#     nSets[i]:= [ S[i], D, #II_[1], #II_[2] ], where
#        nSets:= [S_1,...,S_n], an exaustive list of sets S_i such that 
#            S_i:= [2,p_1,...,p_{n_i}] such that 2*p_1*...*p_{n_i} <= m
#        D:= the values ( (p_1)^(b_1) )* /cdots *( (p_{n_i})^(b_{n_i}) ) where b_i in {0,1} corresponding to S_i
#        #II_[1] #II_[2]:= the number of elements in I, I_
#            Nb. if D == 1 or there are no split primes (ie. II_ == []), these values are both 0
#        
# COMMENTS:
#    The sets S_i are used to generate all solutions of X + Y = Z^2 with good reduction outside the primes of S_i, where
#        X,Y are (S_i)-units, and Z is a rational integer
#    It suffices to take the sets S_i to be maximal; that is, in the sense that S_i is not a subset of S_j for all i,j (i != j)
#        ie. if S_i is a subset of S_j, then the solutions of X + Y = Z^2 corresponding to S_i will be a subset of the solutions of X + Y = Z^2 corresponding to S_j
#    
#    This algorithm is a bash function. It opens Magma to compute all maximal sets S_i whose product is <= m, and outputs them into the text file nSets.out
#       To load the function, enter 'source $ABSOLUTE_PATH/SUnitCode/GenerateSCode/BashGenerateSetsCode/GenerateSets.sh' in a terminal window
#                
# EXAMPLE:
#    $ source $ABSOLUTE_PATH/SUnitCode/GenerateSCode/BashGenerateSetsCode/GenerateSets.sh
#    $ GenerateSetsXYZ2 100
#    $ less $ABSOLUTE_PATH/SUnitOutput/nSets.out
#        [[ 2, 3, 5 ], 1, 0, 0]
#        [[ 2, 3, 5 ], 2, 0, 0]
#        [[ 2, 3, 5 ], 3, 0, 0]
#        [[ 2, 3, 5 ], 5, 0, 0]
#        [[ 2, 3, 5 ], 6, 1, 0]
#        [[ 2, 3, 5 ], 10, 1, 0]
#        [[ 2, 3, 5 ], 15, 0, 0]
#        [[ 2, 3, 5 ], 30, 0, 0]
#        [[ 2, 3, 7 ], 1, 0, 0]
#        [[ 2, 3, 7 ], 2, 0, 0]
#        [[ 2, 3, 7 ], 3, 0, 0]
#        [[ 2, 3, 7 ], 6, 1, 0]
#        [[ 2, 3, 7 ], 7, 1, 0]
#        [[ 2, 3, 7 ], 14, 1, 0]
#        [[ 2, 3, 7 ], 21, 1, 0]
#        [[ 2, 3, 7 ], 42, 0, 0]
#        [[ 2, 3, 11 ], 1, 0, 0]
#        [[ 2, 3, 11 ], 2, 0, 0]
#        [[ 2, 3, 11 ], 3, 0, 0]
#        [[ 2, 3, 11 ], 6, 1, 0]
#        [[ 2, 3, 11 ], 11, 1, 0]
#        [[ 2, 3, 11 ], 22, 1, 0]
#        [[ 2, 3, 11 ], 33, 1, 0]
#        [[ 2, 3, 11 ], 66, 1, 0]
#        [[ 2, 3, 13 ], 1, 0, 0]
#        [[ 2, 3, 13 ], 2, 0, 0]
#        [[ 2, 3, 13 ], 3, 0, 0]
#        [[ 2, 3, 13 ], 6, 1, 0]
#        [[ 2, 3, 13 ], 13, 1, 0]
#        [[ 2, 3, 13 ], 26, 1, 0]
#        [[ 2, 3, 13 ], 39, 1, 0]
#        [[ 2, 3, 13 ], 78, 0, 0]
#        [[ 2, 5, 7 ], 1, 0, 0]
#        [[ 2, 5, 7 ], 2, 0, 0]
#        [[ 2, 5, 7 ], 5, 0, 0]
#        [[ 2, 5, 7 ], 7, 1, 0]
#        [[ 2, 5, 7 ], 10, 1, 0]
#        [[ 2, 5, 7 ], 14, 1, 0]
#        [[ 2, 5, 7 ], 35, 0, 0]
#        [[ 2, 5, 7 ], 70, 1, 0]
#        [[ 2, 17 ], 1, 0, 0]
#        [[ 2, 17 ], 2, 0, 0]
#        [[ 2, 17 ], 17, 1, 0]
#        [[ 2, 17 ], 34, 2, 0]
#        [[ 2, 17 ], 34, 1, 0]
#        [[ 2, 17 ], 34, 1, 0]
#        [[ 2, 17 ], 34, 1, 1]
#        [[ 2, 19 ], 1, 0, 0]
#        [[ 2, 19 ], 2, 0, 0]
#        [[ 2, 19 ], 19, 2, 0]
#        [[ 2, 19 ], 19, 1, 0]
#        [[ 2, 19 ], 19, 1, 0]
#        [[ 2, 19 ], 19, 1, 1]
#        [[ 2, 19 ], 38, 0, 0]
#        [[ 2, 23 ], 1, 0, 0]
#        [[ 2, 23 ], 2, 0, 0]
#        [[ 2, 23 ], 23, 0, 0]
#        [[ 2, 23 ], 46, 2, 0]
#        [[ 2, 23 ], 46, 1, 0]
#        [[ 2, 23 ], 46, 1, 0]
#        [[ 2, 23 ], 46, 1, 1]
#        [[ 2, 29 ], 1, 0, 0]
#        [[ 2, 29 ], 2, 0, 0]
#        [[ 2, 29 ], 29, 1, 0]
#        [[ 2, 29 ], 58, 1, 0]
#        [[ 2, 31 ], 1, 0, 0]
#        [[ 2, 31 ], 2, 0, 0]
#        [[ 2, 31 ], 31, 2, 0]
#        [[ 2, 31 ], 31, 1, 0]
#        [[ 2, 31 ], 31, 1, 0]
#        [[ 2, 31 ], 31, 1, 1]
#        [[ 2, 31 ], 62, 0, 0]
#        [[ 2, 37 ], 1, 0, 0]
#        [[ 2, 37 ], 2, 0, 0]
#        [[ 2, 37 ], 37, 1, 0]
#        [[ 2, 37 ], 74, 1, 0]
#        [[ 2, 41 ], 1, 0, 0]
#        [[ 2, 41 ], 2, 0, 0]
#        [[ 2, 41 ], 41, 2, 0]
#        [[ 2, 41 ], 41, 1, 0]
#        [[ 2, 41 ], 41, 1, 0]
#        [[ 2, 41 ], 41, 1, 1]
#        [[ 2, 41 ], 82, 1, 0]
#        [[ 2, 43 ], 1, 0, 0]
#        [[ 2, 43 ], 2, 0, 0]
#        [[ 2, 43 ], 43, 1, 0]
#        [[ 2, 43 ], 86, 1, 0]
#        [[ 2, 47 ], 1, 0, 0]
#        [[ 2, 47 ], 2, 0, 0]
#        [[ 2, 47 ], 47, 0, 0]
#        [[ 2, 47 ], 94, 2, 0]
#        [[ 2, 47 ], 94, 1, 0]
#        [[ 2, 47 ], 94, 1, 0]
#        [[ 2, 47 ], 94, 1, 1]


GenerateSetsXYZ2 () {
    mkdir $ABSOLUTE_PATH/SUnitOutput    # generates SUnitOutput directory
    m=$1        # selects the first item as the input, m        
    cd $ABSOLUTE_PATH/SUnitCode         # changes the directory to SUnitCode so that all called-upon Magma functions are loaded properly
    echo "Generating all sets of primes whose product is less than m = $m..."
    nohup $MAGMA -b m:=$m $ABSOLUTE_PATH/SUnitCode/GenerateSCode/BashToMagmaGenerateS/GenerateSetsXYZ2.m > $ABSOLUTE_PATH/SUnitOutput/nSets.out 2>&1 &  
                # generates all sets S; prints to nSets.out
    wait        # waits until all processes have finished
    sed '/2,/!d' $ABSOLUTE_PATH/SUnitOutput/nSets.out > temp.txt ; mv temp.txt $ABSOLUTE_PATH/SUnitOutput/nSets.out     # deletes unnecesary lines so that output file contains only the nSets and no additional text

    sed 's/</[/g' $ABSOLUTE_PATH/SUnitOutput/nSets.out > temp.txt ; mv temp.txt $ABSOLUTE_PATH/SUnitOutput/nSets.out    # changes "<" character to "[" so that the set may be input into magma startup without error
    sed 's/>/]/g' $ABSOLUTE_PATH/SUnitOutput/nSets.out > temp.txt ; mv temp.txt $ABSOLUTE_PATH/SUnitOutput/nSets.out    # changes ">" character to "]" so that the set may be input into magma startup without error
    sed 's/ //g' $ABSOLUTE_PATH/SUnitOutput/nSets.out > temp.txt ; mv temp.txt $ABSOLUTE_PATH/SUnitOutput/nSets.out     # removes all whitespace between characters on each line
    
    cd $ABSOLUTE_PATH   # returns to parent directory
    return 0
}
