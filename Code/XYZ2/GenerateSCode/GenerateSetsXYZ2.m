// GenerateSetsXYZ2.m

/*
INPUT:
    N/A
        
OUTPUT:
    nSets:= [S_1,...,S_n], where nSets is an exaustive list of sets S_i such that 
        S_i:= [2,p_1,...,p_{n_i}] such that 2*p_1*...*p_{n_i} <= m
        
COMMENTS:
    This algorithm is an intermediary procedure between bash and magma, to output all maximal sets S_i whose product is <= m in a text file

EXAMPLE:
    N/A
       
        
*/

load "./GenerateSCode/GenerateSFunctions/nSetsXYZ2.m";

m:= StringToInteger(m);         // converts bash input to an integer in magma

Set:= nSetsXYZ2(m);     // generates all sets of primes whose product is <= m
print Set;      // prints output to a text file
exit;   // exits magma