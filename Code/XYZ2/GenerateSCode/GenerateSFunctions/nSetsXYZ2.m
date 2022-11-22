// nSetsXYZ2.m

/*
INPUT:
    m:= the upper bound on the product of the primes of S:= [2,p_1,...,p_{n-1}]
        
OUTPUT:
    nSets[i]:= < S_i, D, Ind >, I, I_, where
        nSets:= [S_1,...,S_n] is an exaustive list of sets S_i such that 
            S_i:= [2,p_1,...,p_{n_i}] such that 2*p_1*...*p_{n_i} <= m
        D:= the values ( (p_1)^(b_1) )* /cdots *( (p_{n_i})^(b_{n_i}) ) where b_i in {0,1} corresponding to S_i
        Ind:= the index of II_ where
            II_:= [I,I_], where, assuming the primes p_i of S_i which split in Q(sqrt(D)) are p_1,...,p_t, for some t in [0,s],
                I:= { i | 1 <= i <= t, a_i > b_i} and 
                I_:= { i | 1 <= i <= t, a_i < b_i}    
            Nb. if D == 1 or there are no split primes (ie. II_ == []), Ind is 0
        
COMMENTS:
    The sets S_i are used to generate all solutions of X + Y = Z^2 with good reduction outside the primes of S_i, where
        X,Y are (S_i)-units, and Z is a rational integer
    It suffices to take the sets S_i to be maximal; that is, in the sense that S_i is not a subset of S_j for all i,j (i != j)
        ie. if S_i is a subset of S_j, then the solutions of X + Y = Z^2 corresponding to S_i will be a subset of the solutions of X + Y = Z^2 corresponding to S_j
    
    This algorithm outputs all maximal sets S_i whose product is <= m

EXAMPLE:
    > m:= 100;
    > nSetsXYZ2(m);
    <[ 2, 3, 5 ], 1, 0>, [], []
    <[ 2, 3, 5 ], 2, 0>, [], []
    <[ 2, 3, 5 ], 3, 0>, [], []
    <[ 2, 3, 5 ], 5, 0>, [], []
    <[ 2, 3, 5 ], 6, 1>, [ 5 ], []
    <[ 2, 3, 5 ], 10, 1>, [ 3 ], []
    <[ 2, 3, 5 ], 15, 0>, [], []
    <[ 2, 3, 5 ], 30, 0>, [], []
    <[ 2, 3, 7 ], 1, 0>, [], []
    <[ 2, 3, 7 ], 2, 1>, [ 7 ], []
    <[ 2, 3, 7 ], 3, 0>, [], []
    <[ 2, 3, 7 ], 6, 0>, [], []
    <[ 2, 3, 7 ], 7, 1>, [ 3 ], []
    <[ 2, 3, 7 ], 14, 0>, [], []
    <[ 2, 3, 7 ], 21, 0>, [], []
    <[ 2, 3, 7 ], 42, 0>, [], []
    <[ 2, 3, 11 ], 1, 0>, [], []
    <[ 2, 3, 11 ], 2, 0>, [], []
    <[ 2, 3, 11 ], 3, 1>, [ 11 ], []
    <[ 2, 3, 11 ], 6, 0>, [], []
    <[ 2, 3, 11 ], 11, 0>, [], []
    <[ 2, 3, 11 ], 22, 1>, [ 3 ], []
    <[ 2, 3, 11 ], 33, 1>, [ 2 ], []
    <[ 2, 3, 11 ], 66, 0>, [], []
    <[ 2, 3, 13 ], 1, 0>, [], []
    <[ 2, 3, 13 ], 2, 0>, [], []
    <[ 2, 3, 13 ], 3, 1>, [ 13 ], []
    <[ 2, 3, 13 ], 6, 0>, [], []
    <[ 2, 3, 13 ], 13, 1>, [ 3 ], []
    <[ 2, 3, 13 ], 26, 0>, [], []
    <[ 2, 3, 13 ], 39, 0>, [], []
    <[ 2, 3, 13 ], 78, 0>, [], []
    <[ 2, 5, 7 ], 1, 0>, [], []
    <[ 2, 5, 7 ], 2, 1>, [ 7 ], []
    <[ 2, 5, 7 ], 5, 0>, [], []
    <[ 2, 5, 7 ], 7, 0>, [], []
    <[ 2, 5, 7 ], 10, 0>, [], []
    <[ 2, 5, 7 ], 14, 1>, [ 5 ], []
    <[ 2, 5, 7 ], 35, 0>, [], []
    <[ 2, 5, 7 ], 70, 0>, [], []

*/

load "./GenerateSCode/GenerateSProcedures/RecursiveSets.m";
load "./XYZ2Code/MagmaXYZ2Functions/ExponentsXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/DecompositionOfPrimesXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/SplitPrimePropertiesXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/IIPrimeXYZ2.m";


function nSetsXYZ2(m)
    nSets:= [];         // stores sets S_i, where S_i:= [2, p_1, ..., p_{n_i}] and 2*p_1*...*p_{n_i} <= m

    S:= [2];    // generates the initial set S; this is the largest set S:= [2, p_1, ..., p_n] such that 2*p_1*...*p_n <= m 
    P:= 3;
    while ( (&*S)*P le m ) do              
        Append(~S,P);   // appends P to S if [ prod_{s \in S} s ]*P <= m
        P:= NextPrime(P);       // generates the next largest prime not already in S
    end while;

    Append(~nSets, S);  // store S in nSets
    n:= #S;     // computes the size of the initial set S

    for i in [n..2 by -1] do
        prod:= &*([ S[j] : j in [1..(i-1)] ]);  // computes the product of the first (i-1) primes of the set of primes S, where i:= n, n-1, ..., 3
        bound:= Floor(m/prod);  // generates the upper bound on the possible ith prime of S; ie. for all primes P with P <= bound, we have [ prod_{s \in S} s ]*P <= m 
        a:= NextPrime(S[i]);    // generates the next largest prime of S[i]; this is the smallest potential ith prime 
        
        if a le bound then      // while the interval [a,b] is not empty 
            b:= PreviousPrime(bound + 1);       // generates the largest prime smaller than bound; this is the largest potential ith prime
            
            Int:= PrimesInInterval(a, b);       // generates the interval of all possible primes to add to S; all of these are valid candidates for the ith prime
            for p in Int do
                tS:= ([ S[j] : j in [1..(i-1)] ]);       // generates a set of primes tS; these are the first (i-1) primes of S
                Append(~tS, p);         // appends p to tS
                RecursiveSets(~tS, ~nSets, i+1, n, m);  // generates all sets S_j of length at most n whose first i primes are [2, p_1, ..., p_{i-2},p] of tS
            end for;
        end if;
    end for;
    
    NnSets:= [nSets[1]];        // generates a copy of nSets, to contain only maximal elements
    for s in nSets do
        if (&or[s subset t: t in NnSets | t ne s] eq false) and (s notin NnSets) then   // for each set s, verifies that s is not already a subset of another set in NnSets
            Append(~NnSets, s);
        end if;     
    end for;
    
    nSets:= NnSets;
    
    for s in nSets do
        if #s gt 2 then
            D0:= ExponentsXYZ2(s,0);    // generates all possible values for D:= ( (p_1)^(b_1) )* /cdots *( (p_n)^(b_n) ), where s:= [p_1, ..., p_n], b_i in {0,1} 
            for D in D0 do
                if D eq 1 then
                    printf "<%o, %o, %o>, %o, %o\n", s, D, 0, [], [];
                else
                    b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(s,D);   // determines prime splitting
                    if IsEmpty(SplitPrimes) then
                        printf "<%o, %o, %o>, %o, %o\n", s, D, 0, [], [];      // outputs the set of primes s, the corresponding D value (ie. the symmetric case), Ind, I, I_, where Ind = 0, I = [], I_ = []
                    else 
                        J:= IIPrimeXYZ2(SplitPrimes,D);         // generates all possible I, I_
                        for II_ in J do
                            Ind:= Index(J,II_);         // computes the index of II_ in J
                            p:= [j[1] : j in II_[1]];   // determines the primes of I
                            p_:= [j[1]: j in II_[2]];   // determines the primes of I_
                            printf "<%o, %o, %o>, %o, %o\n", s, D, Ind, p, p_;  // outputs the set of primes s, the corresponding D value, Ind, and I,I_ in the non-symmetric case 
                        end for;
                    end if;
                end if;
            end for;
        end if;
    end for;
    return "";
end function;