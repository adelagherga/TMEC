//IIPrimeXYZ2.m    

/*
INPUT:
    SplitPrimes:= [[p_1, h_1, P_1],...,[p_n, h_n, P_n]], where
        p_i:= prime of S which splits in K:= Q(Sqrt(D))Sp
        h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
        P_i:= (unique) choice of prime ideal in K lying above p
    D:= squarefree part of x, where
        x:= Du^2 and x + y = z^2

OUTPUT:
    C:= [[Qi_1,(Qi_)_1],...,[Qi_n,(Qi_)_n]];
        Qi_j:= [[p_1, h_1, pi_1, (pi_)_1, P_1],...,[p_n, h_n, pi_n, (pi_)_n, P_n]],
        corresponding to the set I_j of primes that split in K such that a_i > b_i, where
            p_i:= prime of S which splits in K
            h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
            pi_i:= generator of the ideal P^(h_i) = (pi)R, chosen so that Abs(pi) > Abs(pi_)
            (pi_)_i:= conjugate of pi, chosen so that Abs(pi) > Abs(pi_)
            P_i:= (unique) choice of prime ideal in K lying above p
        (Qi_)_j:= [[p_1, h_1, pi_1, (pi_)_1, P_1],...,[p_n, h_n, pi_n, (pi_)_n, P_n]],
        corresponding to the set (I_)_j of primes that split in K such that a_i < b_i, where
            p_i:= prime of S which splits in K
            h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
            pi_i:= generator of the ideal P^(h_i) = (pi)R, chosen so that Abs(pi) < Abs(pi_)
            (pi_)_i:= conjugate of pi, chosen so that Abs(pi) < Abs(pi_)
            P_i:= (unique) choice of prime ideal in K lying above p
         
COMMENTS:
    Given a solution (x,y,z) to x + y = z^2, if x = D*u^2, define X:= z + u*SqrtD(D) and
        (X):= Prod_{i=1}^s ((P_i)^{a_i})*((P_i')^{b_i}) and 
        (X'):= Prod_{i=1}^s ((P_i')^{a_i})*((P_i)^{b_i}), ideals of K:= Q(Sqrt(D)) where
            s:= number of primes of S:= [p_1,...,p_s]
            P_i:= prime ideal of K above the prime p_i of S
            P_i':= the conjugate ideal of P_i
    For each i, exactly one of (a_i > b_i) or (a_i < b_i) or (a_i = b_i = 0) must occur

    Assume the primes p_i of S which split in K are p_1,...,p_t, for some t in [0,s]
    This algorithm computes all possible combinations (without mirroring) of the sets of splitting primes
        I:= { i | 1 <= i <= t, a_i > b_i} and 
        I_:= { i | 1 <= i <= t, a_i < b_i}
        
    This information is used to compute u and z, and hence (x,y,z)  

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,11,23];
    > D:= 253;   
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);  
    > SplitPrimes;
    [
        <3, 1, Principal Prime Ideal
        Generator:
            -2*$.2 - 15>
    ]
    > IIPrimeXYZ2(SplitPrimes,D);
    [
        [
            [
                <3, 1, -sqrtD - 16, sqrtD - 16, Principal Prime Ideal
                Generator:
                    -2*$.2 - 15>
            ],
            []
        ]
    ]

    S:= [2,3,5,7];
    > D:= 3;
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
    > SplitPrimes;  
    []
    > IIPrimeXYZ2(SplitPrimes, D);
    []

*/

            
function IIPrimeXYZ2(SplitPrimes,D)
    n:= #SplitPrimes;   // computes the number of primes of S which split in K = Q(Sqrt(D))
    C0:= Subsets({1..n});       // generates all subsets of the set {1..n}
    C:= {};     // stores all subsets of {1..n} as indexed sets
    for c in C0 do
        Include(~C,SetToIndexedSet(c));         // converts subsets of {1..n} to indexed sets (needed to access indices)
    end for;
    F:= { {@u, v@} : u,v in C | IsDisjoint(u, v) and #u ge #v };        // computes all possible combinations of subsets of {1..n} which are disjoint
                                                                        // sorted so that F[i]:= [F[i][1],F[i][2]] has #F[i][1] >= #F[i][2]
    Exclude(~F,{@ {@@} @});     // removes empty set
    
    C:= [];     // stores all possible combinations of I,I_    
    for d in F do
        I:= [SplitPrimes[i] : i in d[1]];       // stores the set of primes that split in K such that a_i > b_i
        I_:= [SplitPrimes[i] : i in d[2]];      // stores the set of primes that split in K such that a_i < b_i
        Qi, Qi_:= SplitPrimePropertiesXYZ2(I,I_,D);
        Append(~C,[Qi,Qi_]);
    end for;
    
    /*
    C0:= C;
    for i in [1..#C0] do
        for j in [1..#C0] do
            if (i ne j) and ( C0[i][2] eq C0[j][2] ) and ( C0[i][1] subset C0[j][1] ) then
                Exclude(~C, C0[i]);
            end if;
        end for;
    end for; */
    
    return C;
end function;
