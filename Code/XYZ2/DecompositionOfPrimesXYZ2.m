//DecompositionOfPrimesXYZ2.m

/*
INPUT:
    S:= [p_1,...,p_s], p_i primes in S
    D:= squarefree part of x, where
        x:= Du^2 and x + y = z^2

OUTPUT:
    b0:= an integer in {0,1} indicating whether 2 splits in K:= Q(Sqrt(D))
        b0:= 0 indicates that 2 does not split in K
        b0:= 1 indicates that 2 splits in K
        
    SplitPrimes:= [[p_1, h_1, P_1],...,[p_n, h_n, P_n]], where
        p_i:= prime of S which splits in K
        h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
        P_i:= (unique) choice of prime ideal in K lying above p
        
    NonSplitPrimes:= [[p_1, P_1],...,[p_n, P_n]], where
        p_i:= prime of S which does not split in K but contributes to the product ideal (X)
        P_i:= prime ideal in K lying above p 

COMMENTS:
    For all primes p in the set S, computes the decomposition of primes above p in K = Q(Sqrt(D))

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,5,7,11,13];
    > D:= 55;
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D); 
    > b0;
    0
    > SplitPrimes;
    [
        <3, 2, Prime Ideal
        Two element generators:
            3
            $.2 + 2>,
        <13, 2, Prime Ideal
        Two element generators:
            13
            $.2 + 9>
    ]
    > NonSplitPrimes;
    [
        <2, Prime Ideal
        Two element generators:
            2
            $.2 + 1>,
        <5, Prime Ideal
        Two element generators:
            5
            $.2>,
        <11, Prime Ideal
        Two element generators:
            11
            $.2>
    ]

*/


function DecompositionOfPrimesXYZ2(S,D)
    Sort(~S);   // orders primes by size p_1 < ... < p_s
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    R:= RingOfIntegers(K);      // generates ring of integers of K; every element of R is of the form a + Nu*b, where a,b are integers, Nu is defined as below
    
    if (D mod 4 eq 2) or (D mod 4 eq 3) then
        Nu:= sqrtD;
        NumericalNu:= Sqrt(D);  // the numerical approximation of Nu as a real number
    elif (D mod 4 eq 1) then
        Nu:= (1 + sqrtD)/2;
        NumericalNu:= (1 + Sqrt(D))/2;  // the numerical approximation of Nu as a real number
    end if;     // D mod 4 eq 0 is impossible since it would imply that D would be divisible by the square 4
    Divs:= Divisors(ClassNumber(K));     // computes the divisors of the class number of K
    
    SplitPrimes:= [];   // stores < p, h_i, P >, where p splits in K; P is the choice of prime ideal above p in K; h_i is an integer such that P^(h_i) = (pi) is principal
    NonSplitPrimes:= [];         // stores < p, P >, where p ramifies in K; P^2 = (p)R is the prime ideal above p in K
                                // stores < 2, P >, if 2 remains inert in K, where P = (2)R is the prime ideal above p in K
    
    b0:= 0;     // b0 = 1 only if 2 splits in K
    
    for p in S do
        // if p remains prime in K, then it contributes to the product ideal (X) only if p = 2, with exponent 0 or 1:
        if (p eq 2) and (IsInert(2,R) eq true) then      // if 2 remains prime in K 
            Append(~NonSplitPrimes,< p, ideal<R|2> >);  // stores p = 2 and the ideal (2)R in K above 2
        
        // if p ramifies in K, then it contributes to the product ideal (X) with exponent 1
        elif IsRamified(p,R) eq true then       // if p ramifies in K
            Append(~NonSplitPrimes,< p, Decomposition(R,p)[1][1] >);      // stores p and the ideal (p)R in K above p
        
        // if p splits in K, then it contributes to the product ideal (X) with exponent a_i
        elif IsSplit(p,R) eq true then
            Decomp:= Decomposition(R,p);     // factors the ideal (p)R above p into (two - since p splits in a quadratic field) prime ideals of R with multiplicity
            for j in [1..(#Divs)] do
                h_i:= Divs[j];  // h_i necessarily must divide h, the class number of K
                t,pi:= IsPrincipal(Decomp[1][1]^(h_i));     // determines the smallest integer h_i such that the 1st ideal of (p)R above p is principal; if true, a generator is also returned
                t_,pi_:= IsPrincipal(Decomp[2][1]^(h_i));   // determines the smallest integer h_i such that the 2nd ideal of (p)R above p is principal; if true, a generator is also returned  
                if (t eq true) or (t_ eq true) then
                    break j;
                end if;
            end for;
            
            NumericalPi:= pi[1] + NumericalNu*pi[2];    // computes the numerical approximation of pi as a real number
            NumericalPi_:= pi_[1] + NumericalNu*pi_[2];         //  computes the numerical approximation of pi_ as a real number
            
            if Abs(NumericalPi) gt Abs(NumericalPi_) then       // NumericalPi, NumericalPi_, the generators of the 2 primes above p in K 
                Append(~SplitPrimes, <p, h_i, Decomp[1][1]>);       // makes a choice of unique prime lying above p, based on which generator, NumericalPi or NumericalPi_, is larger
            else
                Append(~SplitPrimes, <p, h_i, Decomp[2][1]>);
            end if;
            
            if (p eq 2) then    // if 2 splits in K
                b0:= 1;
            end if;
        end if;
    end for;
    return b0, SplitPrimes, NonSplitPrimes;
    
end function;
