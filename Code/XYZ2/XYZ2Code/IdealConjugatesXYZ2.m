// IdealConjugatesXYZ2.m
   
/*
INPUT:
    Iset:= I or I_ (as output by IIPrimeXYZ2.m), where
        I:= { i | 1 <= i <= t, a_i > b_i }
            I:= [[p_1, h_1, pi_1, (pi_)_1, P_1],...,[p_n, h_n, pi_n, (pi_)_n, P_n]],
                p_i:= prime of S which splits in K
                h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
                pi_i:= generator of the ideal P^(h_i) = (pi)R, chosen so that Abs(pi) > Abs(pi_)
                (pi_)_i:= conjugate of pi, chosen so that Abs(pi) > Abs(pi_)
                P_i:= (unique) choice of prime ideal in K lying above p
        I_:= { i | 1 <= i <= t, a_i < b_i }
            I_:= [[p_1, h_1, pi_1, (pi_)_1, P_1],...,[p_n, h_n, pi_n, (pi_)_n, P_n]],
                p_i:= prime of S which splits in K
                h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
                pi_i:= generator of the ideal P^(h_i) = (pi)R, chosen so that Abs(pi) < Abs(pi_)
                (pi_)_i:= conjugate of pi, chosen so that Abs(pi) < Abs(pi_)
                P_i:= (unique) choice of prime ideal in K lying above p
    D:= squarefree part of x, where
        x:= Du^2 and x + y = z^2

OUTPUT:
    if Iset:= I then:
        Iset_:= I' = [[p_1, h_1, pi_1, (pi_)_1, P_1'],...,[p_n, h_n, pi_n, (pi_)_n, P_n']], where
            p_i:= prime of S which splits in K
            h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
            pi:= generator of the ideal P^(h_i) = (pi)R, chosen so that Abs(pi) > Abs(pi_)
            pi_:= conjugate of pi, chosen so that Abs(pi) > Abs(pi_)
            P_i':= conjugate of the (unique) choice of prime ideal in K lying above p
        corresponds to the set I of primes that split in K such that a_i > b_i
        
    if Iset:= I_ then:
        Iset_:= I_' = [[p_1, h_1, pi_1, (pi_)_1, P_1'],...,[p_n, h_n, pi_n, (pi_)_n, P_n']], where
            p_i:= prime of S which splits in K
            h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
            pi:= generator of the ideal P^(h_i) = (pi)R, chosen so that Abs(pi) < Abs(pi_)
            pi_:= conjugate of pi, chosen so that Abs(pi) < Abs(pi_)
            P_i':= conjugate of the (unique) choice of prime ideal in K lying above p
        corresponds to the set I_ of primes that split in K such that a_i < b_i
         
COMMENTS:
    Given an ideal P in K = Q(Sqrt(D)) above p in S, computes the ideal conjugate ideal P_, where
        (p)R = P*P_ as an ideal in K
         
    This algorithm is used to compute alpha, where
        a(R):= (alpha) = (2)^(b_0)*[prod_{i in I}(P_i)^(d_i)]*[prod_{i in I_}((P_i)_i)^(d_i)]*prod_{i = t + 1 ... s}P_i^(a_i)  

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,11,23];
    > D:= 253;   
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);  
    > II_:= (IIPrimeXYZ2(SplitPrimes,D))[1]; 
    > Iset:= II_[1];
    > Iset;
    [
        <3, 1, -sqrtD - 16, sqrtD - 16, Principal Prime Ideal
        Generator:
            -2*$.2 - 15>
    ]
    > IdealConjugatesXYZ2(Iset,D);
    [
        <3, 1, -sqrtD - 16, sqrtD - 16, Principal Ideal
        Generator:
            2*$.2 - 17>
    ]

    > Iset:= II_[2];
    > IdealConjugatesXYZ2(Iset,D);
    []

*/


function IdealConjugatesXYZ2(Iset,D)
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    R:= RingOfIntegers(K);      // generates ring of integers of K
    Iset_:= Iset;
    for i in [1..#Iset] do
        Iset_[i][5]:= Conjugate(Iset[i][5]);    // computes the conjugate P_ of the ideal P, where (p)R:= P*P_
    end for;
    return Iset_;
end function;