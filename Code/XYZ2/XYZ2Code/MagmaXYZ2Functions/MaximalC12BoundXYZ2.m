// MaximalC12BoundXYZ2.m

/*
INPUT:
    II_:= [I,I_] (as output by IIPrimeXYZ2.m), where
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
    
    b0:= an integer in {0,1} indicating whether 2 splits in K
        
    SplitPrimes:= primes of S which split in K = Q(Sqrt(D)) 
        SplitPrimes is stored in Magma as SplitPrimes[i]:= [p_i, h_i, P_i], where
            p_i:= prime of S which splits in K
            h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
            P_i:= (unique) choice of prime ideal in K lying above p
        
    NonSplitPrimes:= primes of S which do not split in K = Q(Sqrt(D)) but potentially contribute to alpha
        NonSplitPrimes is stored in Magma as SplitPrimes[i]:= [p_i, P_i], where
            p_i:= prime of S which does not split in K but contributes to the product ideal (X)
            P_i:= prime ideal in K lying above p 

    S:= [p_1,...,p_s], p_i primes in S

    D:= squarefree part of x, where
        x:= Du^2 and x + y = z^2

OUTPUT:
    FA:= < D, MaxA, MaxC12star, < MaxC5, MaxC6, MaxC7, MaxC8, MaxC9 >, MaxN, Maxhstar, MaxIUstar >, where
        MaxA:= [ q, q_, alpha, tplus1 ] as output by AlphasXYZ2.m giving the largest C12star value
        MaxC12star:= the maximal positive constant C12star (where the maximum runs over all values of alpha) such that Bstar < C12star, where 
            Bstar, C12star are derived from B, C12 such that B < C12, and B:= Max([M, U, |n|]), where, given 
                G_{alpha}:= prod_{i in IU}(p_i)^(u_i), 
                G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_])
                    Prod[I]:= prod_{i in I} (pi_i)^{m_i}
                    Prod_[I_]:= prod_{i in I_} ((pi_)_i)^{m_i}
                    Prod_[I]:= prod_{i in I} ((pi_)_i)^{m_i}
                    Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}
                M:= Max_{i in I cat I_}(c_i), where c_i > 0 are the exponents m_i in G_{alpha} giving G_{alpha} = u, under x = D*u^2
                    ie. M is an upper bound on the exponents of the split primes appearing in u
                U:= Max_{i in IU}(u_i)
                    ie. U is an upper bound on the exponents of the primes of S appearing in the factorization of u
                n:= the exponent of the fundamental unit, eps, and its conjugate, eps_, appearing in G_{alpha} giving G_{alpha} = u, under x = D*u^2
        MaxC5:= the maximal constant C5 appearing in the formulation of C12star, where the maximum runs over all values of alpha
        MaxC6:= the maximal constant C6 appearing in the formulation of C12star, where the maximum runs over all values of alpha
        MaxC7:= the maximal constant C7 appearing in the formulation of C12star, where the maximum runs over all values of alpha
        MaxC8:= the maximal constant C8 appearing in the formulation of C12star, where the maximum runs over all values of alpha
        MaxC9:= the maximal constant C9 appearing in the formulation of C12star, where the maximum runs over all values of alpha
        MaxN:= the maximal rational integer N = Max(|n_0|,...,|n_t|, |n_{t+1} - v_{t+1}|,...,|n_{s} - v_{s}|), where the maximum runs over all values of alpha, and 
            n_0:= the exponent appearing in alpha^{hstar} of +/-eps
            n_i:= the exponent appearing in alpha^{hstar} of pi_i of I, respectively (pi_)_i of I_, respectively p_i of tplus1
            v_i:= (1/2)*hstar*ord_{p_i}(4*D) for p_i in tplus1
        hstar:= the maximal hstar = LCM(2,h_1,...,h_s), where the maximum runs over all values of alpha and 
            h_i:= smallest exponent such that P^(h_i) is principal, for P an ideal in K, lying above p in S 
        IUstar:= IU cat { i | t+1 <= i <= s, v_i != 0}, giving the largest C12star value, where
            IU:= the set of primes p_i of S such that 
                G_{alpha} = prod_{i in IU}(p_i)^(u_i)
            v_i := (1/2)*hstar*ord_{p_i}(4D) for p_i in tplus1 
    
COMMENTS:
    Given a solution (x,y,z) to x + y = z^2, if x = D*u^2, define u:= G_{alpha}, where
        G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_])
            Prod[I]:= prod_{i in I} (pi_i)^{m_i}
            Prod_[I_]:= prod_{i in I_} ((pi_)_i)^{m_i}
            Prod_[I]:= prod_{i in I} ((pi_)_i)^{m_i}
            Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}, where
                m_i:= an integer > 0
    Solving for x, y, z is equivalent to solving for n, m_i = c_i such that 
        +/- u = G_{alpha} = prod_{i in IU}(p_i)^(u_i)
        
    This algorithm computes the maximal possible c_i, u_i, n appearing in u over all possible values of alpha associated to D 

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,5,7,11,13];
    > D:= 5;
    > K<sqrtD>:= QuadraticField(D);
    > R:= RingOfIntegers(K);
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
    > C:= IIPrimeXYZ2(SplitPrimes,D);
    > II_:= C[1];
    > FA:= MaximalC12BoundXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > FA;
    <5, <<Principal Ideal of R
    Generator:
        1, [
        <11, 1, 1/2*(-3*sqrtD - 1), 1/2*(3*sqrtD - 1), Principal Prime Ideal of R
        Generator:
            -3*$.2 + 1, 0>
    ]>, [], 2, 2>, 3.07682284851718749252059008295E33, 
    <0.720210045206278239508775749793, 0.720210045206278239508775749793, 
    0.312108686183103147624758921650, 3.11269602859711122518787062598, 
    16.4419419393923019107568843370>, 1, 2, [ 2, 3, 5, 7, 13 ]>

*/
    
    
function MaximalC12BoundXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D)
    
    Alphas:= AlphasXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);         // computes every possible alpha associated to this value of D
    C12s:= [];  // for each alpha, stores C12star and associated values 
    for A in Alphas do
        C12star, C5, C6, C7, C8, C9, N, hstar, IUstar:= C12BoundXYZ2(II_,b0,A,S,D);
        Append(~C12s, < D, A, C12star, C5, C6, C7, C8, C9, N, hstar, IUstar >);
    end for;

    MaxC5:= Max([C[4] : C in C12s]);
    MaxC6:= Max([C[5] : C in C12s]);
    MaxC7:= Max([C[6] : C in C12s]);
    MaxC8:= Max([C[7] : C in C12s]);
    MaxC9:= Max([C[8] : C in C12s]);
    MaxN:= Max([C[9] : C in C12s]);
    Maxhstar:= Max([C[10] : C in C12s]);
    
    MaxC12, MaxC12Index:= Max([C[3] : C in C12s]);
    
    FA:= < C12s[MaxC12Index][1], C12s[MaxC12Index][2], MaxC12, < MaxC5, MaxC6, MaxC7, MaxC8, MaxC9 >, MaxN, Maxhstar, C12s[MaxC12Index][11] >;
    return FA;  

end function;

