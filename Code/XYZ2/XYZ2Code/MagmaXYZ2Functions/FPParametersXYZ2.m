// FPParametersXYZ2.m

/*
INPUT:
    FA:= < D, MaxA, MaxC12star, < MaxC5, MaxC6, MaxC7, MaxC8, MaxC9 >, MaxN, Maxhstar, MaxIUstar >, where
        D:= squarefree part of x, where
            x:= Du^2 and x + y = z^2
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
    U0:= [ [u_i, p_i] : p_i in IU], where
        u_i:= ord_{p_i}(u) such that 
            u:= G_{alpha}, or, said differently +/- u = G_{alpha} = prod_{i in IU}(p_i)^(u_i)
            ie. +/- u = G_{alpha} = prod_{i in IU}(p_i)^(u_i), under x = D*u^2 
        p_i:= prime in IU
    M0:= [ [c_i, p_i] : p_i in I], where
        c_i:= the exponents m_i on Prod[I], Prod_[I] in G_{alpha} giving G_{alpha} = u, under x = D*u^2, where
            G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_]) and
                Prod[I]:= prod_{i in I} (pi_i)^{m_i}
                Prod_[I_]:= prod_{i in I_} ((pi_)_i)^{m_i}
                Prod_[I]:= prod_{i in I} ((pi_)_i)^{m_i}
                Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}, where
                    m_i:= an integer > 0
    M0_:= [ [c_i, p_i] : p_i in I_], where
        c_i:= the exponents m_i on Prod[I_], Prod_[I_] in G_{alpha} giving G_{alpha} = u, under x = D*u^2, where
            G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_]) and
                Prod[I]:= prod_{i in I} (pi_i)^{m_i}
                Prod_[I_]:= prod_{i in I_} ((pi_)_i)^{m_i}
                Prod_[I]:= prod_{i in I} ((pi_)_i)^{m_i}
                Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}, where
                    m_i:= an integer > 0
    absn:= |n|, where n is an integer representing the exponent on eps, eps_ in G_{alpha} giving G_{alpha} = u, under x = D*u^2, where
            G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_]) and 
                eps:= the fundamental unit of the rind of integers of K = Q(sqrt(D))
                eps_:= the conjugate of the fundamental unit, eps
    b0:= an integer in {0,1} indicating whether 2 splits in K
    A:= [ q, q_, alpha, tplus1 ] (as output by AlphasXYZ2.m), where 
        q:= < prod_i, [< p, h_i, pi, pi_, P , d_i > : j in [1..n] ] >, such that
            prod_i:= prod_{i in I}(P_i)^(d_i), the product of prime ideals of I which contribute to alpha
            p:= prime of S which splits in K
            h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
            pi:= generator of the ideal P^(h_i) = (pi)R
            (pi_):= conjugate of pi
            P:= (unique) choice of prime ideal in K lying above p
            n:= the length of Iset
            d_i:= the exponent on the prime ideal P_i in prod_i
        q_:= < prod_i, [< p, h_i, pi, pi_, P_ , d_i > : j in [1..n] ] >, such that
            prod_i:= prod_{i in I}((P_)_i)^(d_i), the product of the conjugate of the prime ideals of I_ which contribute to alpha
            p:= prime of S which splits in K
            h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
            pi:= generator of the ideal P^(h_i) = (pi)R
            (pi_):= conjugate of pi
            P_:= conjugate of the (unique) choice of prime ideal in K lying above p
            n:= the length of Iset
            d_i:= the exponent on the prime ideal P_i in prod_i
        alpha:= element of K, generating the ideal (alpha), where
            (alpha):= (2)^(b_0)*[prod_{i in I}(P_i)^(d_i)]*[prod_{i in I_}((P_i)_i)^(d_i)]*prod_{i = t + 1 ... s}P_i^(a_i)
        tplus1:= the product of primes of S from NonSplitPrimes whose prime ideals above contribute to alpha as 
            prod_{i = t + 1 ... s}P_i^(a_i)
    S:= [p_1,...,p_s], p_i primes in S

OUTPUT:
    U0:= [ [u_i, p_i] : p_i in IU], where 
        u_i:= ord_{p_i}(u) updated such that 
            u:= G_{alpha}, or, said differently +/- u = G_{alpha} = prod_{i in IU}(p_i)^(u_i)
            ie. +/- u = G_{alpha} = prod_{i in IU}(p_i)^(u_i), under x = D*u^2 
        p_i:= prime in IU
    M0:= [ [c_i, p_i] : p_i in I], where
        c_i:= the updated exponents m_i on Prod[I], Prod_[I] in G_{alpha} giving G_{alpha} = u, under x = D*u^2, where
            G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_]) and
                Prod[I]:= prod_{i in I} (pi_i)^{m_i}
                Prod_[I_]:= prod_{i in I_} ((pi_)_i)^{m_i}
                Prod_[I]:= prod_{i in I} ((pi_)_i)^{m_i}
                Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}, where
                    m_i:= an integer > 0
    M0_:= [ [c_i, p_i] : p_i in I_], where
        c_i:= the updated exponents m_i on Prod[I_], Prod_[I_] in G_{alpha} giving G_{alpha} = u, under x = D*u^2, where
            G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_]) and
                Prod[I]:= prod_{i in I} (pi_i)^{m_i}
                Prod_[I_]:= prod_{i in I_} ((pi_)_i)^{m_i}
                Prod_[I]:= prod_{i in I} ((pi_)_i)^{m_i}
                Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}, where
                    m_i:= an integer > 0
    absn:= updated |n|, where n is an integer representing the exponent on eps, eps_ in G_{alpha} giving G_{alpha} = u, under x = D*u^2, where
            G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_]) and 
                eps:= the fundamental unit of the rind of integers of K = Q(sqrt(D))
                eps_:= the conjugate of the fundamental unit, eps
    eqns:= [[x_1,y_1,z_1],...,[x_n,y_n,z_n]], all S-unit equations x_i + y_i = (z_i)^2 such that 
        u_i:= ord_p(u) is in the range [m + m0 - lambda_i - Ord_p(hstar), U0[p]], where 
            p:= a prime in U0
        c_j is in the range [(m + m0 - kappa_i - Ord_p(hstar))/h_i, M0[p]], where
            p:= a prime in M0
        (c_j)_ is in the range [(m + m0 - kappa__i - Ord_p(hstar))/h_i, M0_[p]], where
            p:= a prime in M0_
            
COMMENTS:
    This algorithm selects an appropriate memory bound for the Fincke-Pohst algorithm (based on the number of split primes in K), Sbound, and automates the FPRestrictionsXYZ2.m script

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,5,7,11];
    > D:= 15;
    > K<sqrtD>:= QuadraticField(D);
    > R:= RingOfIntegers(K);
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
    > C:= IIPrimeXYZ2(SplitPrimes,D); 
    > II_:= C[4];
    > A0:= AlphasXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > FA:= MaximalC12BoundXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > U0,M0,M0_,absn:= SmallestLatticeBoundXYZ2(II_,FA,S);
    > U0;
    [
        [ 44, 2 ],
        [ 24, 3 ],
        [ 17, 5 ]
    ]
    > M0;
    [
        [ 32, 7 ]
    ]
    > M0_;
    [
        [ 25, 11 ]
    ]
    > absn;
    442
    > A:= A0[2];
    > U0,M0,M0_,absn,eqns:= FPParametersXYZ2(FA,U0,M0,M0_,absn,b0,A,S);
    > U0;
    [
        [ 11, 2 ],
        [ 5, 3 ],
        [ 3, 5 ]
    ]
    > M0;
    [
        [ 7, 7 ]
    ]
    > M0_;
    [
        [ 11, 11 ]
    ]
    > absn;
    119
    > eqns;
    [
        [ 199290375, -686, 14117 ]
    ]

*/

                                                        
function FPParametersXYZ2(FA,U0,M0,M0_,absn,b0,A,S)
    eqns:= [];  // stores S-unit solutions x + y = z^2
    
    if (#M0 + #M0_) eq 2 then   // if there are 2 split primes in K
        Sbound:= 70000;         // sets bound on Fincke-Pohst memory allotment
    elif (#M0 + #M0_) ge 3 then         // if there are 3 or more split primes in K
        Sbound:= 5000000;       // sets bound on Fincke-Pohst memory allotment
    end if;
    StillEnumerating:= true;
    MinExp:= NonEmptyMin([u[1] : u in U0 cat M0 cat M0_]);      // computes the smallest possible fp in U0, M0, M0_; 
                                                                    // to avoid negative m values in the lattice constructions, this algorithm terminates when this value reaches 0 
    
    Mean:= Ceiling( (&+[u[1] : u in U0 cat M0 cat M0_])/(#[u[1] : u in U0 cat M0 cat M0_]) );   // computes the mean value of U0, M0, M0_ values
    rep:= [];   // stores < U0,M0,M0_,absn > to ensure that Fincke-Pohst is no longer applied when these values repeat
    while (StillEnumerating eq true) and (MinExp ne 0) and (Mean gt 10) do
        Append(~rep, < U0,M0,M0_,absn >);       // ensures the while loop terminates if there is no change to U0, M0, M0_, absn after the Fincke-Pohst algorithm is applied 
        U0,M0,M0_,absn,soln,StillEnumerating:= FPRestrictionsXYZ2(FA,U0,M0,M0_,absn,b0,A,0.25,Sbound,S);
        Mean:= Ceiling( (&+[u[1] : u in U0 cat M0 cat M0_])/(#[u[1] : u in U0 cat M0 cat M0_]) );   // updates the mean value of U0, M0, M0_ values
        
        if < U0,M0,M0_,absn > in rep then
            for s in soln do
                if (s in eqns) eq false then
                    Append(~eqns, s);
                end if;
            end for;
            break;
        else
            MinExp:= NonEmptyMin([u[1] : u in U0 cat M0 cat M0_]);
            for s in soln do
                if (s in eqns) eq false then
                    Append(~eqns, s);
                end if;
            end for;
        end if;
    end while;
    
    return U0,M0,M0_,absn,eqns;
end function;