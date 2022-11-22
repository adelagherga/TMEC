// SmallestLatticeBoundXYZ2.m

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
    
    S:= [p_1,...,p_s], p_i primes in S

OUTPUT:
    U0:= [ [u_i, p_i] : p_i in IU], where
        u_i:= ord_{p_i}(u) such that 
            u:= G_{alpha}, or, said differently +/- u = G_{alpha} = prod_{i in IU}(p_i)^(u_i)
            ie. +/- u = G_{alpha} = prod_{i in IU}(p_i)^(u_i), under x = D*u^2 
            This algorithm computes the upper bound on the exponent u_i
        p_i:= prime in IU
    M0:= [ [c_i, p_i] : p_i in I], where
        c_i:= the exponents m_i on Prod[I], Prod_[I] in G_{alpha} giving G_{alpha} = u, under x = D*u^2, where
            G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_]) and
                Prod[I]:= prod_{i in I} (pi_i)^{m_i}
                Prod_[I_]:= prod_{i in I_} ((pi_)_i)^{m_i}
                Prod_[I]:= prod_{i in I} ((pi_)_i)^{m_i}
                Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}, where
                    m_i:= an integer > 0
            This algorithm computes the upper bound on the exponent c_i
    M0_:= [ [c_i, p_i] : p_i in I_], where
        c_i:= the exponents m_i on Prod[I_], Prod_[I_] in G_{alpha} giving G_{alpha} = u, under x = D*u^2, where
            G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_]) and
                Prod[I]:= prod_{i in I} (pi_i)^{m_i}
                Prod_[I_]:= prod_{i in I_} ((pi_)_i)^{m_i}
                Prod_[I]:= prod_{i in I} ((pi_)_i)^{m_i}
                Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}, where
                    m_i:= an integer > 0
            This algorithm computes the upper bound on the exponent c_i
    absn:= |n|, where n is an integer representing the exponent on eps, eps_ in G_{alpha} giving G_{alpha} = u, under x = D*u^2, where
                G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_]) and 
                    eps:= the fundamental unit of the rind of integers of K = Q(sqrt(D))
                    eps_:= the conjugate of the fundamental unit, eps
           This algorithm computes the upper bound on the absolute value of the exponent, n 
    
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
        
    This algorithm computes the maximal possible c_i, u_i, |n| appearing in u over all possible values of alpha associated to D 

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,5,7];
    > D:= 2;
    > K<sqrtD>:= QuadraticField(D);
    > R:= RingOfIntegers(K);
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
    > C:= IIPrimeXYZ2(SplitPrimes,D); 
    > II_:= C[1];
    > A0:= AlphasXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > FA:= MaximalC12BoundXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > U0,M0,M0_,absn:= SmallestLatticeBoundXYZ2(II_,FA,S);
    > U0;
    [
        [ 19, 2 ],
        [ 15, 3 ],
        [ 12, 5 ]
    ]
    > M0;
    [
        [ 21, 7 ]
    ]
    > M0_;
    []
    > absn;
    75

*/


function SmallestLatticeBoundXYZ2(II_,FA,S)
    D:= FA[1];
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    
    A:= FA[2];  
    X0:= FA[3];

    C5:= FA[4][1];
    C6:= FA[4][2];
    C7:= FA[4][3];
    C8:= FA[4][4];
    C9:= FA[4][5];
    
    N:= FA[5];
    hstar:= FA[6];
    IUstar:= FA[7];
    
    q:= A[1];
    q_:= A[2];
    a:= A[3];
    a_:= Conjugate(a);
    
    if #q eq 0 then
        I:= [];
    else
        I:= q[2];
    end if;
    if #q_ eq 0 then
        I_:= [];
    else
        I_:= q_[2];
    end if;
    
    X:= [];     // stores values X0 such that l > sqrt(#terms)*X0
    IU:= IUFactorsXYZ2(II_,a,S,D);
    
    fixX0:= true;
    while fixX0 eq true do
        U0:= [];        // stores u = mu + m0 - 1 - lambda_i values
        M0:= [];        // stores c = mu + m0 - 1 - kappa_i values
        M0_:= [];       // stores c_ = mu + m0 - 1 - (kappa_)_i values
       
        Lterms:= 1 + #[i : i in I] + #[i_:i_ in I_];    // computes the number of terms in Lambda_i*, including eps/eps_
        Kterms:= 1 + #(IUstar) + #[i : i in I] + #[i_:i_ in I_];         // computes the number of terms in Kappa_i*, including eps or eps_, respectively
        
        for p in IU do  // runs through each prime in IU
            lambda_i:= Ordp(K!(2*sqrtD/a_),p,D);    // computes ord_p(2*Sqrt(D)/a_)
            
            mu:= Ceiling((Lterms*Log(X0))/Log(p)) + 5;  // estimates initial mu value for the given prime, p 
            while p^mu lt X0^(Lterms) do     // ensures that mu is chosen so that p^{mu} is of the size of X0^{# or terms in Lambda_i*}
                mu:= mu + 3;
            end while;

            changeCL:= true;
            while changeCL eq true do
                B,m0:= LambdaLatticeXYZ2(A,p,mu,D);     // computes the lattice associated to Lambda_i*, outputs lattice matrix, m0 value
                s:= NumberOfColumns(B);         // computes the number of columns of B

                if s le 2 then  // if B is a 2 x 2 matrix only 
                    l:= F3Approx(B);    // computes short lattice length explicitly via the F3Approx algorithm
                else
                    C:= LLL(Transpose(B));      // for s > 2, computes the lower bound for lattice lenght l(Gamma_{mu}) via LLL algorithm
                    x:= C[1];
                    const:= 2^(-(s-1)/2);
                    l:= const*Sqrt(Norm(x));    // computes the lower bound for l(Gamma_{mu}) using |c_1|, the norm of the first column of the lattice (cf. Lemma 3.4)
                end if;
                
                if l le Sqrt(s)*X0 then
                    mu:= mu + 3;        // changeCL still True, returns to start of changeC with updated mu value, repeats until l = l(Gamma_{mu}) > sqrt(s)*X0
                else
                    u:= mu + Ceiling(m0) - 1 - Ceiling(lambda_i);       // application of Lemma 3.14: no solutions with mu <= ord_p(Lambda_i*), hence ord_p(Lambda_i) <= mu + m0 - 1
                                                                        // by Lemma 7.6(i), u_i = ord_p(Lambda_i) - ord_p(hstar) - lambda_i <= mu + m0 - 1 - ord_p(hstar) - lambda_i
                    
                    Append(~U0,[u,p]);  // stores u_i = mu + m0 - 1 - lambda_i values; corresponding p_i in IU
                    changeCL:= false;   // exits changeCL algorithm, continues to next value of p in IU
                end if;
            end while;
        end for;
        
        for j in I do   // runs through each prime in I
            p:= j[1];
            P:= j[5];
            kappa_i:= Ordp(K!(a/a_),p,D);     // computes ord_p(a/a_)
            
            mu:= Ceiling((Kterms*Log(X0))/Log(p)) + 5;  // estimates initial mu value for the given prime, p 
            while p^mu lt X0^(Kterms) do        // ensures that mu is chosen so that p^{mu} is of the size of X0^{# or terms in K_i*}
                mu:= mu + 3;
            end while;
            
            changeCK:= true;
            while changeCK eq true do
                B,m0:= KappaLatticeXYZ2(A, IUstar, p, P, mu, D);        // computes the lattice associated to K_i*, outputs lattice matrix, m0 value
                s:= NumberOfColumns(B);         // computes the number of columns of B

                if s le 2 then  // if B is a 2 x 2 matrix only 
                    l:= F3Approx(B);    // computes short lattice length explicitly via the F3Approx algorithm
                else
                    C:= LLL(Transpose(B));      // for s > 2, computes the lower bound for lattice lenght l(Gamma_{mu}) via LLL algorithm
                    x:= C[1];
                    const:= 2^(-(s-1)/2);
                    l:= const*Sqrt(Norm(x));    // computes the lower bound for l(Gamma_{mu}) using |c_1|, the norm of the first column of the lattice (cf. Lemma 3.4)
                end if;
                
                if l le Sqrt(s)*X0 then
                    mu:= mu + 3;        // changeCK still True, returns to start of changeC with updated mu value, repeats until l = l(Gamma_{mu}) > sqrt(s)*X0
                else
                    c:= mu + Ceiling(m0) - 1 - Ceiling(kappa_i);        // application of Lemma 3.14: no solutions with mu <= ord_p(K_i*), hence ord_p(K_i) <= mu + m0 - 1
                                                                        // by Lemma 7.6(ii), c_i <= (h_i)*(c_i) = ord_p(K_i) - ord_p(hstar) - kappa_i <= mu + m0 - 1 - ord_p(hstar) - kappa_i
                    
                    Append(~M0,[c,p]);  // stores c_i = mu + m0 - 1 - kappa_i values; corresponding p_i in I
                    changeCK:= false;   // exits changeCK algorithm, continues to next value of p in I
                end if;
            end while;
        end for;
            
        for j in I_ do  // runs through each prime in I_
            p:= j[1];
            P:= j[5];
            kappa__i:= Ordp(K!(a_/a),p,D);      // computes ord_p(a_/a)
            
            mu:= Ceiling((Kterms*Log(X0))/Log(p)) + 5;  // estimates initial mu value for the given prime, p 
            while p^mu lt X0^(Kterms) do        // ensures that mu is chosen so that p^{mu} is of the size of X0^{# or terms in (K_)_i*}
                mu:= mu + 3;
            end while;
            
            changeCK_:= true;
            while changeCK_ eq true do
                B,m0:= Kappa_LatticeXYZ2(A, IUstar, p, P, mu, D);       // computes the lattice associated to (K_)_i*, outputs lattice matrix, m0 value
                s:= NumberOfColumns(B);         // computes the number of columns of B

                if s le 2 then  // if B is a 2 x 2 matrix only 
                    l:= F3Approx(B);    // computes short lattice length explicitly via the F3Approx algorithm
                else
                    C:= LLL(Transpose(B));      // for s > 2, computes the lower bound for lattice lenght l(Gamma_{mu}) via LLL algorithm
                    x:= C[1];
                    const:= 2^(-(s-1)/2);
                    l:= const*Sqrt(Norm(x));    // computes the lower bound for l(Gamma_{mu}) using |c_1|, the norm of the first column of the lattice (cf. Lemma 3.4)
                end if;
                
                if l le Sqrt(s)*X0 then
                    mu:= mu + 3;        // changeCK_ still True, returns to start of changeC with updated mu value, repeats until l = l(Gamma_{mu}) > sqrt(s)*X0
                else
                    c:= mu + Ceiling(m0) - 1 - Ceiling(kappa__i);       // application of Lemma 3.14: no solutions with mu <= ord_p((K_)_i*), hence ord_p((K_)_i) <= mu + m0 - 1
                                                                        // by Lemma 7.6(ii), c_i <= (h_i)*(c_i) = ord_p((K_)_i) - ord_p(hstar) - (kappa_)_i <= mu + m0 - 1 - ord_p(hstar) - (kappa_)_i
                    
                    Append(~M0_,[c,p]); // stores c_i = mu + m0 - 1 - (kappa_)_i values; corresponding p_i in I_
                    changeCK_:= false;  // exits changeCK_ algorithm, continues to next value of p in I_
                end if;
            end while;
        end for;
        
        U:= NonEmptyMax([u0[1] : u0 in U0]);    // computes largest u value
        M:= NonEmptyMax([c0[1] : c0 in (M0 cat M0_)]);  // computes largest c value  
                  
        absn:= Floor(Max([C5, C6 + C7*M, C8 + C9*U]));  // updates |n|
        B:= Max([M,U,absn]);    // updates B
        X0:= hstar*B + N;       // computes B* so that max{all exponents appearing in x + y = z^2} <= B* 
        Append(~X,X0);  // stores all X0 values
        if &+[1: x in X | x eq X0] eq 2 then    // computes number of times X0 appears in X; terminates algorithm if X0 twice - ie. lowest upper bound cannot be further improved
            fixX0 := false;
        end if;
    end while;
    
    return U0,M0,M0_,absn;
end function;