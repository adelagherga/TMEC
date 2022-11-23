// FPRestrictionsXYZ2.m

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
    C:= a constant 0 < C < 1 determining m, where
            m:= the precision on the p-adic field Q_p), a lattice invariant
    Sbound:= a constant > 0, determining the maximum number of short vectors to be computed in the Fincke-Pohst algorithm
            This is necessary to avoid maxing out the memory; Fincke-Pohst solutions are only used when this bound is not reached
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
    StillEnumerating:= boolean variable indicator whether the preset memory-maximum, Sbound, on Fincke-Pohst has been reached; ie.
        if "false":= indicates that the number of short vector solutions in Fincke-Pohst is larger than SBound in all three lattices;
                short vectors search in all three lattices is terminated early; not enough memory, resort to FinalSearchXYZ2.m algorithm
        if "true": indicates that the number of short vector solutions in Fincke-Pohst is smaller than SBound in at least one of the three lattices;
                short vectors search in at least one lattice has completed, with memory below preset memory-maximum; may re-run this algorithm

COMMENTS:
    For each p in U0 (respectively M0, M0_), this algorithm computes the lattice Lambda_m* (respectively Kappa_m*, (Kappa_)_m*) for an appropriate choice of m. The relevant short
    vectors are found within an ellipsoid in this lattice; the volume of this ellipsoid is
        (( Pi(RealField())^(n/2) )/( Gamma(n/2+1) ))*(bound^(n/2)/Det(B)), where B is an n x n matrix representing the lattice Lambda* (respectively Kappa_m*, (Kappa_)_m*).
    Heuristically, the smaller the volume of the ellipsoid, the less short vectors to be found.
    This algorithm selects the smallest-volume lattice among all Lambda_m*, Kappa_m*, and (Kappa_)_m*, and uses the Fincke-Pohst algorithm to find all relevant short vectors. These vectors
    are then sorted to give solutions (x,y,z) of x + y = z^2.

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
    > C:= 0.5;
    > Sbound:= 1000000;
    > U0,M0,M0_,absn,eqns,StillEnumerating:= FPRestrictionsXYZ2(FA,U0,M0,M0_,absn,b0,A,C,Sbound,S);
    > U0;
    [
        [ 44, 2 ],
        [ 24, 3 ],
        [ 17, 5 ]
    ]
    > M0;
    [
        [ 23, 7 ]
    ]
    > M0_;
    [
        [ 25, 11 ]
    ]
    > absn;
    442
    > eqns;
    []
    > StillEnumerating;
    true

*/


function FPRestrictionsXYZ2(FA,U0,M0,M0_,absn,b0,A,C,Sbound,S)
    D:= FA[1];
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))

    C5:= FA[4][1];
    C6:= FA[4][2];
    C7:= FA[4][3];
    C8:= FA[4][4];
    C9:= FA[4][5];

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

    IU:= IUFactorsXYZ2([I,I_],a,S,D);   // computes prime factors (belonging to S) of G_{alpha}
    IUstar:= IU;

    n, nM, nM_, nEps, N:= nExponentsXYZ2(b0,A,S,D);
    for j in n do
        if (j[1] in IUstar eq false) and (j[3] ne 0) then       // computes IU cat { i | t+1 <= i <= s, v_i != 0}; v_i := (1/2)*hstar*ord_{p_i}(4D) for i = t+1, ..., s, (p_i != 2)
                                                                // Nb. v_i > 0 when p_i|4D, where p_i != 2, hence v_i > 0 when p_i|D, precisely when p_i ramifies in K, thus n_i != 0
                                                                    // that is, do not have to include primes p_i where n_i = v_i
            Append(~IUstar,j[1]);
        end if;
    end for;
    Sort(~IUstar);

    H:= [h[2] : h in I cat I_];         // computes set of h_i for i = 1, ..., s
    Append(~H,2);
    hstar:= Lcm(H);     // Nb. prime ideals above non-split primes are either ramified in K, so h_i = 2, or remain prime in K, so h_i = 1
                            // hence, hstar = lcm(2, h_i of split primes)

    Us:= [u[2] : u in U0];      // generates the set of primes p of U0 (hence of IU)
    Ms:= [m[2] : m in M0];      // generates the set of primes p of M0 (hence of I)
    M_s:= [m[2] : m in M0_];    // generates the set of primes p of M0_ (hence of I_)
    ns:= [en[1] : en in n];     // generates the set of primes p of n

    U:= NonEmptyMax([u[1] : u in U0]);
    M:= NonEmptyMax([m[1] : m in M0 cat M0_]);

    LambdaVols:= [];    // stores the volume of each lattice Lambda* corresponding to each p in IU
                            // the relevant short vectors are found within an ellipsoid in this lattice; the volume of this ellipsoid is
                                // (( Pi(RealField())^(n/2) )/( Gamma(n/2+1) ))*(bound^(n/2)/Det(B)), where B is an n x n matrix representing the lattice Lambda*
                            // hence (heuristically) the smaller the volume of the ellipsoid, the less short vectors to be found
    KappaVols:= [];     // stores the volume of each lattice Kappa* corresponding to each p in I
                            // the relevant short vectors are found within an ellipsoid in this lattice; the volume of this ellipsoid is
                                // (( Pi(RealField())^(n/2) )/( Gamma(n/2+1) ))*(bound^(n/2)/Det(B)), where B is an n x n matrix representing the lattice Kappa*
                            // hence (heuristically) the smaller the volume of the ellipsoid, the less short vectors to be found

    Kappa_Vols:= [];    // stores the volume of each lattice Kappa_* corresponding to each p in I_
                            // the relevant short vectors are found within an ellipsoid in this lattice; the volume of this ellipsoid is
                                // (( Pi(RealField())^(n/2) )/( Gamma(n/2+1) ))*(bound^(n/2)/Det(B)), where B is an n x n matrix representing the lattice Kappa_*
                            // hence (heuristically) the smaller the volume of the ellipsoid, the less short vectors to be found

    LBound:= (hstar*absn + nEps)^2 + NonEmptySum([ (hstar*M0[ Index(Ms,j[1]) ][1] + nM[Index(Ms,j[1]) ][2])^2 : j in I ])       // computes the bound corresponding to the lattice Lambda*
            + NonEmptySum([ (hstar*M0_[ Index(M_s,j[1]) ][1] + nM_[Index(M_s,j[1]) ][2])^2 : j in I_ ]);                                // ie. (n*)^2 + sum_{i \in I}(c_j*)^2 + sum_{i \in I_}(c_j*)^2
    KBound:= (hstar*absn + nEps)^2 + NonEmptySum([ (hstar*M0[ Index(Ms,j[1]) ][1] + nM[Index(Ms,j[1]) ][2])^2 : j in I ])
            + NonEmptySum([ (hstar*M0_[ Index(M_s,j[1]) ][1] + nM_[Index(M_s,j[1]) ][2])^2 : j in I_ ]);        // computes the bound corresponding to the lattices Kappa*, Kappa_*
                                                                                                                    // ie. (n*)^2 + sum_{i \in I}(c_j*)^2 + sum_{i \in I_}(c_j*)^2 + sum_{i \in IU}(u_j*)^2

    for i in Us do
        if i in ns then
            KBound:= KBound + ( hstar*U0[ Index(Us,i) ][1] - (n[ Index(ns, i) ][2] - n[ Index(ns, i) ][3]) )^2;         // if p_i in IU also appears in the factorization of alpha (ie. in ns), then
                                                                                                                            // n_j, v_j may be non-zero, so u_j*:= (hstar)*)(u_j) - (n_j - v_j)
        else
            KBound:= KBound + ( hstar*U0[ Index(Us,i) ][1] )^2;         // if p_i in IU does not appear in the factorization of alpha (ie. in ns), then
                                                                            // necessarily, n_j:= v_j:= 0, so u_j*:= (hstar)*)(u_j)
        end if;
    end for;

    for p in Us do
        Indexp:= Index(Us, p);  // computes the index of p in U0
        fp:= U0[Indexp][1];     // computes the fp value associated to p, where ord_p(u) <= fp

        // if a solution exists with ord_{p_i}(u):= u_i <= fp, for a given p_i in IU, then pick m such that m < fp - m0 + lambda_i + ord_{p_i}(hstar)
            // then m + m0 - lambda_i - ord_{p_i}(hstar) < fp
            // if a solution exists with ord_{p_i} = u_i >= m + m0 - lambda_i - ord_{p_i}(hstar); ie. if a solution exists with u_i in [m + m0 - lambda_i - ord_{p_i}(hstar), fp], then
            // since Lambda_m*:= LambdaLattice(A,p_i,m,D) is the set of all solutions (x_1, ..., x_n)^T such that u_i >= m + m0 - lambda_i - ord_{p_i}(hstar),
            // said solution would necessarily have to be in Lambda_m*
        // this computation is carried out below

        if fp ge 2 then
            terms, p0, m0:= LambdaInitialSortXYZ2(A,p,100,D);
            lambda_i:= Ordp(K!(2*sqrtD/a_),p,D);        // computes ord_p(2*Sqrt(D)/a_)

            if Floor( fp - m0 + lambda_i + Valuation(hstar,p) ) gt 0 then       // fail-safe ensuring m is chosen so that final fp value remains positive, despite C value
                if fp gt 8 then
                    m:= Floor( fp - m0 + lambda_i + Valuation(hstar,p) ) - Ceiling(C*fp);       // for "large" (i.e. > 8) values of fp, chooses m < Floor(fp - m0 + lambda_i + ord_p(hstar))
                                                                                                    // according to input value, C
                else
                    m:= Floor( fp - m0 + lambda_i + Valuation(hstar,p) ) - 1;   // for "small" (i.e. <= 8) values of fp, decreases m only by 1
                end if;
            elif Floor( fp - m0 + lambda_i + Valuation(hstar,p) ) lt 0 then
                m:= 1;
            else
                m:= Floor( fp - m0 + lambda_i + Valuation(hstar,p) );   // if Floor( fp - m0 + lambda_i + ord_p(hstar) ) == 0, chooses m so that resulting fp = 0
            end if;

            B,m0:= LambdaLatticeXYZ2(A,p,m,D);
            cols:= NumberOfColumns(B);
            VolumeOfLattice:= (( Pi(RealField())^(cols/2) )/( Gamma(cols/2+1) ))*(LBound^(cols/2)/Determinant(B));      // computes the volume of the ellipsoid
            Append(~LambdaVols,< A, p, m, VolumeOfLattice >);
        end if;
    end for;

    for p in Ms do
        Indexp:= Index(Ms,p);   // computes the index of p in M0
        P:= I[Indexp][5];
        h_i:= I[Indexp][2];
        fp:= M0[Indexp][1];     // computes the fp value associated to p, where c_i <= fp

        // if a solution exists with c_i <= fp, for a given p_i in I, then pick m such that m < (h_i)*fp - m0 + kappa_i + ord_{p_i}(hstar)
            // then (m + m0 - kappa_i - ord_{p_i}(hstar))/h_i < fp
            // if a solution exists with c_i >= (m + m0 - kappa_i - ord_{p_i}(hstar))/h_i; ie. if a solution exists with c_i in [(m + m0 - kappa_i - ord_{p_i}(hstar))/h_i, fp], then
            // since Kappa_m*:= KappaLatticeXYZ2(A,IUstar,p,P,m,D) is the set of all solutions (x_1, ..., x_n)^T such that c_i >= (m + m0 - kappa_i - ord_{p_i}(hstar))/h_i,
            // said solution would necessarily have to be in Kappa_m*
        // this computation is carried out below

        terms, p0, m0:= KappaInitialSortXYZ2(A,IUstar,p,P,100,D);
        kappa_i:= Ordp(K!(a/a_),p,D);     // computes ord_p(a/a_)

        if Floor( (h_i)*fp - m0 + kappa_i + Valuation(hstar,p) ) gt 0 then      // fail-safe ensuring m is chosen so that final fp value remains positive, despite C value
            if fp gt 8 then
                m:= Floor( (h_i)*fp - m0 + kappa_i + Valuation(hstar,p) ) - Ceiling(C*fp);      // for "large" (i.e. > 8) values of fp, chooses m < Floor((h_i)*fp - m0 + kappa_i + ord_p(hstar))
                                                                                                    // according to input value, C
            else
                m:= Floor( (h_i)*fp - m0 + kappa_i + Valuation(hstar,p) ) - 1;  // for "small" (i.e. <= 8) values of fp, decreases m only by 1
            end if;
        elif Floor( (h_i)*fp - m0 + kappa_i + Valuation(hstar,p) ) lt 0 then
            m:= 1;
        else
            m:= Floor( (h_i)*fp - m0 + kappa_i + Valuation(hstar,p) );  // if Floor( (h_i)*fp - m0 + kappa_i + ord_p(hstar) ) == 0, chooses m so that resulting fp = 0
        end if;

        B,m0:= KappaLatticeXYZ2(A,IUstar,p,P,m,D);
        cols:= NumberOfColumns(B);
        VolumeOfLattice:= (( Pi(RealField())^(cols/2) )/( Gamma(cols/2+1) ))*(KBound^(cols/2)/Determinant(B));  // computes the volume of the ellipsoid
        Append(~KappaVols,< A, IUstar, p, P, h_i, m , VolumeOfLattice >);
    end for;

    for p in M_s do
        Indexp:= Index(M_s,p);  // computes the index of p in M0_
        P:= I_[Indexp][5];
        h_i:= I_[Indexp][2];
        fp:= M0_[Indexp][1];    // computes the fp value associated to p, where c_i <= fp

        // if a solution exists with c_i <= fp, for a given p_i in I_, then pick m such that m < (h_i)*fp - m0 + kappa__i + ord_{p_i}(hstar)
            // then (m + m0 - kappa__i - ord_{p_i}(hstar))/h_i < fp
            // if a solution exists with c_i >= (m + m0 - kappa__i - ord_{p_i}(hstar))/h_i; ie. if a solution exists with c_i in [(m + m0 - kappa__i - ord_{p_i}(hstar))/h_i, fp], then
            // since Kappa__m*:= Kappa_LatticeXYZ2(A,IUstar,p,P,m,D) is the set of all solutions (x_1, ..., x_n)^T such that c_i >= (m + m0 - kappa__i - ord_{p_i}(hstar))/h_i,
            // said solution would necessarily have to be in Kappa__m*
        // this computation is carried out below

        terms, p0, m0:= Kappa_InitialSortXYZ2(A,IUstar,p,P,100,D);
        kappa__i:= Ordp(K!(a_/a),p,D);      // computes ord_p(a_/a)

         if Floor( (h_i)*fp - m0 + kappa__i + Valuation(hstar,p) ) gt 0 then    // fail-safe ensuring m is chosen so that final fp value remains positive, despite C value
            if fp gt 8 then
                m:= Floor( (h_i)*fp - m0 + kappa__i + Valuation(hstar,p) ) - Ceiling(C*fp);     // for "large" (i.e. > 8) values of fp, chooses m < Floor((h_i)*fp - m0 + kappa__i + ord_p(hstar))
                                                                                                    // according to input value, C
            else
                m:= Floor( (h_i)*fp - m0 + kappa__i + Valuation(hstar,p) ) - 1; // for "small" (i.e. <= 8) values of fp, decreases m only by 1
            end if;
        elif Floor( (h_i)*fp - m0 + kappa__i + Valuation(hstar,p) ) lt 0 then
            m:= 1;
        else
            m:= Floor( (h_i)*fp - m0 + kappa__i + Valuation(hstar,p) );         // if Floor( (h_i)*fp - m0 + kappa__i + ord_p(hstar) ) == 0, chooses m so that resulting fp = 0
        end if;

        B,m0:= Kappa_LatticeXYZ2(A, IUstar,p,P,m,D);
        cols:= NumberOfColumns(B);
        VolumeOfLattice:= (( Pi(RealField())^(cols/2) )/( Gamma(cols/2+1) ))*(KBound^(cols/2)/Determinant(B));  // computes the volume of the ellipsoid
        Append(~Kappa_Vols,< A, IUstar, p, P, h_i, m , VolumeOfLattice >);
    end for;

    Vols:= [];
    if LambdaVols ne [] then
        MinLambdaVol, IndexMinLambdaVol:= Min([d[4] : d in LambdaVols]);        // selects the prime p of U0 such that the corresponding matrix B gives the smallest-volume ellipsoid corresponding to Lambda*
        Append(~Vols, [MinLambdaVol, IndexMinLambdaVol]);
    else
        MinLambdaVol:= 0;
        IndexMinLambdaVol:= 0;
    end if;
    if KappaVols ne [] then
        MinKappaVol, IndexMinKappaVol:= Min([d[7] : d in KappaVols]);   // selects the prime p of M0 such that the corresponding matrix B gives the smallest-volume ellipsoid corresponding to Kappa*
        Append(~Vols, [MinKappaVol, IndexMinKappaVol]);
    else
        MinKappaVol:= 0;
        IndexMinKappaVol:= 0;
    end if;
    if Kappa_Vols ne [] then
        MinKappa_Vol, IndexMinKappa_Vol:= Min([d[7] : d in Kappa_Vols]);        // selects the prime p of M0_ such that the corresponding matrix B gives the smallest-volume ellipsoid corresponding to Kappa_*
        Append(~Vols, [MinKappa_Vol, IndexMinKappa_Vol]);
    else
        MinKappa_Vol:= 0;
        IndexMinKappa_Vol:= 0;
    end if;

    Sort(~Vols);        // sorts volumes of lattices in order of increasing volume
    soln:= {};  // stores S-unit solutions x + y = z^2

    i:= 1;
    StillEnumerating:= false;
    while i le #Vols do     // runs through Vols sequentially, beginning with the lattice of smallest volume
        if Vols[i][1] eq MinLambdaVol then
            A:= LambdaVols[IndexMinLambdaVol][1];
            p:= LambdaVols[IndexMinLambdaVol][2];
            m:= LambdaVols[IndexMinLambdaVol][3];
            lambda_i:= Ordp(K!(2*sqrtD/a_),p,D);    // computes ord_p(2*Sqrt(D)/a_)
            terms, p0, m0:= LambdaInitialSortXYZ2(A,p,100,D);
            B,m0:= LambdaLatticeXYZ2(A,p,m,D);
            W,w,maxReached:= FinckePohst(B,LBound,2*Sbound);    // computes the vectors of B, w, whose norm is <= LBound, enumerated using the Fincke-Pohst algorithm

            if maxReached eq false then         // if all short vectors are found before reaching the memory-maximum, 2*Sbound
                if #w ne 0 then
                    LambdaSol:= FPRestrictionsLambdaXYZ2(A,U0,M0,M0_,absn,nM,nM_,nEps,hstar,terms,p0,w,W,S,D);
                    soln:=soln join LambdaSol;
                end if;

                U0[IndexMinLambdaVol][1]:= Ceiling ( m + m0 - lambda_i - Valuation(hstar,p) ) - 1;      // updates U0
                U:= Max([u[1] : u in U0]);      // updates U
                absn:= Floor(Max([C5, C6 + C7*M, C8 + C9*U]));
                i:= 5;  // sets i to 5, exits loop
                StillEnumerating:= true;        // sets StillEnumerating == true, indicating that solutions were found
            else
                i:= i + 1;
            end if;

        elif Vols[i][1] eq MinKappaVol then
            A:= KappaVols[IndexMinKappaVol][1];
            IUstar:= KappaVols[IndexMinKappaVol][2];
            p:= KappaVols[IndexMinKappaVol][3];
            P:= KappaVols[IndexMinKappaVol][4];
            h_i:= KappaVols[IndexMinKappaVol][5];
            m:= KappaVols[IndexMinKappaVol][6];
            kappa_i:= Ordp(K!(a/a_),p,D);     // computes ord_p(a/a_)
            terms, p0, m0:= KappaInitialSortXYZ2(A, IUstar, p, P, 100, D);
            B,m0:= KappaLatticeXYZ2(A, IUstar, p, P, m, D);
            W,w,maxReached:= FinckePohst(B,KBound,2*Sbound);    // computes the vectors of B whose norm is <= bound, enumerated using the Fincke-Pohst algorithm

            if maxReached eq false then         // if all short vectors are found before reaching the memory-maximum, 2*Sbound
                if #w ne 0 then
                    KappaSol:= FPRestrictionsKappaXYZ2(A,U0,M0,M0_,absn,n,nM,nM_,nEps,hstar,terms,p0,w,W,S,D);
                    if IsEmpty(KappaSol) eq false then
                        soln:=soln join KappaSol;
                    end if;
                end if;

                M0[IndexMinKappaVol][1]:= Ceiling( (m + m0 - kappa_i - Valuation(hstar,p))/h_i ) - 1;   // updates M0
                M:= NonEmptyMax([m[1] : m in M0 cat M0_]);      // updates M
                absn:= Floor(Max([C5, C6 + C7*M, C8 + C9*U]));
                i:= 5;  // sets i to 5, exits loop
                StillEnumerating:= true;        // sets StillEnumerating == true, indicating that solutions were found
            else
                i:= i + 1;
            end if;

        elif Vols[i][1] eq MinKappa_Vol then
            A:= Kappa_Vols[IndexMinKappa_Vol][1];
            IUstar:= Kappa_Vols[IndexMinKappa_Vol][2];
            p:= Kappa_Vols[IndexMinKappa_Vol][3];
            P:= Kappa_Vols[IndexMinKappa_Vol][4];
            h_i:= Kappa_Vols[IndexMinKappa_Vol][5];
            m:= Kappa_Vols[IndexMinKappa_Vol][6];
            kappa__i:= Ordp(K!(a_/a),p,D);      // computes ord_p(a_/a)
            terms, p0, m0:= Kappa_InitialSortXYZ2(A, IUstar, p, P, 100, D);
            B,m0:= Kappa_LatticeXYZ2(A, IUstar, p, P, m, D);
            W,w,maxReached:= FinckePohst(B,KBound,2*Sbound);    // computes the vectors of B whose norm is <= bound, enumerated using the Fincke-Pohst algorithm

            if maxReached eq false then         // if all short vectors are found before reaching the memory-maximum, 2*Sbound
                if #w ne 0 then
                    Kappa_Sol:= FPRestrictionsKappa_XYZ2(A,U0,M0,M0_,absn,n,nM,nM_,nEps,hstar,terms,p0,w,W,S,D);
                    if IsEmpty(Kappa_Sol) eq false then
                        soln:=soln join Kappa_Sol;
                    end if;
                end if;

                M0_[IndexMinKappa_Vol][1]:= Ceiling( (m + m0 - kappa__i - Valuation(hstar,p))/h_i ) - 1;        // updates M0_
                M:= NonEmptyMax([m[1] : m in M0 cat M0_]);      // updates M
                absn:= Floor(Max([C5, C6 + C7*M, C8 + C9*U]));
                i:= 5;  // sets i to 5, exits loop
                StillEnumerating:= true;        // sets StillEnumerating == true, indicating that solutions were found
            else
                i:= i + 1;
            end if;
        end if;
    end while;

    return U0,M0,M0_,absn,soln,StillEnumerating;
end function;
