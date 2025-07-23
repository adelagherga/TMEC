// C12BoundXYZ2.m

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

    D:= squarefree part of x, where
        x:= Du^2 and x + y = z^2

OUTPUT:
    C12star:= a large positive constant such that Bstar < C12star, where Bstar, C12star are derived from B, C12 such that B < C12,
        B:= Max([M, U, |n|]), where, given
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
    C5:= constant appearing in the formulation of C12star
    C6:= constant appearing in the formulation of C12star
    C7:= constant appearing in the formulation of C12star
    C8:= constant appearing in the formulation of C12star
    C9:= constant appearing in the formulation of C12star
    N:= Max(|n_0|,...,|n_t|, |n_{t+1} - v_{t+1}|,...,|n_{s} - v_{s}|), a rational integer, where
        n_0:= the exponent appearing in alpha^{hstar} of +/-eps
        n_i:= the exponent appearing in alpha^{hstar} of pi_i of I, respectively (pi_)_i of I_, respectively p_i of tplus1
        v_i:= (1/2)*hstar*ord_{p_i}(4*D) for p_i in tplus1
    hstar:= LCM(2,h_1,...,h_s), where
        h_i:= smallest exponent such that P^(h_i) is principal, for P an ideal in K, lying above p in S
    IUstar:= IU cat { i | t+1 <= i <= s, v_i != 0}, where
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

    For a given value of alpha, this algorithm computes the maximal possible c_i, u_i, n appearing in u using p-adic linear forms in logarithms

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
    > A0:= AlphasXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > A:= A0[2];
    > A;
    <<Principal Ideal of R
    Generator:
        1, [
        <11, 1, 1/2*(-3*sqrtD - 1), 1/2*(3*sqrtD - 1), Principal Prime Ideal of R
        Generator:
            -3*$.2 + 1, 0>
    ]>, [], 2*sqrtD, 10>
    > C12star, C5, C6, C7, C8, C9, N, hstar, IUstar:= C12BoundXYZ2(II_,b0,A,S,D);
    > C12star;
    3.07682284851718749252059008295E33
    > C5;
    0.720210045206278239508775749793
    > C6;
    0.720210045206278239508775749793
    > C7;
    0.312108686183103147624758921650
    > C8;
    1.44042009041255647901755149959
    > C9;
    13.0973900630231924184162460842
    > N;
    0
    > hstar;
    2
    > IUstar;
    [ 2, 3, 5, 7, 13 ]

*/


function C12BoundXYZ2(II_,b0,A,S,D)
    Sort(~S);   // orders primes by size p_1 < ... < p_s
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    pls:=InfinitePlaces(K);
    assert Evaluate(sqrtD,pls[1]) eq Sqrt(D);

    eps:= FundamentalUnitXYZ2(D);
    //    Numericaleps:= QsqrtDPrecision(eps[1], D, eps[2]);  // computes numerical approximation of the fundamental unit, eps, as a real number
    Numericaleps:=Evaluate(eps,pls[1]);
    eps_:= Conjugate(eps);
//    Numericaleps_:= QsqrtDPrecision(eps_[1], D, eps_[2]);       // computes the numerical approximation of the conjugate of the fundamental unit, eps_, as a real number
    Numericaleps_:=Evaluate(eps,pls[1]);

    //    Numericala:= QsqrtDPrecision(A[3][1], D, A[3][2]);  // computes the numerical approximation of alpha as a real number
    Numericala:=Evaluate(A[3],pls[1]);
    a_:= Conjugate(A[3]);
//    Numericala_:= QsqrtDPrecision(a_[1], D, a_[2]);     // computes the numerical approximation of the conjugate of alpha, a_, as a real number
    Numericala_:=Evaluate(a_,pls[1]);

    I:= II_[1];
    I_:= II_[2];

    NumericalI:=[Evaluate(q[3]/q[4],pls[1]) : q in I];
    NumericalI_:=[Evaluate(q[4]/q[4],pls[1]) : q in I_];

    IU:= IUFactorsXYZ2(II_,A[3],S,D);   // computes prime factors (belonging to S) of G_{alpha}
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

    c1:= [];    // stores all potential values of C1 for each p in IU
    c2:= [];    // stores all potential values of C2 for each p in IU

    t2:= [];    // stores all potential values of (hstar*(gamma_i - lambda_i) + N) for each p in IU
    t3:= [];    // stores all potential values of (hstar*((gamma_i - kappa_i)/h_i) + N) for each p in I cat I_
    t5:= [];    // stores all potential values of ((4/3)*t*( p^(f_P/2) -1 )) for each p in I cat I_ cat IU

    for p in IU do
        C1_i, C2_i, e_P, f_P, t:= LambdaBoundXYZ2(II_,p,D);     // computes C1_i, C2_i such that ord_p(Lambda_i*) < C1_i + C2_i*log(B*)

        lambda_i:= Ordp(K!(2*sqrtD/a_),p,D);    // computes ord_p(2*Sqrt(D)/a_)
        if p eq 2 then
            gamma_i:= 3/2;
        elif p eq 3 then
            gamma_i:= 1;
        else
            gamma_i:= 1/2;
        end if;

        Append(~c1, -( lambda_i + Valuation(hstar,p) ) + C1_i/e_P);
        Append(~c2, C2_i/e_P);

        Append(~t2, hstar*( gamma_i - lambda_i) + N);
        Append(~t5, (4/3)*t*( p^(f_P/2) -1 ));
    end for;

    C1:= NonEmptyMax(c1);       // computes C1 such that U < C1 + C2*log(B*)
    C2:= NonEmptyMax(c2);       // computes C2 such that U < C1 + C2*log(B*)

    c3:= [];    // stores all potential values of C3 for each p in I cat I_
    c4:= [];    // stores all potential values of C4 for each p in I cat I_

    for k in I do
        C3_i, C4_i, e_P, f_P, t:= KappaBoundXYZ2(II_,IUstar,k[1],D);    // computes C3_i, C4_i such that ord_p(Kappa_i*) < C3_i + C4_i*log(B*)

        kappa_i:= Ordp(K!(A[3]/a_),k[1],D);     // computes ord_p(a/a_)
        if k[1] eq 2 then
            gamma_i:= 3/2;
        elif k[1] eq 3 then
            gamma_i:= 1;
        else
            gamma_i:= 1/2;
        end if;

        Append(~c3, (kappa_i + Ordp(K!(hstar),k[1],D))/k[2] + C3_i/(k[2]*e_P) );
        Append(~c4, C4_i/(k[2]*e_P) );

        Append(~t3, hstar*((gamma_i - kappa_i)/k[2]) + N);
        Append(~t5, (4/3)*t*( k[1]^(f_P/2) -1 ));
    end for;

    for k in I_ do
        C3_i_, C4_i_, e_P, f_P, t:= Kappa_BoundXYZ2(II_,IUstar,k[1],D);         // computes C3_i, C4_i such that ord_p((Kappa_)_i*) < C3_i + C4_i*log(B*)

        kappa_i:= Ordp(K!(a_/A[3]),k[1],D);     // computes ord_p(a_/a)
        if k[1] eq 2 then
            gamma_i:= 3/2;
        elif k[1] eq 3 then
            gamma_i:= 1;
        else
            gamma_i:= 1/2;
        end if;

        Append(~c3, (kappa_i + Ordp(K!(hstar),k[1],D))/k[2] + C3_i_/(k[2]*e_P) );
        Append(~c4, C4_i_/(k[2]*e_P) );

        Append(~t3, hstar*((gamma_i - kappa_i)/k[2]) + N);
        Append(~t5, (4/3)*t*( k[1]^(f_P/2) -1 ));
    end for;

    C3:= NonEmptyMax(c3);       // computes C3 such that M < C3 + C4*log(B*)
    C4:= NonEmptyMax(c4);       // computes C4 such that M < C3 + C4*log(B*)

    C5:= ( Log( 2*Abs(Numericala_/Numericala) ) )/(2*Log(Numericaleps));
    C6:= ( Log( 2*Abs(Numericala/Numericala_) ) )/(2*Log(Numericaleps));

//    C7:= (NonEmptySum([ Log(Abs( QsqrtDPrecision(q[3][1], D, q[3][2])/QsqrtDPrecision(q[4][1], D, q[4][2]) )) : q in I]) + NonEmptySum([ Log(Abs((QsqrtDPrecision(q[4][1], D, q[4][2]))/(QsqrtDPrecision(q[3][1], D, q[3][2])))) : q in I_]))/(2*Log(Numericaleps));
    C7:=(NonEmptySum([Log(Abs(q)) : q in NumericalI]) +
	 NonEmptySum([Log(Abs(q)) : q in NumericalI_]))/(2*Log(Numericaleps));

    //GammaXYZ2:= NonEmptyProduct([ Min( 1, Abs(QsqrtDPrecision(q[4][1], D, q[4][2])) ) : q in I])*NonEmptyProduct( [ Min( 1, Abs(QsqrtDPrecision(q[3][1], D, q[3][2])) ) : q in I_] );
    GammaXYZ2:= NonEmptyProduct([Min(1,Abs(Evaluate(q[4],pls[1]))) : q in I])*
		NonEmptyProduct([Min(1,Abs(Evaluate(q[3],pls[1]))) : q in I_]);

    PXYZ2:= NonEmptyProduct([p : p in IU]);

    C8:= NonEmptyMax([ Log( 4*Sqrt(D)/Abs(Numericala) )/Log(Numericaleps) , Log( (4*Sqrt(D)/Abs(Numericala_))*GammaXYZ2^(-C6/C7) )/Log( Numericaleps*GammaXYZ2^(1/C7) ), -C5,-C6 ]);

    C9:= Log(PXYZ2)/Log( Numericaleps*RealField(40)!(GammaXYZ2^(1/C7)) );
    C10:= NonEmptyMax([C1, C3, Abs(C5), Abs(C6)+C3*C7, C8+C1*C9]);
    C11:= NonEmptyMax([C2, C4, C4*C7, C2*C9]);

    T1:= 2*(N + hstar*C10 + hstar*C11*Log(hstar*C11));
    T2:= NonEmptyMax(t2);
    T3:= NonEmptyMax(t3);
    T4:= 2;
    T5:= NonEmptyMax(t5);

    C12star:= NonEmptyMax([T1,T2,T3,T4,T5]);    // computes C12* such that B* < C12*

    return C12star, C5, C6, C7, C8, C9, N, hstar, IUstar;
end function;
