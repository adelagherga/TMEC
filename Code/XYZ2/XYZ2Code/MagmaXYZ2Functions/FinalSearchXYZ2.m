// FinalSearchXYZ2.m

/*
INPUT:
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
    eqns:= [[x_1,y_1,z_1],...,[x_n,y_n,z_n]], all S-unit equations x_i + y_i = (z_i)^2 such that
        u_i:= ord_p(u) is in the range [0, ..., U0[p]], where
            p:= a prime in U0
        c_j is in the range [0, ..., M0[p]], where
            p:= a prime in M0
        (c_j)_ is in the range [0, ..., M0_[p]], where
            p:= a prime in M0_

COMMENTS:
    This algorithm computes all solutions (x,y,z) of x + y = z^2 below the bounds set out by U0, M0, M0_ via brute force. This is only used after an initial bound reduction via the Fincke-Pohst algorithm.

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,5,7,11];
    > D:= 15;
    > K<sqrtD>:= QuadraticField(D);
    > R:= RingOfIntegers(K);
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
    > C:= IIPrimeXYZ2(SplitPrimes,D);
    > II_:= C[2];
    > A0:= AlphasXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > FA:= MaximalC12BoundXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > U0,M0,M0_,absn:= SmallestLatticeBoundXYZ2(II_,FA,S);
    > U0;
    [
        [ 20, 2 ],
        [ 14, 3 ],
        [ 10, 5 ],
        [ 9, 11 ]
    ]
    > M0;
    [
        [ 26, 7 ]
    ]
    > M0_;
    []
    > absn;
    57
    > A:= A0[1];
    > A;
    <<Principal Ideal of R
    Generator:
        1, [
        <7, 2, sqrtD + 8, -sqrtD + 8, Prime Ideal of R
        Two element generators:
            7
            $.2 + 1, 0>
    ]>, [], 1, 1>
    > FinalSearchXYZ2(U0,M0,M0_,absn,A,S,D);
    [
        [ 1070530560, 2401, 32719 ],
        [ 960, 1, 31 ],
        [ 16335, 49, 128 ],
        [ 240, 49, 17 ],
        [ 3375, 2401, 76 ],
        [ 117649, 29040, 383 ],
        [ 3840, 2401, 79 ],
        [ 2160, 49, 47 ],
        [ 15, 1, 4 ],
        [ 49, 15, 8 ]
    ]

*/


function FinalSearchXYZ2(U0,M0,M0_,absn,A,S,D)
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    Z:= IntegerRing();

    eps:= FundamentalUnitXYZ2(D);
    eps_:= Conjugate(eps);
    q:= A[1];
    q_:= A[2];
    a:= A[3];
    a_:= Conjugate(a);

    if #q eq 0 then
        Prod1:= [1];    // I == [], hence Prod[I]:= prod_{i in I} (pi_i)^{m_i}: = 1
        Prod3:= [1];    // I == [], hence Prod_[I]:= prod_{i in I} ((pi_)_i)^{m_i}:= 1
    else
        I:= q[2];
        Prod1:= IdealExponentsXYZ2(I,<"FinalSearch", M0>);      // computes every possible Prod[I]:= prod_{i in I} (pi_i)^{m_i} value; m_i bounded by M0[i]
        Prod3:= IdealExponentsXYZ2(I,<"FinalSearch_", M0>);     // computes every possible Prod_[I]:= prod_{i in I} ((pi_)_i)^{m_i} values; m_i bounded by M0[i]
    end if;

    if #q_ eq 0 then
        I_:= [];
        Prod2:= [1];    // I_ == [], hence Prod_[I_]:= prod_{i in I_} ((pi_)_i)^{m_i}: = 1
        Prod4:= [1];    // I_ == [], hence Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}: = 1
    else
        I_:= q_[2];
        Prod2:= IdealExponentsXYZ2(I_,<"FinalSearch_", M0_>);        // computes every possible Prod_[I_]:= prod_{i in I_} ((pi_)_i)^{m_i}; m_i bounded by M0_[i]
        Prod4:= IdealExponentsXYZ2(I_,<"FinalSearch", M0_>);         // computes every possible Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}; m_i bounded by M0_[i]
    end if;

    P1:= [];
    for p1 in Prod1 do
        for p2 in Prod2 do
            Append(~P1, p1*p2);         // computes possible (Prod[I])*(Prod_[I_]) values
        end for;
    end for;

    P2:= [];
    for p3 in Prod3 do
        for p4 in Prod4 do
            Append(~P2, p3*p4);         // computes possible (Prod_[I])*(Prod[I_]) values
        end for;
    end for;

    G_alphas:= [];      // stores [ G_alpha, n, m ] where G_alpha is a non-zero integer divisible only by primes of S; n is the exponent on eps in G_alpha; m is the exponent on the splitting primes in G_alpha
    Ones:= [];  // stores [ n, m ] where G_alpha == 1; n is the exponent on eps in G_alpha; m is the exponent on the splitting primes in G_alpha

    n:= -absn;  // minimap possible exponent on the fundamental unit, eps
    m:= 1;      // sets the exponent on the splitting primes; Nb. the case m = 0 is precisely the symmetric case, hence may be omitted
    while n le absn do
        while m le #P1 do
            G_alpha:= (a/(2*sqrtD))*((eps)^(n))*(P1[m]) - (a_/(2*sqrtD))*((eps_)^(n))*(P2[m]);

            if (G_alpha in Z) and (G_alpha ne 0) and &or([Z! (G_alpha) mod p eq 0 : p in S]) then
                rem, factors:= SFactors(Z !(G_alpha), S);
                    if (rem eq 1) or (rem eq -1) then
                        Append(~G_alphas,< G_alpha, n, m >);
                    end if;
                m:= m + 1;
            elif (G_alpha eq 1) or (G_alpha eq -1) then
                Append(~Ones, < n,m >);
                m:= m + 1;
            else
                m:= m + 1;
            end if;
        end while;
        n:= n + 1;
        m:= 1;  // resets m value
    end while;

    soln:={};         // stores [x,y,z] where x + y = z^2
    for g in G_alphas do
        u:= Z! g[1];    // computes u, where u:= +/-G_alpha
        n:= g[2];
        m:= g[3];
        z:= (a/2)*((eps)^(n))*(P1[m]) + (a_/2)*((eps_)^(n))*(P2[m]);    // computes z, where z:= +/-H_alpha
        if z in Z then
            x:= Z! ((sqrtD*u)^2);
            y:= Z!(z^2) - x;
            z:= Abs(Z! z);
            if (x ne 0) and (y ne 0) and (z ne 0) then
                if (x ge y) then
                    soln:=soln join {[x,y,z]};
                elif (x lt y) then
                    soln:=soln join {[y,x,z]};
                end if;
            end if;
        end if;
    end for;

    for g in Ones do
        u:= 1;  // computes u, where u:= +/-G_alpha = +/-1
        n:= g[1];
        m:= g[2];
        z:= (a/2)*((eps)^(n))*(P1[m]) + (a_/2)*((eps_)^(n))*(P2[m]);    // computes z, where z:= +/-H_alpha
        if z in Z then
            x:= Z! ((sqrtD*u)^2);
            y:= Z!(z^2) - x;
            z:= Abs(Z! z);
            if (x ne 0) and (y ne 0) and (z ne 0) then
                if (x ge y) then
                    soln:=soln join {[x,y,z]};
                elif (x lt y) then
                    soln:=soln join {[y,x,z]};
                end if;
            end if;
        end if;
    end for;

    return soln;
end function;
