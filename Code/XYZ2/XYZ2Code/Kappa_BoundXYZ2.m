// Kappa_BoundXYZ2.m

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
    
    IUstar:= IU Cat { i | t+1 <= i <= s, v_i != 0 } Cat { i | t+1 <= i <= s, v_i - n_i != 0 }, where
        IU:= [p_1,...,p_n], p_i primes in S such that
            G_{alpha} mod p_i == 0, where
                G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_])
                    Prod[I]:= prod_{i in I} (pi_i)^{m_i}
                    Prod_[I_]:= prod_{i in I_} ((pi_)_i)^{m_i}
                    Prod_[I]:= prod_{i in I} ((pi_)_i)^{m_i}
                    Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}, where
                        m_i:= an integer > 0
        n_i:= the exponent appearing in alpha^{hstar} of pi_i of I, respectively (pi_)_i of I_, respectively p_i of tplus1
            alpha^{hstar} = +/- [(eps)^(n_0)]*[prod_{i in I} (pi_i)^(n_i)]*[prod_{i in I_} ((pi_)_i)^(n_i)]*[prod_{i = t+1,...,s}(p_i)^(n_i)]*[2^(hstar*b0)]
        v_i:= (1/2)*hstar*ord_{p_i}(4*D) for p_i in tplus1
    
    p:= prime of I_
    
    D:= squarefree part of x, where
        x:= Du^2 and x + y = z^2
        
OUTPUT:
    C3_i_:= the real constant such that ord_P((Kappa_)_i*) < C3_i_ + C4_i_*Log(B)
    C4_i_:= the real constant such that ord_P((Kappa_)_i*) < C3_i_ + C4_i_*Log(B)
    e_P:= ramification index of the ideal P of K:= Q(Sqrt(D)) above p
    f_P:= residue class degree of the ideal P of K:= Q(Sqrt(D)) above p
    n:= the number of terms in (Kappa_)_i*
        ie. the number of elements alpha_1,...,alpha_n:= eps, p_1,...,p_{n_1}, pi_1,...,pi_{n_2}, (pi_)_1,...,(pi_)_{n_3}
        
COMMENTS:
    Lemma 2.6 (Yu) of the Reference, in this instance:
        Let alpha_1,...,alpha_n:= eps, p_1,...,p_{n_1}, pi_1,...,pi_{n_2}, (pi_)_1,...,(pi_)_{n_3} (n >= 2) be nonzero algebraic numbers. 
        Put L:= Q(alpha_1,...,alpha_n), d:= [L:Q]. Let b_1,...,b_n be rational integers.
        Let P be a prime ideal of L, lying above the rational prime p of I. Then (given some lengthy conditions), 
            ord_P((Kappa_)_i*) < C3_i_ + C4_i_*Log(B), 
        where C3_i_, C4_i_, B are explicitly given, and 
            (Kappa_)_i*:= n*Log_p(eps) - sum_{j in IUstar} u_j*Log_p(p_j) + sum_{j in I} c_j*Log_p(pi_j) + sum_{j in I_} c_j*Log_p((pi_)_j)
    This result is used to compute an upper bound for a given linear form in logarithms, (Kappa_)_i*

    This algorithm computes C3_i_, C4_i_ 

EXAMPLE:
    > S:= [2,3,5,7,11,13];
    > D:=10;
    > K<sqrtD>:= QuadraticField(D);
    > R:= RingOfIntegers(K);
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
    > C:= IIPrimeXYZ2(SplitPrimes,D); 
    > II_:= C[4];
    > II_;
    [
        [
            <3, 2, sqrtD + 1, -sqrtD + 1, Prime Ideal of R
            Two element generators:
                3
                $.2 + 4>
        ],
        [
            <13, 2, 58973*sqrtD - 186489, -58973*sqrtD - 186489, Prime Ideal of R
            Two element generators:
                13
                $.2 + 7>
        ]
    ]
    > A0:= AlphasXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > A:= A0[1];
    > IUstar:= [ 2, 5, 7, 11 ];
    > p:= 13;
    > C3_i_, C4_i_, e_P, f_P, n:= Kappa_BoundXYZ2(II_,IUstar,p,D);
    > C3_i_;                                                      
    7.14481146325023399195655516030E36
    > C4_i_;
    8.18078100485534190662722167859E34
    > e_P;
    1
    > f_P;
    1
    > n;
    7

*/


function Kappa_BoundXYZ2(II_,IUstar,p,D)
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    R:= RingOfIntegers(K);      // generates ring of integers of K; every element of R is of the form a + Nu*b, where a,b are integers, Nu is defined as below
    
    eps:= FundamentalUnitXYZ2(D);
    
    e_P:= 1;
    f_P:= 1;  // # p splits in K, where [K:Q] = 2 = (r_P)*(e_P)*(f_P), and r_P = 2, the ramification and inertial degrees of P above p are both 1

    terms:= [eps];
    for q in II_[1] do
        Append(~terms,q[3]);    // appends pi_j for j in I
    end for;
    for q in II_[2] do
        Append(~terms,q[4]);    // appends (pi_)_j for j in I_
    end for;
    for k in IUstar do
        Append(~terms,k);       // appends p_j for j in IUstar
    end for;
    
    n:= #terms;         // computes the number of terms in (Kappa_)_i*
    c1np,a1:= C1pnXYZ2(p,n);    // computes constants C1np, a1, as per Lemma 2.6

    h:= ClassNumber(K);
    
    Q:= 2;
    while (h*p*(p^(f_P) - 1)) mod Q eq 0 do
        Q:= NextPrime(Q);       // computes smallest prime q such that q does not divide h*p*(p^(f_P) - 1)
    end while;
    q:= Q;
    
    Vi:= [];    // stores V_i = Max( h(alpha_i), (f_P)*Log(p)/d ), where h(alpha_i) is the logarithmic height of some element alpha_i of terms, d:= [K:Q] = 2
    for t in terms do
        if t eq 2 then 
            Append(~Vi,Max(1, f_P*Log(p)/2));    // de Weger defines h(2) = 1
        else
            Append(~Vi,Max(AbsoluteLogarithmicHeight(t), f_P*Log(p)/2));
        end if;
    end for;
    
    Vplus:= Max(1,Max(Vi));     // computes Max(1,V_1,...,V_n) as per Yu's Lemma, de Weger
    V1Vn:= &*[v : v in Vi];     // computes (V_1)*...*(V_n) to be used in bound
    
    U:= (c1np)*(a1^n)*(n^(n+5/2))*(q^(2*n))*(q-1)*((Log(n*q))^2)*((p^(f_P))-1)*((2+(1/(p-1)))^n)*((f_P*(Log(p)/2))^(-(n+2)))*V1Vn*(Log(4*2*Vplus) + f_P*Log(p)/(8*n));
    
    C3_i_:= U*Log(4*2);
    C4_i_:= U/(6*n);

    return C3_i_, C4_i_, e_P, f_P, n;
end function;