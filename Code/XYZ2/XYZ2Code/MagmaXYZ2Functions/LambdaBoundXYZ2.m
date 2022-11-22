// LambdaBoundXYZ2.m

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
    p:= prime of IU
    D:= squarefree part of x, where
        x:= Du^2 and x + y = z^2

OUTPUT:
    C1_i:= the real constant such that ord_P(Lambda_i*) < C1_i + C2_i*Log(B)
    C2_i:= the real constant such that ord_P(Lambda_i*) < C1_i + C2_i*Log(B)
    e_P:= ramification index of the ideal P of K:= Q(Sqrt(D)) above p
    f_P:= residue class degree of the ideal P of K:= Q(Sqrt(D)) above p
    n:= the number of terms in Lambda_i*
        ie. the number of elements alpha_1,...,alpha_n:= eps/eps_, pi_1/(pi_)_1,...,pi_n/(pi_)_n

COMMENTS:
    Lemma 2.6 (Yu) of the Reference, in this instance:
        Let alpha_1,...,alpha_n:= eps/eps_, pi_1/(pi_)_1,...,pi_n/(pi_)_n (n >= 2) be nonzero algebraic numbers. 
        Put L:= Q(alpha_1,...,alpha_n), d:= [L:Q]. Let b_1,...,b_n be rational integers.
        Let P be a prime ideal of L, lying above the rational prime p of IU. Then (given some lengthy conditions), 
            ord_P(Lambda_i*) < C1_i + C2_i*Log(B), 
        where C1_i, C2_i, B are explicitly given, and 
            Lambda_i*:= n*Log_p(eps/eps_) + sum_{j in I} c_j*Log_p(pi_j/(pi_)_j) - sum_{j in I_} c_j*Log_p(pi_j/(pi_)_j)
    This result is used to compute an upper bound for a given linear form in logarithms, Lambda_i*

    This algorithm computes C1_i, C2_i 

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,5,7];
    > D:= 105;
    > K<sqrtD>:= QuadraticField(D);
    > R:= RingOfIntegers(K);
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
    > C:= IIPrimeXYZ2(SplitPrimes,D); 
    > II_:= C[1];
    > A0:= AlphasXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > A:= A0[3];
    > IUFactorsXYZ2(II_,A[3],S,D);
    [ 3, 7 ]
    > p:= 3;
    > C1_i, C2_i, e_P, f_P, n:= LambdaBoundXYZ2(II_,p,D);
    > C1_i;
    55977858513857767.9432537031855
    > C2_i;
    2243304968820184.46995538583390
    > e_P;
    2
    > f_P;
    1
    > n;
    2
    
*/


function LambdaBoundXYZ2(II_,p,D)
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    R:= RingOfIntegers(K);      // generates ring of integers of K; every element of R is of the form a + Nu*b, where a,b are integers, Nu is defined as below
    
    eps:= FundamentalUnitXYZ2(D);
    eps_:= Conjugate(eps);
    
    if IsRamified(p,R) eq true then     // if p is ramified in K = Q(Sqrt(D)), where [K:Q] = 2 = (r_P)*(e_P)*(f_P):
        e_P:= 2;                                // r_P:= 1, e_P:= 2, f_P:= 1
        f_P:= 1;
    else        // if p does not ramify in K = Q(Sqrt(D)), where [K:Q] = 2 = (r_P)*(e_P)*(f_P), since p is in IU (disjunct from I,I_, hence does not split):
        e_P:= 1;        // r_P:= 1, e_P:= 1, f_P:= 2
        f_P:= 2;
    end if;
    
    terms:= [eps/eps_];         // stores the terms of Lambda_i*:= n*Log_p(eps/eps_) + sum_{j in I} c_j*Log_p(pi_j/(pi_)_j) - sum_{j in I_} c_j*Log_p(pi_j/(pi_)_j)
    
    for q in II_[1] do
        Append(~terms,q[3]/q[4]);       // stores pi_j/(pi_)_j for j in I
    end for;
    for q in II_[2] do
        Append(~terms,q[3]/q[4]);       // stores pi_j/(pi_)_j for j in I_
    end for;
    
    n:= #terms;         // computes the number of terms in Lambda_i*
    c1np,a1:= C1pnXYZ2(p,n);    // computes constants C1np, a1, as per Lemma 2.6
    h:= ClassNumber(K);
    
    Q:= 3;
    while (h*p*(p^(f_P) - 1)) mod Q eq 0 do
        Q:= NextPrime(Q);       // computes smallest prime odd q such that q does not divide h*p*(p^(f_P) - 1)
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

    C1_i:= U*Log(4*2);
    C2_i:= U/(6*n);

    return C1_i, C2_i, e_P, f_P, n;
end function;