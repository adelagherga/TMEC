// LambdaInitialSortXYZ2.m

/*
INPUT:
    A:= [ q, q_, alpha, tplus1 ] as output by MaximalC12BoundXYZ2.m giving the largest C12star value, where 
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
    
    p:= prime in IU
    mu:= the precision on the p-adic field Q_p

OUTPUT:
    terms:= the terms {eps/eps_, pi_j/(pi_)_j} appearing Lambda*, excluding p0, where
        Lambda*:= nstar*log_p(eps/eps_) + sum_{j in I} cstar_j*log_p( pi_j/(pi_)_j ) - sum_{j in I_} cstar_j*log_p( pi_j/(pi_)_j )
    p0:= the term of Lambda* giving minimal ord_p(log_p(p0))
    m0:= ord_p(log_p(p0))
    
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
    > a:= A[3];
    > a;
    2*sqrtD
    > IU:= IUFactorsXYZ2(II_,a,S,D);
    > IU;
    [ 2, 3, 7, 13 ]
    > p:= 2;
    > mu:= 100;
    > terms, p0, m0:= LambdaInitialSortXYZ2(A,p,mu,D);
    > terms;
    [ 1/22*(-3*sqrtD - 23) ]
    > p0;
    1/2*(-sqrtD - 3)
    > m0;
    2
        
*/


function LambdaInitialSortXYZ2(A,p,mu,D)
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    
    eps:= FundamentalUnitXYZ2(D);
    eps_:= Conjugate(eps);
    q:= A[1];
    q_:= A[2];
    
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

    Qp:= pAdicField(p,mu+20);
    terms:= [ eps/eps_ ];       // stores all terms appearing in the linear form in logarithm, Lambda*
    
    for q in I do
        Append(~terms, q[3]/q[4]);
    end for;
    for q_ in I_ do
        Append(~terms, q_[3]/q_[4]);
    end for;
    
    P0:= [Ordp( LambdaLogpXYZ2(i,p,D,mu) , p, D ) : i in terms];        // computes ord_p(log_p(p_i)) for each term Lambda* as a p-adic integer
    m0, p0index:= Min(P0);      // finds the minimal ord_p(log_p(p_i)), m0, and it's location, p0index, in P0
    p0:= terms[p0index];        // computes p0 such that m0:= ord_p(log_p(p0))
    Exclude(~terms, p0);        
    
    return terms, p0, m0;
end function;

