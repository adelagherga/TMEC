// Kappa_InitialSortXYZ2.m

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
    IUstar:= IU cat { i | t+1 <= i <= s, v_i != 0}, giving the largest C12star value, where
        IU:= the set of primes p_i of S such that 
            G_{alpha} = prod_{i in IU}(p_i)^(u_i)
        v_i := (1/2)*hstar*ord_{p_i}(4D) for p_i in tplus1
    p:= prime in I_
    P:= (unique) choice of prime ideal in K lying above p
    mu:= the precision on the p-adic field Q_p
    D:= squarefree part of x, where
        x:= Du^2 and x + y = z^2 
    
OUTPUT:
    terms:= the terms {eps_, p_j, pi_j, (pi_)_j} appearing Kappa_*, excluding p0, where
        Kappa_*:= nstar*log_p(eps) - sum_{j in IUstar} ustar_j*log_p( p_j ) + sum_{j in I} cstar_j*log_p( pi_j ) + sum_{j in I_} cstar_j*log_p( (pi_)_j )
    p0:= the term of Lambda* giving minimal ord_p(log_p(p0))
    m0:= ord_p(log_p(p0))
    
REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,5,7,11,13];
    > D:= 10;
    > K<sqrtD>:= QuadraticField(D);
    > R:= RingOfIntegers(K);
    > mu:= 100;
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
    > C:= IIPrimeXYZ2(SplitPrimes,D); 
    > II_:= C[4];
    > I:= II_[1];
    > I_:= II_[2];
    > I;
    [
        <3, 2, sqrtD + 1, -sqrtD + 1, Prime Ideal of R
        Two element generators:
            3
            $.2 + 4>
    ]
    > I_;
    [
        <13, 2, 58973*sqrtD - 186489, -58973*sqrtD - 186489, Prime Ideal of R
        Two element generators:
            13
            $.2 + 7>
    ]
    > p:= I_[1][1];
    > P:= I_[1][5];
    > p;
    13
    > P; 
    Prime Ideal of R
    Two element generators:
        13
        $.2 + 7
    > A0:= AlphasXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > A:= A0[2];
    > FA:= MaximalC12BoundXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > IUstar:= FA[7];
    > IUstar;        
    [ 2, 5, 7, 11 ]
    > terms, p0, m0:= Kappa_InitialSortXYZ2(A, IUstar, p, P, mu, D);
    > terms;
    [ sqrtD + 1, -58973*sqrtD - 186489, 2, 5, 7, 11 ]
    > p0;
    sqrtD + 3
    > m0;
    1
        
*/


function Kappa_InitialSortXYZ2(A, IUstar, p, P, mu, D)
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
    SqrtDModp:= HenselLiftXYZ2(p,P,D,mu);       // computes Sqrt(D) in Qp with precision mu 
    
    terms:= [ eps ];    // stores all terms appearing in the linear form in logarithm, Kappa_*
    P0:= [Valuation(pAdicLog( (Qp!eps_[1] + Qp!eps_[2]*SqrtDModp) ,p) )];        // computes ord_p(log_p(p_i)) for eps as a p-adic integer
    
    for q in I do
        Append(~terms, q[3]);
        Append(~P0, Valuation(pAdicLog( (Qp!q[3][1] + Qp!q[3][2]*SqrtDModp) ,p) ) );    // computes ord_p(log_p(p_i)) for each term pi_j, j in I, as a p-adic integer
    end for;
    for q_ in I_ do
        Append(~terms, q_[4]);
        Append(~P0, Valuation(pAdicLog( (Qp!q_[4][1] + Qp!q_[4][2]*SqrtDModp) ,p) ) );  // computes ord_p(log_p(p_i)) for each term (pi_)_j, j in I_, as a p-adic integer
    end for;
    
    for j in IUstar do
        Append(~terms, K!j);
        Append(~P0, Valuation( pAdicLog(Qp!j, p) ) );   // computes ord_p(log_p(p_i)) for each term p_j, j in IUstar, as a p-adic integer
    end for;
    
    m0, p0index:= Min(P0);      // finds the minimal ord_p(log_p(p_i)), m0, and it's location, p0index, in P0
    p0:= terms[p0index];        // computes p0 such that m0:= ord_p(log_p(p0))
    Exclude(~terms, p0);        
    
    return terms, p0, m0;
end function;
