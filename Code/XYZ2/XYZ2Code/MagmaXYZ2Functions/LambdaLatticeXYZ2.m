// LambdaLatticeXYZ2.m

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
    B:= [b_1,...,b_{#terms}, b_0], a matrix over ZZ with columns b_i,
        b_0,b_1,...,b_{#terms} form a basis for the p-adic approximation lattice Gamma_{Lambda,mu} associated to the linear form in logarithms, Lambda*
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
    > B, m0:= LambdaLatticeXYZ2(A,p,mu,D);
    > B;
    [                              1                               0]
    [1122320629957454298493490877891 1267650600228229401496703205376]
    > m0;
    2
   
*/


function LambdaLatticeXYZ2(A,p,mu,D)
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    Qp:= pAdicField(p,mu+20);
    
    terms, p0, m0:= LambdaInitialSortXYZ2(A,p,mu,D);
    
    t:= [-LambdaLogpXYZ2(t,p,D,mu)/LambdaLogpXYZ2(p0,p,D,mu) : t in terms];     // computes theta_i = -log_p(t_i)/log_p(p0) for each term t_i in Lambda*, excluding p0
                                                                                    // Nb. log_p is computed via the Taylor expansion, hence theta_i:= theta_i[1] + SqrtD*(theta_i[2])
    tm:= [];    // stores theta_i^{mu} = theta_i (mod p^mu) for for each term in Lambda*, excluding p0
    
    for j in t do
        if Ordp( K!j[1] , p, D ) lt mu then     // if ord_p(theta_i[1]) < mu
            tm0:= ConvertpAdic(Qp!j[1], p, mu);         // computes theta_i[1]^{mu} = theta_i[1] (mod p^mu) 
        else
            tm0:= 0;
        end if;
        
        if Ordp( K!j[2] , p, D ) lt mu then     // if ord_p(theta_i[2]) < mu
            tm1:= Sqrt(D)*ConvertpAdic(Qp!j[2], p, mu);         // computes SqrtD*(theta_i[2]^{mu}) = Sqrt(D)*(theta_i[2]) (mod p^mu) 
        else
            tm1:= 0;
        end if;
        Append(~tm, tm0+tm1);
    end for;
   
    Append(~tm,p^mu);   // adjoins p^mu to tm
    
    row_tm := Matrix(IntegerRing(),1,#tm,tm);   // writes theta_i^{mu} values as a row matrix
    I:= ScalarMatrix(#t,1);     // #(terms)-identity matrix
    Z:= ZeroMatrix(IntegerRing(),#t,1);        // creates (#terms)x1 zero matrix
    
    B := VerticalJoin(HorizontalJoin(I,Z),row_tm);      // creates matrix associated to basis of Gamma_{Lambda,mu}

    return B,m0;
end function;