// HenselLiftXYZ2.m

/*
INPUT:
    p:= prime of S which splits in K = Q(Sqrt(D))
    P:= (unique) choice of prime ideal in K lying above p
    D:= squarefree part of x, where
        x:= Du^2 and x + y = z^2 
    mu:= the precision on the p-adic field Q_p

OUTPUT:
    Sqrt(D) in Qp with precision mu
    
COMMENTS:
    A prime p splits in K:= Q(Sqrt(D)) iff D = n^2 mod p for some integer n; ie. Sqrt(D) = n mod p

    MAGMA built-in functions are able to compute Sqrt(D) in Qp, so this algorithm is not a direct implementation of Hensel's Lemma.
    This algorithm ensures that MAGAMA's built-in function makes a choice of Sqrt(D) mod p that is consistent with that used to compute the ideal P above p
    
EXAMPLE:
    > S:= [2,3,5,7,11,13];
    > D:= 5;
    > K<sqrtD>:= QuadraticField(D);
    > R:= RingOfIntegers(K);
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
    > C:= IIPrimeXYZ2(SplitPrimes,D); 
    > II_:= C[1];
    > I:= II_[1];
    > I_:= II_[2];
    > I;
    [
        <11, 1, 1/2*(-3*sqrtD - 1), 1/2*(3*sqrtD - 1), Principal Prime Ideal of R
        Generator:
            -3*$.2 + 1>
    ]
    > I_;
    []
    > p:= I[1][1];
    > P:= I[1][5];
    > p;
    11
    > P;
    Principal Prime Ideal of R
    Generator:
        -3*$.2 + 1
    > mu:= 100;
    > x:= HenselLiftXYZ2(p,P,D,mu);
    > x;
    1708782003238830266078387337961468271155513215387929503265756220983592669840974\
    4944538057329069788358770658886982864610133609584986976391093482880346371451921\
    8272598872920414753609320568923887571107577915756850052427504039588347191665862\
    393403591257246299301710

*/


function HenselLiftXYZ2(p,P,D,mu)
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    R:= RingOfIntegers(K);      // generates ring of integers of K
    Qp:= pAdicField(p,2*mu+50);
    
    x:= IntegerRing()! (Sqrt(Qp!D));    // computes Sqrt(D) in Qp, ignoring the O(p^mu) part
    x:= x mod p^(2*mu+50);      // rewrites x = Sqrt(D) with desired precision
    choicex:= Basis(P)[2][1];   // computes the choice made from the two possibilities for Sqrt(D) mod p; this choice determines the unique choice of P above p
    
    if choicex ne x mod p then  // if the user-defined choice of Sqrt(D) mod p differs from MAGMA's choice of Sqrt(D) mod p, shifts x = Sqrt(D) accordingly to match 
        x:= -x mod p^(2*mu+50);
    end if;
    
    return x;
end function;
