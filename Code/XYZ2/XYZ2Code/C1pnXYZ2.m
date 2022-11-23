// C1pnXYZ2.m

/*
INPUT:
    p:= a prime of S
    n:= an integer >= 2,
        represents the number of terms in the linear form in logarithms being evaluated by Lemma 2.6

OUTPUT:
    c1pn:= a constant used in Yu's Lemma (Lemma 2.6) to bound a given linear form in logarithms
    a1:= a constant depending on n; used in Yu's Lemma (Lemma 2.6) to bound a given linear form in logarithms

COMMENTS:
    Lemma 2.6 (Yu) of the Reference, in this instance:
        Let alpha_1,...,alpha_n (n >= 2) be nonzero algebraic numbers. Put L:= Q(alpha_1,...,alpha_n), d:= [L:Q]. 
        Let b_1,...,b_n be rational integers. Let P be a prime ideal of L, lying above the rational prime p. 
        Then (given some lengthy conditions), 
            ord_P((alpha_1)^(b_1)*...*(alpha_n)^(b_n) - 1) < c1pn*(a1)^n*B, 
        where c1pn, a1, B are explicitly given.
    This result is used to compute an upper bound for a given linear form in logarithms

    This algorithm computes c1pn and a1 
    
REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

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
    > n:= 7;
    > c1pn,a1:= C1pnXYZ2(p,n);                 
    > c1pn;
    1920625/18
    > a1;
    10.1482521595804355453450732264
    
*/


function C1pnXYZ2(p,n)
    C12n:= [768523, 476217, 373024, 318871, 284931, 261379, 2770008];   // C1(2,n) as in Lemma 2.6
    C13n:= [167881, 104028, 81486, 69657, 62243, 57098, 116055];        // C1(3,n) as in Lemma 2.6
    C1pn:= [87055, 53944, 42255, 36121, 32276, 24584, 311077];  // C1(p,n) as in Lemma 2.6

    if p eq 2 then
        if n lt 8 then
            c1pn:= C12n[n-1];
            a1:= 56*Exp(1)/15;
        else
            c1pn:= C12n[7];
            a1:= 8*Exp(1)/3;
        end if;
    elif p eq 3 then
        if n lt 8 then
            c1pn:= C13n[n-1];
            a1:= 56*Exp(1)/15;
        else
            c1pn:= C13n[7];
            a1:= 8*Exp(1)/3;
        end if;
    else
        if n lt 8 then
            c1pn:= (C1pn[n-1])*(2 + 1/(p-1))^2;
            a1:= 56*Exp(1)/15;
        else
            c1pn:= (C1pn[7])*(2 + 1/(p-1))^2;
            a1:= 8*Exp(1)/3;
        end if;
    end if;
    
    return c1pn,a1;
end function;