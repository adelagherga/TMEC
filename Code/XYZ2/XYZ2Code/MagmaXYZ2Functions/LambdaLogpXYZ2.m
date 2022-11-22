// LambdaLogpXYZ2.m

/*
INPUT:
    x:= an element of the quadratic field K:= Q(Sqrt(D))
    p:= a rational prime
    D:= a nonsquare rational integer defining K:= Q(Sqrt(D))
    mu:= the precision on the p-adic logarithm of x 
    
OUTPUT:
    The p-adic logarithm of x, Log_p(x)
    
COMMENTS:
    The built-in MAGMA function Log does not handle non-rational FldPadElt as the input for x
     
    Let Omega_p denote the completion of the algebraic closure of Qp and let z in Omega_p with ord_p(z-1) > 0. Define the p-adic logarithm by
        Log_p(z) = sum_{n = 1..infty}[(-1)^(n+1)(z-1)]/n.
    This function satisfies the usual logarithm properties, Log_p(zw) = Log_p(z) + Log_p(w), and can be extended to all of Omega_p* (the set of nonzero elements of Omega_p) by setting 
        Log_p(p):= 0,
        ie. by choosing the branch of the logarithm, p_branch:= 0, for which to send p to. 
    Specifically, every element w of Omega_p* can be written as w = (p^r)*(u)*(z) with r a rational number, u a root of unity, and ord_p(z-1) > 0, in which case
        Log_p(w) = Log_p(z)
    (In fact, there is an extension of the logarithm from ord_p(z-1) > 0 to all of Omega_p* for each choice of Log_p(p) Omega_p.)
    This implementation uses p_branch = 0 for consistency with other algorithms.     
    Alternatively, if ord_p(z) = 0 and ord_p(z-1) = 0, there exists an integer k such that ord_p(z^k-1) > 0. Hence Log_p(z) may be extended to z with ord_p(z) = 0 by setting
        Log_p(z) = (1/k)*Log_p(1+(z^k-1)), where k is a positive integer such that ord_p(z^k-1) > 0
    
    In computing the above mentioned Taylor series, there will be factors p in the denominators p of the terms. 
    To find the first mu p-adic digits of Log_p(1+z), the first k terms must be taken into account, where k is the smallest integer satisfying
        k*ord_p(z)- Log(k)/Log(p) >= mu.
        
    Based on the algorithm found in Chapter 2 of the Reference

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.
    
    "P-adic Exponential Function." Wikipedia. Wikimedia Foundation, n.d. Web. <https://en.wikipedia.org/wiki/P-adic_exponential_function>.
    
EXAMPLE:
    > D:= 5;
    > K<sqrtD>:= QuadraticField(D);
    > mu:= 100;
    > x:= sqrtD;
    > p:= 5;
    > Log_p:= LambdaLogpXYZ2(x,p,D,mu);
    > Log_p;
    0
    
    > D:= 5;
    > K<sqrtD>:= QuadraticField(D);
    > mu:= 100;
    > x:= sqrtD;
    > p:= 2;
    > Log_p:= LambdaLogpXYZ2(x,p,D,mu);
    > Log_p;
    3126847970517100631362383567968884243527390708523334497635635414314956998/6414924694381721303722858446525
    
    > D:= 5;
    > K<sqrtD>:= QuadraticField(D);
    > mu:= 100;
    > x:= D;
    > p:= 2;
    > Q:= pAdicField(p,mu); 
    > Q;
    2-adic field mod 2^100
    > Log_p:= LambdaLogpXYZ2(x,p,D,mu);
    > Log_p;
    6253695941034201262724767135937768487054781417046668995271270828629913996/6414924694381721303722858446525
    > Q! 6253695941034201262724767135937768487054781417046668995271270828629913996\6414924694381721303722858446525;
    -8671155750987115445059024481*2^2 + O(2^102)
    > x:= 5;
    > Log(Q!x);
    -8671155750987115445059024481*2^2 + O(2^100)
    
*/


function LambdaLogpXYZ2(x,p,D,mu)
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    
    ord_p:= Ordp(x, p, D);      // computes the valuation of x at p, ord_p(x), where ord_p(x) is a rational number with denominator in {1,2}
                                    // ie. x = p^(ord_p(x))*u, where u is not divisible by p
                                    // since log_p(p^r) = r*log_p(p) = r*0 = 0, need only the p-free part of x, u
    if Denominator(ord_p) ne 1 then     // if ord_p(x) has denominator != 1, D must be divisible by p
        x:= x/sqrtD;    // updates x so that ord_p(x) is an integer
        ord_p:= Ordp(x, p, D);  // updates ord_p(x)
    end if;

    if ord_p gt 0 then
        x:= x/p^(ord_p);        // updates x to give the p-free part
    end if;
    
    k:= 1;
    while Ordp(x^k-1, p, D) le 0 do
        k:= k + 1;      // computes the integer k such that ord_p(x^k-1) > 0; if ord_p(x) = 0 and ord_p(x-1) = 0, such an integer k exists  
    end while;
    ordpk:= Ordp(x^k-1, p, D);
    
    t:= 1;
    while t*(ordpk) - Log(t)/Log(p) lt mu do    // computes the number of terms required for the desired precision
        t:= t + 1;
    end while;
    
    t:= t + 20;         // increases precision by 20 to avoid round-off errors

    Log_p:= &+[ ( (-1)^(n+1) )*( (x^k-1)^n )/(n*k) : n in [1..t]  ];    // computes the Taylor expansion defining the p-adic logarithm
    
    return Log_p;
end function;
