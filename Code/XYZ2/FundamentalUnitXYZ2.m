//FundamentalUnitXYZ2.m

/*
INPUT:
    D:= squarefree part of x, where
        x:= Du^2 and x + y = z^2

OUTPUT:
    u:= the unique fundamental unit of K = Q(Sqrt(D)) such that u > 1

COMMENTS:
    The unit group, U, of the number field K = Q(Sqrt(D)) is U = {+/- u^k, k in Z}, so that U is generated by -1 and powers of u
    Fundamental unit:= unique u when normalized so that u is the smallest unit among those u with the property that u > 1 

    This algorithm normalizes the MAGMA built-in FundamentalUnit function so that the fundamental unit u is such that
        u > 1
    This is necessary in the computation of G_{alpha}

REFERENCE:
    
    D.A. Marcus. Number Fields. Springer, 1977

EXAMPLE:
    > D:= 5;
    > FundamentalUnitXYZ2(D);
    1/2*(sqrtD + 1)
    
    > D:= 17;
    > FundamentalUnitXYZ2(D);
    sqrtD + 4

*/


function FundamentalUnitXYZ2(D)
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    
    u:= FundamentalUnit(K);     // computes the fundamental unit of K, potentially such that < 1
    u_:= Conjugate(u);  
    Numericalu:= u[1] + u[2]*Sqrt(D);           // the numerical approximation of u as a real number
    Numericalu_:= u_[1] + u_[2]*Sqrt(D);        // the numerical approximation of the conjugate of u as a real number
    
    if (Numericalu lt 0) then
        u:= -u;
    end if;
    if (Numericalu_ lt 0) then
        u_:= -u_;
    end if; 
    
    if (Abs(Numericalu) lt 1) and (Abs(Numericalu_) gt 1) then    // if u < 1 and u_ > 1
        Newu:= u_;
    else
        Newu:= u;
    end if;
    
    Numericalu:= Newu[1] + Newu[2]*Sqrt(D);
    if Numericalu lt 1 then     // verification that indeed u > 1; if false, there is an error in FundamentalUnitXYZ2.m
        print "Something is wrong in FundamentalUnitXYZ2: eps < 1.";
    else
        return Newu;
    end if;
end function;