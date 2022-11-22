// Ordp.m

/*
INPUT:
    x:= an element of the quadratic field K:= Q(Sqrt(D))
    p:= a rational prime
    D:= a nonsquare rational integer defining K:= Q(Sqrt(D))

OUTPUT:
    The valuation of x at p, ord_p(x)
    
COMMENTS:
    The built-in MAGMA function Valuation does not handle FldQuadElt as the input for x
    This function uses the multiplicativity of valuations to compute ord_p(x) for a FldQuadElt element x

EXAMPLE:
    > D:= 2;
    > K<SqrtD>:= QuadraticField(D);
    > x:= 1/2 + SqrtD;
    > p:= 2;
    > Valuation(x,2);

    >> Valuation(x,2);
                ^
    Runtime error in 'Valuation': Bad argument types
    Argument types given: FldQuadElt[FldRat], RngIntElt
    > Ordp(x,p,D);
    -1

*/


intrinsic Ordp(x::FldQuadElt[FldRat], p::RngIntElt, D::RngIntElt) -> RngIntElt
    {Returns ord_p(x) for a FldQuadElt element x at the prime p.}
    
    valx1:= Valuation(x[1],p);  // if x[1]:= 0, then Valuation returns Infinity, so it does not contribute to Min
    valx2:= Valuation(x[2],p);  // if x[2]:= 0, then Valuation returns Infinity, so it does not contribute to Min
    valD:= (1/2)*Valuation(D,p);        // ord_p(SqrtD)
    ord_p:= Min(valx1, valx2 + valD);   // computes ord_p(x) using additive/multiplicative properties of ord_p
    
    return ord_p;    
end intrinsic;
