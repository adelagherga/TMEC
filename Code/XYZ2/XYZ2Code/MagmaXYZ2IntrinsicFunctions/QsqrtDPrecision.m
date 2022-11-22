// QsqrtDPrecision.m

/*
INPUT:
    a:= a rational field element, making up part of x, where x is an element of K:= Q(Sqrt(D)) and 
        x:= a + Sqrt(D))*b, a and b not both 0
    D:= a nonsquare rational integer defining K:= Q(Sqrt(D))
    b:= a rational field element, making up part of x, where x is an element of K:= Q(Sqrt(D)) and 
        x:= a + Sqrt(D))*b, a and b not both 0

OUTPUT:
    The numerical approximation of x:= a + Sqrt(D))*b, a and b not both 0, with large enough precision so that x != 0
    
COMMENTS:
    For x:= a + Sqrt(D))*b, a and b not both 0, the preset RealField precision in MAGMA sometimes returns x == 0.
    This function increases the precision until this value returns x != 0. 

EXAMPLE:
    > K<SqrtD>:= QuadraticField(D);
    > x:= -1221880517890354122164659057717 + 216000000000000000000000000000*SqrtD;
    > a:= x[1];
    > b:= x[2];
    > Numericalx:= a + b*Sqrt(D));
    > Numericalx;
    0.000000000000000000000000000000
    > QsqrtDPrecision(a,D,b);
    -0.2717200645711272954940795898437500000000

*/


intrinsic QsqrtDPrecision(a::FldRatElt, D::RngIntElt, b::FldRatElt) -> FldRealElt
    {Returns the numerical approximation of x:= a + Sqrt(D))*b, a and b not both 0, with large enough precision so that x != 0.}
    prc:= 30;   // initializes the precision of the RealField to 30
    x:= a + b*Sqrt(D);
    while x eq 0 do
        x:= a + b*(RealField(prc) ! (Sqrt(D)));
        prc:= prc + 10;         // increases the precision on the RealField until MAGMA returns x != 0
    end while;
    
    return x;
end intrinsic;