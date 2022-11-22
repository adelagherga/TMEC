// pAdicLog.m

/*
INPUT:
    x:= a p-adic element belonging to a finite extension of Q_p
    p:= a rational prime

OUTPUT:
    The p-adic logarithm of x
    
COMMENTS:
    MAGMA built-in Log function gives an error when unit is not a 1-unit

EXAMPLE:
    > p:= 2;
    > Qp:= pAdicRing(p,20);
    > x:= Qp!22;
    > pAdicLog(x,p);
    -23221*2^2 + O(2^19)
   
    > Log(x);  
    >> Log(x);
          ^
    Runtime error in 'Log': Log is only defined on units
    
*/


intrinsic pAdicLog(x::FldPadElt, p::RngIntElt) -> RngPadElt
    {Returns the p-adic logarithm of x.}
    
    r:= Valuation(x);   // computes the p-adic order of o
    if r gt 0 then
        x:= x/(p^r);    // updates o to be the p-free-part of o, since log_p(p^r) = r*log_p(p) = r*0 = 0
    end if;
    return Log(x);
end intrinsic;    