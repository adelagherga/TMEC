// ConvertpAdic.m

/*
INPUT:
    theta:= a p-adic element belonging to a finite extension of Q_p
    p:= a rational prime
    mu:= the precision on the p-adic field Q_p

OUTPUT:
    The rational integer theta^(mu) such that 0 <= theta^(mu) < p^mu with ord_p(theta - theta^(mu)) >= mu
    
COMMENTS:
    Used in the construction of lattices for SUnitXYZ.m and SUnitXYZ2.m 

EXAMPLE:
    > p:= 5;
    > mu:= 25;
    > Qp:= pAdicField(p, mu);
    > theta:= pAdicLog(Qp!3,p);
    > theta;
    -12997808195050206*5 + O(5^25)
    > ConvertpAdic(theta, p, mu);  
    233034182901702095
    
*/


intrinsic ConvertpAdic(theta::FldPadElt, p::RngIntElt, mu::RngIntElt) -> RngIntElt
    {Returns the rational integer theta^(mu) such that 0 <= theta^(mu) < p^mu with ord_p(theta - theta^(mu)) >= mu.}
    x:= IntegerRing()! theta;  // computes x = IntegerPart(theta), where theta = IntegerPart(theta) + O(p^precision)
    x:= x mod p^mu;
    if x lt 0 then
        x:= x + p^mu;
    end if;
    return x;
end intrinsic;     
