//SFactors.m

/*
INPUT:
    n:= a rational integer
    S:= {p_1,...,p_s}, p_i primes in S

OUTPUT:
    rem:= the remainder of n after dividing by all elements of S
        ie. rem = n/(ord_{p_1}(n)*...*ord_{p_s}(n))
    Factors:= [[p_1,count_1],...,[p_s,count_s]], where
        p_i:= a prime in S
        count_i:= ord_{p_i}(n)
    
COMMENTS:
    Used to avoid the expensive built-in MAGMA Factorization function

EXAMPLE:
    > n:= 34;
    > S:= [2,5,11];
    > rem, Factors:= SFactors(n,S);
    > rem;
    17
    > Factors;
    [
        [ 1, 2 ],
        [ 0, 5 ],
        [ 0, 11 ]
    ]
    
*/


intrinsic SFactors(n::RngIntElt, S::SeqEnum[RngIntElt]) -> RngIntElt, SeqEnum[RngIntElt]
    {Returns the remainder of n after dividing by all elements of S, along with the set [ [p, ord_p(n)] : p in S].}
    
    Sort(~S);
    rem:= IntegerRing() ! n;
    Factors:= []; 
    for p in S do
        count:= 0;
        while IntegerRing()! rem mod p eq 0 do  // checks if rem is divisible by p
            rem:= IntegerRing()! rem/p;         // updates rem by dividing by p
            count:= count + 1;
        end while;
        Append(~Factors,[count,p]);
    end for;
    
    return rem, Factors;        
end intrinsic;