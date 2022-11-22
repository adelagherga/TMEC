// SymmetricCaseZerosXYZ2.m

/*
INPUT:
    b:= an integer partially determining the recurrence relation
    c:= an integer partially determining the recurrence relation
    u_0:= an integer; the 0th term of the binary recurrence sequence
    u_1:= an integer; the 1st term of the binary recurrence sequence
    p:= a rational prime
    n:= an integer > 0
            
OUTPUT:
    Is:= [i_1,...,i_k], where i_j is the index of G(i) such that G(i_j):= 0 mod (p^n);
        G(i):= the binary recurrence sequence defined by initial conditions u_0 and u_1 and recurrence relation
            u_{n+2} = b*( u_{n+1} )+c*( u_n )
    
COMMENTS:
    Used in the Symmetric Case for SUnitXYZ2.m 

EXAMPLE:
    > b:= 3;
    > c:= 1;
    > u_0:= 2;
    > u_1:= 1;
    > p:= 37;
    > n:= 20;
    > SymmetricCaseZerosXYZ2(b,c,u_0,u_1,p,n);
    [ 7977515069008045655059677577607, 19851222897834046082276943035694, 
    31724930726660046509494208493781, 43598638555486046936711473951868 ]
    
    > // verification:
    > i:= SymmetricCaseZerosXYZ2(b,c,u_0,u_1,p,n)[1];
    > BinaryRecurrenceSequence(b,c,u_0,u_1,i,p^n);
    0
    
*/


function SymmetricCaseZerosXYZ2(b,c,u_0,u_1,p,n)
    periodn0:= BinaryRecurrenceSequencePeriod(b,c,u_0,u_1,p);   // computes the period mod p of the binary recurrence sequence, G(i)
    Is0:= [];   // stores indices i for which the binary recurrence sequence, G(i), is divisible by p in its period mod p; ie. where G(i) = 0 mod p
    
    for i in [0..(periodn0 - 1)] do   // runs through all elements in the period mod p
        if BinaryRecurrenceSequence(b,c,u_0,u_1,i,p) eq 0 then  // if G(i) mod p is 0 for all i in the period mod p
            Append(~Is0,i);
        end if;
    end for;
    
    if IsEmpty(Is0) eq true then
        Is:= [];        // if no elements in the period mod p are divisible by p, returns the empty set, Is; ie there are no i such that G(i) = 0 mod p

    elif n eq 1 then  
        Is:= Is0;       // if n == 1, returns Is, the set of indices for which the binary recurrene sequence, G(i), is 0 mod p in its period mod p, [1, ..., periodn0]
    else        // if n > 1 and at least one index, i, in the range [1, ... , periodn], is divisible by p, it may also be divisible by p^n, where
                    // periodn:= the period mod (p^n)
                    
        Isc:= Is0;      // generates the set of all indices in the range [1,...,periodn] which are divisible by p, where
                            // periodn:= the period mod (p^n0)
                            // everything not divisible by p cannot be divisible by p^n for n > 1, so we may ignore those indices
        n0:= 2;         // computes SymmetricCaseZerosXYZ2 successively, for each exponent n0 <= n
        while n0 le n do
            periodn:= BinaryRecurrenceSequencePeriod(b,c,u_0,u_1,p^n0);         // computes the period mod (p^n0), n0 > 1, of the binary recurrence sequence, G(i)
            t:= #(Is0)*periodn/periodn0 - #(Is0);       // computes the multiplication factor needed to increase the period of G(i) mod ( p^(n0-1) ) to the period of G(i) mod (p^n0)
                                                        // ie. this is ( the number of 0s in the period mod (p^(n0-1)) )*( the scale difference between the period mod p^n0 and the period mod p^(n0-1) )
                                                                // - ( the number of 0s in the period mod (p^(n0-1)) )
            for j in [1..t] do  
                Append(~Isc, Isc[j] + periodn0);        // computes the indices in the period mod p^n0 which are divisible by p^(n0-1) by extending the 0s from the period mod (p^(n0-1)) to that of p^n0
            end for;
            
            A:= Isc;    // generates a copy of Isc, the elements in the period mod p^n0 which are divisible by p^(n0-1)
            for i in A do
                if  BinaryRecurrenceSequence(b,c,u_0,u_1,i,p^n0) ne 0 then      // removes the index i if G(i) - which is 0 mod (p^(n0-1)) - is not also 0 mod (p^n0)
                    Exclude(~Isc, i);   // stores only those indices i for which G(i):= 0 mod (p^n0)
                end if;
            end for;
            
            if IsEmpty(Isc) eq true then
                break;  // terminates algorithm if none of the indices in the range of the period mod (p^(n0)) are 0 mod p^(n0); ie. since n0 <= n, there are no indices which are 0 mod (p^n)
            else
                periodn0:= periodn;     // updates the period; periodn0 represents the period mod (p^(n0-1))
                n0:= n0 + 1;
                Is0:= Isc;      // updates the set of indices which are 0 mod (p^(n0-1))
            end if;
        end while;

        Is:= Isc;
    end if;
    return Is;
end function;
