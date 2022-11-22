// BinaryRecurrenceSequence.m

/*
INPUT:
    b:= an integer partially determining the recurrence relation
    c:= an integer partially determining the recurrence relation
    u_0:= an integer; the 0th term of the binary recurrence sequence
    u_1:= an integer; the 1st term of the binary recurrence sequence
    n:= an integer
    m:= an integer such that 
        m:= 0 if:
            the nth term of the binary recurrence sequence is to be returned
        m > 0 if:
            the nth term of the binary recurrence sequence modulo m is to be returned
            
OUTPUT:
    if m:= 0:
        The nth term in the binary recurrence sequence defined by initial conditions u_0 and u_1 and recurrence relation
            u_{n+2} = b*( u_{n+1} )+c*( u_n )
    if m > 0:
        The nth term modulo m in the binary recurrence sequence defined by initial conditions u_0 and u_1 and recurrence relation
            u_{n+2} = b*( u_{n+1} )+c*( u_n )
    
COMMENTS:
    This algorithm generates the nth term in the binary recurrence sequence defined by initial conditions u_0 and u_1 and recurrence relation
        u_{n+2} = b*( u_{n+1} )+c*( u_n )
    Used in the Symmetric Case for SUnitXYZ2.m 

EXAMPLE:
    > b:= 3;
    > c:= 3;
    > u_0:= 2;
    > u_1:= 1;
    > n:= 101;
    > BinaryRecurrenceSequence(b,c,u_0,u_1,n);  
    16158686318788579168659644539538474790082623100896663971001
    
*/


intrinsic BinaryRecurrenceSequence(b:: RngIntElt, c:: RngIntElt, u_0:: RngIntElt, u_1:: RngIntElt, n:: RngIntElt, m:: RngIntElt) -> RngIntElt
    {Returns the nth term of the binary recurrence sequence u_(n+2) = b*( u_(n+1) ) + c*( u_n ), possibly modulo m.}
    if (m eq 0) then
        R:= IntegerRing();      // generates Z if m == 0
    else
        R:= IntegerRing(m);     // generates Z/mZ if m > 0
    end if;
    F:= Matrix(R, [[0,1],[c,b]]);
    v:= Matrix(R, [[u_0], [u_1]]);
    
    return ((F^n)*v)[1,1];
end intrinsic;     

