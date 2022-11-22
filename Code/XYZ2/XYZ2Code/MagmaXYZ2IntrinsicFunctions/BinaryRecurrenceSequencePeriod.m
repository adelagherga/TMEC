// BinaryRecurrenceSequencePeriodXYZ2.m

/*
INPUT:
    b:= an integer partially determining the recurrence relation
    c:= an integer partially determining the recurrence relation
    u_0:= an integer; the 0th term of the binary recurrence sequence
    u_1:= an integer; the 1st term of the binary recurrence sequence
    m:= an integer
    
OUTPUT:
    The period of the binary recurrence sequence modulo m, where the binary recurrence sequence is defined by initial conditions u_0 and u_1 and recurrence relation
        u_{n+2} = b*( u_{n+1} )+c*( u_n )
    
COMMENTS:
    This algorithm is an improvement on the period function used in SAGE; given a binary recurrence sequence, a prime p, and an integer e > 0, 
    if the period mod p is 1, then equally, the period mod p may be j, for any integer j > 0. 
    This algorithm takes this into account, and therefore does not run into the time-out errors encountered in SAGE.
    
    NB. If d:= gcd(b,c) != 1, each term in the binary recurrence sequence is divisible by d, 
    and hence BinaryRecurrenceSequencePeriod(b,c,u_0,u_1,m) does not terminate for m such that d:= gcd(b,c) mod m != 1.
    
    Used in the Symmetric Case for SUnitXYZ2.m; in this algorithm, b:= eps + eps_ and c:= -(eps)*(eps_) = +/-1, and therefore gcd(b,c) = 1 and above issue is not encountered. 

EXAMPLE:
    > b:= 42590;
    > c:= -1;
    > u_0:= 0;
    > u_1:= 1716;
    > m:= 11^3;
    > BinaryRecurrenceSequence(b,c,u_0,u_1,m);  
    242
    
*/


intrinsic BinaryRecurrenceSequencePeriod(b:: RngIntElt, c:: RngIntElt, u_0:: RngIntElt, u_1:: RngIntElt, m:: RngIntElt) -> RngIntElt
    {Returns the period of the binary recurrence sequence u_(n+2) = b*( u_(n+1) ) + c*( u_n ) modulo m.}
    R:= IntegerRing(m);
    A:= Matrix(R, [[0,1],[c,b]]);
    w:= Matrix(R, [[u_0], [u_1]]);
    Fac:= Factorization(m);     // generates all prime factors of m:= ( (p_1)^(e_1) )* \cdots *( (p_r)^(e_r) )
    Periods:= [];       // stores the period mod (p^e), where p^e is a factor of m

    for i in Fac do
        p:= i[1];
        e:= i[2];
        
        // generates the setup to compute the period mod p
        F:= ChangeRing(A, IntegerRing(p));      // writes A as a matrix, F, in GL_2(F_p); the period mod p divides the order of the matrix A over GL_2(F_p), ie. it divides the order of F
                                                    // GL_2(F_p):= the general linear group:= the group of invertible 2 x 2 matrices with entires in F_p under matrix multiplication
        v:= ChangeRing(w, IntegerRing(p));      // writes w as a vector, v, in GL_2(F_p) 
        p1fac:= Factorization(p-1);     // generates all prime factors of p-1; the order of F divides either p(p-1) or (p-1)*(p+1) 
        FF:= F^(p-1);   
        
        // determines the order of F (the period mod p divides this)
        if FF*v eq v then       // checks if the order of F divides (p-1)
            M:= p-1;    
            Mfac:= p1fac;       // computes the prime factors of M:= (p-1); the order of F divides one of these
        
        elif Trace(FF) eq 2 then        // checks if Trace(F^(p-1)) = 2; if so, it must be a triangular matrix of the form [[1,a],[0,1]]
                                            // the order of the subgroup of matrices of triangular matrices of the form [[1,a],[0,1]] is p
                                            // hence the order of F divides p*(p-1), and therefore the order of F is a multiple of p
            M:= p-1;    
            Mfac:= p1fac;       // computes the prime factors of M:= (p-1)
            F:= F^p;    // replaces F by F^p; F has order dividing p*(p-1), hence F^p has order dividing (p-1)
        
        else    // checks if the order of F divides (p-1)*(p+1)
            M:= (p+1)*(p-1);
            p2fac:= Factorization(p+1);         // factors (p+1)
            Mfac:= p1fac;       // stores the prime factors of M:= (p-1)*(p+1)
            Mfac_dic:= [m[1] : m in Mfac];      // lists prime factors of (p-1) 
            for i in p2fac do   // combines the prime factors of the (p+1) and (p-1) terms for speed to determine the factorization of (p-1)*(p+1)
                if (i[1] in Mfac_dic) eq true then
                    Mfac[Index(Mfac_dic,i[1])][2]:= Mfac[Index(Mfac_dic,i[1])][2] + i[2] ;
                else
                    Append(~Mfac, i);
                    Append(~Mfac_dic, i[1]);
                end if;
            end for;
        end if;
        
        C:= [];         // stores the prime factors of M in such a way that every factor has multiplicity 1; ie. if M:= (i_1)* \cdots *(i_k), where 
                            // the i_j denote not necessarily distinct prime factors, then C:= [i_1,...,i_k] 
        
        for i in [1..#Mfac] do  // runs through all factors of M
            for j in [1.. Mfac[i][2]] do
                Append(~C, Mfac[i][1]);         // appends each prime factor of M with multiplicity 1
            end for;
        end for;

        // determines the period mod p
        Mfac:= C;       // rewrites the prime factors of M in such a way that every factor has multiplicity 1
        M0:= M;         // generates the initial value that the period divides, M0
        for i in Mfac do
            b:= IntegerRing()! (M0/i);  // computes the integer M/i where i is in Mfac:= [i_1,...,i_k]; the period divides M:= (i_1)* \cdots *(i_k), therefore 
                                            // since (F^M)*v == v, if ( F^(M/i) )*v == v, then the period divides (M/i)
            if (F^b)*v eq v then
                M0:= b;         // if ( F^(M/i) )*v == v, updates the value that the period must divide, M0 
            end if;
        end for;
        perp:= M0;      // stores the period mod p           
        
        // determines the period mod p^e by steping up by multiples of p
        F:= ChangeRing(A, IntegerRing(p^e));       // writes A as a matrix, F, in GL_2( F_{p^e} ); the period mod (p^e) divides the order of the matrix A over GL_2( F_{p^e} ), ie. it divides the order of F
        v:= ChangeRing(w, IntegerRing(p^e));       // writes w as a vector, v, in GL_2( F_{p^e} )
        FF:= F^(perp);
        
        perpe:= 0;      // initializes the period mod (p^e), perpe
        if perp eq 1 then
            PeriodCheck:= true;         // if the period mod p is 1, the period mod p may equally be j, for any integer j > 0
        end if;
        
        if FF*v eq v then
            perpe:= perp;       // if ( F^(perp) )*v == v, the period mod (p^e) is the same as the period mod p
        else
            tries:= 0;
            while tries le (4*e-4) do   // the period mod (p^e) divides period(p)*( p^(4e-4) ):= perp*( p^(4e-4) ); runs through all divisors of perp*( p^(4e-4) ) by successively multiplying perp by powers of p
                tries:= tries + 1;
                FF:= FF^p;
                if FF*v eq v then
                    perpe:= perp*(p^tries);       // if ( F^( perp*(p^tries) ) )*v == v, the period mod (p^e) is perp*( p^tries )
                    break;
                elif (tries eq (4*e-4)+1) and (perpe eq 0) and (PeriodCheck eq true) then       // if the period has not been found in the alloted time and perp = 1 (ie. period(p) = 1)
                    perp:= perp + 1;    // resets perp:= perp + 1; ie. the period mod p, perp, is 1, hence the period mod p, perp, may equally be j, for any integer j > 0
                    FF:= F^perp;        // computes F^(perp) with this new value of the period mod p, perp
                    tries:= 0;  // resets tries so that the loop is re-entered
                end if;
            end while;
        end if;
        Append(~Periods, perpe);
    end for;
           
    period:= 1;
    for perpe in Periods do
        period:= Lcm(perpe,period);     // computes the lcm of the periods mod (p_i)^(e_i) for all distinct prime factors p_i dividing m
    end for;
    return period;
end intrinsic;