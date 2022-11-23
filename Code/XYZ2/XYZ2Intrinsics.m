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

// NonEmptyMaxMinProductSum.m

/*
INPUT:
    C:= [c_1,...,c_n], an arbitrary sequence of elements c_i, potentially empty
        ie. C:= [ ]

OUTPUT:
    if C != [ ]:
        The maximal element of C, Max(C)
    if C = [ ]:
        Returns 0

COMMENTS:
    Built-in MAGMA function Max(C) returns error if C == [ ]

EXAMPLE:
    > C:= [ ];
    > Max(C);

    >> Max(C);
          ^
    Runtime error in 'Max': Argument 1 is not non-empty

    > NonEmptyMax(C);
    0

*/


intrinsic NonEmptyMax(C::SeqEnum) -> .  // any type
    {Returns Max(C) if C is nonempty, otherwise returns 0.}
    if IsEmpty(C) eq true then
        MaxC:= 0;       // consistent with SAGE output
    else
        MaxC:= Max(C);
    end if;
    return MaxC;
end intrinsic;


/*
INPUT:
    C:= [c_1,...,c_n], an arbitrary sequence of elements c_i, potentially empty
        ie. C:= [ ]

OUTPUT:
    if C != [ ]:
        The minimal element of C, Min(C)
    if C = [ ]:
        Returns 0

COMMENTS:
    Built-in MAGMA function Min(C) returns error if C == [ ]

EXAMPLE:
    > C:= [ ];
    > Min(C);

    >> Min(C);
          ^
    Runtime error in 'Min': Argument 1 is not non-empty

    > NonEmptyMin(C);
    0

*/


intrinsic NonEmptyMin(C::SeqEnum) -> .  // any type
    {Returns Min(C) if C is nonempty, otherwise returns 0.}
    if IsEmpty(C) eq true then
        MinC:= 0;       // consistent with SAGE output
    else
        MinC:= Min(C);
    end if;
    return MinC;
end intrinsic;


/*
INPUT:
    C:= [c_1,...,c_n], an arbitrary sequence of elements c_i, potentially empty
        ie. C:= [ ]

OUTPUT:
    if C != [ ]:
        The product of all element of C, &*(C)
    if C = [ ]:
        Returns 1

COMMENTS:
    Built-in MAGMA function &*(C) returns error if C == [ ]

EXAMPLE:
    > C:= [ ];
    > &*(C);

    >> &*(C);
       ^
    Runtime error in '&*': Illegal null sequence

    > NonEmptyProduct(C);
    1

*/


intrinsic NonEmptyProduct(C::SeqEnum) -> .      // any type
    {Returns the product of elements of C, &*(C), if C is nonempty, otherwise returns 1.}
    if IsEmpty(C) eq true then
        ProductC:= 1;   // consistent with SAGE output
    else
        ProductC:= &*C;
    end if;
    return ProductC;
end intrinsic;


/*
INPUT:
    C:= [c_1,...,c_n], an arbitrary sequence of elements c_i, potentially empty
        ie. C:= [ ]

OUTPUT:
    if C != [ ]:
        The sum of all element of C, &*(C)
    if C = [ ]:
        Returns 0

COMMENTS:
    Built-in MAGMA function &+(C) returns error if C == [ ]

EXAMPLE:
    > C:= [ ];
    > &+(C);

    >> &+(C);
       ^
    Runtime error in '&+': Illegal null sequence

    > NonEmptySum(C);
    0

*/


intrinsic NonEmptySum(C::SeqEnum) -> .
    {Returns the sum of elements of C, &+(C), if C is nonempty, otherwise returns 0.}
    if IsEmpty(C) eq true then
        SumC:= 0;       // consistent with SAGE output
    else
        SumC:= &+C;
    end if;
    return SumC;
end intrinsic;
