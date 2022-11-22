//SymmetricCaseXYZ2.m

/*
INPUT:
    A:= [ q, q_, alpha, tplus1 ] (as output by AlphasXYZ2.m), where 
        q:= < prod_i, [< p, h_i, pi, pi_, P , d_i > : j in [1..n] ] >, such that
            prod_i:= prod_{i in I}(P_i)^(d_i), the product of prime ideals of I which contribute to alpha
            p:= prime of S which splits in K
            h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
            pi:= generator of the ideal P^(h_i) = (pi)R
            (pi_):= conjugate of pi
            P:= (unique) choice of prime ideal in K lying above p
            n:= the length of Iset
            d_i:= the exponent on the prime ideal P_i in prod_i
        q_:= < prod_i, [< p, h_i, pi, pi_, P_ , d_i > : j in [1..n] ] >, such that
            prod_i:= prod_{i in I}((P_)_i)^(d_i), the product of the conjugate of the prime ideals of I_ which contribute to alpha
            p:= prime of S which splits in K
            h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
            pi:= generator of the ideal P^(h_i) = (pi)R
            (pi_):= conjugate of pi
            P_:= conjugate of the (unique) choice of prime ideal in K lying above p
            n:= the length of Iset
            d_i:= the exponent on the prime ideal P_i in prod_i
        alpha:= element of K, generating the ideal (alpha), where
            (alpha):= (2)^(b_0)*[prod_{i in I}(P_i)^(d_i)]*[prod_{i in I_}((P_i)_i)^(d_i)]*prod_{i = t + 1 ... s}P_i^(a_i)
        tplus1:= the product of primes of S from NonSplitPrimes whose prime ideals above contribute to alpha as 
            prod_{i = t + 1 ... s}P_i^(a_i)
    S:= [p_1,...,p_s], p_i primes in S
    D:= squarefree part of x, where
        x:= Du^2 and x + y = z^2
            
OUTPUT:
    eqns:= [[x_1,y_1,z_1],...,[x_n,y_n,z_n]], all S-unit equations x_i + y_i = (z_i)^2 corresponding to D, a, in the symmetric case
        ie. in the case where I == I_ == []
    
COMMENTS:
    In the case where I == I_ == [] (the "symmetric" case), none of the primes of S split in K, hence (since t == 0) G_{alpha} becomes
        G_a(n):= ( alpha/(2*sqrtD) )*( (eps)^n ) - ( alpha_/(2*sqrtD) )*( (eps_)^n ).
    That is, G_a(n) is a linear binary recurrence sequence, with G_a(n+2) = A*G_a(n+1) - B*G_a(n), where 
        A:= eps + eps_ 
        B:= (eps)*(eps_)
    In this situation, elementary congruence arguments are applied to solve
        G_a(n):= ( alpha/(2*sqrtD) )*( (eps)^n ) - ( alpha_/(2*sqrtD) )*( (eps_)^n ):= +/- prod_{ i \in IU} (p_i)^(u_i):= +/- u,
    hence to compute all solutions to x + y = z^2 in the symmetric case
    
EXAMPLE:
    > a:= sqrtD;
    > S:= [2,3,5,7];
    > D:= 3;
    > SymmetricCaseXYZ2(a,S,D);
    [
        [ 12, -3, 3 ],
        [ 147, -3, 12 ]
    ]
    
*/


function SymmetricCaseXYZ2(A,S,D)
    Z:= IntegerRing();
    Sort(~S);   // orders primes by size p_1 < ... < p_s
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    R:= RingOfIntegers(K);      // generates ring of integers of K
    
    a:= A[3];
    a_:= Conjugate(a);
    eps:= FundamentalUnitXYZ2(D);
    eps_:= Conjugate(eps);
    
    if ideal<R|a> ne ideal<R|a_> then   // verification that indeed (a)O_K == (a_)O_K: if false, there is an problem in SymmetricCaseXYZ2.m
        print "There is a problem in SymmetricCaseXYZ2: the ideals (a)O_K and (a_)O_K are not the same";
    end if;
    
    IU:= IUFactorsXYZ2([],a,S,D);       // generates the set of primes [p_1,...,p_n] of S such that G_a:= G_{alpha} mod p_i == 0, where G_a:= G_{alpha}:= prod_{i in IU}(p_i)^(u_i); in this case, II_:= []
    
    b:= Z! (eps+eps_);  // equivalently, 2*Re(eps), an integer
    c:= Z! -(eps*eps_);         // equivalently, Norm(eps):= +/- 1
    u_0:= Z! ( (a/(2*sqrtD))*(eps)^(0) - (a_/(2*sqrtD))*(eps_)^(0) );   // generates G_a(0), an integer; this is the 0th term of the binary recurrence sequence G_a(n)
    u_1:= Z! ( (a/(2*sqrtD))*(eps)^(1) - (a_/(2*sqrtD))*(eps_)^(1) );   // generates G_a(1), an integer; this is the 1st term of the binary recurrence sequence G_a(n)
    
    l:= [];     // stores [[p_1,m_1],...,[p_n,m_n]], where p_i is a prime of IU and m_i is an integer such that G_a(n) != 0 mod (p_i)^( (m_i)+1 ) for all integers, n, but 
                    // G_a(n):= 0 mod (p_i)^(m_i) for at least one integer n; hence ord_{p_i}( G_a(n) ) <= m_i
    for p in IU do      // runs through all primes p in IU; ie. these are the only primes of S for which G_a(n):= 0 mod p (so p|G_a(n))
        
        m:= 1;  // sets initial exponent, m:= 1
        bound:= 3;      // sets initial bound on m to 3
        r:= [];         // for each m, stores the number of times G_a(n) is divisible by p^m
        while (m le bound) do   // note that after m:= 5, computation time takes noticable longer since G_a(n) increases rapidly and the periods are harder to compute

            I:= SymmetricCaseZerosXYZ2(b,c,u_0,u_1,p,m);        // determines the indices i of G_a(i) such that G_a(i):= 0 mod (p^m)
            
            if IsEmpty(I) eq true then  // if I == [], then G_a(n) != 0 mod (p^m) for all n; otherwise if G_a(n):= 0 mod (p^m), then I would contain 'n'
                Append(~l, [p,m-1]);    // if G_a(n) != 0 mod (p^m) for all n in the period mod (p^m), hence for all n, then ord_p(G_a(n)) <= m-1
                break;  // exits the while loop
            
            elif u_0 eq 0 then  // if G_a(0) == 0
                break;  // exits the while loop; if G_a(0) == 0, then I will contain at least the index 0 for every m

            elif m eq bound then        // if, for all m in [1,...,bound], G_a(n) == 0 mod (p^m) for at least one integer, n
                Append(~r, #I);         // stores the number of elements in the period mod (p^m) for which G_a(n) == 0 mod (p^m); n is an index in the period mod (p^m)
                rep:= &+[1: x in r | x eq #I];  // counts the number of times r has contained the last previously added value, #I

                if rep ge 3 then
                    break;      // exits the while loop if the previously added value, #I, has been repeated at least 3 times
                else
                    bound:= bound + 1;  // increases the bound until I == [] (ie. G_a(n) != 0 mod (p^m) for all n), or rep >= 3 (ie. G_a(n) == 0 mod (p^m) for at least one n, possibly indefinitely)
                    m:= m + 1;
                end if;
            else
                Append(~r, #I);         // stores the number of elements in the period mod (p^m) for which G_a(n) == 0 mod (p^m); n is an index in the period mod (p^m)
                m:= m + 1;
            end if;
        end while;
        
        if (p in [ell[1] : ell in l] eq false) then     // if l does not contain the prime p
            I:= SymmetricCaseZerosXYZ2(b,c,u_0,u_1,p,1);        // determines the indices i of G_a(i) such that G_a(i):= 0 mod p; Nb. I != [] as otherwise l would contain p
            
            if (#I eq 1) and (u_0 eq 0) then    // if the only element that is divisible by p is the first, u_0
                j:= BinaryRecurrenceSequencePeriod(b,c,u_0,u_1,p);      // computes the period of G_a mod p, j; 
                                                                            // j corresponds to the first element of the period mod p such that G_a(j) != 0 and G_a(j) == 0 mod p 
            else
                for k in I do
                    if (BinaryRecurrenceSequence(b,c,u_0,u_1,k,0) ne 0) then
                        j:= k;  // computes the index j, the first element of the period mod p such that G_a(j) != 0 and G_a(j) == 0 mod p 
                        break k;
                    end if;
                end for;
            end if;
            
            m:= 1;
            while Abs(BinaryRecurrenceSequence(b,c,u_0,u_1,j,0)) le &*( [q^(Valuation(BinaryRecurrenceSequence(b,c,u_0,u_1,j,0),q)) : q in IU] ) do
            // if Abs(G_a(j)) <= prod( (p_i)^(h_{p_i}) ) over all primes, p_i, of IU, where h_{p_i} = ord_{p_i}( G_a(j) )
                m:= m + 1;      // increases the exponent on p until the index j is reached such that G_a(j) > prod( (p_i)^(h_{p_i}) ) over all primes, p_i, of IU, where h_{p_i} = ord_{p_i}( G_a(j) )
                I:= SymmetricCaseZerosXYZ2(b,c,u_0,u_1,p,m);        // determines the indices i of G_a(i) such that G_a(i):= 0 mod (p^m)
                
                if (#I eq 1) and (u_0 eq 0) then    // if the only element that is divisible by p^m is the first, u_0
                    j:= BinaryRecurrenceSequencePeriod(b,c,u_0,u_1,p^m);        // computes the period of G_a mod (p^m), j; 
                                                                                    // j corresponds to the first element of the period mod (p^m) such that G_a(j) != 0 and G_a(j) == 0 mod (p^m) 
                elif IsEmpty(I) eq true then    // if I == [], then G_a(n) != 0 mod (p^m) for all n; otherwise if G_a(n):= 0 mod (p^m), then I would contain 'n'
                    break;      // exits the while loop       
            
                else
                    for k in I do
                        if (BinaryRecurrenceSequence(b,c,u_0,u_1,k,0) ne 0) then
                            j:= k;  // computes the index j, the first element of the period mod p such that G_a(j) != 0 and G_a(j) == 0 mod p 
                            break k;
                        end if;
                    end for;
                end if;
           
            end while;
            Append(~l, [p,m-1]);    // if G_a(n) != 0 mod (p^m) for all n in the period mod (p^m), hence for all n, then ord_p(G_a(n)) <= m-1
        end if;
    end for;

    i:= 0;
    if BinaryRecurrenceSequence(b,c,u_0,u_1,i,0) eq 0 then
        i:= i + 1;      // computes the ith term of G_a(n) such that G_a(i) != 0
    end if;
    
    soln:= [];  // stores S-unit solutions x + y = z^2
    LargestG_alpha:= NonEmptyProduct([ell[1]^ell[2] : ell in l]);       // computes the largest G_{alpha} possible, based on the upper bounds on each exponent m_i on p_i in l
    G_alpha:= Z! BinaryRecurrenceSequence(b,c,u_0,u_1,i,0);
    
    while Abs(G_alpha) le LargestG_alpha do
        rem, factors:= SFactors( G_alpha, IU );    
        if (rem eq 1) or (rem eq -1) then       // if G_a:= G_alpha is a nonzero integer <= LargestG_alpha, divisible only by primes of IU
    
            u:= G_alpha;        // u:= G_a(i)
            z:= Z! ( (a/(2))*(eps)^(i) + (a_/2)*(eps_)^(i) );
            x:= Z! ((sqrtD*u)^2);
            y:= z^2 - x;
            if (x eq 0) or (y eq 0) or (z eq 0) then
                i:= i + 1;      // omits G_alpha if any one of x,y,z is identically 0
                G_alpha:= Z! BinaryRecurrenceSequence(b,c,u_0,u_1,i,0);         // updates G_alpha
            else
                if (x ge y) and ([x,y,z] in soln eq false) then
                    Append(~soln,[x,y,z]);
                elif (x lt y) and ([y,x,z] in soln eq false) then
                    Append(~soln,[y,x,z]);
                end if;
                
                i:= i + 1;
                G_alpha:= Z! BinaryRecurrenceSequence(b,c,u_0,u_1,i,0);         // updates G_alpha
            end if;
        else
            i:= i + 1;
            G_alpha:= Z! BinaryRecurrenceSequence(b,c,u_0,u_1,i,0);         // updates G_alpha
        end if;
    end while;
    return soln;
end function;