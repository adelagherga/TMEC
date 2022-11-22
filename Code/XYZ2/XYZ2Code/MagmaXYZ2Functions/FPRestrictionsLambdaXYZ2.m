// FPRestrictionsLambdaXYZ2.m

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
    U0:= [ [u_i, p_i] : p_i in IU], where
        u_i:= ord_{p_i}(u) such that 
            u:= G_{alpha}, or, said differently +/- u = G_{alpha} = prod_{i in IU}(p_i)^(u_i)
            ie. +/- u = G_{alpha} = prod_{i in IU}(p_i)^(u_i), under x = D*u^2 
        p_i:= prime in IU
    M0:= [ [c_i, p_i] : p_i in I], where
        c_i:= the exponents m_i on Prod[I], Prod_[I] in G_{alpha} giving G_{alpha} = u, under x = D*u^2, where
            G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_]) and
                Prod[I]:= prod_{i in I} (pi_i)^{m_i}
                Prod_[I_]:= prod_{i in I_} ((pi_)_i)^{m_i}
                Prod_[I]:= prod_{i in I} ((pi_)_i)^{m_i}
                Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}, where
                    m_i:= an integer > 0
    M0_:= [ [c_i, p_i] : p_i in I_], where
        c_i:= the exponents m_i on Prod[I_], Prod_[I_] in G_{alpha} giving G_{alpha} = u, under x = D*u^2, where
            G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_]) and
                Prod[I]:= prod_{i in I} (pi_i)^{m_i}
                Prod_[I_]:= prod_{i in I_} ((pi_)_i)^{m_i}
                Prod_[I]:= prod_{i in I} ((pi_)_i)^{m_i}
                Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}, where
                    m_i:= an integer > 0
    absn:= |n|, where n is an integer representing the exponent on eps, eps_ in G_{alpha} giving G_{alpha} = u, under x = D*u^2, where
            G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_]) and 
                eps:= the fundamental unit of the rind of integers of K = Q(sqrt(D))
                eps_:= the conjugate of the fundamental unit, eps
    nM:= [[p_i, n_i] | i in I], where
        p_i:= a prime of I
        n_i:= the exponent of p_i appearing in alpha^{hstar}  
    nM_:= [[p_i, n_i] | i in I_], where
        p_i:= a prime of I_
        n_i:= the exponent of p_i appearing in alpha^{hstar}
    nEps:= the exponent n_0 of the fundamental unit, eps, appearing in alpha^{hstar}
    hstar:= LCM(2,h_1,...,h_s), where
        h_i:= smallest exponent such that P^(h_i) is principal, for P an ideal in K, lying above p in S 
    terms:= the terms {eps/eps_, pi_j/(pi_)_j} appearing Lambda*, excluding p0, where
        Lambda*:= nstar*log_p(eps/eps_) + sum_{j in I} cstar_j*log_p( pi_j/(pi_)_j ) - sum_{j in I_} cstar_j*log_p( pi_j/(pi_)_j )
    p0:= the term of Lambda* giving minimal ord_p(log_p(p0))
    w:= the set of all lattice coefficient vectors x_i of the lattice L such that y_i:= Transpose(W*Transpose(x_i)) has norm <= C,
            ie. Q(y_i) <= C, where Q is the binary quadratic form Q associated to the lattice L
            Equivalently, (y_i)*A*Transpose(y_i) = [z_i], where z_i <= C
        y_i:= row vector whose entries give nstar, cstar_j of Lambda*, where
            Lambda*:= nstar*log_p(eps/eps_) + sum_{j in I} cstar_j*log_p( pi_j/(pi_)_j ) - sum_{j in I_} cstar_j*log_p( pi_j/(pi_)_j )
    W:= associated lattice matrix such that y_i:= Transpose(M*Transpose(x_i)) is a vector if the lattice with norm <= C, where
            x_i:= a lattice coefficient vector
    S:= [p_1,...,p_s], p_i primes in S
    D:= squarefree part of x, where
        x:= Du^2 and x + y = z^2

OUTPUT:
    soln:= [[x_1,y_1,z_1],...,[x_n,y_n,z_n]], all S-unit equations x_i + y_i = (z_i)^2 such that 
        u_i:= ord_p(u) is in the range [m + m0 - lambda_i - Ord_p(hstar), U0[p]], where
            p:= a given prime in U0
            u:= G_{alpha}, or, said differently +/- u = G_{alpha} = prod_{i in IU}(p_i)^(u_i)
                ie. +/- u = G_{alpha} = prod_{i in IU}(p_i)^(u_i), under x = D*u^2 
        
COMMENTS:
    Computes all S-unit equations x + y = z^2 involving ord_p(u) in the range [m + m0 - lambda_i - Ord_p(hstar), U0[p]], for given prime p in U0:
        For a given prime p, to find all solutions (x,y,z) coprime with ord_p(u) <= U0[p]:
            Choose m < U0[p] - m0 + lambda_i + ord_{p_i}(hstar) and consider the lattice Lambda_{m}* = LambdaLatticeXYZ2(A,p,m,D);  
            If a solution (x,y,z) exists with ord_p(u) in [m + m0 - lambda_i - Ord_p(hstar), U0[p]], then the vector (x_1*,...,x_s*,x_0*)^T is in the lattice, where x_i*:= n* or c_j*, and
                n:= (n* - n_0)/hstar, c_j:= (c_j* - n_j)/hstar such that 
                    n:= an integer representing the exponent on eps
                    c_j:= the exponents m_j on Prod[I], Prod_[I] in G_{alpha} giving G_{alpha} = u, under x = D*u^2, where
                G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_]) and
                    Prod[I]:= prod_{i in I} (pi_i)^{m_i}
                    Prod_[I_]:= prod_{i in I_} ((pi_)_i)^{m_i}
                    Prod_[I]:= prod_{i in I} ((pi_)_i)^{m_i}
                    Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}, where
                        m_i:= an integer > 0
            The length of this vector is bounded by ((absn*)^2 + [(c_j*)^2 | j in I] + [(c_j*)^2 | j in I_]) 
    This algorithm uses the Fincke-Pohst algorithm to find all vectors of length at most ((absn*)^2 + [(c_j*)^2 | j in I] + [(c_j*)^2 | j in I_]) and sorts them to give (x,y,z) for the specified p in U0

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,5,7,11];
    > D:= 15;
    > K<sqrtD>:= QuadraticField(D);
    > R:= RingOfIntegers(K);
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
    > C:= IIPrimeXYZ2(SplitPrimes,D); 
    > II_:= C[4];
    > A0:= AlphasXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > FA:= MaximalC12BoundXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > U0,M0,M0_,absn:= SmallestLatticeBoundXYZ2(II_,FA,S);
    > U0;
    [
        [ 44, 2 ],
        [ 24, 3 ],
        [ 17, 5 ]
    ]
    > M0;
    [
        [ 32, 7 ]
    ]
    > M0_;
    [
        [ 25, 11 ]
    ]
    > absn;
    442
    > A:= A0[1];
    > A;
    <<Principal Ideal of R
    Generator:
        1, [
        <7, 2, sqrtD + 8, -sqrtD + 8, Prime Ideal of R
        Two element generators:
            7
            $.2 + 1, 0>
    ]>, <Principal Ideal of R
    Generator:
        1, [
        <11, 1, -7314*sqrtD + 28327, 7314*sqrtD + 28327, Principal Ideal of R
        Generator:
            -$.2 + 2, 0>
    ]>, 1, 1>
    > p:= 5;
    > m:= 1;
    > Bound:= 4000;
    > B,m0:= LambdaLatticeXYZ2(A,p,m,D);                
    > Sbound:= 10000000;                           
    > W,w,maxReached:= FinckePohst(B,Bound,Sbound);
    > #w;
    105921
    > maxReached;
    false
    > soln:= FPRestrictionsLambdaXYZ2(A,U0,M0,M0_,absn,nM,nM_,nEps,hstar,terms,p0,w,W,S,D);
    > soln;
    [
        [ 75937500, -7891499, 8249 ],
        [ 138240, 14641, 391 ],
        [ 135, 121, 16 ],
        [ 290521, 2160, 541 ],
        [ 540, -539, 1 ],
        [ 3375, 2401, 76 ]
    ]
    
*/


function FPRestrictionsLambdaXYZ2(A,U0,M0,M0_,absn,nM,nM_,nEps,hstar,terms,p0,w,W,S,D);
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    Z:= IntegerRing();
    
    eps:= FundamentalUnitXYZ2(D);
    eps_:= Conjugate(eps);
    q:= A[1];
    q_:= A[2];
    a:= A[3];
    a_:= Conjugate(a);
    
    if #q eq 0 then
        I:= [];
    else
        I:= q[2];
    end if;
    if #q_ eq 0 then
        I_:= [];
    else
        I_:= q_[2];
    end if;
    
    // LargestG_alpha:= NonEmptyProduct([u[2]^u[1] : u in U0]);    // computes the largest G_{alpha} possible, based on the upper bounds on each u_i of U0 
    U:= [u[2] : u in U0];       // generates the set of primes p_i of U0 (hence of IU)
    Ms:= [m[2] : m in M0];      // generates the set of primes p_i of M0 (hence of I)
    M_s:= [m[2] : m in M0_];    // generates the set of primes p_i of M0_ (hence of I_)
    
    Append(~terms, p0);         // computes the terms {eps/eps_, pi_j/(pi_)_j} appearing Lambda*, including p0, where
                                    // Lambda*:= nstar*log_p(eps/eps_) + sum_{j in I} cstar_j*log_p( pi_j/(pi_)_j ) - sum_{j in I_} cstar_j*log_p( pi_j/(pi_)_j )
                                    // NB. the order of terms corresponds to the ordering of the columns of B, and hence corresponds to the ordering of the indices 
                                        // of the short vectors m_i of w. (In particular, p0 corresponds to the last column of B and therefore to the last index of m_i)
    
    epsIndex:= Index(terms, eps/eps_);  // computes the index of eps/eps_ in terms
    IIndices:= [];      // stores [ < pi_j, (pi_)_j, Index of (pi_j/(pi_)_j), M0[j] , nM[j] >: j in I ]; Index of (pi_j/(pi_)_j) is the index of pi_j/(pi_)_j as it appear in terms; 
                            // M0[j] is the upper bound on the exponent of this term, as in M0; nM[j] is the exponent n_j corresponding to c_j
    I_Indices:= [];     // stores [ < pi_j, (pi_)_j, Index of (pi_j/(pi_)_j), M0_[j] , nM_[j] >: j in I_ ]; Index of (pi_j/(pi_)_j) is the index of pi_j/(pi_)_j as it appear in terms;
                            // M0_[j] is the upper bound on the exponent of this term, as in M0_; nM_[j] is the exponent n_j corresponding to c_j
    for j in I do
        Append(~IIndices, < j[3], j[4], Index(terms, j[3]/j[4]), M0[ Index(Ms,j[1]) ][1], nM[Index(Ms,j[1]) ][2] > );
    end for;
    for j_ in I_ do
        Append(~I_Indices, < j_[3], j_[4], Index(terms, j_[3]/j_[4]), M0_[ Index(M_s,j_[1]) ][1], nM_[Index(M_s,j_[1]) ][2] > );
    end for;
    
    VectorSort1:= [];   // stores only the vectors m_i = w[i] found in Fincke-Pohst such that each vector element m_i[j] correspond to elements of the lattice Lambda
                            // ie. the vector element m_i[j] of m_i corresponding to n*, c_j* of Lambda* should necessarily also correspond to n, c_j of Lambda; 
    for m in w do       // checks each vector w[j] found in Fincke-Pohst, weeding them out accordingly
        m_i:= Transpose(W*Transpose(m));        // computes the short vector m_i of the lattice Lambda* 
        if ( Denominator( (m_i[1,epsIndex] - nEps)/hstar ) eq 1) and (&and[ Denominator( (m_i[1,i[3]] - i[5])/hstar ) eq 1 : i in IIndices]) and (&and[ Denominator( (m_i[1,i[3]] - i[5])/hstar ) eq 1 : i in I_Indices]) then
            Append(~VectorSort1, m_i);
        end if;
        if ( Denominator( (-m_i[1,epsIndex] - nEps)/hstar ) eq 1) and (&and[ Denominator( (-m_i[1,i[3]] - i[5])/hstar ) eq 1 : i in IIndices]) and (&and[ Denominator( (-m_i[1,i[3]] - i[5])/hstar ) eq 1 : i in I_Indices]) then
            Append(~VectorSort1, -m_i);
        end if;
    end for;
    
    VectorSort2:= [];   // stores only the vectors m_i found in VectorSort1 such that each vector element m_i[j] has the correct sign
                            // ie. the vector elements m_i[j] of m_i corresponding to c_j should necessarily all be positive integers; 
                            // the vector element m_i[j] of m_i corresponding to n should be an integer (this is automatically true)
    for m_i in VectorSort1 do   // checks each vector w[j] found in VectorSort1, weeding them out accordingly
        if (&and[ ( (m_i[1,i[3]] - i[5])/hstar ) ge 0 : i in IIndices]) and (&and[ ( (m_i[1,i[3]] - i[5])/hstar ) ge 0 : i in I_Indices]) then
            Append(~VectorSort2, m_i);
        end if;
    end for;
    
    VectorSort3:= [];   // stores only the vectors m_i found in VectorSort2 such that each vector element m_i[j] is the correct size
                            // ie. the vector elements m_i[j] of m_i corresponding to c_j should necessarily have m_i[j] <= M0[j] (or resp. m_i[j] <= M0_[j])
                            // the vector element m_i[j] of m_i corresponding to n should necessarily have |m_i[j]| <= absn
    for m_i in VectorSort2 do   // checks each vector w[j] found in VectorSort2, weeding them out accordingly
        if (Abs( (m_i[1,epsIndex] - nEps)/hstar ) le absn) and (&and[ ( (m_i[1,i[3]] - i[5])/hstar ) le i[4] : i in IIndices]) and (&and[ ( (m_i[1,i[3]] - i[5])/hstar ) le i[4] : i in I_Indices]) then
            Append(~VectorSort3, m_i);
        end if;
    end for;
   
    G_alphas:= [];      // stores [ G_alpha, G_alpha1, G_alpha2 ] where G_alpha:= G_alpha1 - G_alpha2 and G_alpha is a non-zero integer divisible only by primes of S 
    Ones:= [];  // stores [ G_alpha1, G_alpha2 ] where G_alpha:= G_alpha1 - G_alpha2 whenever G_alpha = 1 
    
    for m_i in VectorSort3 do   // checks each vector w[j] found in VectorSort3, weeding them out accordingly
        G_alpha1:= ( a/(2*sqrtD) )*( (eps)^(Z!((m_i[1,epsIndex] - nEps)/hstar)) )*( NonEmptyProduct([ i[1]^(Z!((m_i[1,i[3]] - i[5])/hstar)) : i in IIndices ]) )*( NonEmptyProduct([ i[2]^(Z!((m_i[1,i[3]] - i[5])/hstar)) : i in I_Indices ]) ); 
        G_alpha2:= ( a_/(2*sqrtD) )*( (eps_)^(Z!((m_i[1,epsIndex] - nEps)/hstar)) )*( NonEmptyProduct([ i[2]^(Z!((m_i[1,i[3]] - i[5])/hstar)) : i in IIndices ]) )*( NonEmptyProduct([ i[1]^(Z!((m_i[1,i[3]] - i[5])/hstar)) : i in I_Indices ]) );
        
        G_alpha:= G_alpha1 - G_alpha2;  // computes G_alpha
        if (G_alpha in Z) and (G_alpha ne 0) and &or([Z!(G_alpha) mod p eq 0 : p in S]) then    // if G_alpha is a nonzero integer <= LargestG_alpha, divisible by some prime of S
            rem, factors:= SFactors(Z!(G_alpha), S);
            if (rem eq 1) or (rem eq -1) then   // if G_alpha is a nonzero integer <= LargestG_alpha, divisible only by primes of S
                Append(~G_alphas,< G_alpha, G_alpha1, G_alpha2 >);
            end if;
        elif (G_alpha eq 1) or (G_alpha eq -1) then     // if G_alpha = 1 or G_alpha = - 1
            Append(~Ones, < G_alpha1, G_alpha2 >);
        end if;
    end for;
    
    soln:= [];  // stores [x,y,z] where x + y = z^2
    for g in G_alphas do
        u:= Z! g[1];    // computes u, where u:= +/-G_alpha
        z:= sqrtD*g[2] + sqrtD*g[3];    // computes z, where z:= +/-H_alpha, where H_alpha:= sqrtD*G_alpha1 + sqrtD*G_alpha2
        if z in Z then
            x:= Z! ((sqrtD*u)^2);
            y:= Z!(z^2) - x;
            z:= Abs(Z! z);
            if (x ne 0) and (y ne 0) and (z ne 0) then
                if (x ge y) and ([x,y,z] in soln eq false) then
                    Append(~soln,[x,y,z]);
                elif (x lt y) and ([y,x,z] in soln eq false) then
                    Append(~soln,[y,x,z]);
                end if;
            end if;
        end if;
    end for;
    
    for g in Ones do
        u:= 1;  // computes u, where u:= +/-G_alpha = +/-1
        z:= sqrtD*g[1] + sqrtD*g[2];    // computes z, where z:= +/-H_alpha, where H_alpha:= sqrtD*G_alpha1 + sqrtD*G_alpha2
        if z in Z then
            x:= Z! ((sqrtD*u)^2);
            y:= Z!(z^2) - x;
            z:= Abs(Z! z);
            if (x ne 0) and (y ne 0) and (z ne 0) then
                if (x ge y) and ([x,y,z] in soln eq false) then
                    Append(~soln,[x,y,z]);
                elif (x lt y) and ([y,x,z] in soln eq false) then
                    Append(~soln,[y,x,z]);
                end if;
            end if;
        end if;
    end for;
    
    return soln;
end function;