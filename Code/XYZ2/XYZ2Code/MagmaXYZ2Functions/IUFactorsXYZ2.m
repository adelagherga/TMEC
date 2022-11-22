// IUFactorsXYZ2.m

/*
INPUT:
    II_:= [I,I_] (as output by IIPrimeXYZ2.m), where
        I:= { i | 1 <= i <= t, a_i > b_i }
            I:= [[p_1, h_1, pi_1, (pi_)_1, P_1],...,[p_n, h_n, pi_n, (pi_)_n, P_n]],
                p_i:= prime of S which splits in K
                h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
                pi_i:= generator of the ideal P^(h_i) = (pi)R, chosen so that Abs(pi) > Abs(pi_)
                (pi_)_i:= conjugate of pi, chosen so that Abs(pi) > Abs(pi_)
                P_i:= (unique) choice of prime ideal in K lying above p
        I_:= { i | 1 <= i <= t, a_i < b_i }
            I_:= [[p_1, h_1, pi_1, (pi_)_1, P_1],...,[p_n, h_n, pi_n, (pi_)_n, P_n]],
                p_i:= prime of S which splits in K
                h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
                pi_i:= generator of the ideal P^(h_i) = (pi)R, chosen so that Abs(pi) < Abs(pi_)
                (pi_)_i:= conjugate of pi, chosen so that Abs(pi) < Abs(pi_)
                P_i:= (unique) choice of prime ideal in K lying above p
    
    a:= element of K generating the principal ideal (a)R, where
        a(R):= (alpha) = (2)^(b_0)*[prod_{i in I}(P_i)^(d_i)]*[prod_{i in I_}(P_i')^(d_i)]*prod_{i = t + 1 ... s}P_i^(a_i)
    
    S:= [p_1,...,p_s], p_i primes in S
    
    D:= squarefree part of x, where
        x:= Du^2 and x + y = z^2

OUTPUT:
    IU:= [p_1,...,p_n], p_i primes in S such that
        G_{alpha} mod p_i == 0, where
            G_{alpha}:= prod_{i in IU}(p_i)^(u_i)
                
COMMENTS:
    Given a solution (x,y,z) to x + y = z^2, if x = D*u^2, define u:= G_{alpha}, where
        G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_])
            Prod[I]:= prod_{i in I} (pi_i)^{m_i}
            Prod_[I_]:= prod_{i in I_} ((pi_)_i)^{m_i}
            Prod_[I]:= prod_{i in I} ((pi_)_i)^{m_i}
            Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}, where
                m_i:= an integer > 0
    IU is the set of primes p_i of S such that 
        G_{alpha} = prod_{i in IU}(p_i)^(u_i)
    disjunct from I, I_
        
    This algorithm computes IU by checking G_{alpha} mod p_i for all p_i in S

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > D:= 10;
    > S:= [2,3,5,7,11,13];
    > K<sqrtD>:= QuadraticField(D);
    > R:= RingOfIntegers(K);
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
    > C:= IIPrimeXYZ2(SplitPrimes,D); 
    > II_:= C[4];
    > II_;
    [
        [
            <3, 2, -sqrtD - 1, sqrtD - 1, Prime Ideal of R
            Two element generators:
                3
                $.2 + 4>
        ],
        [
            <13, 2, -58973*sqrtD + 186489, 58973*sqrtD + 186489, Prime Ideal of R
            Two element generators:
                13
                $.2 + 7>
        ]
    ]
    > A0:= AlphasXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > A:= A0[2];
    > A;
    <<Prime Ideal of R
    Two element generators:
        3
        $.2 + 4, [
        <3, 2, -sqrtD - 1, sqrtD - 1, Prime Ideal of R
        Two element generators:
            3
            $.2 + 4, 1>
    ]>, <Ideal of R
    Two element generators:
        13
        $.2 + 6, [
        <13, 2, -58973*sqrtD + 186489, 58973*sqrtD + 186489, Ideal of R
        Two element generators:
            13
            $.2 + 6, 1>
    ]>, 2*sqrtD - 1, 1>
    > a:= A[3];
    > a;
    2*sqrtD - 1
    > IUFactorsXYZ2(II_,a,S,D);
    [ 2, 5, 7, 11 ]
    
*/


function IUFactorsXYZ2(II_,a,S,D)
    Sort(~S);   // orders primes by size p_1 < ... < p_s
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    
    a_:= Conjugate(a);
    eps:= FundamentalUnitXYZ2(D);
    eps_:= Conjugate(eps);
    
    IU:= [];    // stores primes p of S which divide G_{alpha} for some (n,m_i)
    
    if IsEmpty(II_) eq true then
        Prod1:= [1];
        Prod2:= [1];
        Prod3:= [1];
        Prod4:= [1];
    else
        if II_[1] eq [] then
            Prod1:= [1];
            Prod3:= [1];
        else
            Prod1:= IdealExponentsXYZ2(II_[1], <"IUFactors">);         // computes possible Prod[I]:= prod_{i in I} (pi_i)^{m_i} values
            Prod3:= IdealExponentsXYZ2(II_[1], <"IUFactors_">);        // computes possible Prod_[I]:= prod_{i in I} ((pi_)_i)^{m_i} values
        end if;
        
        if II_[2] eq [] then
            Prod2:= [1];
            Prod4:= [1];
        else
            Prod2:= IdealExponentsXYZ2(II_[2], <"IUFactors_">);        // computes possible Prod_[I_]:= prod_{i in I_} ((pi_)_i)^{m_i} values
            Prod4:= IdealExponentsXYZ2(II_[2], <"IUFactors">);         // computes possible Prod[I_]:= prod_{i in I_} (pi_i)^{m_i} values
        end if;
    end if;
        
    P1:= [];    
    for p1 in Prod1 do
        for p2 in Prod2 do
            Append(~P1, p1*p2);         // computes possible (Prod[I])*(Prod_[I_]) values
        end for;
    end for;

    P2:= [];    
    for p3 in Prod3 do
        for p4 in Prod4 do
            Append(~P2, p3*p4);         // computes possible (Prod_[I])*(Prod[I_]) values
        end for;
    end for;
    
    P:= [];
    for i in [1..#P1] do
        Append(~P,[P1[i],P2[i]]);       // comptutes possible tuples [(Prod[I])*(Prod_[I_]), (Prod_[I])*(Prod[I_])]
    end for;
    
    for p in S do
        n:= -50;        // exponent on the fundamental unit, eps
        m:= 1;          // exponent on the splitting primes; nb. the case m = 0 is precisely the symmetric case, hence may be omitted
        while n lt 50 do
            while m le #P do
                G_alpha:= (a/(2*sqrtD))*((eps)^(n))*(P[m][1]) - (a_/(2*sqrtD))*((eps_)^(n))*(P[m][2]);
                
                if (G_alpha in IntegerRing()) and (IntegerRing()! G_alpha ne 0) and (IntegerRing()! G_alpha mod p eq 0) then      // verifies if p|G_{alpha}
                    Append(~IU,p);
                    n:= 102;
                    m:= 2*#P;   // exit loop for this value of p
                else
                    m:= m + 1;
                end if;
            end while;
            n:= n + 1;
            m:= 1;
        end while;
    end for;
    
    if IsEmpty(II_) eq false then
        for i in II_[1] do
            if i[1] in IU then
                Exclude(~IU,i[1]);  // IU, I are disjoint, where I is predetermined, hence removes the primes from IU which also appear in I
            end if;
        end for; 
        
        for i in II_[2] do
            if i[1] in IU then
                Exclude(~IU,i[1]);  // IU, I_ are disjoint, where I_ is predetermined, hence removes the primes from IU which also appear in I_
            end if;
        end for; 
    end if; 
    
    return IU;
end function;