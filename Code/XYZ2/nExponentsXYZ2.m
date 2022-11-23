// nExponentsXYZ2.m

/*
INPUT:
    b0:= an integer in {0,1} indicating whether 2 splits in K
    
    A:= < q, (q_), alpha, tplus1 >, (as output by AlphasXYZ2.m) where
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
    n:= [[p_{t+1}, n_{t+1}, v_{t+1}],...,[p_s, n_s, v_s]], where
        p_i:= a nonsplit prime of S
        n_i:= the exponent of p_i appearing in alpha^{hstar}
        v_i:= (1/2)*hstar*ord_{p_i}(4*D)
    nM:= [[p_i, n_i] | i in I], where
        p_i:= a prime of I
        n_i:= the exponent of p_i appearing in alpha^{hstar}  
    nM_:= [[p_i, n_i] | i in I_], where
        p_i:= a prime of I_
        n_i:= the exponent of p_i appearing in alpha^{hstar}
    nEps:= the exponent n_0 of the fundamental unit, eps, appearing in alpha^{hstar}
    N:= Max(|n_0|,...,|n_t|, |n_{t+1} - v_{t+1}|,...,|n_{s} - v_{s}|), a rational integer, where
        n_0:= the exponent appearing in alpha^{hstar} of +/-eps
        n_i:= the exponent appearing in alpha^{hstar} of pi_i of I, respectively (pi_)_i of I_, respectively p_i of tplus1
        v_i:= (1/2)*hstar*ord_{p_i}(4*D) for p_i in tplus1
        
COMMENTS:
    Given alpha, where
        (alpha):= (2)^(b_0)*[prod_{i in I}(P_i)^(d_i)]*[prod_{i in I_}((P_i)_i)^(d_i)]*prod_{i = t + 1 ... s}P_i^(a_i), and
        h^star:= h*:= LCM(2,h_1,...,h_s), where
            h_i:= smallest exponent such that P^(h_i) is principal, for P an ideal in K, lying above p in S 
    This algorithm computes the exponents n_i of alpha^{hstar}, where
        alpha^{hstar} = +/- [(eps)^(n_0)]*[prod_{i in I} (pi_i)^(n_i)]*[prod_{i in I_} ((pi_)_i)^(n_i)]*[prod_{i = t+1,...,s}(p_i)^(n_i)]*[2^(hstar*b0)]
    
    Based on the algorithm found in Chapter 7 of the Reference

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,5,7,11,13];
    > D:= 5;
    > K<sqrtD>:= QuadraticField(D);
    > R:= RingOfIntegers(K);
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
    > C:= IIPrimeXYZ2(SplitPrimes,D); 
    > II_:= C[1];
    > A0:= AlphasXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > A:= A0[2];
    > A;
    <<Principal Ideal of R
    Generator:
        1, [
        <11, 1, 1/2*(-3*sqrtD - 1), 1/2*(3*sqrtD - 1), Principal Prime Ideal of R
        Generator:
            -3*$.2 + 1, 0>
    ]>, [], 2*sqrtD, 10>
    > b0;
    0
    > IU:= IUFactorsXYZ2([I,I_],A[3],S,D); 
    > IU;
    [ 2, 3, 7, 13 ]
    > n, nM, nM_, nEps, N:= nExponentsXYZ2(b0,A,S,D);
    > n;
    [ <2, 2, 2>, <3, 0, 0>, <5, 1, 1>, <7, 0, 0>, <13, 0, 0> ]
    > nM;
    [ <11, 0> ]
    > nM_;
    []
    > nEps;
    0
    > N;
    0

*/


function nExponentsXYZ2(b0,A,S,D)
    Sort(~S);   // orders primes by size p_1 < ... < p_s
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    R:= RingOfIntegers(K);      // generates ring of integers of K
    
    eps:= FundamentalUnitXYZ2(D);
    eps_:= Conjugate(eps);
    q:= A[1];
    q_:= A[2];
    
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
    
    IU:= IUFactorsXYZ2([I,I_],A[3],S,D);        // computes prime factors (belonging to S) of G_{alpha}
    H:= [h[2] : h in I cat I_];         // stores h_i values corresponding to each elemet of I, I_ 
    Append(~H,2);       // appends h_i:= 2, corresponding to the primes which ramify in K
    hstar:= Lcm(H);     // computes h*:= LCM(2,h_1,...,h_s), corresponding to primes which split (h_i such that h_i|ClassNumber(K)), ramify (h_i:= 2), or remain inert (h_i:= 1)
   
    astar:= A[3]^(hstar);       // computes alpha^{hstar} = +/- [(eps)^(n_0)]*[prod_{i in I} (pi_i)^(n_i)]*[prod_{i in I_} ((pi_)_i)^(n_i)]*[prod_{i = t+1,...,s}(p_i)^(n_i)]*[2^(hstar*b0)]
                                // used to verify correct n_i values
   
    n:= [];     // stores [p, n_p, v_p, IndexIU], where p is a non-split prime of S, n_p is the exponent of p appearing in alpha^{hstar}, v_p = (1/2)*hstar*ord_p(4*D)
    nM:= [];    // stores [p, n_p] where p is a prime of I, n_p is the exponent of p appearing in alpha^{hstar}
    nM_:= [];   // stores [p, n_p] where p is a prime of I_, n_p is the exponent of p appearing in alpha^{hstar}
    N:= [];     // stores |n_0|,...,|n_t|, |n_{t+1} - v_{t+1}|,...,|n_{s} - v_{s}|

    rem, tplus1:= SFactors(A[4],S);     // returns tplus1:= [[p_1,count_1],...,[p_s,count_s]], where p_i is a prime in S, count_i:= ord_{p_i}(tplus1product); necessarily, rem == 1
    
    if b0 eq 1 then     // if 2 splits (hence it is not in the set tplus1)
        astar:= astar/(2^(IntegerRing()!(hstar*b0)));   // updates astar verification by dividing astar by 2^(hstar*b0)
    end if;
    
    for i in I do
        ni:= IntegerRing() ! ((hstar*i[6])/i[2]);       // computes n_i:= (hstar*d_i)/h_i; since P^(h_i) = (pi)K, necessarily (P^(d_i))^(hstar) = P^(h_i*(d_i*hstar/h_i)) = (pi)^((d_i*hstar)/h_i)K
        Append(~N, Abs(ni));
        Append(~nM, < i[1], ni >);
        astar:= astar/(i[3]^ni);        // updates astar verification
        Exclude(~tplus1, tplus1[Index(tplus1,[0,i[1]])]);       // removes p (appearing in I) from tplus1; tplus1 includes all primes of S, including split primes
    end for;
    
    for i in I_ do    
        ni:= IntegerRing() ! ((hstar*i[6])/i[2]);       // computes n_i:= (hstar*d_i)/h_i; since P^(h_i) = (pi)K, necessarily (P^(d_i))^(hstar) = P^(h_i*(d_i*hstar/h_i)) = (pi)^((d_i*hstar)/h_i)K
        Append(~N, Abs(ni));
        Append(~nM_, < i[1], ni >);
        astar:= astar/(i[4]^ni);        // updates astar verification
        Exclude(~tplus1,tplus1[Index(tplus1,[0,i[1]])]);        // removes p (appearing in I_) from tplus1; tplus1 includes all primes of S, including split primes
    end for;
    
    for t in tplus1 do  // runs through all nonsplit primes of S corresponding to alpha
        if t[1] ne 0 then       // ord_{p_i}(tplus1product) != 0, so n_i != 0
            vi:= (1/2)*hstar*Valuation(4*D,t[2]);
            if IsInert(t[2],R) eq true then     // if 2 is inert in K, h_i = 1
                Append(~n, < t[2], hstar, vi >);        // computes and stores n_i:= hstar
                Append(~N, Abs(hstar - vi));    // stores |n_i - v_i|
                astar:= astar/(t[2]^(IntegerRing()!hstar));     // updates astar verification
            else
                Append(~n, < t[2], hstar/2, vi >);      // if 2 ramifies in K, h_i = 2, so n_i:= hstar/2
                Append(~N, Abs(hstar/2 - vi));  // stores |n_i - v_i|
                astar:= astar/(t[2]^(IntegerRing()!(hstar/2)));         // updates astar verification
            end if;
        else    // ord_{p_i}(tplus1product) = 0, so n_i:= 0 
            vi:= (1/2)*hstar*Valuation(4*D,t[2]);
            Append(~n, < t[2], 0 , vi >);
            Append(~N, Abs(-vi));       // stores |n_i - v_i| = |-v_i|
        end if;
    end for;
    
    if (astar eq 1) or (astar eq -1) then
        nEps:= 0;       // stores the exponent n_0 of the fundamental unit eps in alpha^{hstar}
    else
        i:= 0;
        while (astar ne 1) and (astar ne -1) do     // computes n_0, the exponent eps appearing in astar
            if (astar/(eps^i) eq 1) or (astar/(eps^i) eq -1) then   // if n_0 > 0
                Append(~N, Abs(i));
                nEps:= i;       // stores the exponent n_0 of the fundamental unit eps in alpha^{hstar}
                astar:= astar/(eps^i);      // updates astar verification
            elif (astar*(eps^i) eq 1) or (astar*(eps^i) eq -1) then         // if n_0 <= 0
                Append(~N, Abs(-i));
                nEps:= -i;      // stores the exponent n_0 of the fundamental unit eps in alpha^{hstar}
                astar:= astar*(eps^i);      // updates astar verification
            else
                i:= i + 1;
            end if;
        end while;
    end if;
   
    return n, nM, nM_, nEps, Max(N);
end function;
