// AlphasXYZ2.m

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
    
    b0:= an integer in {0,1} indicating whether 2 splits in K
        
    SplitPrimes:= primes of S which split in K = Q(Sqrt(D)) 
        SplitPrimes is stored in Magma as SplitPrimes[i]:= [p_i, h_i, P_i], where
            p_i:= prime of S which splits in K
            h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
            P_i:= (unique) choice of prime ideal in K lying above p
        
    NonSplitPrimes:= primes of S which do not split in K = Q(Sqrt(D)) but potentially contribute to alpha
        NonSplitPrimes is stored in Magma as SplitPrimes[i]:= [p_i, P_i], where
            p_i:= prime of S which does not split in K but contributes to the product ideal (X)
            P_i:= prime ideal in K lying above p 

    S:= [p_1,...,p_s], p_i primes in S

    D:= squarefree part of x, where
        x:= Du^2 and x + y = z^2

OUTPUT:
    if SplitPrimes = []:
        ie. if no primes of S split in K = Q(Sqrt(D))
        Alphas:= [< "-","-",alpha_i,tplus1_i >, ..., < "-","-",alpha_n, tplus1_n >], where
            alpha_i:= element of K, generating the ideal (alpha), where
                (alpha):= prod_{i = t + 1 ... s}P_i^(a_i)
            tplus1_i:= the product of primes of S from NonSplitPrimes whose prime ideals above contribute to alpha as 
                prod_{i = t + 1 ... s}P_i^(a_i)
                
    if SplitPrimes != []:
        ie. if there are primes of S which split in K = Q(Sqrt(D))
        Alphas:= [< q_i, (q_)_i, alpha_i, tplus1_i >, ..., < q_n, (q_)_n, alpha_n, tplus1_n >], where
            q_i:= < prod_i, [< p, h_i, pi, pi_, P , d_i > : j in [1..n] ] >, such that
                prod_i:= prod_{i in I}(P_i)^(d_i), the product of prime ideals of I which contribute to alpha
                p:= prime of S which splits in K
                h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
                pi:= generator of the ideal P^(h_i) = (pi)R
                (pi_):= conjugate of pi
                P:= (unique) choice of prime ideal in K lying above p
                n:= the length of Iset
                d_i:= the exponent on the prime ideal P_i in prod_i
            (q_)_i:= < prod_i, [< p, h_i, pi, pi_, P_ , d_i > : j in [1..n] ] >, such that
                prod_i:= prod_{i in I}((P_)_i)^(d_i), the product of the conjugate of the prime ideals of I_ which contribute to alpha
                p:= prime of S which splits in K
                h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
                pi:= generator of the ideal P^(h_i) = (pi)R
                (pi_):= conjugate of pi
                P_:= conjugate of the (unique) choice of prime ideal in K lying above p
                n:= the length of Iset
                d_i:= the exponent on the prime ideal P_i in prod_i
            alpha_i:= element of K, generating the ideal (alpha), where
                (alpha):= (2)^(b_0)*[prod_{i in I}(P_i)^(d_i)]*[prod_{i in I_}((P_i)_i)^(d_i)]*prod_{i = t + 1 ... s}P_i^(a_i)
            tplus1_i:= the product of primes of S from NonSplitPrimes whose prime ideals above contribute to alpha as 
                prod_{i = t + 1 ... s}P_i^(a_i)

COMMENTS:
    Given a solution (x,y,z) to x + y = z^2, if x = D*u^2, define X:= z + u*SqrtD(D) and
        (X):= Prod_{i=1}^s ((P_i)^{a_i})*((P_i')^{b_i}), where
            s:= number of primes of S:= [p_1,...,p_s]
            P_i:= prime ideal of K above the prime p_i of S
            P_i':= the conjugate ideal of P_i
    Alternatively, (X) may be written as 
        (X):= (alpha)*(beta) for some ideals alpha, beta, of K = Q(Sqrt(D)) 
    The ideal alpha is principal, with
        (alpha):= (2)^(b_0)*[prod_{i in I}(P_i)^(d_i)]*[prod_{i in I_}((P_i)_i)^(d_i)]*prod_{i = t + 1 ... s}P_i^(a_i)
            d_i:= an integer such that 0 <= d_i <= h_i - 1 
            
    There are only finitely many possibilities for the generators of alpha; this algorithm computes every possible such generator

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,11,23];
    > D:= 253;   
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
    > II_:= (IIPrimeXYZ2(SplitPrimes,D))[1]; 
    > Alphas:= AlphasXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > Alphas;
    [
        <<Principal Ideal
        Generator:
            1, [
            <3, 1, -sqrtD - 16, sqrtD - 16, Principal Prime Ideal
            Generator:
                -2*$.2 - 15, 0>
        ]>, [], 2, 2>,
        <<Principal Ideal
        Generator:
            1, [
            <3, 1, -sqrtD - 16, sqrtD - 16, Principal Prime Ideal
            Generator:
                -2*$.2 - 15, 0>
        ]>, [], 9*sqrtD + 143, 22>,
        <<Principal Ideal
        Generator:
            1, [
            <3, 1, -sqrtD - 16, sqrtD - 16, Principal Prime Ideal
            Generator:
                -2*$.2 - 15, 0>
        ]>, [], 13*sqrtD + 207, 46>,
        <<Principal Ideal
        Generator:
            1, [
            <3, 1, -sqrtD - 16, sqrtD - 16, Principal Prime Ideal
            Generator:
                -2*$.2 - 15, 0>
        ]>, [], 2*sqrtD, 506>
    ]

*/


function AlphasXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D)
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    R:= RingOfIntegers(K);      // generates ring of integers of K; every element of R is of the form a + Nu*b, where a,b are integers, Nu is defined as below
    
    if (D mod 4 eq 2) or (D mod 4 eq 3) then
        Nu:= sqrtD;
        NumericalNu:= Sqrt(D);  // the numerical approximation of Nu as a real number
    elif (D mod 4 eq 1) then
        Nu:= (1 + sqrtD)/2;
        NumericalNu:= (1 + Sqrt(D))/2;  // the numerical approximation of Nu as a real number
    end if;     // D mod 4 eq 0 is impossible since it would imply that D would be divisible by the square 4
    
    NonSplitProduct:= ExponentsXYZ2(NonSplitPrimes,1);  // computes the possiblilities for prod_{i = t + 1 ... s}P_i^(a_i), the product of non-split primes contributing to alpha 
    Alphas:= [];        // stores values of alpha and their invariants as:
                            // if SplitPrimes = []: Alpha[i]:= < alpha_i,tplus1_i >
                            // if SplitPrimes != []: Alpha[i]:= < q, q_, alpha_i,tplus1_i >
     
    if IsEmpty(SplitPrimes) eq true then        // if no primes of S split in K = Q(Sqrt(D)); in this case (alpha):= prod_{i = t + 1 ... s}P_i^(a_i)
        for PossibleAlpha in NonSplitProduct do
            Pr:= PossibleAlpha[2];      // computes the ideal prod_{i = t + 1 ... s}P_i^(a_i) in K
            t, gen:= IsPrincipal(Pr);   // Pr must be principal; if it is, (alpha):= (Pr)
            if t eq true then   // if Pr is principal, necessarily, (alpha):= (Pr), hence gen = alpha 
                
                a1:= gen[1] + Nu*gen[2]; // writes a1 in the form a + Nu*b
                Numericala1:= gen[1] + NumericalNu*gen[2];      // computes the numerical approximation of a1 as a real number
                if Numericala1 lt 0 then
                    a1:= -1*a1;
                end if;
                    
                a2:= Conjugate(a1);
                Numericala2:= a2[1] + Sqrt(D)*a2[2];
                if Numericala2 lt 0 then
                    a2:= -1*a2;
                end if;
                
                if (Abs(Numericala1) ge Abs(Numericala2)) and (ideal<R|a1> eq ideal<R|a2>) then    
                    alpha:= a1;         // chooses alpha to be the larger value of a1, a2 (the conjugate of a1) if both a1, a2 generate the ideal (alpha)
                elif (Abs(Numericala1) lt Abs(Numericala2)) and (ideal<R|a1> eq ideal<R|a2>) then
                    alpha:= a2;
                else
                    alpha:= a1;         // chooses alpha to be a1 if (a1)R != (a2)R generate different ideals (alpha) as G_(a1) = +/- G_(a2)
                end if;
               
                Append(~Alphas,< "-","-",alpha,PossibleAlpha[1] >);
            end if;
        end for;
       
    else        // if primes of S split in K = Q(Sqrt(D)); in this case (alpha):= (2)^(b_0)*[prod_{i in I}(P_i)^(d_i)]*[prod_{i in I_}((P_i)_i)^(d_i)]*prod_{i = t + 1 ... s}P_i^(a_i)
        I:= II_[1];
        I_:= II_[2];
        
        if IsEmpty(I) eq false then
            Q:= IdealExponentsXYZ2(I,<"Alphas">);         // computes the possiblilities for [prod_{i in I}(P_i)^(d_i)], the product of split primes of I contributing to alpha 
        else
            Q:= [[1]];    // if I is empty, sets [prod_{i in I}(P_i)^(d_i)] = 1 
        end if;
        
        I_prime:= IdealConjugatesXYZ2(I_,D); 
        if IsEmpty(I_prime) eq false then
            Q_:= IdealExponentsXYZ2(I_prime,<"Alphas">);  // computes the possiblilities for [prod_{i in I_}((P_i)_i)^(d_i)], the product of split (conjugate) primes of I_ contributing to alpha 
        else
            Q_:= [[1]];   // if I_ is empty, sets [prod_{i in I_}((P_i)_i)^(d_i)] = 1
        end if;

        for PossibleAlpha in NonSplitProduct do
            for q in Q do
                for q_ in Q_ do
                    Pr:= PossibleAlpha[2]*(q[1])*(q_[1]);     // computes the ideal [prod_{i in I}(P_i)^(d_i)]*[prod_{i in I_}((P_i)_i)^(d_i)]*prod_{i = t + 1 ... s}P_i^(a_i)
                    
                    if b0 eq 1 then     // if p = 2 splits in K
                        Pr:= ideal<R|2>*Pr;     // computes the ideal (2)^(b_0)*[prod_{i in I}(P_i)^(d_i)]*[prod_{i in I_}((P_i)_i)^(d_i)]*prod_{i = t + 1 ... s}P_i^(a_i)
                    end if;
                    
                    t, gen:= IsPrincipal(Pr);   // Pr must be principal; if it is, (alpha):= (Pr)
                    if t eq true then   // if Pr is principal, necessarily, (alpha):= (Pr), hence gen = alpha // if Pr is principal, necessarily, (alpha):= (Pr), hence gen = alpha
                        
                        a1:= gen[1] + Nu*gen[2]; // writes a1 in the form a + Nu*b
                        Numericala1:= gen[1] + NumericalNu*gen[2];      // computes the numerical approximation of a1 as a real number
                        
                        if Numericala1 lt 0 then
                            a1:= -1*a1;
                        end if;
                            
                        a2:= Conjugate(a1);
                        Numericala2:= a2[1] + Sqrt(D)*a2[2];
                        if Numericala2 lt 0 then
                            a2:= -1*a2;
                        end if;
                        
                        if (Abs(Numericala1) ge Abs(Numericala2)) and (ideal<R|a1> eq ideal<R|a2>) then    
                            alpha:= a1;         // chooses alpha to be the larger value of a1, a2 (the conjugate of a1) if both a1, a2 generate the ideal (alpha)
                        elif (Abs(Numericala1) lt Abs(Numericala2)) and (ideal<R|a1> eq ideal<R|a2>) then
                            alpha:= a2;
                        else
                            alpha:= a1;         // chooses alpha to be a1 if (a1)R != (a2)R generate different ideals (alpha) as G_(a1) = +/- G_(a2)
                        end if;
                        
                        if #q eq 1 then         // #q = 2 for each q in Q if I is not empty
                            Append(~Alphas,< [],q_,alpha,PossibleAlpha[1]>);
                        elif #q_ eq 1 then      // #q_ = 2 for each q in Q_ if I_ is not empty
                            Append(~Alphas,< q,[],alpha,PossibleAlpha[1]>);
                        else
                            Append(~Alphas,< q,q_,alpha,PossibleAlpha[1]>);
                        end if;
                    end if;
                end for;
            end for;
        end for;
    end if;
    
    if (D mod 8 eq 1) or (D mod 8 eq 5) then    // if D = 1,5 mod 8, may reduce the set of alphas
        A0:= Alphas;
        for a in A0 do
            for b in A0 do
                if (a[3] eq 2*b[3]) and (a[1] eq b[1]) and (a[2] eq b[2]) then
                    Exclude(~Alphas,b);         // if alpha_i = 2*(alpha_j), discards alpha_i since both give the same solutions (x,y,z)
                elif (a[3] eq -2*b[3]) and (a[1] eq b[1]) and (a[2] eq b[2]) then
                    Exclude(~Alphas,b);         // if alpha_i = 2*(alpha_j), discards alpha_i since both give the same solutions (x,y,z)
                elif (-a[3] eq 2*b[3]) and (a[1] eq b[1]) and (a[2] eq b[2]) then
                    Exclude(~Alphas,b);         // if alpha_i = 2*(alpha_j), discards alpha_i since both give the same solutions (x,y,z)
                end if;
            end for;
        end for;
    end if;

    return Alphas;
end function;