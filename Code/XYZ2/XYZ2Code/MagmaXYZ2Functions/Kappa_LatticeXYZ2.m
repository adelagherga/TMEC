// Kappa_LatticeXYZ2.m

/*
INPUT:
    A:= [ q, q_, alpha, tplus1 ] as output by MaximalC12BoundXYZ2.m giving the largest C12star value, where 
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
    IUstar:= IU cat { i | t+1 <= i <= s, v_i != 0}, giving the largest C12star value, where
        IU:= the set of primes p_i of S such that 
            G_{alpha} = prod_{i in IU}(p_i)^(u_i)
        v_i := (1/2)*hstar*ord_{p_i}(4D) for p_i in tplus1
    p:= prime in I_
    P:= (unique) choice of prime ideal in K lying above p
    mu:= the precision on the p-adic field Q_p
    D:= squarefree part of x, where
        x:= Du^2 and x + y = z^2 
    
OUTPUT:
    B:= [b_1,...,b_{#terms}, b_0], a matrix over ZZ with columns b_i,
        b_0,b_1,...,b_{#terms} form a basis for the p-adic approximation lattice Gamma_{Kappa_,mu} associated to the linear form in logarithms, Kappa_*
    m0:= ord_p(log_p(p0))
    
REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,5,7,11,13];
    > D:= 10;
    > K<sqrtD>:= QuadraticField(D);
    > R:= RingOfIntegers(K);
    > mu:= 100;
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
    > C:= IIPrimeXYZ2(SplitPrimes,D); 
    > II_:= C[4];
    > I:= II_[1];
    > I_:= II_[2];
    > I;
    [
        <3, 2, sqrtD + 1, -sqrtD + 1, Prime Ideal of R
        Two element generators:
            3
            $.2 + 4>
    ]
    > I_;
    [
        <13, 2, 58973*sqrtD - 186489, -58973*sqrtD - 186489, Prime Ideal of R
        Two element generators:
            13
            $.2 + 7>
    ]
    > p:= I_[1][1];
    > P:= I_[1][5];
    > p;
    13
    > P; 
    Prime Ideal of R
    Two element generators:
        13
        $.2 + 7
    > A0:= AlphasXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > A:= A0[2];
    > FA:= MaximalC12BoundXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > IUstar:= FA[7];
    > IUstar;        
    [ 2, 5, 7, 11 ]
    > B, m0:= Kappa_LatticeXYZ2(A, IUstar, p, P, mu, D);
    > B;
    [1 0 0 0 0 0 0]
    [0 1 0 0 0 0 0]
    [0 0 1 0 0 0 0]
    [0 0 0 1 0 0 0]
    [0 0 0 0 1 0 0]
    [0 0 0 0 0 1 0]
    [122116377240422946867991911175862412122968424738130137386958138895213077762274\
        8142584321681651790195828390430342 5285679804821030356029055271699728310530\
        12819250563414973246057408269069940559360781761261681776796241260097742 
        131785484578212331871329215561177231948411491626379194250982918549540329671\
        251657314054070814968580551043978082 15569517504961746538708219892459828848\
        3324123067838731939754721949300507289206466823005010974594344772458982857 
        109778418406967690768941967821126102846473791523797141729654168926260230213\
        4813453536811627545659449160352113699 9054216401221582184846383081729618042\
        48448392478290306495067978596640172360332397177897560429424566838066586876 
        247933511096597253351107288473486513623877446787494114981218909940615869983\
        7975560158285662939982180192171594001]
    > m0;
    1
        
*/


function Kappa_LatticeXYZ2(A, IUstar, p, P, mu, D)
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    Qp:= pAdicField(p,mu+20);
    
    SqrtDModp:= HenselLiftXYZ2(p,P,D,mu);       // writes Sqrt(D) in Qp with precision mu
    terms, p0, m0:= Kappa_InitialSortXYZ2(A, IUstar, p, P, mu, D);
    
    p0:= Qp!p0[1] + (Qp!p0[2])*SqrtDModp;       // writes p0 as an element of Qp
    
    t:= [-pAdicLog( (Qp!t[1] + Qp!t[2]*SqrtDModp) ,p)/pAdicLog(p0 ,p) : t in terms];    // computes theta_i = -log_p(t_i)/log_p(p0) for each term t_i in Kappa_*, excluding p0
    tm:= [ConvertpAdic(Qp!i,p,mu) : i in t];    // computes theta_i^{mu} = theta_i (mod p^mu) for each term in Kappa_*, excluding p0
                   
    Append(~tm,p^mu);   // adjoins p^mu to tm
    
    row_tm := Matrix(IntegerRing(),1,#tm,tm);   // writes theta_i^{mu} values as a row matrix
    I:= ScalarMatrix(#t,1);     // #(terms)-identity matrix
    Z:= ZeroMatrix(IntegerRing(),#t,1);        // creates (#terms)x1 zero matrix
    
    B := VerticalJoin(HorizontalJoin(I,Z),row_tm);      // creates matrix associated to basis of Gamma_{Kappa_,mu}
    
    return B,m0;
end function;
