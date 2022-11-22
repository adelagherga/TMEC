// KappaLatticeXYZ2.m

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
    p:= prime in I
    P:= (unique) choice of prime ideal in K lying above p
    mu:= the precision on the p-adic field Q_p
    D:= squarefree part of x, where
        x:= Du^2 and x + y = z^2 
    
OUTPUT:
    B:= [b_1,...,b_{#terms}, b_0], a matrix over ZZ with columns b_i,
        b_0,b_1,...,b_{#terms} form a basis for the p-adic approximation lattice Gamma_{Kappa,mu} associated to the linear form in logarithms, Kappa*
    m0:= ord_p(log_p(p0))
    
REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,5,7,11,13];
    > D:= 5;
    > K<sqrtD>:= QuadraticField(D);
    > R:= RingOfIntegers(K);
    > mu:= 100;
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
    > C:= IIPrimeXYZ2(SplitPrimes,D); 
    > II_:= C[1];
    > I:= II_[1];
    > I_:= II_[2];
    > I;
    [
        <11, 1, 1/2*(-3*sqrtD - 1), 1/2*(3*sqrtD - 1), Principal Prime Ideal of R
        Generator:
            -3*$.2 + 1>
    ]
    > I_;
    []
    > p:= I[1][1];
    > P:= I[1][5];
    > p;
    11
    > P;
    Principal Prime Ideal of R
    Generator:
        -3*$.2 + 1
    > A0:= AlphasXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > A:= A0[2];
    > FA:= MaximalC12BoundXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);
    > IUstar:= FA[7];         
    > IUstar;
    [ 2, 3, 5, 7, 13 ]
    > B, m0:= KappaLatticeXYZ2(A, IUstar, p, P, mu, D);               
    > B;
    [1 0 0 0 0 0 0]
    [0 1 0 0 0 0 0]
    [0 0 1 0 0 0 0]
    [0 0 0 1 0 0 0]
    [0 0 0 0 1 0 0]
    [0 0 0 0 0 1 0]
    [600472106347985799285064285395435897714474793748768730531964215372073706004280\
        12777233697606263077445547 940469807474561234462806118349601498126313620244\
        6059218670374497769838244981259667236829613297448058534 
        113538319391720521946690240082115824079140929070860577898285428494058330779\
        03654974006308043366383079250 675733224930055104643767736900527616538005300\
        79879770289480695022583953970595364342341709992453503993090 
        215767766299338292924324060137049867250792295852823557628537487564510612459\
        35886230040171606987031745911 112807037970714697705999584014691319554650886\
        376910628344532620183281704038414319770126238698379546796914 
        137806123398222701841183371720896367762643312000384664331464775521549852095\
        523076769401159497458526446001]
    > m0;
    1

*/


function KappaLatticeXYZ2(A, IUstar, p, P, mu, D)
    K<sqrtD>:= QuadraticField(D);       // generates real quadratic field K = Q(Sqrt(D))
    Qp:= pAdicField(p,mu+20);
    
    SqrtDModp:= HenselLiftXYZ2(p,P,D,mu);       // writes Sqrt(D) in Qp with precision mu
    terms, p0, m0:= KappaInitialSortXYZ2(A, IUstar, p, P, mu, D);
    
    p0:= Qp!p0[1] + (Qp!p0[2])*SqrtDModp;       // writes p0 as an element of Qp
    
    t:= [-pAdicLog( (Qp!t[1] + Qp!t[2]*SqrtDModp) ,p)/pAdicLog(p0 ,p) : t in terms];    // computes theta_i = -log_p(t_i)/log_p(p0) for each term t_i in Kappa*, excluding p0
    tm:= [ConvertpAdic(Qp!i,p,mu) : i in t];    // computes theta_i^{mu} = theta_i (mod p^mu) for each term in Kappa*, excluding p0
                   
    Append(~tm,p^mu);   // adjoins p^mu to tm
    
    row_tm := Matrix(IntegerRing(),1,#tm,tm);   // writes theta_i^{mu} values as a row matrix
    I:= ScalarMatrix(#t,1);     // #(terms)-identity matrix
    Z:= ZeroMatrix(IntegerRing(),#t,1);        // creates (#terms)x1 zero matrix
    
    B := VerticalJoin(HorizontalJoin(I,Z),row_tm);      // creates matrix associated to basis of Gamma_{Kappa,mu}
    
    return B,m0;
end function;
