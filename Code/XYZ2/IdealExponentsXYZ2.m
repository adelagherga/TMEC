// IdealExponentsXYZ2.m

/*
INPUT:
    Iset:= I or I_ (as output by IIPrimeXYZ2.m), where
        I:= { i | 1 <= i <= t, a_i > b_i}
            I:= [[p_1, h_1, pi_1, (pi_)_1, P_1],...,[p_n, h_n, pi_n, (pi_)_n, P_n]],
                p_i:= prime of S which splits in K
                h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
                pi_i:= generator of the ideal P^(h_i) = (pi)R, chosen so that Abs(pi) > Abs(pi_)
                (pi_)_i:= conjugate of pi, chosen so that Abs(pi) > Abs(pi_)
                P_i:= (unique) choice of prime ideal in K lying above p
        I_:= { i | 1 <= i <= t, a_i < b_i}
            I_:= [[p_1, h_1, pi_1, (pi_)_1, P_1],...,[p_n, h_n, pi_n, (pi_)_n, P_n]],
                p_i:= prime of S which splits in K
                h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
                pi_i:= generator of the ideal P^(h_i) = (pi)R, chosen so that Abs(pi) < Abs(pi_)
                (pi_)_i:= conjugate of pi, chosen so that Abs(pi) < Abs(pi_)
                P_i:= (unique) choice of prime ideal in K lying above p

    type:= string designating which subroutine this algorithm is to be used in
        ["Alphas"]:= used to compute every possible D[j], hence every possible alpha, where
            a:= (alpha) = (2)^(b_0)*D[j]*prod_{i = t + 1 ... s} P_i^(a_i)
            D[j]:= [prod_{i in I}(P_i)^(d_i)]*[prod_{i in I_}(P_i')^(d_i)]

        ["IUFactors"]:= used to compute every possible Prod[I], Prod[I_], hence every possible G_{alpha}, where
            G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_])
            Prod[I]:= prod_{i in I} (pi_i)^{m_i}
            Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}, where
                m_i:= an integer in [0,...,5]

        ["IUFactors_"]:= used to compute every possible Prod_[I], Prod_[I_], hence every possible G_{alpha}, where
            G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_])
            Prod_[I]:= prod_{i in I} ((pi_)_i)^{m_i}
            Prod_[I_]:= prod_{i in I_} ((pi_)_i)^{m_i}, where
                m_i:= an integer in [0,...,5]

        ["FinalSearch", M]:= used to compute every possible Prod[I], Prod[I_], hence every possible G_{alpha}, where
            G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_])
            Prod[I]:= prod_{i in I} (pi_i)^{m_i}
            Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}
            M:= M0 or M0_ (as output by SmallestLatticeBoundXYZ2.m) such that
                M0:= [ [c_i, p_i] : p_i in I], where
                    c_i:= the maximal possible exponents m_i on Prod[I] (and Prod_[I]) in G_{alpha} giving G_{alpha} = u
                        m_i:= an integer in [0,...,c_i]
                M0_:= [ [c_i, p_i] : p_i in I_], where
                    c_i:= the maximal possible exponents m_i on Prod[I_] (and Prod_[I_]) in G_{alpha} giving G_{alpha} = u
                        m_i:= an integer in [0,...,c_i]

        ["FinalSearch_", M]:= used to compute every possible Prod_[I], Prod_[I_], hence every possible G_{alpha}, where
            G_{alpha}:= (a/(2*sqrtD))*((eps)^(n))*(Prod[I])*(Prod_[I_]) - (a_/(2*sqrtD))*((eps_)^(n))*(Prod_[I])*(Prod[I_])
            Prod[I]:= prod_{i in I} (pi_i)^{m_i}
            Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}
            M:= M0 or M0_ (as output by SmallestLatticeBoundXYZ2.m) such that
                M0:= [ [c_i, p_i] : p_i in I], where
                    c_i:= the maximal possible exponents m_i on Prod_[I] (and Prod[I]) in G_{alpha} giving G_{alpha} = u
                        m_i:= an integer in [0,...,c_i]
                M0_:= [ [c_i, p_i] : p_i in I_], where
                    c_i:= the maximal possible exponents m_i on Prod_[I_] (and Prod[I_]) in G_{alpha} giving G_{alpha} = u
                        m_i:= an integer in [0,...,c_i]

OUTPUT:
    if type:= ["Alphas"]:
        D:= [c_1,...,c_m], where
            c_i:= prod_i, [< p, h_i, pi, pi_, P , d_i > : j in [1..n] ], such that
                prod_i:= the product of the P_i with exponent d_i, where 0 <= d_i <= (h_i)-1
                p:= prime of S which splits in K
                h_i:= smallest positive integer such that P^(h_i) = (pi) is principal, where P lies above p in K
                pi:= generator of the ideal P^(h_i) = (pi)R
                pi_:= conjugate of pi
                P:= (unique) choice of prime ideal in K lying above p
                d_i:= the exponent on the prime ideal P_i in prod_i
                n:= the length of Iset

    if type:= "IUFactors":
        D:= [c_1,...,c_m], all possible products of the pi_i with exponents m_i, where 1 <= m_i <= 5 and
            if Iset:= I:
                c_i:= Prod[I]:= prod_{i in I} (pi_i)^{m_i}
            if Iset:= I_:
                c_i:= Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}

    if type:= "IUFactors_":
        D:= [c_1,...,c_m], all possible products of the (pi_)_i with exponents m_i, where 1 <= m_i <= 5 and
            if Iset:= I:
                c_i:= Prod_[I] = prod_{i in I} ((pi_)_i)^{m_i}
            if Iset:= I_:
                c_i:= Prod_[I_] = prod_{i in I_} ((pi_)_i)^{m_i}

    if type:= ["FinalSearch", M]:
        D:= [d_1,...,d_m], all possible products of the pi_i with exponents m_i, where 1 <= m_i <= M[i] and
            if Iset:= I, hence M:= M0:
                d_i:= Prod[I]:= prod_{i in I} (pi_i)^{m_i}
            if Iset:= I_, hence M:= M0_:
                d_i:= Prod[I_]:= prod_{i in I_} (pi_i)^{m_i}

    if type:= ["FinalSearch_", M]:
        D:= [d_1,...,d_m], all possible products of the (pi_)_i with exponents m_i, where 1 <= m_i <= M[i] and
            if Iset:= I, hence M:= M0:
                d_i:= Prod_[I] = prod_{i in I} ((pi_)_i)^{m_i}
            if Iset:= I_, hence M:= M0_:
                d_i:= Prod_[I_] = prod_{i in I_} ((pi_)_i)^{m_i}

COMMENTS:
    if type:= ["Alphas"]:
        Computes every possible product of (splitting) prime ideals to appear in alpha:

    if type:= ["IUFactors"]:
        Computes every possible product of the pi_i to appear in G_{alpha}, where
            pi_i:= an ideal of K:= Q(Sqrt(D)) such that
                (pi_i)R:= (P_i)^{h_i}
            P_i:= a splitting prime ideal above some prime p of S

    if type:= ["IUFactors_"]:
        Computes every possible product of the (pi_)_i to appear in G_{alpha}, where
            (pi_)_i:= the conjugate of the element pi_i of K:= Q(Sqrt(D)), where
                (pi_i)R:= (P_i)^{h_i}
            P_i:= a splitting prime ideal above some prime p of S

     if type:= ["FinalSearch", M]:
        Computes every possible product of the pi_i to appear in G_{alpha} with exponents <= M[i], where
            pi_i:= an ideal of K:= Q(Sqrt(D)) such that
                (pi_i)R:= (P_i)^{h_i}
            P_i:= a splitting prime ideal above some prime p of S

    if type:= ["FinalSearch_", M]:
        Computes every possible product of the (pi_)_i to appear in G_{alpha} with exponents <= M[i], where
            (pi_)_i:= the conjugate of the element pi_i of K:= Q(Sqrt(D)), where
                (pi_i)R:= (P_i)^{h_i}
            P_i:= a splitting prime ideal above some prime p of S

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,11,23];
    > D:= 253;
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
    > II_:= (IIPrimeXYZ2(SplitPrimes,D))[1];
    > Iset:= II_[1];
    > Iset;
    [
        <3, 1, -sqrtD - 16, sqrtD - 16, Principal Prime Ideal
        Generator:
            -2*$.2 - 15>
    ]
    > type:= ["Alphas"];
    > IdealExponentsXYZ2(Iset,["Alphas"]);
    [
        <Principal Ideal
        Generator:
            1, [
            <3, 1, -sqrtD - 16, sqrtD - 16, Principal Prime Ideal
            Generator:
                -2*$.2 - 15, 0>
        ]>
    ]

    > D:= 7;
    > S:= [2,3,5,7,11,13];
    > K<sqrtD>:= QuadraticField(D);
    > R:= RingOfIntegers(K);
    > b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
    > C:= IIPrimeXYZ2(SplitPrimes,D);
    > II_:= C[1];
    > Iset:= II_[1];
    > Iset;
    [
        <3, 1, sqrtD + 2, -sqrtD + 2, Principal Prime Ideal of R
        Generator:
            $.2 + 2>
    ]
    > type:= ["IUFactors"];
    > IdealExponentsXYZ2(Iset,type);
    [ 1, sqrtD + 2, 4*sqrtD + 11, 19*sqrtD + 50, 88*sqrtD + 233, 409*sqrtD + 1082 ]

    > type:= ["IUFactors_"];
    > IdealExponentsXYZ2(Iset,type);
    [ 1, -sqrtD + 2, -4*sqrtD + 11, -19*sqrtD + 50, -88*sqrtD + 233, -409*sqrtD +
        1082 ]

    > M0;
    [
        [ 4, 7 ],
        [ 6, 11 ]
    ]
    > type:= <"FinalSearch", M0>;
    > IdealExponentsXYZ2(Iset,type);
    [ 1, -sqrtD - 2, 4*sqrtD + 19, -27*sqrtD - 98, 152*sqrtD + 601, -905*sqrtD -
        3482, 5292*sqrtD + 20539, sqrtD + 8, -10*sqrtD - 31, 51*sqrtD + 212,
        -314*sqrtD - 1189, 1817*sqrtD + 7088, -10722*sqrtD - 41431, 62875*sqrtD +
        243692, 16*sqrtD + 79, -111*sqrtD - 398, 620*sqrtD + 2461, -3701*sqrtD -
        14222, 21624*sqrtD + 83959, -127207*sqrtD - 492278, 746692*sqrtD + 2892661,
        207*sqrtD + 872, -1286*sqrtD - 4849, 7421*sqrtD + 28988, -43830*sqrtD - 169291,
        256951*sqrtD + 996032, -1509934*sqrtD - 5846329, 8866197*sqrtD + 34341668,
        2528*sqrtD + 10081, -15137*sqrtD - 58082, 88356*sqrtD + 343219,
        -519931*sqrtD - 2011778, 3051640*sqrtD + 11822521, -17925801*sqrtD -
        69419642, 105271244*sqrtD + 407726299 ]

    > type:= <"FinalSearch_", M0>;
    > IdealExponentsXYZ2(Iset,type);
    [ 1, sqrtD - 2, -4*sqrtD + 19, 27*sqrtD - 98, -152*sqrtD + 601, 905*sqrtD -
        3482, -5292*sqrtD + 20539, -sqrtD + 8, 10*sqrtD - 31, -51*sqrtD + 212,
        314*sqrtD - 1189, -1817*sqrtD + 7088, 10722*sqrtD - 41431, -62875*sqrtD +
        243692, -16*sqrtD + 79, 111*sqrtD - 398, -620*sqrtD + 2461, 3701*sqrtD -
        14222, -21624*sqrtD + 83959, 127207*sqrtD - 492278, -746692*sqrtD + 2892661,
        -207*sqrtD + 872, 1286*sqrtD - 4849, -7421*sqrtD + 28988, 43830*sqrtD - 169291,
        -256951*sqrtD + 996032, 1509934*sqrtD - 5846329, -8866197*sqrtD + 34341668,
        -2528*sqrtD + 10081, 15137*sqrtD - 58082, -88356*sqrtD + 343219,
        519931*sqrtD - 2011778, -3051640*sqrtD + 11822521, 17925801*sqrtD -
        69419642, -105271244*sqrtD + 407726299 ]

*/


function IdealExponentsXYZ2(Iset,type)
    n:= #Iset;  // computes the length of Iset
    if type[1] eq "Alphas" then
        Hs:= [[0..d_i[2]-1]: d_i in Iset];      // compiles the ranges of all possible d_i values from Iset, where 0 <= d_i <= (h_i) - 1
        ExponentCombos:= CartesianProduct(Hs);  // computes every possible combination of d_i values
        D:= [];         // stores all possible products of elements of Iset raised to the respective exponent d_i of ExponentCombos, along with Iset, d_i values

        for c in ExponentCombos do
            Append(~D,< &*{Iset[j][5]^c[j] : j in [1..n]}, [ < Iset[j][1], Iset[j][2], Iset[j][3], Iset[j][4], Iset[j][5], c[j] > : j in [1..n] ] >);
        end for;

    elif type[1] eq "IUFactors" then
        Hs:= [[0..5]: i in [1..n]];     // compiles the ranges [0..5] only, representing the potential m_i values
        ExponentCombos:= CartesianProduct(Hs);  // computes every possible combination of the potential m_i values
        D:= [];         // stores all possible products of elements of Iset raised to the respective exponent m_i of ExponentCombos

        for c in ExponentCombos do
            Append(~D,&*{Iset[j][3]^c[j] : j in [1..n]});
        end for;

    elif type[1] eq "IUFactors_" then
        Hs:= [[0..5]: i in [1..n]];     // compiles the ranges [0..5] only, representing the potential m_i values
        ExponentCombos:= CartesianProduct(Hs);  // computes every possible combination of the potential m_i values
        D:= [];         // stores all possible products of elements of Iset raised to the respective exponent m_i of ExponentCombos

        for c in ExponentCombos do
            Append(~D,&*{Iset[j][4]^c[j] : j in [1..n]});
        end for;

    elif type[1] eq "FinalSearch" then
        Hs:= [[0..type[2][i][1]]: i in [1..n]];         // compiles the ranges [0..M[i]] only, representing the potential m_i values
        ExponentCombos:= CartesianProduct(Hs);  // computes every possible combination of the potential m_i values
        D:= [];         // stores all possible products of elements of Iset raised to the respective exponent m_i of ExponentCombos

        for c in ExponentCombos do
            Append(~D,&*{Iset[j][3]^c[j] : j in [1..n]});
        end for;

    elif type[1] eq "FinalSearch_" then
        Hs:= [[0..type[2][i][1]]: i in [1..n]];         // compiles the ranges [0..M[i]] only, representing the potential m_i values
        ExponentCombos:= CartesianProduct(Hs);  // computes every possible combination of the potential m_i values
        D:= [];         // stores all possible products of elements of Iset raised to the respective exponent m_i of ExponentCombos

        for c in ExponentCombos do
            Append(~D,&*{Iset[j][4]^c[j] : j in [1..n]});
        end for;

    end if;
    return D;
end function;
