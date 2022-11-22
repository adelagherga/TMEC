// SUnitXYZ2FromCommandLineInput.m

/*
INPUT:
    S:= [ S_i, D, Ind ], where
        S_i:= [2,p_1,...,p_{n_i}] such that 2*p_1*...*p_{n_i} <= m
        D:= the values ( (p_1)^(b_1) )* /cdots *( (p_{n_i})^(b_{n_i}) ) where b_i in {0,1} corresponding to S_i
        Ind:= the index of II_ where
            II_:= [I,I_], where, assuming the primes p_i of S_i which split in Q(sqrt(D)) are p_1,...,p_t, for some t in [0,s],
                I:= { i | 1 <= i <= t, a_i > b_i} and 
                I_:= { i | 1 <= i <= t, a_i < b_i}    
            Nb. if D == 1 or there are no split primes (ie. II_ == []), Ind is 0

OUTPUT:
    FinalSolns:= [[x_1,y_1,z_1],...,[x_n,y_n,z_n]], all S-unit equations x_i + y_i = z_i^2, where
        x_i, y_i:= prod_{i:= 1 to s} p_i^{a_i} for rational integers x_i, y_i such that
            gcd(x_i, y_i) is squarefree
            x_i >= y_i and x_i >= 0
        z_i:= a rational integer, > 0 OUTPUT:
   
COMMENTS:
    Each element of S is input as Type(S[i]) = MonStgElt in Magma, where the symbols "[", "[", and "," are included. 
    This algorithm converts the primes of S into RngIntElts in Magma so that they may be used as a valid input in SUnitXYZ2.m

EXAMPLE:
    $ $MAGMA S:=[[2,3,5,7,11],22,3],[7],[] $ABSOLUTE_PATH/SUnitCodede/XYZ2Code/BashToMagmaRunSUnit/SUnitXYZ2FromCommandLineInput.m
    Magma V2.19-9     Sat Dec  3 2016 12:51:33 on thumper  [Seed = 3754211116]
    Type ? for help.  Type <Ctrl>-D to quit.

    Loading file "/Users/Adela/Documents/SUnitCode/XYZ2Code/BashToMagmaRunSUnit/SUnitXYZ2FromCommandLineInput.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/SUnitDXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZFunctions/SUnitXYZ.m"
    Loading "./FinckePohstCode/F3Approx.m"
    Loading "./FinckePohstCode/MyCholesky.m"
    Loading "./FinckePohstCode/ShortVectors.m"
    Loading "./FinckePohstCode/FinckePohst.m"
    Loading "./XYZ2Code/MagmaXYZFunctions/InitialSortXYZ.m"
    Loading "./XYZ2Code/MagmaXYZFunctions/LatticeXYZ.m"
    Loading "./XYZ2Code/MagmaXYZFunctions/LatticeBoundXYZ.m"
    Loading "./XYZ2Code/MagmaXYZFunctions/SmallestLatticeBoundXYZ.m"
    Loading "./XYZ2Code/MagmaXYZFunctions/FPRestrictionsXYZ.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/SUnitXYZtoSUnitXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/ExponentsXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/DecompositionOfPrimesXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/IdealExponentsXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/IdealConjugatesXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/AlphasXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/FundamentalUnitXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/IUFactorsXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/SymmetricCaseZerosXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/SymmetricCaseXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/SplitPrimePropertiesXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/IIPrimeXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/nExponentsXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/C1pnXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/LambdaBoundXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/KappaBoundXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/Kappa_BoundXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/C12BoundXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/MaximalC12BoundXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/LambdaLogpXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/LambdaInitialSortXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/LambdaLatticeXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/HenselLiftXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/KappaInitialSortXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/Kappa_InitialSortXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/KappaLatticeXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/Kappa_LatticeXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/SmallestLatticeBoundXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/FPRestrictionsKappaXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/FPRestrictionsKappa_XYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/FPRestrictionsLambdaXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/FPRestrictionsXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/FPParametersXYZ2.m"
    Loading "./XYZ2Code/MagmaXYZ2Functions/FinalSearchXYZ2.m"
    S = [ 2, 3, 5, 7, 11 ], D = 22, I = [ 7 ], I_ = []
    [
        [ 126720000, 49, 11257 ],
        [ 117649, 12672, 361 ],
        [ 38808, 1, 197 ],
        [ 5632, -7, 75 ],
        [ 792, 49, 29 ],
        [ 352, -343, 3 ],
        [ 19800, 2401, 149 ],
        [ 88, -7, 9 ],
        [ 28512, 49, 169 ],
        [ 16038, -4802, 106 ],
        [ 198, -2, 14 ],
        [ 198, -98, 10 ],
        [ 66550, 14, 258 ],
        [ 22, 14, 6 ],
        [ 1078, 11, 33 ],
        [ 198, -77, 11 ],
        [ 550, 539, 33 ],
        [ 123750, 154, 352 ],
        [ 7546, 198, 88 ],
        [ 1782, 154, 44 ],
        [ 22, -22, 0 ]
    ]

    Total time: 1.260 seconds, Total memory usage: 32.09MB

*/


load "./XYZ2Code/MagmaXYZ2Functions/SUnitDXYZ2.m";

s:= #S;         // computes the total number of elements of S, including symbols "[", "[", and ","

sb:= [];        // stores indices of all brackets "["; this set should have exactly 3 elements, one for each set S, I, I_
eb:= [];        // stores indices of all brackets "]"; this set should have exactly 3 elements, one for each set S, I, I_
cs:= [];        // stores indices of all commas "," occuring between the end of the set NewS and "]" in S:= [NewS,D,Ind],I,I_

for i in [1..s] do  
    if S[i] eq "[" then
        Append(~sb, i);
    elif S[i] eq "]" then
        Append(~eb, i);
    end if;
end for;

for i in [eb[1]..eb[2]] do
    if S[i] eq "," then
        Append(~cs,i);
    end if;
end for;

NewS:= [];      // stores [p_1,...,p_s], p_i primes in S, as RngIntElts in Magma
j:= sb[2] + 1;  // sets the index of the first element of the set, NewS; ie. S:=[NewS,D,Ind],I,I_, where NewS:= [p_1,...,p_s]
while j le eb[1] do
    i:= 0;
    n:= [];     // stores consecutive symbols of S which are not "[", "[", or "," 
    while (j+i le eb[1]) and (S[j+i] ne "[") and (S[j+i] ne "]") and (S[j+i] ne ",") do         // ignores non-integer members of S
        Append(~n, S[j+i]);
        i:= i + 1;
    end while;
    j:= j+i+1;
    if IsEmpty(n) eq false then
        Append(~NewS, StringToInteger(&cat(n))); 
    end if;
end while;

D0:= [S[i] : i in [(cs[1]+1)..(cs[2]-1)]];      // stores all integer members of S between the symbols "," following NewS and the symbols "," before Ind; ie. those of D 
D:= StringToInteger(&cat(D0));

I0:= [S[i] : i in [(cs[2]+1)..(eb[2]-1)]];      // stores all integer members of S between the symbols "," following Ind and the symbols "]"; ie. those of Ind 
Ind:= StringToInteger(&cat(I0));

I:= [];         // stores [p_1,...,p_{s_I}], p_i primes in I, as RngIntElts in Magma
j:= sb[3] + 1;  // sets the index of the first element of the set, I; ie. S:=[NewS,D,Ind],I,I_, where I:= [p_1,...,p_{s_I}]
while j le eb[3] do
    i:= 0;
    n:= [];     // stores consecutive symbols of S which are not "[", "[", or "," 
    while (j+i le eb[3]) and (S[j+i] ne "[") and (S[j+i] ne "]") and (S[j+i] ne ",") do         // ignores non-integer members of S
        Append(~n, S[j+i]);
        i:= i + 1;
    end while;
    j:= j+i+1;
    if IsEmpty(n) eq false then
        Append(~I, StringToInteger(&cat(n))); 
    end if;
end while;

I_:= [];        // stores [p_1,...,p_{s_I_}], p_i primes in I_, as RngIntElts in Magma
j:= sb[4] + 1;  // sets the index of the first element of the set, I; ie. S:=[NewS,D,Ind],I,I_, where I_:= [p_1,...,p_{s_I_}]
while j le eb[4] do
    i:= 0;
    n:= [];     // stores consecutive symbols of S which are not "[", "[", or "," 
    while (j+i le eb[4]) and (S[j+i] ne "[") and (S[j+i] ne "]") and (S[j+i] ne ",") do         // ignores non-integer members of S
        Append(~n, S[j+i]);
        i:= i + 1;
    end while;
    j:= j+i+1;
    if IsEmpty(n) eq false then
        Append(~I_, StringToInteger(&cat(n))); 
    end if;
end while;

if #NewS ge 3 then      // computes SUnitXYZ2 only for those sets NewS which have at least 3 primes
    Sol:= SUnitDXYZ2(NewS,D,Ind);       // if D == 1 or (I == [] and I_ == []) corresponding to the symmetric case, Ind is set to 0; ie. does not parallelize among I,I_
    printf "S = %o, D = %o, I = %o, I_ = %o\n", NewS, D, I, I_;         // prints NewS, D, I, I_
    print Sol;
end if;
exit;   // exits Magma to enforce garbage disposal
