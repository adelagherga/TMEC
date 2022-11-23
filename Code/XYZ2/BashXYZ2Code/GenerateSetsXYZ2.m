// nSetsXYZ2.m

/*
INPUT:
    m:= the upper bound on the product of the primes of S:= [2,p_1,...,p_{n-1}]

OUTPUT:
    nSets[i]:= < S_i, D, Ind >, I, I_, where
        nSets:= [S_1,...,S_n] is an exaustive list of sets S_i such that
            S_i:= [2,p_1,...,p_{n_i}] such that 2*p_1*...*p_{n_i} <= m
        D:= the values ( (p_1)^(b_1) )* /cdots *( (p_{n_i})^(b_{n_i}) ) where b_i in {0,1} corresponding to S_i
        Ind:= the index of II_ where
            II_:= [I,I_], where, assuming the primes p_i of S_i which split in Q(sqrt(D)) are p_1,...,p_t, for some t in [0,s],
                I:= { i | 1 <= i <= t, a_i > b_i} and
                I_:= { i | 1 <= i <= t, a_i < b_i}
            Nb. if D == 1 or there are no split primes (ie. II_ == []), Ind is 0

COMMENTS:
    The sets S_i are used to generate all solutions of X + Y = Z^2 with good reduction outside the primes of S_i, where
        X,Y are (S_i)-units, and Z is a rational integer
    It suffices to take the sets S_i to be maximal; that is, in the sense that S_i is not a subset of S_j for all i,j (i != j)
        ie. if S_i is a subset of S_j, then the solutions of X + Y = Z^2 corresponding to S_i will be a subset of the solutions of X + Y = Z^2 corresponding to S_j

    This algorithm outputs all maximal sets S_i whose product is <= m

EXAMPLE:
    > m:= 100;
    > nSetsXYZ2(m);
    <[ 2, 3, 5 ], 1, 0>, [], []
    <[ 2, 3, 5 ], 2, 0>, [], []
    <[ 2, 3, 5 ], 3, 0>, [], []
    <[ 2, 3, 5 ], 5, 0>, [], []
    <[ 2, 3, 5 ], 6, 1>, [ 5 ], []
    <[ 2, 3, 5 ], 10, 1>, [ 3 ], []
    <[ 2, 3, 5 ], 15, 0>, [], []
    <[ 2, 3, 5 ], 30, 0>, [], []
    <[ 2, 3, 7 ], 1, 0>, [], []
    <[ 2, 3, 7 ], 2, 1>, [ 7 ], []
    <[ 2, 3, 7 ], 3, 0>, [], []
    <[ 2, 3, 7 ], 6, 0>, [], []
    <[ 2, 3, 7 ], 7, 1>, [ 3 ], []
    <[ 2, 3, 7 ], 14, 0>, [], []
    <[ 2, 3, 7 ], 21, 0>, [], []
    <[ 2, 3, 7 ], 42, 0>, [], []
    <[ 2, 3, 11 ], 1, 0>, [], []
    <[ 2, 3, 11 ], 2, 0>, [], []
    <[ 2, 3, 11 ], 3, 1>, [ 11 ], []
    <[ 2, 3, 11 ], 6, 0>, [], []
    <[ 2, 3, 11 ], 11, 0>, [], []
    <[ 2, 3, 11 ], 22, 1>, [ 3 ], []
    <[ 2, 3, 11 ], 33, 1>, [ 2 ], []
    <[ 2, 3, 11 ], 66, 0>, [], []
    <[ 2, 3, 13 ], 1, 0>, [], []
    <[ 2, 3, 13 ], 2, 0>, [], []
    <[ 2, 3, 13 ], 3, 1>, [ 13 ], []
    <[ 2, 3, 13 ], 6, 0>, [], []
    <[ 2, 3, 13 ], 13, 1>, [ 3 ], []
    <[ 2, 3, 13 ], 26, 0>, [], []
    <[ 2, 3, 13 ], 39, 0>, [], []
    <[ 2, 3, 13 ], 78, 0>, [], []
    <[ 2, 5, 7 ], 1, 0>, [], []
    <[ 2, 5, 7 ], 2, 1>, [ 7 ], []
    <[ 2, 5, 7 ], 5, 0>, [], []
    <[ 2, 5, 7 ], 7, 0>, [], []
    <[ 2, 5, 7 ], 10, 0>, [], []
    <[ 2, 5, 7 ], 14, 1>, [ 5 ], []
    <[ 2, 5, 7 ], 35, 0>, [], []
    <[ 2, 5, 7 ], 70, 0>, [], []

*/

ChangeDirectory("./Code");
load "./parseIO.m";
SetAutoColumns(false);
SetColumns(235);

load "./XYZ2/XYZ2Code/ExponentsXYZ2.m";
load "./XYZ2/XYZ2Code/DecompositionOfPrimesXYZ2.m";
load "./XYZ2/XYZ2Code/SplitPrimePropertiesXYZ2.m";
load "./XYZ2/XYZ2Code/IIPrimeXYZ2.m";

function sortByD(S)
    Sort(~S);
    assert 2 in S;
    s:=#S;
    assert s ge 2;
    Z:= IntegerRing();
    D0:= ExponentsXYZ2(S,0);    // generates all possible values for D:= ( (p_1)^(b_1) )*
    Sort(~D0);

    caseNo:=[];
    for i in [1..#D0] do
        D:=D0[i];
	if D eq 1 then
	    Append(~caseNo,[i,0]);
        else
            b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);   // determines prime splitting
            if IsEmpty(SplitPrimes) then
		Append(~caseNo,[i,0]);
            else
                J:= IIPrimeXYZ2(SplitPrimes,D);         // generates all possible I, I_
                for j in [1..#J] do
		    Append(~caseNo,[i,j]);
		end for;
	    end if;
	end if;
    end for;
    return caseNo;
end function;

extractForm:=function(set)
    /*
      Extracts Nlist,alist,a,primelist,[i,j] from the string set.

      Parameters
          set: MonStgElt
              A string in the format "Nlist,alist,a,primelist,[i,j]".
      Returns
          Nlist: SeqEnum
              A list of conductors.
          alist: SeqEnum
              A list of coefficients a_0, a_1,...,a_3.
          a: RngIntElt
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_v.
          ij: SeqEnum
              The index (i,j) of the corresponding S-unit equation.
   */
    bracketSplit:=Split(set,"[]");
    assert (#bracketSplit eq 3);
    Nlist:=[StringToInteger(N) : N in Split(bracketSplit[1],",")];
    primelist:=[StringToInteger(p) : p in Split(bracketSplit[3],",")];
    return Nlist,primelist;
end function;





// GenerateSetsXYZ2.m

/*
INPUT:
    N/A

OUTPUT:
    nSets:= [S_1,...,S_n], where nSets is an exaustive list of sets S_i such that
        S_i:= [2,p_1,...,p_{n_i}] such that 2*p_1*...*p_{n_i} <= m

COMMENTS:
    This algorithm is an intermediary procedure between bash and magma, to output all maximal sets S_i whose product is <= m in a text file

EXAMPLE:
    N/A


*/

OutFile:="../" cat dir cat "/" cat set cat ".csv";
Nlist,primelist:=extractForm(set);
assert (IsEmpty(primelist) eq false);

caseNo:=sortByD(primelist);
for i in [1..#caseNo] do
    assert i eq caseNo[i][1];
    j:=caseNo[i][2];
    fprintf OutFile, "%o,%o,%o\n",seqEnumToString(Nlist),
	    seqEnumToString(primelist),seqEnumToString([i,j]);
end for;
exit;
