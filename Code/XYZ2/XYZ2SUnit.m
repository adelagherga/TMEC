/*
XYZ2SUnit.m

This function generates all quadratic fields associated to X+Y=Z^2, where X, Y
are made up of {p_1,...,p_s}, and determines the indices to iterate over all
associated S-unit equations.

Parameters
    set: MonStgElt
        A string in the format "Nlist,primelist".
    dir: MonStgElt
        A string in the format "Data/[N1,N2,...]i/XYZ2" which serves as the name
	of the directory to which all output files are printed.
Returns
    OutFile: MonStgElt
        A .csv file named Nlist,primelist.csv containing rows written in the
	format "Nlist,primelist,[i,j]", with i,j iterating through
	the indices of the associated S-unit equations to solve.
Authors
    Adela Gherga <adelagherga@gmail.com>
Created
    23 November 2022
*/

SetOutputFile(dir cat "/" cat set cat "tmp.txt");
ChangeDirectory("./Code/XYZ2");
load "./SUnitXYZ2.m";

XYZ2Sunit:=function(primelist)
    /*
      Generate all quadratic fields associated to X+Y=Z^2, where X, Y are made up
      of {p_1,...,p_s}, and determine the indices to iterate over all associated
      S-unit equations.

      Parameters
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_s.
      Returns
          caseNo_i: SeqEnum
              A list of elements [i,j] corresponding to primelist, where i
              denotes the index of the potential quadratic field generator D
	      associated to primelist, and j denotes the number of S-unit
	      equations corresponding to that field.
   */
    assert &and[IsPrime(p) : p in primelist];
    assert 2 in primelist;
    Sort(~primelist);
    s:=#primelist;
    assert s ge 2;
    C:=CartesianProduct([[0,1] : i in [1..s]]);
    DList:=Sort([&*[primelist[i]^c[i] : i in [1..s]] : c in C]);
    // Generate all possible values for D.
    caseNo:=[];
    for i in [1..#DList] do
        D:=DList[i];
	if D eq 1 then
	    Append(~caseNo,[i,0]);
        else
            b0,SplitPrimes,NonSplitPrimes:=DecompositionOfPrimesXYZ2(S,D);
            if IsEmpty(SplitPrimes) then
		Append(~caseNo,[i,0]);
            else
                J:=IIPrimeXYZ2(SplitPrimes,D);
                for j in [1..#J] do
		    Append(~caseNo,[i,j]);
		end for;
	    end if;
	end if;
    end for;
    return caseNo;
end function;

OutFile:="../../" cat dir cat "/" cat set cat ".csv";
Nlist,primelist,_:=extractForm(set);
assert (IsEmpty(primelist) eq false);
caseNo:=sortByD(primelist);
for i in [1..#caseNo] do
    assert i eq caseNo[i][1];
    j:=caseNo[i][2];
    fprintf OutFile, "%o,%o,%o\n",seqEnumToString(Nlist),
	    seqEnumToString(primelist),seqEnumToString([i,j]);
end for;
exit;
