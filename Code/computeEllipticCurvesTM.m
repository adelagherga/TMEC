/*
computeEllipticCurvesTM.m

This function solves the Thue--Mahler equations defined in the string "set",
computes the corresponding elliptic curves, and outputs the logfile and elliptic
curve data in seperate files. Note that these output files may contain multiple
Thue--Mahler form solutions, one for each integer a in "a values" and for each
monic form arising from "optimal form".

Parameters
    set: MonStgElt
        A string in the format
	"N,\"form\",\"optimal form\",\"min poly\",\"partial obstructions\",
	class number,r,no ideal eq,no Thue eq,\"a values\",\"primelist\"".
Returns
    OutFile: MonStgElt
        A .csv file named "N,[optimal form]Out.csv" containing, for each
	elliptic curve E, the row
	"N,aInvariants(E),optimal form,a,primelist,sol". If there are no
	elliptic curves corresponding to this set, no such file is created.
    LogFile: MonStgElt
        A .txt file named "N,[optimal form]Log.txt" containing
	all output, timings, and solutions.
Authors
    Adela Gherga <adelagherga@gmail.com>
Created
    24 August 2022
*/

ChangeDirectory("./Code");
load "./solveThueMahler.m";
load "./convertTMToEllipticCurves.m";

Nlist,alist,a,primelist,ij:=extractForm(set);
OutFile:="../Data/TMOutfiles/" cat set cat ".csv";
LogFile:="../Data/TMLogfiles/" cat set cat ".txt";
SetLogFile(LogFile);
printf "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
printf "N:=%o; alist:=%o; a:=%o; primelist:=%o; \n",N,alist,a,primelist;
printf "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
sols:=solveTMSUnit(alist,a,primelist,ij);
printf "sols:=%o\n",sols;
for N in Nlist do
    ECs:=convertTMToEllipticCurves(N,alist,sols);
    printf "%o\n"ECs;
    for E in ECs do
	assert E[1] eq N;
	fprintf OutFile, "%o,%o,%o,%o,%o,%o\n",
		N,seqEnumToString(E[2]),seqEnumToString(alist),a,
		seqEnumToString(primelist),seqEnumToString(E[3]);
    end for;
end for;
UnsetLogFile();
