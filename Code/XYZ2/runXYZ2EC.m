/*
runXYZ2EC.m

This function solves the equation X+Y=Z^2 for non-zero integers X, Y whose prime
divisors all lie in the prime list defined in the string "set",
computes the corresponding elliptic curves, and outputs the logfile and elliptic
curve data in seperate files.

Parameters
    set: MonStgElt
        A string in the format "Nlist,primelist,[i,j]", with [i,j] denoting the
	index of the associated S-unit equation.
    dir: MonStgElt
        A string in the format "Data/[N1,N2,...]i/XYZ2" which serves as the name
	of the directory to which all output files are printed.
Returns
    OutFile: MonStgElt
        A .csv file named "N,primelist,[i,j]Out.csv" containing, for each
	elliptic curve E, the row
	"N,aInvariants(E),primelist,sol". If there are no elliptic curves
	corresponding to this set, no such file is created.
    LogFile: MonStgElt
        A .txt file named "N,primelist,[i,j]Log.txt" containing all output,
	timings, and solutions.
Authors
    Adela Gherga <adelagherga@gmail.com>
Created
    24 August 2022
*/

ChangeDirectory("./Code/XYZ2");
LogFile:="../../" cat dir cat "/XYZ2Logfiles/" cat set cat "Log.txt";
SetOutputFile(LogFile);
load "./solveXYZ2.m";
load "./convertXYZ2ToEllipticCurves.m";

Nlist,primelist,ij:=extractForm(set);
OutFile:="../../" cat dir cat "/XYZ2Outfiles/" cat set cat "Out.csv";
printf "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
printf "Nlist:=%o; ij:=%o;\n",Nlist,ij;
printf "primelist:=%o;\n",primelist;
printf "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
sols:=solveXYZ2SUnit(primelist,ij);
printf "sols:=%o\n",sols;
ECs:={};
for N in Nlist do
    ECs:=ECs join convertXYZ2ToEllipticCurves(N,sols);
end for;
printf "%o\n",ECs;
for E in ECs do
    fprintf OutFile, "%o %o\n",E[1],seqEnumToString(E[2]);
end for;
exit;
