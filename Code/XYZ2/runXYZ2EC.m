/*
runXYZ2EC.m
// EDIT STILL!!!!!!!!!!!!!!!!!!!!!!!!!!
This function solves the Thue--Mahler equations defined in the string "set",
computes the corresponding elliptic curves, and outputs the logfile and elliptic
curve data in seperate files.

Parameters
    set: MonStgElt
        A string in the format "Nlist,alist,a,primelist,[i,j]", with [i,j]
	denoting the index of the associated S-unit equation.
    dir: MonStgElt
        A string in the format "Data/[N1,N2,...]i/TM" which serves as the name of
	the directory to which all output files are printed.
Returns
    OutFile: MonStgElt
        A .csv file named "N,alist,a,primelist,[i,j]Out.csv" containing, for
	each elliptic curve E, the row
	"N,aInvariants(E),alist,a,primelist,sol". If there are no
	elliptic curves corresponding to this set, no such file is created.
    LogFile: MonStgElt
        A .txt file named "N,alist,a,primelist,[i,j]Log.txt" containing
	all output, timings, and solutions.
Authors
    Adela Gherga <adelagherga@gmail.com>
Created
    24 August 2022
*/

ChangeDirectory("./Code/TM");
LogFile:="../../" cat dir cat "/TMLogfiles/" cat set cat "Log.txt";
SetOutputFile(LogFile);
load "./solveThueMahler.m";
load "./convertTMToEllipticCurves.m";

Nlist,alist,a,primelist,ij:=extractForm(set);
OutFile:="../../" cat dir cat "/TMOutfiles/" cat set cat "Out.csv";
printf "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
printf "Nlist:=%o; ij:=%o;\n",Nlist,ij;
printf "alist:=%o; a:=%o; primelist:=%o;\n",alist,a,primelist;
printf "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
sols:=solveTMSUnit(alist,a,primelist,ij);
printf "sols:=%o\n",sols;
ECs:={};
for N in Nlist do
    ECs:=ECs join convertTMToEllipticCurves(N,alist,sols);
end for;
printf "%o\n",ECs;
for E in ECs do
    fprintf OutFile, "%o %o\n",E[1],seqEnumToString(E[2]);
end for;
exit;
