#!/usr/bin/python
# gatherRedundancy.py

# These functions remove redundant equations across conductors for the
# Thue--Mahler solver and the XYZ2 solver. In the Thue--Mahler case, for any 2
# lines of the file IF, "N_1,alist_1,a_1,primelist_1" and
# "N_2,alist_2,a_2,primelist_2", where alist_1 = alist_2, a_1 = a_2, and
# primelist_1 is a subset of primelist_2, writes
# "[N_1,N_2],alist_1,a_1,primelist_2" in the file OF. In the XYZ2 case, for any 2
# lines of the file IF, "N_1,primelist_1" and "N_2,primelist_2", if primelist_1
# is a subset of primelist_2, writes "primelist_2,[N_1,N_2]" in the file OF.

# Parameters
#     IF: <class 'str'>
#         The string denoting the input file, wherein each line is in the
#         format "N,alist,a,primelist" in the Thue--Mahler case and
#         "N: primelist" in the XYZ2 case.
#     OF: <class 'str'>
#         A string denoting an empty output file.
# Returns
#     OF: <class 'str'>
#         The string denoting the output file, wherein each line is in the
#         "Nlist,alist,a,primelist" in the Thue--Mahler case and
#         "primelist,Nlist" in the XYZ2 case.
# Authors
#     Adela Gherga <adelagherga@gmail.com>
# Created
#     28 September 2022

import sys
import os
from collections import OrderedDict

def parseXYZ2Form(line):
    """
      Extracts N, primelist from the string line.

      Parameters
          line: <class 'str'>
              A string in the format "N: primelist".
      Returns
          N: <class 'int'>
          primelist: <class 'list'>
          sprimelist: <class 'str'>
    """
    spaceSplit=line.split(" ")
    assert (spaceSplit[0][-1] == ":")
    spaceSplit.insert(1,'2')
    N=int(spaceSplit[0][0:-1])
    primelist0=list(OrderedDict.fromkeys(spaceSplit[1:]))
    primelist=[int(p) for p in primelist0]
    sprimelist="["+''.join([p+"," for p in primelist0[:-1]])+primelist0[-1]+"]"
    return N,primelist,sprimelist

def parseTMForm(line):
    """
      Extracts N,alist,a,primelist from the string line.

      Parameters
          line: <class 'str'>
              A string in the format "N,alist,a,primelist".
      Returns
          N: <class 'int'>
          alist: <class 'str'>
          aprimelist: <class 'str'>
    """
    commaSplit=line.split(",")
    assert (commaSplit[1][0] == "[")
    assert (commaSplit[4][len(commaSplit[4])-1] == "]")
    N=int(commaSplit[0])
    alist="["+line.split("[")[1].split("]")[0]+"]"
    aprimelist=commaSplit[5]+",["+line.split("[")[2].split("]")[0]+"]"
    return N,alist,aprimelist

def extractTMPrimelist(line):
    """
      Extracts a,primelist from the string line.

      Parameters
          line: <class 'str'>
              A string in the format "a,primelist".
      Returns
          a: <class 'int'>
          aprimelist: <class 'list'>
    """
    commaSplit=line.split(",")
    bracketSplit=line.split("[")
    assert (commaSplit[1][0] == "[")
    a=int(commaSplit[0])
    if bracketSplit[1] == "]":
        primelist=[]
    else:
        primelist=[int(i) for i in bracketSplit[1].split("]")[0].split(",")]
    return a,primelist

def prime_divisors(n):
    """
      Returns all the prime divisors of the positive integer n.

      Parameters
          n: <class 'int'>
      Returns
          factors: <class 'set'>
              The list of prime divisors of n.
    """
    factors=set()
    d=2
    while n>1:
        while n%d == 0:
            factors.add(d)
            n/=d
        d=d+1
        if d*d > n:
            if n > 1:
                factors.add(int(n))
            break
    return sorted(factors)

def determineOF(IF):
    ind=[pos for pos, char in enumerate(IF) if char == "/"][-1]
    assert IF[ind] == "/"
    OF=IF[:ind+1]+"tmp"+IF[ind+1:]
    print(OF)
    return OF

def gatherXYZ2Redundancy(IF):
    """
      Removes redundant conductor factorizations and appends 2 to each list.
      That is, for any 2 lines of the file IF, "N_1,primelist_1" and
      "N_2,primelist_2", if primelist_1 is a subset of primelist_2, writes
      "primelist_2,[N_1,N_2]" in the file OF.

      Parameters
          IF: <class 'str'>
              The string denoting the input file, wherein each line is in the
              format "N: primelist".
          OF: <class 'str'>
              A string denoting an empty output file.
      Returns
          OF: <class 'str'>
              The string denoting the output file, wherein each line is in the
              format "primelist,Nlist", with Nlist now a list of
              conductors.
    """
    OF=determineOF(IF)
    primesN={}
    for line in open(IF):
        # Sort data by primelist.
        N,primelist,sprimelist=parseXYZ2Form(line.rstrip())
        if sprimelist in primesN:
            primesN[sprimelist].append([N,primelist])
        else:
            primesN[sprimelist]=[[N,primelist]]

    for sprimelist in sorted(primesN):
        # Given 2 sets of conductor factorizations, (Nlist1,primelist1)
        # and (Nlist2,primelist2), if primelist1 is contained in primelist2,
        # collapse both conductor lists into
        # (Nlist1+Nlist2,primelist2).
        primelist=primesN[sprimelist][0][1]
        for sprimelist2 in sorted(primesN):
            if sprimelist != sprimelist2:
                primelist2=primesN[sprimelist2][0][1]
                if set(primelist) <= set(primelist2):
                    primesN[sprimelist2]=primesN[sprimelist2]+primesN[sprimelist]
                    del primesN[sprimelist]
                    break

    OutFile=open(OF,"w")
    # Output all data to the file OF in the format "primelist,Nlist".
    for sprimelist in sorted(primesN):
        Nlist=str(sorted([S[0] for S in primesN[sprimelist]])).replace(" ","")
        out=Nlist+","+sprimelist
        OutFile.write("%s\n" % out)
    OutFile.close()
    os.rename(OF,IF)

def gatherTMFormRedundancy(IF):
    """
      Removes redundant Thue--Mahler equations across conductors. That is, for
      any 2 lines of the file IF, "N_1,alist_1,a_1,primelist_1" and
      "N_2,alist_2,a_2,primelist_2", where alist_1 = alist_2, a_1 = a_2,
      and primelist_1 is a subset of primelist_2, writes
      "[N_1,N_2],alist_1,a_1,primelist_2" in the file OF.

      Parameters
          IF: <class 'str'>
              The string denoting the input file, wherein each line is in the
              format "N,alist,a,primelist".
          OF: <class 'str'>
              A string denoting an empty output file.
      Returns
          OF: <class 'str'>
              The string denoting the output file, wherein each line is in the
              format "Nlist,alist,a,primelist", with Nlist now a list of
              conductors.
    """
    OF=determineOF(IF)
    forms={}
    for line in open(IF):
        # Sort data by alist.
        N,alist,aprimelist=parseTMForm(line)
        if alist in forms:
            forms[alist].append([N,aprimelist])
        else:
            forms[alist]=[[N,aprimelist]]

    formsRHS={}
    for alist in sorted(forms):
        formsRHS[alist]={}
        for rhsN in forms[alist]:
            # Sort data by alist,aprimelist.
            N=rhsN[0]
            aprimelist=rhsN[1]
            if aprimelist in formsRHS[alist]:
                formsRHS[alist][aprimelist].append(N)
            else:
                formsRHS[alist][aprimelist]=[N]

    for alist in sorted(formsRHS):
        for aprimelist in sorted(formsRHS[alist]):
            # For a given form with 2 rhs possibilities, (Nlist1,a1,primelist1)
            # and (Nlist2,a2,primelist2), if a1 = a2 and primelist1 is contained
            # in primelist2, collapse both rhs possibilities into
            # (Nlist1+Nlist2,a1,primelist2).
            Nlist=formsRHS[alist][aprimelist]
            a,primelist=extractTMPrimelist(aprimelist)
            for aprimelist2 in sorted(formsRHS[alist]):
                if aprimelist != aprimelist2:
                    a2,primelist2=extractTMPrimelist(aprimelist2)
                    if (a == a2) and (set(primelist) <= set(primelist2)):
                        formsRHS[alist][aprimelist2]=(
                            formsRHS[alist][aprimelist2]+Nlist)
                        formsRHS[alist][aprimelist2].sort()
                        del formsRHS[alist][aprimelist]
                        break

    for alist in sorted(formsRHS):
        for aprimelist in sorted(formsRHS[alist]):
            # For a given form with 2 rhs possibilities, (Nlist1,a1,primelist1)
            # and (Nlist2,a2,primelist2), if a1 != 1 and a2 = 1, primelist1 is
            # contained in primelist2, and the prime divisors of a1 are contained
            # in primelist2, collapse both rhs possibilities into
            # (Nlist1+Nlist2,a2,primelist2).
            Nlist=formsRHS[alist][aprimelist]
            a,primelist=extractTMPrimelist(aprimelist)
            if (a != 1):
                afacs=set(prime_divisors(a))
                for aprimelist2 in sorted(formsRHS[alist]):
                    if aprimelist != aprimelist2:
                        a2,primelist2=extractTMPrimelist(aprimelist2)
                        if ((a2 == 1) and (afacs <= set(primelist2)) and
                            (set(primelist) <= set(primelist2))):
                            if not (set(Nlist) <=
                                    set(formsRHS[alist][aprimelist2])):
                                formsRHS[alist][aprimelist2]=(
                                    formsRHS[alist][aprimelist2]+Nlist)
                                formsRHS[alist][aprimelist2].sort()
                            del formsRHS[alist][aprimelist]
                            break

    OutFile=open(OF,"w")
    # Output all data to the file OF in the format "Nlist,alist,a,primelist".
    for alist in sorted(formsRHS):
        for aprimelist in sorted(formsRHS[alist]):
            Nlist=str(formsRHS[alist][aprimelist]).replace(" ","")
            out=Nlist+","+alist+","+aprimelist
            OutFile.write("%s\n" % out)
    OutFile.close()
    os.rename(OF,IF)

if __name__ == '__main__':
    # Map command line arguments to function arguments.
    args=sys.argv
    if args[2] == "TM":
        gatherTMFormRedundancy(args[1])
    elif args[2] == "XYZ2":
        gatherXYZ2Redundancy(args[1])
