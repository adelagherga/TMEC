#!/usr/bin/python
# gatherFormRedundancy.py

# These functions remove redundant Thue--Mahler equations across conductors.
# That is, for any 2 lines of the file IF, "N_1,alist_1,a_1,primelist_1" and
# "N_2,alist_2,a_2,primelist_2", where alist_1 = alist_2, a_1 = a_2, and
# primelist_1 is a subset of primelist_2, writes
# "[N_1,N_2],alist_1,a_1,primelist_1" in the file OF.

# Parameters
#     IF: <class 'str'>
#         The string denoting the input file, wherein each line is in the
#         format "N,alist,a,primelist".
#     OF: <class 'str'>
#         A string denoting an empty output file.
# Returns
#     OF: <class 'str'>
#         The string denoting the output file, wherein each line is in the
#         "Nlist,alist,a,primelist".
# Authors
#     Adela Gherga <adelagherga@gmail.com>
# Created
#     28 September 2022

import sys

def parceForm(line):
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

def extractPrimelist(line):
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

def gatherFormRedundancy(IF,OF):
    """
      Removes redundant Thue--Mahler equations across conductors. That is, for
      any 2 lines of the file IF, "N_1,alist_1,a_1,primelist_1" and
      "N_2,alist_2,a_2,primelist_2", where alist_1 = alist_2, a_1 = a_2,
      and primelist_1 is a subset of primelist_2, writes
      "[N_1,N_2],alist_1,a_1,primelist_1" in the file OF.

      Parameters
          IF: <class 'str'>
              The string denoting the input file, wherein each line is in the
              format "N,alist,a,primelist".
          OF: <class 'str'>
              A string denoting an empty output file.
      Returns
          OF: <class 'str'>
              The string denoting the output file, wherein each line is in the
              "Nlist,alist,a,primelist".
    """
    forms={}
    for line in open(IF):
        # Sort data by alist.
        N,alist,aprimelist=parceForm(line)
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
            # For a given form with 2 rhs possibilities, (Nlist1,a1,primelist1) and
            # (Nlist2,a2,primelist2), if a1 = a2 and primelist1 is contained in
            # primelist2, collapse both rhs possibilities into
            # (Nlist1+Nlist2,a1,primelist2).
            Nlist=formsRHS[alist][aprimelist]
            a,primelist=extractPrimelist(aprimelist)
            for aprimelist2 in sorted(formsRHS[alist]):
                if aprimelist != aprimelist2:
                    a2,primelist2=extractPrimelist(aprimelist2)
                    if (a == a2) and (set(primelist) <= set(primelist2)):
                        formsRHS[alist][aprimelist2]=(formsRHS[alist][aprimelist2]+Nlist)
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

if __name__ == '__main__':
    # Map command line arguments to function arguments.
    #gatherFormRedundancy(*sys.argv[1:])
    args = sys.argv
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])
