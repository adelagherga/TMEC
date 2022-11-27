#!/usr/bin/python
# compareLMFDB.py

# These functions compare the computed elliptic curve file
# "ID/EllipticCurves/AllCurves.csv" with the LMFDB, for those curves whose
# conductor is less than 500000. If there are elliptic curves in the LMFDB which
# are not in the directory ID, the missing curves are printed to the output file
# "ID/EllipticCurves/MissingCurves.csv".

# Parameters
#     ID: <class 'str'>
#         The string denoting the directory, in the format Data/${name}.
# Returns
#     OF: <class 'str'>
#         The string denoting the output file, wherein each line is in the
#         format "N aInvs" to indicate the elliptic curve that is present in the
#         LMFDB, but not in the directory, ID.
# Authors
#     Adela Gherga <adelagherga@gmail.com>
# Created
#     27 November 2022

import sys
import os
from urllib import request

def extractAInv(strX):
    """
      Given a string in the format "[x_1,...,x_n]", return the tuple of integers
      (x_1,...,x_n).

      Parameters
          strX: <class 'str'>
              A string in the format "[x_1,...,x_n]".
      Returns
          tupX: <class 'tuple'>
              A tuple of integers (x_1,...,x_n).
    """
    strX=strX.lstrip("[").rstrip("]")
    tupX=tuple([int(char) for char in strX.split(",")]):
    return tupX

def determineNlist(ID):
    """
      Extract the list of conductors from the directory name ID.

      Parameters
          ID: <class 'str'>
              The string denoting the directory, in the format Data/${name}.
      Returns
          Nlist: <class 'set'>
              The list of conductors pertaining to the directory.
    """
    name=ID.split("/")[1].split("]")[0][1:]
    if ".." in name:
        Ns=[int(name.split("..")[0]),int(name.split("..")[1])]
        Nlist={N for N in range(Ns[0],Ns[1]+1)}
    else:
        Nlist={int(N) for N in name.split(",")}
    return Nlist

def compareLMFDB(ID):
    """
      Compare the computed elliptic curve file "ID/EllipticCurves/AllCurves.csv"
      with the LMFDB, for those curves whose conductor is less than 500000. If
      there are elliptic curves in the LMFDB which are not in the directory,
      print the missing curves to the output file
      "ID/EllipticCurves/MissingCurves.csv".

      Parameters
          ID: <class 'str'>
              The string denoting the directory, in the format Data/${name}.
      Returns
          OF: <class 'str'>
              The string denoting the output file, wherein each line is in the
              format "N aInvs" to indicate the elliptic curve that is present in
              the LMFDB, but not in the directory, ID.
    """
    Nlist=determineNlist(ID)
    if Nlist.intersection(range(0,500000)) == set():
        # The conductor list is outside the range of the LMFDB.
        exit()

    IF=ID+"/EllipticCurves/AllCurves.csv"
    OF=ID+"/EllipticCurves/MissingCurves.csv"
    tempF=ID+"/tmp"

    rangeLMFDB=[]
    for i in range(0,50):
        inter=Nlist.intersection(range(i*10000,(i+1)*10000))
        if inter != set():
            # Determine the LMFDB file to retrieve, as well as the list of
            # conductors expected to appear in that filename.
            rangeLMFDB.append([i,inter])

    url0="https://raw.githubusercontent.com/JohnCremona/"
    url0+="ecdata/master/allcurves/allcurves."
    LMFDB={}
    for iinter in rangeLMFDB:
        i=iinter[0]
        if i == 0:
            Nrange="00000-09999"
        else:
            Nrange=str(i*10000)+"-"+str(((i+1)*10000)-1)
        url=url0+Nrange
        LMFDBfile=request.urlretrieve(url,tempF)
        # Download the necessary LMFDB file into the temporary file tempF.
        for line in open(tempF):
            line=line.rstrip()
            N=int(line.split(" ")[0])
            if (N in iinter[1]) and (N in LMFDB):
                LMFDB[N].add(extractAInv(line.split(" ")[3]))
            elif (N in iinter[1]):
                LMFDB[N]={extractAInv(line.split(" ")[3])}
        os.remove(tempF)

    TMEC={}
    for line in open(IF):
        line=line.rstrip()
        N=int(line.split(" ")[0])
        if N in TMEC:
            TMEC[N].add(extractAInv(line.split(" ")[1]))
        else:
            TMEC[N]={extractAInv(line.split(" ")[1])}

    OutFile=open(OF,"w")
    if TMEC != LMFDB:
        # Output any missing curves to the file OF in the format "N aInvs"
        for N in sorted(LMFDB):
            if N not in sorted(TMEC):
                for ainv in LMFDB[N]:
                    out=str(N)+" "+str(ainv).replace(" ","")
                    OutFile.write("%s\n" %out)
            else:
                for ainv in LMFDB[N]:
                    if ainv not in TMEC[N]:
                        out=str(N)+" "+str(ainv).replace(" ","")
                        OutFile.write("%s\n" %out)
    OutFile.close()
    exit()

if __name__ == '__main__':
    # Map command line arguments to function arguments.
    args=sys.argv
    compareLMFDB(args[1])
