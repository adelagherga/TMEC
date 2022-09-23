These files generate all elliptic curves of conductor 500000-500999 by solving
a large list of corresponding Thue and Thue--Mahler forms. To run this program
using GNU parallel, running at most 20 jobs in parallel, run

$ cd LMFDB
$ chmod u+x Code/computeEllipticCurves.sh
$ nohup Code/computeEllipticCurves.sh &
$ rm nohup

This program proceeds as follows:
1. Solve all relevant Thue equations and output elliptic curves, in parallel.
   At most 20 jobs are run in parallel, and GNU parallel's progress is
   recorded in the logfile Data/ThueTest1Log.
   a. Each line, "set", of Data/Forms/ThueTestForms.csv is taken as input
      to the magma program Code/computeEllipticCurvesThue.m. This magma
      program solves the Thue equations listed in set, outputs a logfile in
      Data/ThueLogfiles, and if there are relevant elliptic curves, outputs
      these in Data/ThueOutfiles.
      Note that multiple lines of Data/Forms/ThueTestForms.csv may correspond
      to the same conductor, and thus elliptic curves of the same conductor
      may be spread across multiple such outfiles.
2. Solve all relevant Thue--Mahler equations and output elliptic curves, in
   parallel. At most 20 jobs are run in parallel, and GNU parallel's progress
   is recorded in the logfile Data/TMTest1Log.
   a. Each line, "set", of Data/Forms/TMTestForms.csv is taken as input
      to the magma program Code/computeEllipticCurvesTM.m. This magma
      program solves the Thue--Mahler equations listed in set, outputs a
      logfile in Data/TMLogfiles, and if there are relevant elliptic curves,
      outputs these in Data/TMOutfiles.
      Note that multiple lines of Data/Forms/TMTestForms.csv may correspond
      to the same conductor, and thus elliptic curves of the same conductor
      may be spread across multiple such outfiles.
3. For each conductor N in the range 500000-500999, generate a file
   Data/EllipticCurves/N.csv.
4. Iterate through all files in Data/TMOutfiles and, for each file containing
   elliptic curves of conductor N, append the output of these files to
   Data/EllipticCurves/N.csv. Repeat with Data/ThueOutfiles.
