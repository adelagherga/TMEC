// Below is an example of how to use the Thue--Mahler solver,
// for instance, to solve
// 3 X^3 + 2 X^2 Y + 7 X Y^2 + 2 Y^3 = 2^{z_1} 3^{z_2} 7^{z_3} 41^{z_4}
// under the assumptions that gcd(X,Y) = gcd(a_0,Y) = 1:

load "solveThueMahler.m";

alist:=[3,2,7,2];
a:=1;
primelist:=[2,3,7,41];
time sols:=solveThueMahler(alist,a,primelist);
sols;

// To solve the same Thue--Mahler equation without the assumption that
// gcd(a_0,Y) = 1:

alist:=[3,2,7,2];
a:=1;
primelist:=[2,3,7,41];
time sols:=solveThueMahlerAll(alist,a,primelist : coprime:=false);
sols;
