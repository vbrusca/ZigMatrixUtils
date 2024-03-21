#Difference of Vectors 2

The second method for finding a difference of vectors requires 3 vectors to work with and will store the differences from all three vectors at the given index in the first argument, line 12.

//XMTX: ELA - Larson, Edwards: 2.1 Problem 1 test
<pre>
a = .{ 1, -1, 2, -1 };
b = .{ 2, -1, -1, 8 };
aPb = .{ 0, 0, 0, 0 };
exp = .{ -1, 0, 3, -9 };
diff2Xvec(&aPb, &a, &b);

prntXmtx(&aPb, 2);
prntNl();

//0: x: -1.0e+00 y: 0.0e+00
//1: x: 3.0e+00 y: -9.0e+00

prntXmtx(&exp, 2);
prntNl();

//0: x: -1.0e+00 y: 0.0e+00
//1: x: 3.0e+00 y: -9.0e+00

const b: bool = equXmtx(&aPb, &exp);
</pre>