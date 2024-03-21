# Difference of Vectors 2

The second method for finding a difference of vectors requires 3 vectors to work with and will store the differences from all three vectors at the given index in the first argument, line 12.

<!-- //XMTX: ELA - Larson, Edwards: 2.1 Problem 1 test -->
<pre>
01 var a: [4]f32 = .{ 1, -1, 2, -1 };
02 var b: [4]f32 = .{ 2, -1, -1, 8 };
03 var aPb: [4]f32 = .{ 0, 0, 0, 0 };
04 var exp: [4]f32 = .{ -1, 0, 3, -9 };
05 diff2Xvec(&aPb, &a, &b);
06
07 prntXmtx(&aPb, 2);
08 prntNl();
09
10 //0: x: -1.0e+00 y: 0.0e+00
11 //1: x: 3.0e+00 y: -9.0e+00
12 
13 prntXmtx(&exp, 2);
14 prntNl();
15 
16 //0: x: -1.0e+00 y: 0.0e+00
17 //1: x: 3.0e+00 y: -9.0e+00
18 
19 const b: bool = equXmtx(&aPb, &exp);
</pre>