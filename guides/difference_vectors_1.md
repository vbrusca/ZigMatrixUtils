# Difference of Vectors

To find the difference between two vectors use the <b>diff1Xvec</b> function
which takes the difference of <b>m2 - m1</b> and stores the result in <b>m1</b>.

<!-- //"XMTX: ELA - Larson, Edwards: 2.1 Example 3 test" -->
<pre>
01 var m1: [9]f32 = .{ 3, 6, 12, -9, 0, -3, 6, 3, 6 };
02 var m2: [9]f32 = .{ 2, 0, 0, 1, -4, 3, -1, 3, 2 };
03 diff1Xvec(&m1, &m2);
04 var exp1: [9]f32 = .{ 1, 6, 12, -10, 4, -6, 7, 0, 4 };
05 clnXmtx(&m1);
06 const b: bool = equXmtx(&exp1, &m1);
</pre>

Notice that we clean the resulting matrix, <b>m1</b>, before using it as an argument in the <b>equXmtx</b> function. The second function for finding a difference of vectors requires 3 vectors to work with and will store the differences from all three vectors at the given index in the first argument, line 5.

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