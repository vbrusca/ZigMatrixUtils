# Inverse of Matrices

To find the inverse of a matrix we need to reduce it to row eschelon form
while tracking the matrix row manipulations on an identity matrix. This can all be handled with one call to the <b>rdcXmtx</b> function.
<br>
<br>
Arguments 5 and 6 of the function call are a Boolean value, that indicates if an identity matrix has been provided, and the identity matrix for the given matrix being reduced. 

<!-- //"XMTX: ELA - Larson, Edwards: 2.3 Example 2 test" -->
<pre>
01 //Find the inverse of A and verify it is correct.
02 //A = 1  4
03 //   -1  3
04 //
05 //A^-1 = -3 -4
06 //        1  1
07 //
08 var A: [4]f32 = .{ 1, 4, -1, -3 };
09 var I: [4]f32 = .{ 1, 0, 0, 1 };
10 var B: [4]f32 = .{ 0, 0, 0, 0 };
11 var exp: [4]f32 = .{ -3, -4, 1, 1 };
12 var res: bool = false;
13 
14 var sclr: f32 = 0.0;
15 res = rdcXmtx(&A, 2, false, &B, true, &I, 2, false, &sclr);
16 
17 prntXmtx(&B, 2);
18 prntNl();
19 
20 //0: x: 1.0e+00 y: 0.0e+00 
21 //1: x: 0.0e+00 y: 1.0e+00 
22
23 prntXmtx(&I, 2);
24 prntNl();
25 
26 //0: x: -3.0e+00 y: -4.0e+00 
27 //1: x: 1.0e+00 y: 1.0e+00 
28 
29 const b1: bool = equXmtx(&exp, &I);
</pre>