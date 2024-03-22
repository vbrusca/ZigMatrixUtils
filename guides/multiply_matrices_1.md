# Multiplication of Matrices

Multiplying to matrices together using matrix multiplication can be
accomplished using the <b>tmsXmtx</b> function which takes two input matrices, with column counts, as arguments along with an output matrix, with column count, as the last two arguments to the function call.

<!-- //"XMTX: ELA - Larson, Edwards: 2.1 Example 4, E5 test" -->
<pre>
//Example 1
01 var m3: [4]f32 = .{ 3, 4, -2, 5 };
02 var m4: [4]f32 = .{ 1, 0, 0, 1 };
03 var ret2: [4]f32 = .{ 0, 0, 0, 0 };
04 var exp2: [4]f32 = .{ 3, 4, -2, 5 };
05 const b1: bool = tmsXmtx(&m3, 2, &m4, 2, &ret2, 2);
06 const b2: bool = equXmtx(&exp2, &ret2);
</pre>

In the second example, listed subsequently, we have in-line print outs
of the matrix values.

<!-- //XMTX: ELA - Larson, Edwards: 2.1 Example 4, E5 test -->
<pre>
//Example 2
01 var m5: [3]f32 = .{ 1, -2, -3 };
02 var m6: [3]f32 = .{ 2, -1, 1 };
03 var ret3: [1]f32 = .{0};
04 var exp3: [1]f32 = .{1};
05 b = tmsXmtx(&m5, 3, &m6, 1, &ret3, 1);
06 
07 prntXmtx(&m5, 3);
08 prntNl();
09 
10 //0: x: 1.0e+00 y: -2.0e+00 z: -3.0e+00
11 
12 prntXmtx(&m6, 1);
13 prntNl();
14 
15 //0: x: 2.0e+00 
16 //1: x: -1.0e+00
17 //2: x: 1.0e+00
18 
19 prntXmtx(&ret3, 1);
20 prntNl();
21 
22 //0: x: 1.0e+00
</pre>

Notice the use of the column count on lines 7, 12, and 19 to adjust how a matrix is drawn. In some instances as a 1 row 3 column matrix, in others a 3 row 1 column matrix, and still others with just one row and one column.
<br>
<br>
In order to multiply two vectors by a scalar value we use the following 
example code. Note that if the provided argument is a matrix the function will process it as a large vector. So you can use this function with matrix arguments.

<!-- //"XMTX: ELA - Larson, Edwards: 2.1 Example 3 test" -->
<pre>
01 var m1: [9]f32 = .{ 1, 2, 4, -3, 0, -1, 2, 1, 2 };
02 var exp1: [9]f32 = .{ 3, 6, 12, -9, 0, -3, 6, 3, 6 };
03 mulXvec(&m1, 3);
04 const b: bool = equXmtx(&exp1, &m1);
</pre>

The results of the function call are stored in the <b>m1</b> argument. Tests for equality are performed with the <b>equXmtx</b> function.