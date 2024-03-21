# Scalar Multiplication of Matrices

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