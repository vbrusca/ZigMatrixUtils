# Miscellaneous Functions for Matrices

TODO

There are a number of miscellaneous matrix functions that you can use for various calculations. The functions are as follows.

<ul>
    <li>Adjoint Functions</li>
    <li>adjXmtx3</li>
    <li>adjXmtx4</li>
</ul>

<ul>
    <li>Cofactor Functions</li>
    <li>cofXmtx</li>
    <li>cofXmtx3</li>
    <li>cofXmtx4</li>
</ul>

<ul>
    <li>Cofactor Sign Functions</li>    
    <li>cofXmtxSign</li>
    <li>cofXmtxSign2, cofXmtxSign2Ret</li>
    <li>cofXmtxSign3, cofXmtxSign3Ret</li>
    <li>cofXmtxSign4, cofXmtxSign4Ret</li>
</ul>

<ul>
    <li>Basis Functions</li>
    <li>getBasisCnvXmtx</li>
    <li>getBasisHndXmtx3, getBasisHndXmtx3Ret</li>
</ul>

<ul>
    <li>Other Functions</li>
    <li>fndLgstRowAbsXmtx</li>
    <li>getCramerSupportXmtx</li>
</ul>

<ul>
    <li>Matrix Processing</li>
    <li>idnfXmtx</li>
    <li>procXmtx</li>
</ul>

<ul>
    <li>Minor Functions</li>
    <li>mnrXmtx2, mnrXmtx2Ret</li>
    <li>mnrXmtx3, mnrXmtx3Ret</li>
    <li>mnrXmtx4, mnrXmtx4Ret</li>
</ul>

<ul>
    <li>Remove Row/Column Functions</li>    
    <li>rmvRowColXmtx3</li>
    <li>rmvRowColXmtx4</li>
</ul>

<!--
NEEDED:
adjXmtx
adjXmtx1
adjXmtx2
cofXmtx1
cofXmtx2
cofXmtxSign1
cofXmtxSign1Ret
mnrXmtx1
mnrXmtx1Ret
rmvRowColXmtx2
-->

### Adjoint Functions for Matrices

There are a few functions that help with getting the adjoint matrix of a given matrix. Use <b>adjXmtx3</b> for 3x3 matrices and <b>adjXmtx4</b> for 4x4 matrices.

<!-- //"XMTX: adjXmtx3 test" -->
<pre>
01 var A: [9]f32 = .{ -1, 3, 2, 0, -2, 1, 1, 0, -2 };
02 var retA: [9]f32 = std.mem.zeroes([9]f32);
03 retA = adjXmtx3(&A);
04 var expA: [9]f32 = .{ 4, 6, 7, 1, 0, 1, 2, 3, 2 };
05 try std.testing.expectEqual(true, equXmtx(&expA, &retA));
06 prntNl();
</pre>

### Cofactor Functions for Matrices

To generate a cofactor matrix you can use the generix <b>cofXmtx</b> function or the more specific <b>cofXmtx3</b>, <b>cofXmtx4</b> functions. The example below shows the <b>cofXmtx</b> function in action on lines 24 and 33.

<!-- //"XMTX: cofXmtx test" -->
<pre>
01 //A = 1 2 3
02 //    4 5 6
03 //    7 2 9
04 
05 //detA = (1)(33)          - (2)(-6)          + (3)(-27)
06 //detA = (1)(det 5 6 2 9) - (2)(det 4 6 7 9) + (3)(det 4 5 7 2)
07 //detA = -36;
08 
09 //cofA = 33   6    -27
10 //       -12  -12  12
11 //       -3   6    -3
12 
13 //A^-1 = -11/12  1/3  1/12
14 //       -1/6    1/3  -1/6
15 //       3/4     -1/3 1/12
16 
17 var A: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 2, 9 };
18 const cofA: [9]f32 = .{ 33, 6, -27, -12, -12, 12, -3, 6, -3 };
19 var res: [4]f32 = .{ 0, 0, 0, 0 };
20 const resCols: f32 = 2;
21 const cols: f32 = 3;
22 var b: bool = false;
23 
24 b = cofXmtx(&A, cols, 0, 0, &res, resCols);
25 prntXmtxNl(&res, resCols);
26 prntNl();
27 //0: x: 5.0e+00 y: 6.0e+00 
28 //1: x: 2.0e+00 y: 9.0e+00 
29 
30 try std.testing.expectEqual(true, b);
31 try std.testing.expectEqual(cofA[0], (detXmtx2(&res) * cofXmtxSign(0, 0, true)));
32 
33 b = cofXmtx(&A, cols, 0, 1, &res, resCols);
34 prntXmtxNl(&res, resCols);
35 prntNl();
36 //0: x: 4.0e+00 y: 6.0e+00 
37 //1: x: 7.0e+00 y: 9.0e+00 
38 
39 try std.testing.expectEqual(true, b);
40 try std.testing.expectEqual(cofA[1], (detXmtx2(&res) * cofXmtxSign(0, 1, true)));
41 prntNl();
</pre>

A demonstration of using the direct, 3x3 matrix version of the function, <b>cofXmtx3</b> is as follows, line 7.

<!-- //"XMTX: cofXmtx3 test" -->
<pre>
01 //A = 0  2  1
02 //    3  -1 2
03 //    4  0  1
04 var A: [9]f32 = .{ 0, 2, 1, 3, -1, 2, 4, 0, 1 };
05 var exp: [9]f32 = .{ -1, 5, 4, -2, -4, 8, 5, 3, -6 };
06 var cof: [9]f32 = undefined;
07 cof = cofXmtx3(&A);
08 
09 std.debug.print("\nAAA cofXmtx3 test", .{});
10 prntXmtxNl(&cof, 3);
11 //AAA cofXmtx3 test
12 //0: x: -1.0e+00 y: 5.0e+00 z: 4.0e+00 
13 //1: x: -2.0e+00 y: -4.0e+00 z: 8.0e+00 
14 //2: x: 5.0e+00 y: 3.0e+00 z: -6.0e+00 
15 
16 std.debug.print("\nBBB cofXmtx3 test", .{});
17 prntXmtxNl(&exp, 3);
18 //BBB cofXmtx3 test
19 //0: x: -1.0e+00 y: 5.0e+00 z: 4.0e+00 
20 //1: x: -2.0e+00 y: -4.0e+00 z: 8.0e+00 
21 //2: x: 5.0e+00 y: 3.0e+00 z: -6.0e+00 
22 
23 const b: bool = equXmtx(&cof, &exp);
24 std.debug.print("\nCCC cofXmtx3: {}", .{b});
25 //CCC cofXmtx3: true
26
27 try std.testing.expectEqual(true, b);
28 prntNl();
</pre>

### Cofactor Sign Functions for Matrices

<!-- //"XMTX: cofXmtxSign test" -->
<pre>
01 var c: usize = 0;
02 var r: usize = 0;
03 var ret: f32 = -1.0;
04 
05 ret = cofXmtxSign(r, c, true);
06 try std.testing.expectEqual(@as(f32, 1.0), ret);
07 c = 1;
08 r = 0;
09 
10 ret = cofXmtxSign(r, c, true);
11 try std.testing.expectEqual(@as(f32, -1.0), ret);
12 prntNl();
</pre>

<!-- //"XMTX: cofXmtxSign3 test" -->
<pre>
01 var cofSgn: [9]f32 = undefined;
02 cofSgn = cofXmtxSign3();
03 
04 try std.testing.expectEqual(@as(f32, 1.0), cofSgn[0]);
05 try std.testing.expectEqual(@as(f32, -1.0), cofSgn[1]);
06 try std.testing.expectEqual(@as(f32, 1.0), cofSgn[2]);
07 
08 try std.testing.expectEqual(@as(f32, -1.0), cofSgn[3]);
09 try std.testing.expectEqual(@as(f32, 1.0), cofSgn[4]);
10 try std.testing.expectEqual(@as(f32, -1.0), cofSgn[5]);
11 
12 try std.testing.expectEqual(@as(f32, 1.0), cofSgn[6]);
13 try std.testing.expectEqual(@as(f32, -1.0), cofSgn[7]);
14 try std.testing.expectEqual(@as(f32, 1.0), cofSgn[8]);
15 prntNl();
</pre>

## Vector Basis Functions

To get the basis conversion, tranformation, matrix for the given current basis and change of basis use the <b>getBasisCnvXmtx</b> function. The fnction has a few arguments listed as follows, and returns a boolean indicating if the fnction was successful or not.

<ul>
    <li>basis = The matrix that constitutes the current basis and holds the resulting transformation matrix.</li>
    <li>cols = The number of columns in the basis matrix.</li>
    <li>chgBasis = The matrix that constitues the new basis to find coordinates in.</li>
    <li>ret = The result, left-hand side, of the tranformation matrix calculation, should result in an identity matrix.</li>
    <li>verbose = A Boolean indicating if the function is verbose.</li>
</ul>

An example of this function in use is as follows, line 19.

<!-- //"XMTX: getBasisCnvXmtx test" -->
<pre>
01 //B = {(1,0,0), (0, 1, 0), (0, 0, 1)};
02 //B' = {(1,0,1), (0, -1, 2), (2, 3, -5)};
03 
04 //B = | 1  0  0|
05 //    | 0  1  0|
06 //    | 0  0  1|
07 const vbose: bool = false;
08 var ret: [9]f32 = std.mem.zeroes([9]f32); //.{};
09 var B: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
10 const cols: usize = 3;
11
12 //B' = | 1  0  2|
13 //     | 0 -1  3|
14 //     | 1  2 -5|
15 var Bp: [9]f32 = .{ 1, 0, 2, 0, -1, 3, 1, 2, -5 };
16 const colsp: usize = 3;
17 var b: bool = false;
18 
19 b = getBasisCnvXmtx(&B, cols, &Bp, colsp, &ret, vbose);
20 try std.testing.expectEqual(true, b);
21 
22 std.debug.print("\ngetBasisCnvMtx test:getBasisCnvMtx ret {}", .{b});
23 prntXmtxNl(&ret, cols);
24 prntNl();
25 //getBasisCnvMtx test:getBasisCnvMtx ret true
26 //0: x: 1.0e+00 y: 0.0e+00 z: 0.0e+00 
27 //1: x: -0.0e+00 y: 1.0e+00 z: 0.0e+00 
28 //2: x: -0.0e+00 y: -0.0e+00 z: 1.0e+00 
29 
30 std.debug.print("\ngetBasisCnvMtx test:getBasisCnvMtx B", .{});
31 prntXmtxNl(&B, cols);
32 prntNl();
33 //getBasisCnvMtx test:getBasisCnvMtx B
34 //0: x: -1.0e+00 y: 4.0e+00 z: 2.0e+00 
35 //1: x: 3.0e+00 y: -7.0e+00 z: -3.0e+00 
36 //2: x: 1.0e+00 y: -2.0e+00 z: -1.0e+00 
37
38 var exp: [9]f32 = .{ -1, 4, 2, 3, -7, -3, 1, -2, -1 };
39 try std.testing.expectEqual(true, equXmtx(&exp, &B));
40 prntNl();
</pre>

To get the handedness of a vector basis use the <b>getBasisHndXmtx3</b> or the <b>getBasisHndXmtx3Ret</b> function.

<!-- //"XMTX: getBasisHndXmtx3 test" -->
<pre>
01 var m1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
02 var expResHnd: BASIS_HAND = BASIS_HAND.RIGHT;
03 var resHnd: BASIS_HAND = BASIS_HAND.ERROR_ZERO;
04 
05 resHnd = getBasisHndXmtx3(&m1, 3);
06 try std.testing.expectEqual(expResHnd, resHnd);
07 m1 = .{ 1, 0, 0, 0, 1, 0, 0, 0, -1 };
08 expResHnd = BASIS_HAND.LEFT;
09 resHnd = BASIS_HAND.ERROR_ZERO;
10 
11 resHnd = getBasisHndXmtx3(&m1, 3);
12 try std.testing.expectEqual(expResHnd, resHnd);
13 m1 = .{ 0, 0, 0, 0, 0, 0, 0, 0, -1 };
14 expResHnd = BASIS_HAND.ERROR_ZERO;
15 resHnd = BASIS_HAND.LEFT;
16 
17 resHnd = getBasisHndXmtx3(&m1, 3);
18 try std.testing.expectEqual(expResHnd, resHnd);
19 prntNl();
</pre>

### Other Functions

The next method, <b>fndLgstRowAbsXmtx</b>, is a utility function that finds the largest row, by absolute value, in the given matrix, at the specified column, starting on the given row. The function takes the following arguments and returns the index of the row found or an ollegal row index, the row count.

<ul>
    <li>mtx = The matrix to search through.</li>
    <li>cols = The number of columns in the matrix.</li>
    <li>targetCol = The target column to look for values.</li>
    <li>startingRow = The row to start the search on.</li>
</ul>

An example of the method is listed here.

<!-- //"XMTX: fndLgstRowAbsXmtx test" -->
<pre>
01 var mtx: [12]f32 = .{ 1, -3, 21, -10, 3, 6, 7, -13, 5, 0, 0, 0 };
02 const exp: [9]usize = .{ 1, 1, 2, 2, 0, 1, 12, 12, 12 };
03 
04 const v0: usize = fndLgstRowAbsXmtx(&mtx, 3, 0, 0);
05 if (v0 >= mtx.len) {
06     std.debug.print("\nNo matching row found, v0", .{});
07 }
08 
09 const v1: usize = fndLgstRowAbsXmtx(&mtx, 3, 0, 1);
10 if (v1 >= mtx.len) {
11     std.debug.print("\nNo matching row found, v1", .{});
12 }
13 
14 const v2: usize = fndLgstRowAbsXmtx(&mtx, 3, 1, 0);
15 if (v2 >= mtx.len) {
16     std.debug.print("\nNo matching row found, v2", .{});
17 }
18 
19 const v3: usize = fndLgstRowAbsXmtx(&mtx, 3, 1, 1);
20 if (v3 >= mtx.len) {
21     std.debug.print("\nNo matching row found, v3", .{});
22 }
23 
24 const v4: usize = fndLgstRowAbsXmtx(&mtx, 3, 2, 0);
25 if (v4 >= mtx.len) {
26     std.debug.print("\nNo matching row found, v4", .{});
27 }
28 
29 const v5: usize = fndLgstRowAbsXmtx(&mtx, 3, 2, 1);
30 if (v5 >= mtx.len) {
31    std.debug.print("\nNo matching row found, v5", .{});
32 }
33 
34 const v6: usize = fndLgstRowAbsXmtx(&mtx, 3, 3, 0);
35 if (v6 >= mtx.len) {
36     std.debug.print("\nNo matching row found, v6", .{});
37 }
38 
39 const v7: usize = fndLgstRowAbsXmtx(&mtx, 3, 3, 1);
40 if (v7 >= mtx.len) {
41     std.debug.print("\nNo matching row found, v7", .{});
42 }
43 
44 const v8: usize = fndLgstRowAbsXmtx(&mtx, 3, 3, 2);
45 if (v8 >= mtx.len) {
46     std.debug.print("\nNo matching row found, v8", .{});
47 }
48 
49 std.debug.print("\nv0: {} v1: {} v2: {} v3: {} v4: {} v5: {} v6: {} v7: {} v8: {}", .{ v0, v1, v2, v3, v4, v5, v6, v7, v8 });
50 try std.testing.expectEqual(exp[0], v0);
51 try std.testing.expectEqual(exp[1], v1);
52 try std.testing.expectEqual(exp[2], v2);
53 try std.testing.expectEqual(exp[3], v3);
54 try std.testing.expectEqual(exp[4], v4);
55 try std.testing.expectEqual(exp[5], v5);
56 try std.testing.expectEqual(exp[6], v6);
57 try std.testing.expectEqual(exp[7], v7);
58 try std.testing.expectEqual(exp[8], v8);
59 prntNl();
</pre>

### Cramer Suppot Matrix

To create a Cramer’s Rule support matrix based on the provided augmented matrix, mtx, and other arguments, listed here, use the <b>getCramerSupportXmtx</b> function.

<ul>
    <li>mtx = The augmented matrix to create a support matrix for.</li>
    <li>cols = The number of columns in the mtx matrix.</li>
    <li>srcCol = The last column in the mtx matrix that contains constant values.</li>
    <li>dstCol = The target column in the return matrix, ret, to store data.</li>
    <li>ret = The return matrix that has it’s dstCol, destination column, replaced with the mtx matrix’s constant values.</li>
</ul>

An example of this function in use is shown here.

<!-- //"XMTX: getCramerSupportMtx test" -->
<pre>
01 //A = -1  2 -3  1
02 //     2  0  1  0
03 //     3 -4  4  2
04 var A: [12]f32 = .{ -1, 2, -3, 1, 2, 0, 1, 0, 3, -4, 4, 2 };
05 var A1: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
06 var b: bool = false;
07 const cols: usize = 4;
08 const srcCol: usize = 3;
09 var dstCol: usize = 0;
10 var expA1: [9]f32 = .{ 1, 2, -3, 0, 0, 1, 2, -4, 4 };
11 
12 b = getCramerSupportXmtx(&A, cols, srcCol, dstCol, &A1);
13 std.debug.print("\nMatrix A{}:", .{dstCol});
14 prntXmtxNl(&A1, 3);
15 std.debug.print("\nMatrix Exp A{}:", .{dstCol});
16 prntXmtxNl(&expA1, 3);
17 try std.testing.expectEqual(true, b);
18 try std.testing.expectEqual(true, equXmtx(&A1, &expA1));
19 dstCol = 1;
20 var expA2: [9]f32 = .{ -1, 1, -3, 2, 0, 1, 3, 2, 4 };
21 
22 b = getCramerSupportXmtx(&A, cols, srcCol, dstCol, &A1);
23 std.debug.print("\nMatrix A{}:", .{dstCol});
24 prntXmtxNl(&A1, 3);
25 std.debug.print("\nMatrix Exp A{}:", .{dstCol});
26 prntXmtxNl(&expA1, 3);
27 try std.testing.expectEqual(true, b);
28 try std.testing.expectEqual(true, equXmtx(&A1, &expA2));
29 dstCol = 2;
30 var expA3: [9]f32 = .{ -1, 2, 1, 2, 0, 0, 3, -4, 2 };
31 
32 b = getCramerSupportXmtx(&A, cols, srcCol, dstCol, &A1);
33 std.debug.print("\nMatrix A{}:", .{dstCol});
34 prntXmtxNl(&A1, 3);
35 std.debug.print("\nMatrix Exp A{}:", .{dstCol});
36 prntXmtxNl(&expA1, 3);
37 try std.testing.expectEqual(true, b);
38 try std.testing.expectEqual(true, equXmtx(&A1, &expA3));
39 prntNl();
</pre>