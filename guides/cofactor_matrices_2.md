# Cofactor of Matrices 2

Another approach to working with cofactors of matrices is to calculate them in one step by using the <b>cofXmtx</b> function call. In the example shown below two cofactors are calculated by hand and compared to the function output on lines 31 and 41.

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
22 
23 var b: bool = cofXmtx(&A, cols, 0, 0, &res, resCols);
24 prntXmtx(&res, resCols);
25 prntNl();
26 
27 //0: x: 5.0e+00 y: 6.0e+00 
28 //1: x: 2.0e+00 y: 9.0e+00 
29 
30 try std.testing.expectEqual(true, b);
31 try std.testing.expectEqual(cofA[0], (detXmtx2(&res) * cofXmtxSign(0, 0, true)));
32 
33 b = cofXmtx(&A, cols, 0, 1, &res, resCols);
34 prntXmtx(&res, resCols);
35 prntNl();
36 
37 //0: x: 4.0e+00 y: 6.0e+00 
38 //1: x: 7.0e+00 y: 9.0e+00 
39 
40 try std.testing.expectEqual(true, b);
41 try std.testing.expectEqual(cofA[1], (detXmtx2(&res) * cofXmtxSign(0, 1, true)));
</pre>

A more direct way to handle calculating the cofactor matrix is to use a function specific for the matrix size. For a 3x3 matrix you can use <b>cofXmtx3</b> as shown below on line 6.

<!-- //"XMTX: cofXmtx3 test" -->
<pre>
01 //A = 0  2  1
02 //    3  -1 2
03 //    4  0  1
04 var A: [9]f32 = .{ 0, 2, 1, 3, -1, 2, 4, 0, 1 };
05 var exp: [9]f32 = .{ -1, 5, 4, -2, -4, 8, 5, 3, -6 };
06 var cof: [9]f32 = cofXmtx3(&A);
07 const b: bool = equXmtx(&cof, &exp);
08 try std.testing.expectEqual(true, b);
</pre>

Similarly for 4x4 matrices you can use the <b>cofXmtx4</b> function as shown below on line 7.

<!-- //"XMTX: cofXmtx4 test" -->
<pre>
01 //A = 1  2  3  4
02 //    0  1  2  3
03 //    -1 3  5  6
04 //    1  1  1  1
05 var A: [16]f32 = .{ 1, 3, -2, 1, 5, 1, 0, -1, 0, 1, 0, -2, 2, -1, 0, 3 };
06 var exp: [16]f32 = .{ 0, 0, 3, 0, -2, 8, 13, 4, 4, -34, -56, -14, 2, -20, -34, -10 };
07 var cof: [16]f32 = cofXmtx4(&A);
08 clnXmtx(&cof);
09 const b: bool = equXmtx(&cof, &exp);
10 try std.testing.expectEqual(true, b);
</pre>

<!-- //"XMTX: cofXmtxSign test" -->
<pre>
var c: usize = 0;
var r: usize = 0;
try std.testing.expectEqual(@as(f32, 1.0), cofXmtxSign(r, c, true));
c = 1;
r = 0;
try std.testing.expectEqual(@as(f32, -1.0), cofXmtxSign(r, c, true));
</pre>

<!-- //"XMTX: cofXmtx3 test" -->
<pre>
//A = 0  2  1
//    3  -1 2
//    4  0  1
var A: [9]f32 = .{ 0, 2, 1, 3, -1, 2, 4, 0, 1 };
var exp: [9]f32 = .{ -1, 5, 4, -2, -4, 8, 5, 3, -6 };
var cof: [9]f32 = cofXmtx3(&A);
const b: bool = equXmtx(&cof, &exp);
try std.testing.expectEqual(true, b);
</pre>

<!-- //"XMTX: cofXmtxSign4 test" -->
<pre>
01 const cofSgn: [16]f32 = cofXmtxSign4();
02 try std.testing.expectEqual(@as(f32, 1.0), cofSgn[0]);
03 try std.testing.expectEqual(@as(f32, -1.0), cofSgn[1]);
04 try std.testing.expectEqual(@as(f32, 1.0), cofSgn[2]);
05 try std.testing.expectEqual(@as(f32, -1.0), cofSgn[3]);
06 
07 try std.testing.expectEqual(@as(f32, -1.0), cofSgn[4]);
08 try std.testing.expectEqual(@as(f32, 1.0), cofSgn[5]);
09 try std.testing.expectEqual(@as(f32, -1.0), cofSgn[6]);
10 try std.testing.expectEqual(@as(f32, 1.0), cofSgn[7]);
11 
12 try std.testing.expectEqual(@as(f32, 1.0), cofSgn[8]);
13 try std.testing.expectEqual(@as(f32, -1.0), cofSgn[9]);
14 try std.testing.expectEqual(@as(f32, 1.0), cofSgn[10]);
15 try std.testing.expectEqual(@as(f32, -1.0), cofSgn[11]);
16 
17 try std.testing.expectEqual(@as(f32, -1.0), cofSgn[12]);
18 try std.testing.expectEqual(@as(f32, 1.0), cofSgn[13]);
19 try std.testing.expectEqual(@as(f32, -1.0), cofSgn[14]);
20 try std.testing.expectEqual(@as(f32, 1.0), cofSgn[15]);
</pre>