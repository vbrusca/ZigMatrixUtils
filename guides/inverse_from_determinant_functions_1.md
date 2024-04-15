# Inverse from Matrix Determinants

To calculate the inverse of a 2x2 matrix from its determinant you can use the <b>getInvFromDet2</b> function which takes a reference to the original matrix, the determinant of the original matrix, and an empty matrix to store the results as argument, line 5 in the following example.

<!-- "XMTX: getInvFromDet2 test" -->
<pre>
01 var m1: [4]f32 = .{ 2, 1, 7, 4 };
02 var invM1: [4]f32 = .{ 4, -1, -7, 2 };
03 var res: [4]f32 = std.mem.zeroes([4]f32); //.{};
04 const detM1: f32 = detXmtx2(&m1);
05 const b: bool = getInvFromDet2(&m1, detM1, &res);
06 try std.testing.expectEqual(true, b);
07 try std.testing.expectEqual(true, equXmtx(&res, &invM1));
</pre>

To calculate the inverse of a 3x3 matrix from its determinant you can use the <b>getInvFromDet3</b> function as shown on line 17 in the subsequent example.

<!-- "XMTX: getInvFromDet3 test" -->
<pre>
01 //A = 1 2 3
02 //    4 5 6
03 //    7 2 9
04 
05 //detA = (1)(33)          - (2)(-6)          + (3)(-27)
06 //detA = (1)(det 5 6 2 9) - (2)(det 4 6 7 9) + (3)(det 4 5 7 2)
07 //detA = -36;
08 
09 //A^-1 = -11/12  1/3  1/12
10 //       -1/6    1/3  -1/6
11 //       3/4     -1/3 1/12
12 
13 var A: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 2, 9 };
14 var res: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
15 const detA: f32 = -36;
16 var invA: [9]f32 = .{ (-11.0 / 12.0), (1.0 / 3.0), (1.0 / 12.0), (-1.0 / 6.0), (1.0 / 3.0), (-1.0 / 6.0), (3.0 / 4.0), (-1.0 / 3.0), (1.0 / 12.0) };
17 const b: bool = getInvFromDet3(&A, detA, &res);
18 try std.testing.expectEqual(true, b);
19 try std.testing.expectEqual(true, equXmtx(&res, &invA));
</pre>

Lastly, to calculate the inverse of a 4x4 matrix from its determinant you can use the <b>getInvFromDet4</b> function as shown on line 6 in the following example.

<!-- "XMTX: getInvFromDet4 test" -->
<pre>
01 var ret: [16]f32 = std.mem.zeroes([16]f32); //.{};
02 var mtx: [16]f32 = .{ 1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1 };
03 var expMtx: [16]f32 = .{ (1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0), (-1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0), (-1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0), (-1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0), (-1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0) };
04 const det: f32 = detXmtx4(&mtx);
05 const expDet: f32 = -16;
06 const b: bool = getInvFromDet4(&mtx, det, &ret);
07 const cols: usize = 4;
08 
09 std.debug.print("getInvFromDet4 det: {} expDet: {}\n", .{ det, expDet });
10 prntXmtx(&ret, cols);
11 prntNl();
12 
13 //getInvFromDet4 det: -1.6e+01 expDet: -1.6e+01
14 //0: x: 2.5e-01 y: 2.5e-01 z: 2.5e-01 w: -2.5e-01 
15 //1: x: 2.5e-01 y: 2.5e-01 z: -2.5e-01 w: 2.5e-01 
16 //2: x: 2.5e-01 y: -2.5e-01 z: 2.5e-01 w: 2.5e-01 
17 //3: x: -2.5e-01 y: 2.5e-01 z: 2.5e-01 w: 2.5e-01 
18 
19 std.debug.print("getInvFromDet4 det:\n", .{});
20 prntXmtx(&expMtx, cols);
21 prntNl();
22 
23 //getInvFromDet4 det:
24 //0: x: 2.5e-01 y: 2.5e-01 z: 2.5e-01 w: -2.5e-01 
25 //1: x: 2.5e-01 y: 2.5e-01 z: -2.5e-01 w: 2.5e-01 
26 //2: x: 2.5e-01 y: -2.5e-01 z: 2.5e-01 w: 2.5e-01 
27 //3: x: -2.5e-01 y: 2.5e-01 z: 2.5e-01 w: 2.5e-01 
28
29 try std.testing.expectEqual(true, b);
30 try std.testing.expectEqual(true, isEquF32(expDet, det, true));
31 try std.testing.expectEqual(true, equXmtx(&expMtx, &ret));
</pre>
