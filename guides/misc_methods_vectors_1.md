# Miscellaneous Methods for Vectors

There are a few vector based miscellaneous function that you can use for various calculations. The functions are as follows.

<ul>
    <li>aglBtwnXvec</li>
    <li>chgXvecBasis</li>
    <li>getBasisHndXvec3</li>
    <li>getBasisHndXvec3Ret</li>
    <li>magXvec</li>
    <li>nrmXvec</li>
    <li>perpXvec_VecP_To_VecQ</li>
    <li>projXvec_VecP_Onto_VecQ</li>
</ul>

To find the angle between to vectors use the <b>aglBtwnXvec</b> function and pass it a pointer to each vector. The function returns an angle in radians.

<!-- //"XMTX: aglBtwnXvec test" -->
<pre>
01 var v1: [3]f32 = .{ 1, 0, 0 };
02 var v2: [3]f32 = .{ 0, 0, 1 };
03 var angle: f32 = -1.0;
04 angle = aglBtwnXvec(&v1, &v2);
05 
06 const exp: f32 = 1.57079625e+00;
07 prntXvecNl(&v1);
08 //x: 1.0e+00 y: 0.0e+00 z: 0.0e+00 
09 
10 std.debug.print("Angle RAD: {}\n", .{angle});
11 //Angle RAD: 1.57079625e+00
12 
13 std.debug.print("Angle DEG: {}\n", .{rad2Deg(angle)});
14 //Angle DEG: 9.0e+01
15 
16 try std.testing.expectEqual(exp, angle);
17 prntNl();
</pre>

To change the basis of a given vector you can use the <b>chgXvecBasis</b> function. This function is responsible for changing the basis of the provided vector, vec, from the basis, <b>basis</b>, to <b>chgBasis</b>.
</br>
</br>
The arguments for the function are are follows.

<ul>
<li>vec = The vector to convert from one basis to another.</li>
<li>basis = The matrix that constitutes the current basis for vector vec, and holds the resulting transformation matrix.</li>
<li>cols = The number of columns in the matrix, basis.</li>
<li>isStd = A Boolean value indicating if the basis is considered the standard basis.</li>
<li>chgBasis = The mtrix representing the basis to change to.</li>
<li>chgCols = The number of columns in the chgCols matrix.</li>
<li>chgIsStd = A Boolean value indicating if the change basis is considered the standard basis.</li>
<li>idtMtx = The identity associated with the basis and chgBasis matrices.</li>
<li>ret = The result, left-hand side, of the tranformation matrix calculation, should result in an identity matrix.</li>
<li>verbose = A Boolean indicating if the function is verbose.</li>
<li>returns = A Boolean value indicating if the operation was a success or not.</li>
</ul>

An example of the function in use follows, with function call on line 25.

<!-- //"XMTX: chgXvecBasis test" -->
<pre>
01 //B = {(1,0), (1,2)}
02 //Bp = {(1,0), (0,1)}
03 
04 //B = | 1  1|
05 //    | 0  2|
06 //vec = [3 2]
07 
08 //Bp = | 1  0|
09 //     | 0  1|
10 //nvec = [? ?]
11 
12 var B: [4]f32 = .{ 1, 1, 0, 2 };
13 var vec: [2]f32 = .{ 3, 2 };
14 var exp: [2]f32 = .{ 5, 4 };
15 var nvec: [2]f32 = .{ 0, 0 };
16 
17 const cols: usize = 2;
18 var Bp: [4]f32 = .{ 1, 0, 0, 1 };
19 const colsp: usize = 2;
20 var idtMtx: [4]f32 = .{ 1, 0, 0, 1 };
21 
22 var b: bool = false;
23 const vbose: bool = true;
24 var ret: [4]f32 = .{ 0, 0, 0, 0 };
25 b = chgXvecBasis(&vec, &B, cols, false, &Bp, colsp, true, &idtMtx, &ret, &nvec, vbose);
26 
27 try std.testing.expectEqual(true, b);
28 std.debug.print("\nchgXvecBasis ret {}", .{b});
29 
30 //chgXvecBasis ret true
31 //0: x: 0.0e+00 y: 0.0e+00 
32 //1: x: 0.0e+00 y: 0.0e+00 
33 
34 prntXmtxNl(&ret, cols);
35 prntNl();
36 
37 std.debug.print("\nchgXvecBasis B", .{});
38 prntXmtxNl(&B, cols);
39 prntNl();
40 
41 //chgXvecBasis B
42 //0: x: 1.0e+00 y: 1.0e+00 
43 //1: x: 0.0e+00 y: 2.0e+00 
44 
45 std.debug.print("\nchgXvecBasis nvec", .{});
46 prntXmtxNl(&nvec, 1);
47 prntNl();
48 
49 //chgXvecBasis nvec
50 //0: x: 5.0e+00 
51 //1: x: 4.0e+00 
52 
53 try std.testing.expectEqual(true, equXmtx(&exp, &nvec));
54 prntNl();
</pre>

To find the handedness of a basis of 3, 1x3, vectors you can use the <b>getBasisHndXvec3</b> or the <b>getBasisHndXvec3Ret</b> function. The function takes 3 arguments that are ferences to the vectors that form a basis. The <b>getBasisHndXvec3Ret</b> takes one more argument, a reference to the variables that holds the function's return value.

<!-- //"XMTX: getBasisHndXvec3 test" -->
<pre>
01 var v1: [3]f32 = .{ 1, 0, 0 };
02 var v2: [3]f32 = .{ 0, 1, 0 };
03 var v3: [3]f32 = .{ 0, 0, 1 };
04 var expResHnd: BASIS_HAND = BASIS_HAND.RIGHT;
05 var resHnd: BASIS_HAND = BASIS_HAND.ERROR_ZERO;
06 
07 resHnd = getBasisHndXvec3(&v1, &v2, &v3);
08 try std.testing.expectEqual(expResHnd, resHnd);
09 v1 = .{ 1, 0, 0 };
10 v2 = .{ 0, 1, 0 };
11 v3 = .{ 0, 0, -1 };
12 expResHnd = BASIS_HAND.LEFT;
13 
14 resHnd = getBasisHndXvec3(&v1, &v2, &v3);
15 try std.testing.expectEqual(expResHnd, resHnd);
16 v1 = .{ 0, 0, 0 };
17 v2 = .{ 0, 0, 0 };
18 v3 = .{ 0, 0, -1 };
19 expResHnd = BASIS_HAND.ERROR_ZERO;
20 
21 resHnd = getBasisHndXvec3(&v1, &v2, &v3);
22 try std.testing.expectEqual(expResHnd, resHnd);
23 prntNl();
</pre>

To find the magnitude of a vector you can use the <b>magXvec</b> function which takes a pointer to a vector as an argument and returns its magnitude.

<!-- //"XMTX: magXvec test" -->
<pre>
01 var v1: [3]f32 = .{ 1, 2, 3 };
02 const exp: f32 = 3.74165749e+00;
03 const val: f32 = magXvec(&v1);
04 try std.testing.expectEqual(exp, val);
05 prntNl();
</pre>

To normalize a vector you can use the <b>nrmXvec</b> function which takes a vector pointer as an argument.

<!-- //"XMTX: nrmXvec test" -->
<pre>
01 var v1: [3]f32 = .{ 5, 10, 15 };
02 nrmXvec(&v1);
03 
04 prntXvecNl(&v1);
05 //x: 2.67261236e-01 y: 5.34522473e-01 z: 8.01783740e-01
06 
07 try std.testing.expectEqual(@as(f32, 2.67261236e-01), v1[0]);
08 try std.testing.expectEqual(@as(f32, 5.34522473e-01), v1[1]);
09 try std.testing.expectEqual(@as(f32, 8.01783740e-01), v1[2]);
10 prntNl();
</pre>

The last two functions go hand-in-hand. One is used to find the perpendicular vector to a pair of vectors <b>vecP</b> and <b>vecQ</b>, the <b>perpXvec_VecP_To_VecQ</b> function.

<!-- //"XMTX: perpXvec_VecP_To_VecQ test" -->
<pre>
01 //Q = 2 2 1
02 //P = 1 -2 0
03 var Q: [3]f32 = .{ 2, 2, 1 };
04 var P: [3]f32 = .{ 1, -2, 0 };
05 var exp: [3]f32 = .{ 1.44444441e+00, -1.55555558e+00, 2.22222223e-01 };
06 
07 const ret: []f32 = perpXvec_VecP_To_VecQ(&P, &Q);
08 clnXmtx(ret);
09 clnXmtx(&exp);
10 
11 prntXvecNl(ret);
12 prntNl();
13 //x: 1.44444441e+00 y: -1.55555558e+00 z: 2.22222223e-01 
14 
15 prntXvecNl(&exp);
16 prntNl();
17 //x: 1.44444441e+00 y: -1.55555558e+00 z: 2.22222223e-01 
18 
19 try std.testing.expectEqual(true, equXmtx(ret, &exp));
20 prntNl();
</pre>

The next fuction, <b>projXvec_VecP_Onto_VecQ</b>, finds a vector projected from <b>vecP</b> onto <b>vecQ</b>.

<!-- //"XMTX: projXvec_VecP_Onto_VecQ test" -->
<pre>
01 //Q = 2 2 1
02 //P = 1 -2 0
03 var Q: [3]f32 = .{ 2, 2, 1 };
04 var P: [3]f32 = .{ 1, -2, 0 };
05 var exp: [3]f32 = .{ (-4.0 / 9.0), (-4.0 / 9.0), (-2.0 / 9.0) };
06 
07 const ret: []f32 = projXvec_VecP_Onto_VecQ(&P, &Q);
08 clnXmtx(ret);
09 clnXmtx(&exp);
10 
11 prntXvecNl(ret);
12 prntNl();
13 //x: -4.44444447e-01 y: -4.44444447e-01 z: -2.22222223e-01
14 
15 prntXvecNl(&exp);
16 prntNl();
17 //x: -4.44444447e-01 y: -4.44444447e-01 z: -2.22222223e-01
18 
19 try std.testing.expectEqual(true, equXmtx(ret, &exp));
20 prntNl();
</pre>