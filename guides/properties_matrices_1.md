# Properties of Matrices

There are a number of functions available for determining certain properties of matrices. The functions are as follows.

<ul>
    <li>hasInvXmtx</li>
    <li>isDiagXmtx</li>
    <li>isIdtXmtx</li>
    <li>isOrthXmtx</li>
    <li>isRdcFrmXmtx</li>
    <li>isRghtHandedXmtx</li>
    <li>isSqrXmtx</li>
    <li>isSymXmtx</li>
</ul>

To check if a matrix has an inverse you can use the <b>hasInvXmtx</b> function, example shown on line 46, which returns a Boolean indicating if the given matrix has an inverse.

<!-- //"XMTX: hasInvXmtx test" -->
<pre>
01 var m1: [9]f32 = .{ 3, 2, -3, 4, -3, 6, 1, 0, -1 };
02 var m2: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
03 var idtM1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
04 var ret: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
05 var origM1: [9]f32 = .{ 3, 2, -3, 4, -3, 6, 1, 0, -1 };
06 
07 const cols: f32 = 3;
08 var b: bool = false;
09 
10 std.debug.print("\nhasInvXmtx test: initial matrix", .{});
11 prntXmtxNl(&m1, 4);
12 //hasInvXmtx test: initial matrix
13 //0: x: 3.0e+00 y: 2.0e+00 z: -3.0e+00 w: 4.0e+00 
14 //1: x: -3.0e+00 y: 6.0e+00 z: 1.0e+00 w: 0.0e+00 
15 
16 var sclr: f32 = 0.0;
17 b = rdcXmtxInl(&m1, 3, false, true, &idtM1, 3, false, &sclr);
18 try std.testing.expectEqual(true, b);
19 cpyXmtx(&m1, &m2);
20 cpyXmtx(&origM1, &m1);
21 
22 clnXmtx(&m2);
23 try std.testing.expectEqual(true, isDiagXmtx(&m2, 3));
24 try std.testing.expectEqual(true, isIdtXmtx(&m2, 3));
25 
26 std.debug.print("\nClean Short Answer:", .{});
27 prntXmtxNl(&m2, 3);
28 prntNl();
26 //Clean Short Answer:
27 //0: x: 1.0e+00 y: 0.0e+00 z: 0.0e+00 
28 //1: x: 0.0e+00 y: 1.0e+00 z: 0.0e+00 
29 //2: x: 0.0e+00 y: 0.0e+00 z: 1.0e+00
30 
31 try std.testing.expectEqual(true, isRdcFrmXmtx(&m2, 3)); 
32 b = tmsXmtx(&m2, cols, &m1, cols, &ret, cols);
33 try std.testing.expectEqual(true, b);
34 try std.testing.expectEqual(true, equXmtx(&m1, &ret));
35 
36 origM1 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 }; 
37 cpyXmtx(&origM1, &m1);
38 var trnM1: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
39 trnXmtx(&m1, cols, &trnM1);
40 
41 std.debug.print("\nOrig M1:", .{});
42 prntXmtxNl(&origM1, 3);
43 prntNl();
41 //Orig M1:
42 //0: x: 1.0e+00 y: 0.0e+00 z: 0.0e+00 
43 //1: x: 0.0e+00 y: 1.0e+00 z: 0.0e+00 
44 //2: x: 0.0e+00 y: 0.0e+00 z: 1.0e+00 
45 
46 std.debug.print("\nTrn M1:", .{});
47 prntXmtxNl(&trnM1, 3);
48 prntNl();
49 //Trn M1:
50 //0: x: 1.0e+00 y: 0.0e+00 z: 0.0e+00 
51 //1: x: 0.0e+00 y: 1.0e+00 z: 0.0e+00 
52 //2: x: 0.0e+00 y: 0.0e+00 z: 1.0e+00 
53 
54 const b1: bool = hasInvXmtx(&origM1, cols, &trnM1);
55 try std.testing.expectEqual(true, b1);
56 prntNl();
</pre>

To determine if you have a diagonal matrix or not you can use the <b>isDiagXmtx</b> function.

<!-- //"XMTX: isDiagXmtx test" -->
<pre>
01 var mtx1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
02 var mtx2: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
03
04 prntXmtxNl(&mtx1, 3);
05 //0: x: 1.0e+00 y: 2.0e+00 z: 3.0e+00 
06 //1: x: 4.0e+00 y: 5.0e+00 z: 6.0e+00 
07 //2: x: 7.0e+00 y: 8.0e+00 z: 9.0e+00 
08 
09 prntXmtxNl(&mtx2, 3);
10 //0: x: 1.0e+00 y: 0.0e+00 z: 0.0e+00 
11 //1: x: 0.0e+00 y: 1.0e+00 z: 0.0e+00 
12 //2: x: 0.0e+00 y: 0.0e+00 z: 1.0e+00
13
14 const b1: bool = isDiagXmtx(&mtx1, 3);
15 const b2: bool = isDiagXmtx(&mtx2, 3);
16 try std.testing.expectEqual(false, b1);
17 try std.testing.expectEqual(true, b2);
18 prntNl();
</pre>

To determine if a matrix is the identity matrix use the <b>isIdtXmtx</b> function with the given matrix as a pointer argument, lines 13 and 14.

<!-- //"XMTX: isIdtXmtx test" -->
<pre>
01 var mtx1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
02 var mtx2: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
03 prntXmtxNl(&mtx1, 3);
04 //0: x: 1.0e+00 y: 2.0e+00 z: 3.0e+00
05 //1: x: 4.0e+00 y: 5.0e+00 z: 6.0e+00
06 //2: x: 7.0e+00 y: 8.0e+00 z: 9.0e+00
07 
08 prntXmtxNl(&mtx2, 3);
09 //0: x: 1.0e+00 y: 0.0e+00 z: 0.0e+00
10 //1: x: 0.0e+00 y: 1.0e+00 z: 0.0e+00
11 //2: x: 0.0e+00 y: 0.0e+00 z: 1.0e+00
12 
13 const b1: bool = isIdtXmtx(&mtx1, 3);
14 const b2: bool = isIdtXmtx(&mtx2, 3);
15 
16 try std.testing.expectEqual(false, b1);
17 try std.testing.expectEqual(true, b2);
18 prntNl();
</pre>

To determine if a pair of vectors is linealy independent you can use the <b>isLinIndXmtx</b> function which returns a Boolean indicating if the vectors of the matrix, when compared in order 2 at a time, are linearly independent. The version of the function shown below will allocate the memory necessary for the vecL and vecR comparison vectors.

<!-- //"XMTX: isLinIndXmtx test" -->
<pre>
01 const alloc: std.mem.Allocator = std.testing.allocator;
02 var mtx: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
03 var b: bool = false;
04 
05 b = isLinIndXmtx(&mtx, 3, &alloc);
06 try std.testing.expectEqual(true, b);
07 prntNl();
</pre>

To determine if a matrix is orthogonal use the <b>isOrthXmtx</b> function. This function returns a Boolean value indicating if the <b>mtx</b> argument is an orthogonal matrix. The function takes the following arguments.

<ul>
    <li>mtx = The matrix to determine othogonality for.</li>
    <li>cols = The number of columns in the matrix mtx.</li>
    <li>ret = An empty matrix the same size as mtx used to hold return information from the reduce matrix call.</li>
    <li>idtMtx = An identity matrix with the same dimensions as the mtx matrix.</li>
    <li>trnMtx = An empty matrix the same size as mtx used to hold transpose matrix information.</li>
</ul>

An example of this function is use is shown subsequently on line 9.

<!-- //"XMTX: isOrthXmtx test" -->
<pre>
01 var mtx: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
02 const cols: usize = 3;
03 var ret: [9]f32 = std.mem.zeroes([9]f32); //.{};
04 var idtMtx: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
05 var trnMtx: [9]f32 = std.mem.zeroes([9]f32); //.{};
06 const expB: bool = true;
07 var b: bool = false;
08 
09 b = isOrthXmtx(&mtx, cols, &ret, &idtMtx, &trnMtx);
10 try std.testing.expectEqual(expB, b);
11 prntNl();
</pre>

To determine if a function is in reduced form or not se the <b>isRdcFrmXmtx</b> with a pointer to the matrix and the matrix column count as arguments, lines 8, 16, and 24 below.

<!-- //"XMTX: isRdcFrmXmtx test" -->
<pre>
01 var mtx: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
02 var b: bool = false;
03 prntXmtxNl(&mtx, 3);
04 //0: x: 1.0e+00 y: 0.0e+00 z: 0.0e+00 
05 //1: x: 0.0e+00 y: 1.0e+00 z: 0.0e+00 
06 //2: x: 0.0e+00 y: 0.0e+00 z: 1.0e+00 
07 
08 b = isRdcFrmXmtx(&mtx, 3);
09 try std.testing.expectEqual(true, b);
10 mtx = .{ 1, 0, 0, 0, 1, 0, 0, 0, 0 };
11 prntXmtxNl(&mtx, 3);
12 //0: x: 1.0e+00 y: 0.0e+00 z: 0.0e+00 
13 //1: x: 0.0e+00 y: 1.0e+00 z: 0.0e+00 
14 //2: x: 0.0e+00 y: 0.0e+00 z: 0.0e+00 
15 
16 b = isRdcFrmXmtx(&mtx, 3);
17 try std.testing.expectEqual(true, b);
18 mtx = .{ 1, 0, 0, 0, 0, 0, 0, 0, 1 };
19 prntXmtxNl(&mtx, 3);
20 //0: x: 1.0e+00 y: 0.0e+00 z: 0.0e+00 
21 //1: x: 0.0e+00 y: 0.0e+00 z: 0.0e+00 
22 //2: x: 0.0e+00 y: 0.0e+00 z: 1.0e+00 
23 
24 b = isRdcFrmXmtx(&mtx, 3);
25 try std.testing.expectEqual(false, b);
26 prntNl();
</pre>

To determine if a matrix is a square matrix or not use the <b>isSqrXmtx</b> function as shown below on lines 14, and 15.

<!-- //"XMTX: isSqrXmtx test" -->
<pre>
01 var mtx1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
02 var mtx2: [12]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1 };
03 prntXmtxNl(&mtx1, 3);
04 //0: x: 1.0e+00 y: 2.0e+00 z: 3.0e+00 
05 //1: x: 4.0e+00 y: 5.0e+00 z: 6.0e+00 
06 //2: x: 7.0e+00 y: 8.0e+00 z: 9.0e+00 
07
08 prntXmtxNl(&mtx2, 3);
09 //0: x: 1.0e+00 y: 0.0e+00 z: 0.0e+00 
10 //1: x: 0.0e+00 y: 1.0e+00 z: 0.0e+00 
11 //2: x: 0.0e+00 y: 0.0e+00 z: 1.0e+00 
12 //3: x: 1.0e+00 y: 1.0e+00 z: 1.0e+00 
13 
14 const b1: bool = isSqrXmtx(&mtx1, 3);
15 const b2: bool = isSqrXmtx(&mtx2, 3);
16 try std.testing.expectEqual(true, b1);
17 try std.testing.expectEqual(false, b2);
18 prntNl();
</pre>

To check if a matrix is symmatrical or not use the <b>isSymXmtx</b> function, line 3 below.

<!-- //"XMTX: isSymXmtx test" -->
<pre>
01 var m1: [4]f32 = .{ 1, 0, 0, 1 };
02 const cols: usize = 2;
03 const b: bool = isSymXmtx(&m1, cols);
04 try std.testing.expectEqual(true, b);
05 prntNl();
</pre>