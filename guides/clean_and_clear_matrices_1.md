# Clean and Clear Matrices

To clear a matrix, zero out all matrix values, you can use the <b>clrXmtx</b> function as shown below on line 5.

<!-- //"XMTX: clrXmtx test" -->
<pre>
01 var mtx: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
02 var exp: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
03 prntXmtx(&mtx, 3);
04 prntXmtx(&mtx, 3);
05 clrXmtx(&mtx);
06 try std.testing.expectEqual(true, equXvec(&mtx, &exp));
</pre>

When cleaning a matrix to normalize the floating point value of matrix entries use the <b>clnXmtx</b> function, line 18, of the following example.

<!-- //"XMTX: hasInvXmtx test" -->
<pre>
01 var m1: [9]f32 = .{ 3, 2, -3, 4, -3, 6, 1, 0, -1 };
02 var m2: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
03 var idtM1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
04 var ret: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
05 var origM1: [9]f32 = .{ 3, 2, -3, 4, -3, 6, 1, 0, -1 };
06 const cols: f32 = 3;
07 var b: bool = false;
08 
09 std.debug.print("hasInvXmtx test: initial matrix\n", .{});
10 prntXmtx(&m1, 4);
11 
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
</pre>

Note that the <b>clrXmtx</b> and <b>clnXmtx</b> function work on the provided matrix and make changes directly to the matrix.