# Clean and Clear Matrices

To clear a matrix to have only zero values use the following function.

<!-- "XMTX: clrXmtx test" -->
<pre>
01 var mtx: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
02 var exp: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
03 prntXmtx(&mtx, 3);
04 prntXmtx(&mtx, 3);
05 clrXmtx(&mtx);
06 try std.testing.expectEqual(true, equXvec(&mtx, &exp));
</pre>

When cleaning a matrix to normalize the floating point value of matrix entries use the method, line 18, of the following example.

<!-- "XMTX: hasInvXmtx test" -->
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
12 var sclr: f32 = 0.0;
13 b = rdcXmtxInl(&m1, 3, false, true, &idtM1, 3, false, &sclr);
14 try std.testing.expectEqual(true, b);
15 cpyXmtx(&m1, &m2);
16 cpyXmtx(&origM1, &m1);
17 
18 clnXmtx(&m2);
19 try std.testing.expectEqual(true, isDiagXmtx(&m2, 3));
20 try std.testing.expectEqual(true, isIdtXmtx(&m2, 3));
</pre>
