# Adjoint Matrices

There are two functions for finding the adjoint matrix of a given matrix.
For a 3x3 matrix use the <b>adjXmtx3</b> function shown below on line 3.

<!-- //"XMTX: adjXmtx3 test" -->
<pre>
01 var A: [9]f32 = .{ -1, 3, 2, 0, -2, 1, 1, 0, -2 };
02 var retA: [9]f32 = std.mem.zeroes([9]f32);
03 retA = adjXmtx3(&A);
04 var expA: [9]f32 = .{ 4, 6, 7, 1, 0, 1, 2, 3, 2 };
05 try std.testing.expectEqual(true, equXmtx(&expA, &retA));
</pre>

For a 4x4 matrix you can calculate the adjoint matrix using the <b>adjXmtx4</b> function as shown on line 3 below.

<!-- //"XMTX: adjXmtx4 test" -->
<pre>
01 var A: [16]f32 = .{ 1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1 };
02 var retA: [16]f32 = std.mem.zeroes([16]f32);
03 retA = adjXmtx4(&A);
04 var expA: [16]f32 = .{ -4, -4, -4, 4, -4, -4, 4, -4, -4, 4, -4, -4, 4, -4, -4, -4 };
05 try std.testing.expectEqual(true, equXmtx(&expA, &retA));
</pre>

Note that both the <b>adjXmtx3</b> and <b>adjXmtx4</b> functions create a matrix to hold the results and returns it.