# Determinant of Matrices

To calculate the determinant of a diagonal matrix use the <b>detDiagXmtx</b>
function with a reference to the matrix and its column size as shown below on lines 5 and 6.

<!-- //"XMTX: ELA - Larson, Edwards: 3.1 Example 6 test" -->
<pre>
01 var A: [16]f32 = .{ 2, 0, 0, 0, 4, -2, 0, 0, -5, 6, 1, 0, 1, 5, 3, 3 };
02 var B: [25]f32 = .{ -1, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, -2 };
03 var detA: f32 = 0;
04 var detB: f32 = 0;
05 detA = detDiagXmtx(&A, 4);
06 detB = detDiagXmtx(&B, 5);
07 std.debug.print("detA: {} detB: {}\n", .{ detA, detB });
08 try std.testing.expectEqual(true, isEquF32(detA, -12.0, true));
09 try std.testing.expectEqual(true, isEquF32(detB, 48.0, true));
</pre>

To calculate the determinant of a triangular matrix use the <b>detTriangXmtx</b>
function with a reference to the matrix and its column size as shown below on lines 5 and 6.

<!-- //"XMTX: detTriangXmtx test" -->
<pre>
01 var A: [9]f32 = .{ 2, 0, 0, 0, 2, 0, 0, 0, 2 };
02 var B: [16]f32 = .{ 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 3 };
03 var detA: f32 = 0;
04 var detB: f32 = 0;
05 detA = detTriangXmtx(&A, 3);
06 detB = detTriangXmtx(&B, 4);
07 std.debug.print("detA: {} detB: {}\n", .{ detA, detB });
08 try std.testing.expectEqual(true, isEquF32(detA, 8.0, true));
09 try std.testing.expectEqual(true, isEquF32(detB, 81.0, true));
</pre>

For a given matrix of a specific size you can use a more direct calculation method. To find the determinant of a 1x1 matrix you can use the <b>detXmtx1</b> function which takes a reference to the matrix as its only argument.

<!-- //"XMTX: detXmtx1 test" -->
<pre>
01 var m1: [1]f32 = .{2};
02 const detM1: f32 = detXmtx1(&m1);
03 try std.testing.expectEqual(m1[0], detM1);
</pre>

To find the determinant of a 2x2 matrx you can use the <b>detXmtx2</b> function shown below on line 2.

<!-- //"XMTX: ELA - Larson, Edwards: 3.1 Example 1, 2, 3, 4, 5 test" -->
<pre>
01 var A: [4]f32 = .{ 2, -3, 1, 2 };
02 var detA: f32 = detXmtx2(&A);
03 try std.testing.expectEqual(true, isEquF32(detA, 7.0, true));
</pe>

To find the determinant of a 3x3 matrix use the <b>detXmtx3</b> function.

<!-- //"XMTX: detXmtx3 test" -->
<pre>
01 var A: [9]f32 = .{ 6, 1, 1, 4, -2, 5, 2, 8, 7 };
02 const exp: f32 = -306;
03 const detA: f32 = detXmtx3(&A);
04 try std.testing.expectEqual(exp, detA);
</pre>

For a 4x4 matrix use the <b>detXmtx4</b> function as shown below on line 3.

<!-- //"XMTX: detXmtx4 test" -->
<pre>
01 var A: [16]f32 = .{ 4, 3, 2, 2, 0, 1, -3, 3, 0, -1, 3, 3, 0, 3, 1, 1 };
02 const expA: f32 = -240.0;
03 const detA: f32 = detXmtx4(&A);
04 try std.testing.expectEqual(true, isEquF32(expA, detA, true));
</pre>

To find the determinant of an arbitrarily sized matrix you can pass a reference to it along with its column count and a new memory allocator to the <b>detXmtx</b> function. The last argument to the method is for specifying the row to use for the determinant calculation.

<!-- //"XMTX: ELA - Larson, Edwards: 3.1 Example 1, 2, 3, 4, 5 test" Example 3 -->
<pre>
01 const alloc: std.mem.Allocator = std.testing.allocator;
02 var C: [9]f32 = .{ 0, 2, 1, 3, -1, 2, 4, 0, 1 };
03 const cols: usize = 3;
04 const detC: f32 = try detXmtx(&C, cols, &alloc, 0);
05 try std.testing.expectEqual(true, isEquF32(detC, 14.0, true));
</pre>

Note the use of the <b>isEquF32</b> function to determine equality of the value in <b>detC</b>.

