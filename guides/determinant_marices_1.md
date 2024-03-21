# Determinant of Marices 1

To find the determinant of a matrix it is faster to use the methods of known matrix sizes, <b>detXmtx1</b>, <b>detXmtx2</b>, <b>detXmtx3</b>, <b>detXmtx14</b>. Simply pass in a reference to the given 2x2 matrix to the <b>defXmtx2</b> method and it will return the determinant of the matrix.

<!-- //"XMTX: ELA - Larson, Edwards: 3.1 Example 1, 2, 3, 4, 5 test" -->
<pre>
01 var A: [4]f32 = .{ 2, -3, 1, 2 };
02 var detA: f32 = detXmtx2(&A);
03 try std.testing.expectEqual(true, isEquF32(detA, 7.0, true));
</pe>