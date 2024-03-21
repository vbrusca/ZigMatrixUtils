# Determinant of Marices 3

To calculate the determinant of a diagonal matrix use the <b>detDiagXmtx</b>
function with a reference to the matrix and its column ize.

<!-- "XMTX: ELA - Larson, Edwards: 3.1 Example 6 test" -->
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