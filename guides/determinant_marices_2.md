# Determinant of Marices 2

To find the determinant of an arbitrarily sized matrix you can pass a reference to it along with its column count and a new memory allocator to the <b>detXmtx</b> function. The last argument to the method is for specifying the row to use for the determinant calculation.

<!-- "XMTX: ELA - Larson, Edwards: 3.1 Example 1, 2, 3, 4, 5 test" Example 3 -->
<pre>
01 const alloc: std.mem.Allocator = std.testing.allocator;
02 var C: [9]f32 = .{ 0, 2, 1, 3, -1, 2, 4, 0, 1 };
03 const cols: usize = 3;
04 const detC: f32 = try detXmtx(&C, cols, &alloc, 0);
05 try std.testing.expectEqual(true, isEquF32(detC, 14.0, true));
</pre>

Note the use of the isEquF32 function to determine equality of the value in detC.