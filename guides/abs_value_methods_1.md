# Absolute Value Methods

To get the absolue value of an f32 variable you can use the <b>absF32</b> function.

<!-- //"XMTX: absF32 test" -->
<pre>
01 const v: f32 = 7;
02 try std.testing.expectEqual(v, absF32(-7));
</pre>

You can get the absolute value of a vector in a similar way by providing a reference to the given vector to the <b>absXvec</b> function.

<!-- //"XMTX: absXvec test" -->
<pre>
01 std.debug.print("absXvec test:\n", .{});
02 var v1: [9]f32 = .{ -1, -2, -3, 4, 5, 6, -7, -8, -9 };
03 var v2: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
04 absXvec(&v1);
05 try std.testing.expectEqual(true, equXvec(&v1, &v2));
</pre>

Lastly, the same operation can be performed on a matrix and explicitly the absolute value calculated by a call to the <b>absXmtx</b> function.

<!-- //"XMTX: absXmtx test" -->
<pre>
01 std.debug.print("absXmtx test:\n", .{});
02 var v1: [9]f32 = .{ -1, -2, -3, 4, 5, 6, -7, -8, -9 };
03 var v2: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
04 absXmtx(&v1);
05 try std.testing.expectEqual(true, equXmtx(&v1, &v2));
</pre>