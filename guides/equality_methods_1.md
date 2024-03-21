# Equality Methods

<!-- "XMTX: equXmtx test" -->
<pre>
std.debug.print("equXmtx test:\n", .{});
var m1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
var m2: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
prntXmtx(&m1, 3);
try std.testing.expectEqual(true, equXmtx(&m1, &m2));
</pre>

<!-- "XMTX: equXvec test" -->
<pre>
std.debug.print("equXvec test:\n", .{});
var v1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
var v2: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
prntXvec(&v1);
try std.testing.expectEqual(true, equXvec(&v1, &v2));
</pre>

