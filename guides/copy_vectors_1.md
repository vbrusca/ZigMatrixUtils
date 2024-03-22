
# Copy of Vectors

To copy a vector use the <b>cpyXvec</b> function shown on line
4 below.

<!-- //"XMTX: ELA - Larson, Edwards: 2.1 Problem 1 test" -->
<pre>
01 var a: [4]f32 = .{ 1, -1, 2, -1 };
02 var aPb: [4]f32 = .{ 0, 0, 0, 0 };
03 var exp: [4]f32 = .{ 2, -2, 4, -2 };
04 cpyXvec(&a, &aPb);
05 mulXvec(&aPb, 2);
06 const b1: bool = equXmtx(&aPb, &exp);
</pre>

Note that you can use the method with larger matrices it will simply view them as large vectors.

<!-- //"XMTX: cpyXvecNew test" -->
<pre>
const alloc: std.mem.Allocator = std.testing.allocator;
var v1: [3]f32 = .{ 1, 2, 3 };
const v2: []f32 = try cpyXvecNew(&v1, &alloc);
try std.testing.expectEqual(true, equXvec(&v1, v2));
alloc.free(v2);
</pre>