
# Copy of Vectors

To copy a vector use the <b>cpyXvec</b> function shown on line
4 below. It simply takes a reference to the source vector and the destination vector.

<!-- //"XMTX: ELA - Larson, Edwards: 2.1 Problem 1 test" -->
<pre>
01 var a: [4]f32 = .{ 1, -1, 2, -1 };
02 var aPb: [4]f32 = .{ 0, 0, 0, 0 };
03 var exp: [4]f32 = .{ 2, -2, 4, -2 };
04 cpyXvec(&a, &aPb);
05 mulXvec(&aPb, 2);
06 const b1: bool = equXmtx(&aPb, &exp);
</pre>

If you want to copy a vector into a new, undeclared, vector you can use the <b>cpyXvecNew</b> function as shown on line 3 below. Notice that this function requires an allocator to create the new vector to store the copied data.

<!-- //"XMTX: cpyXvecNew test" -->
<pre>
01 const alloc: std.mem.Allocator = std.testing.allocator;
02 var v1: [3]f32 = .{ 1, 2, 3 };
03 const v2: []f32 = try cpyXvecNew(&v1, &alloc);
04 try std.testing.expectEqual(true, equXvec(&v1, v2));
05 alloc.free(v2);
</pre>