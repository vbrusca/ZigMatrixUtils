# Copy of Matrices

To copy a segment of a larger matrix you can use the <b>cpyLessXmtx</b> function as shown below on line 4. The first argument is a reference to the matrix to copy, the second argument is a reference to the new matrix to be copied into. While the next two arguments are the column count of the source matrix followed by the number of columns to copy.

<!-- //"XMTX: cpyLessXmtx test" -->
<pre>
01 var m1: [4]f32 = .{ 1, 2, 3, 4 };
02 var m2: [3]f32 = .{ 0, 0, 0 };
03 var m3: [3]f32 = .{ 1, 2, 3 };
04 cpyLessXmtx(&m1, &m2, 4, 3);
05 try std.testing.expectEqual(true, equXvec(&m2, &m3));
</pre>

To copy the entire contents of a matrix the <b>cpyXmtx</b> function provides the necessary support. An example of its use is shown below on lines 4 and 5.

<!-- //"XMTX: cpyXmtx test" -->
<pre>
01 var m1: [3]f32 = .{ 1, 2, 3 };
02 var m2: [3]f32 = .{ 0, 0, 0 };
03 var m3: [3]f32 = .{ 1, 1, 1 };
04 cpyXmtx(&m1, &m2);
05 cpyXmtx(&m1, &m3);
06 try std.testing.expectEqual(true, equXvec(&m1, &m2));
07 try std.testing.expectEqual(true, equXvec(&m1, &m3));
</pre>

If you want to copy a matrix into a new, undeclared, matrix you can use the <b>cpyXmtxNew</b> function as shown on line 3 below. Notice that this function requires an allocator to create the new matrix to store the copied data.

<!-- //"XMTX: cpyXmtxNew test" -->
<pre>
01 const alloc: std.mem.Allocator = std.testing.allocator;
02 var m1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
03 const m2: []f32 = try cpyXmtxNew(&m1, &alloc);
04 try std.testing.expectEqual(true, equXvec(&m1, m2));
05 alloc.free(m2);
</pre>