# Creation of Vectors and Matrices

There are a few function calls you can use to create a new matrix or vector. These functions require an allocator to create the new variable. To create a new vector you can use the <b>crtXvec</b> function as shown on line 3 in the subsequent example.

<!-- //"XMTX: crtXvec test" -->
<pre>
01 const alloc: std.mem.Allocator = std.testing.allocator;
02 const v1: []f32 = try crtXmtx(9, 3, &alloc, true);
03 const v2: []f32 = try crtXvec(9, &alloc);
04 try std.testing.expectEqual(true, isIdtXmtx(v1, 3));
05 try std.testing.expectEqual(false, isIdtXmtx(v2, 3));
06 alloc.free(v1);
07 alloc.free(v2);
</pre>

To create a new matrix you can use the long form <b>crtXmtx</b> function which takes the total length of the new matrix, the number of columns in the new matrix, an allocator, and a Boolean flag indicating if the new matrix should be an identity matrix as arguments.

<!-- //"XMTX: crtXmtx test" -->
<pre>
01 const alloc: std.mem.Allocator = std.testing.allocator;
02 const l: usize = 9;
03 const cols: usize = 3;
04 const m1: []f32 = try crtXmtx(l, cols, &alloc, false);
05 const m2: []f32 = try crtXmtx(l, cols, &alloc, true);
06 try std.testing.expectEqual(true, isSqrXmtx(m1, cols));
07 try std.testing.expectEqual(true, isIdtXmtx(m2, cols));
08 alloc.free(m1);
09 alloc.free(m2);
</pre>

A slightly more compact version of the matrix creation function that works only with square matrices and doesn't support creating an identity matrix.

<!-- //"XMTX: crtXmtxEz test" -->
<pre>
01 const alloc: std.mem.Allocator = std.testing.allocator;
02 const l: usize = 9;
03 const cols: usize = 3;
04 const m1: []f32 = try crtXmtxEz(l, &alloc);
05 const m2: []f32 = try crtXmtxEz(l, &alloc);
06 try std.testing.expectEqual(true, isSqrXmtx(m1, cols));
07 try std.testing.expectEqual(false, isIdtXmtx(m2, cols));
08 alloc.free(m1);
09 alloc.free(m2);
</pre>