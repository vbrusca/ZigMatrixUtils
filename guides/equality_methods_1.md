# Equality Methods of Vectors and Matrices

To test two arbitrarily vectors for equality you can use the <b>equXvec</b> function as shown on line 4 below.

<!-- "XMTX: equXvec test" -->
<pre>
01 var v1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
02 var v2: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
03 prntXvec(&v1);
04 try std.testing.expectEqual(true, equXvec(&v1, &v2));
</pre>

Similarly for matrices of an arbitrary size you can use the <b>equXmtx</b> function, shown below on line 4.

<!-- "XMTX: equXmtx test" -->
<pre>
01 var m1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
02 var m2: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
03 prntXmtx(&m1, 3);
04 try std.testing.expectEqual(true, equXmtx(&m1, &m2));
</pre>



