# Division of Vectors

To divide a vector by a scalar value use the <b>divXvec</b> function as shown below on line 4.

<!-- //"XMTX: divXvec test" -->
<pre>
01 var v1: [3]f32 = .{ 1, 1, 1 };
02 var v2: [3]f32 = .{ 4.76190485e-02, 4.76190485e-02, 4.76190485e-02 };
03 const a: f32 = 21;
04 divXvec(&v1, a);
05 prntXvec(&v1);
06 try std.testing.expectEqual(true, equXvec(&v1, &v2));
</pre>