# Cross Product of Vectors

In order to calculate the cross product of a given set of vectors you can use the <b>crsPrdXvec</b> function which takes a reference to the two vectors as arguments.

<!-- //"XMTX: crsPrdXvec test" -->
<pre>
01 var v1: [3]f32 = .{ 1, 0, 0 };
02 var v2: [3]f32 = .{ 0, 1, 0 };
03 var v3: [3]f32 = .{ 0, 0, 1 };
04 var v4: [3]f32 = crsPrdXvec(&v1, &v3);
05 var v5: [3]f32 = crsPrdXvec(&v1, &v2);
06 var exp1: [3]f32 = .{ 0, -1, 0 };
07 var exp2: [3]f32 = .{ 0, 0, 1 };
08 prntXvec(&v4);
09 prntXvec(&v5);
10 try std.testing.expectEqual(true, equXvec(&exp1, &v4));
11 try std.testing.expectEqual(true, equXvec(&exp2, &v5));
</pre>