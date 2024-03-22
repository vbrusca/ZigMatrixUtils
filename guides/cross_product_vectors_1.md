# Cross Product of Vectors

In order to calculate the cross product of a given set of vectors you can use the <b>crsPrdXvec</b> function which takes a reference to the two vectors as arguments and returns a vector with the calculated cross-product.

<!-- //"XMTX: crsPrdXvec test" -->
<pre>
01 var v1: [3]f32 = .{ 1, 0, 0 };
02 var v2: [3]f32 = .{ 0, 1, 0 };
03 var v3: [3]f32 = .{ 0, 0, 1 };
04 var v4: [3]f32 = crsPrdXvec(&v1, &v3);
05 var v5: [3]f32 = crsPrdXvec(&v1, &v2);
06 var exp1: [3]f32 = .{ 0, -1, 0 };
07 var exp2: [3]f32 = .{ 0, 0, 1 };
08
09 prntXvec(&v4);
10 //x: 0.0e+00 y: -1.0e+00 z: 0.0e+00 
11 
12 prntXvec(&v5);
13 //x: 0.0e+00 y: 0.0e+00 z: 1.0e+00
14 
15 try std.testing.expectEqual(true, equXvec(&exp1, &v4));
16 try std.testing.expectEqual(true, equXvec(&exp2, &v5));
</pre>