## Multiplication of Vectors

There are a few different ways you can "multiply" two vectors. The functions are listed as follows.

<ul>
    <li>crsPrdXvec3</li>
    <li>dotPrdXvec3</li>
    <li>mulXvec</li>
    <li>tmsXvec</li>
</ul>

To find the corss product of two vectors you can use the <b>crsPrdXvec3</b> function.

<!-- //"XMTX: crsPrdXvec3 test" -->
<pre>
01 var v1: [3]f32 = .{ 1, 0, 0 };
02 var v2: [3]f32 = .{ 0, 1, 0 };
03 var v3: [3]f32 = .{ 0, 0, 1 };
04 var v4: [3]f32 = .{ 0, 0, 0 };
05 var v5: [3]f32 = .{ 0, 0, 0 };
06 
07 v4 = crsPrdXvec3(&v1, &v3);
08 v5 = crsPrdXvec3(&v1, &v2);
09 var exp1: [3]f32 = .{ 0, -1, 0 };
10 var exp2: [3]f32 = .{ 0, 0, 1 };
11 
12 prntXvecNl(&v4);
13 //x: 0.0e+00 y: -1.0e+00 z: 0.0e+00
14 
15 prntXvecNl(&v5);
16 //x: 0.0e+00 y: 0.0e+00 z: 1.0e+00
17 
18 try std.testing.expectEqual(true, equXvec(&exp1, &v4));
19 try std.testing.expectEqual(true, equXvec(&exp2, &v5));
20 prntNl();
</pre>
 
To find the dot product of two vectors you can use the <b>dotPrdXvec3</b> function.

<!-- //"XMTX: dotPrdXvec3 test" -->
<pre>
01 var v1: [3]f32 = .{ 1, 0, 0 };
02 var v2: [3]f32 = .{ 0, 1, 0 };
03 var v3: [3]f32 = .{ 0, 3, 0 };
04 var v4: f32 = -1.0;
05 
06 v4 = dotPrdXvec3(&v1, &v2);
07 try std.testing.expectEqual(true, (v4 == 0));
08 
09 v4 = dotPrdXvec3(&v2, &v3);
10 try std.testing.expectEqual(true, (v4 > 0));
11 prntNl();
</pre>

To multiply a vector by a scalar value use the <b>mulXvec</b> function.

<!-- //"XMTX: mulXvec test" -->
<pre>
01 var v1: [3]f32 = .{ 1, 1, 1 };
02 var v2: [3]f32 = .{ 21, 21, 21 };
03 const a: f32 = 21;
04 
05 mulXvec(&v1, a);
06 prntXvecNl(&v1);
07 //x: 2.1e+01 y: 2.1e+01 z: 2.1e+01
08 
09 try std.testing.expectEqual(true, equXvec(&v1, &v2));
10 prntNl();
</pre>

To multiply two vectors together without altering them use the <b>tmsXvec</b> function. The function performs vector multiplication on the provided vectors storing the result in the <b>ret</b> vector.

<!-- //"XMTX: tmsXvec test" -->
<pre>
01 var v1: [3]f32 = .{ 3, 3, 3 };
02 var v2: [3]f32 = .{ 2, 2, 2 };
03 var v3: [1]f32 = .{0};
04 var exp: [1]f32 = .{18};
05 
06 _ = tmsXvec(&v1, &v2, &v3);
07 
08 prntXvecNl(&v1);
09 prntXvecNl(&v2);
10 prntXvecNl(&v3);
11 try std.testing.expectEqual(true, equXvec(&exp, &v3));
12 prntNl();
</pre>