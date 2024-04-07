# Manipulation of Vectors

There are a number of function that perform different actions on vectors. The list is as follows.

<ul>
    <li>addXvec</li>
    <li>diff1Xvec</li>
    <li>diff2Xvec</li>    
    <li>divXvec</li>
    <li>subXvec</li>
    <li>sum1Xvec</li>
    <li>sum2Xvec</li>
</ul>

To add a scalar value to the components of a vector you can use the <b>addXvec</b> function. It takes a slice, []f32, representing a vector and a scalar as arguments.

<!-- //"XMTX: addXvec test" -->
<pre>
01 var v1: [3]f32 = .{ 1, 1, 1 };
02 var v2: [3]f32 = .{ 3, 3, 3 };
03 const a: f32 = 2;
04
05 addXvec(&v1, a);
06 prntXvecNl(&v1);
07 //x: 3.0e+00 y: 3.0e+00 z: 3.0e+00
08
09 try std.testing.expectEqual(true, equXvec(&v1, &v2));
10 prntNl();
</pre>

To divide a vector by a scalar value use the <b>divXvec</b> function. It takes a slice, []f32, and a scalar as arguments.

<!-- //"XMTX: divXvec test" -->
<pre>
01 var v1: [3]f32 = .{ 1, 1, 1 };
02 var v2: [3]f32 = .{ 4.76190485e-02, 4.76190485e-02, 4.76190485e-02 };
03 const a: f32 = 21;
04 
05 divXvec(&v1, a);
06 prntXvecNl(&v1);
07 //x: 4.76190485e-02 y: 4.76190485e-02 z: 4.76190485e-02
08 
09 try std.testing.expectEqual(true, equXvec(&v1, &v2));
10 prntNl();
</pre>

To subtract a scalar value from a vector use the <b>subXvec</b> function. It takes a slice, []f32, and a scalar as arguments.

<!-- //"XMTX: subXvec test" -->
<pre>
01 var v1: [3]f32 = .{ 3, 3, 3 };
02 var v2: [3]f32 = .{ 1, 1, 1 };
03 const a: f32 = 2;
04 
05 subXvec(&v1, a);
06 prntXvec(&v1);
07 //x: 1.0e+00 y: 1.0e+00 z: 1.0e+00
08 
09 try std.testing.expectEqual(true, equXvec(&v1, &v2));
10 prntNl();
</pre>

To sum two vectors use the <b>sum1Xvec</b> function which takes two slices, []f32, as arguments.

<!-- //"XMTX: sum1Xvec test" -->
<pre>
01 var v1: [3]f32 = .{ 1, 1, 1 };
02 var v2: [3]f32 = .{ 2, 2, 2 };
03 var v3: [3]f32 = .{ 3, 3, 3 };
04  
05 sum1Xvec(&v1, &v2);
06 prntXvec(&v1);
07 //x: 3.0e+00 y: 3.0e+00 z: 3.0e+00
08 
09 try std.testing.expectEqual(true, equXvec(&v1, &v3));
10 prntNl();
</pre>

To add three vectors together use the <b>sum2Xvec</b> function. This function takes 3 arguments. Each one is a slice, []f32, representing a vector. The 3 vectors are added together and the result is stored in the <b>vecL</b> argument. 

<!--//"XMTX: sum2Xvec test" -->
<pre>
01 var v1: [3]f32 = .{ 1, 1, 1 };
02 var v2: [3]f32 = .{ 2, 2, 2 };
03 var v3: [3]f32 = .{ 5, 5, 5 };
04 
05 sum2Xvec(&v1, &v2, &v2);
06 prntXvec(&v1);
07 //x: 5.0e+00 y: 5.0e+00 z: 5.0e+00
08 
09 try std.testing.expectEqual(true, equXvec(&v1, &v3));
10 prntNl();
</pre>

To find the difference between two vectors use the <b>diff1Vec</b> function. It takes two slices, []f32, as vector arguments and stores the results in the first vector argument.

<!-- //"XMTX: diff1Xvec test" -->
<pre>
01 var v1: [3]f32 = .{ 3, 3, 3 };
02 var v2: [3]f32 = .{ 1, 1, 1 };
03 var v3: [3]f32 = .{ 2, 2, 2 };
04 
05 diff1Xvec(&v1, &v2);
06 prntXvec(&v1);
07 //x: 2.0e+00 y: 2.0e+00 z: 2.0e+00
08 
09 try std.testing.expectEqual(true, equXvec(&v1, &v3));
10 prntNl();
</pre>

To find the difference between three vectors with the result stored in the first vector argument use the <b>diff2Xvec</b> function.

<!-- //"XMTX: diff2Xvec test" -->
<pre>
01 var v1: [3]f32 = .{ 3, 3, 3 };
02 var v2: [3]f32 = .{ 1, 1, 1 };
03 var v3: [3]f32 = .{ 3, 3, 3 };
04 
05 diff2Xvec(&v1, &v2, &v2);
06 prntXvec(&v1);
07 //x: 3.0e+00 y: 3.0e+00 z: 3.0e+00
08 
09 try std.testing.expectEqual(true, equXvec(&v1, &v3));
10 prntNl();
</pre>>