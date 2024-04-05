# Properties of Vectors

There are a couple of functions that deterime if certain vector properties exist.
The <b>isZeroXvec</b> function is used to determine if a vector is the zero vector.
The function takes an argument that coerces to a slice, []f32.

<!-- //"XMTX: isZeroXvec test" -->
<pre>
01 var v1: [3]f32 = .{ 0, 0, 0 };
02 var v2: [3]f32 = .{ 1, 2, 3 };
03 var b: bool = false;
04 
05 b = isZeroXvec(&v1);
06 try std.testing.expectEqual(true, b);
07 
08 b = isZeroXvec(&v2);
09 try std.testing.expectEqual(false, b);
10 prntNl();
</pre>

A second property of vectors we can check for is linear independence. We can use the <b>isLinIndXvec</b> function for this. This function takes two arguments a left and right vector to test for linear independence.

<!-- //"XMTX: isLinIndXvec test" -->
<pre>
01 var v1: [3]f32 = .{ 1, 0, 0 };
02 var v2: [3]f32 = .{ 0, 1, 0 };
03 var v3: [3]f32 = .{ 0, 3, 0 };
04 var b: bool = false;
05 
06 b = isLinIndXvec(&v1, &v2);
07 try std.testing.expectEqual(true, b);
08 
09 b = isLinIndXvec(&v2, &v3);
10 try std.testing.expectEqual(false, b);
11 prntNl();
</pre>