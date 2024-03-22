# Cramer's Rule and Marices

To apply Cramer’s rule to the given matrix, <b>B</b>, with <b>BA</b> representing the un-augmented version of matrix <b>B</b>, and <b>BAi</b> representing Cramer rule’s supporting matrix.

<!-- "XMTX: ELA - Larson, Edwards: 3.4 Example 1, 3, 5 test" -->
<pre>
01 //Use Cramer's rule to solve the following system of equations: B = |4,-2,10 ,3,-5,11|
02 var B: [6]f32 = .{ 4, -2, 10, 3, -5, 11 };
03 const cols: usize = 3;
04 var BA: [4]f32 = .{ 4, -2, 3, -5 };
05 const colsA: usize = 2;
06 var BAi: [4]f32 = .{ 0, 0, 0, 0 };
07 var Bret: [2]f32 = .{ 0, 0 };
08 const expX1: f32 = 2;
09 const expX2: f32 = -1;
10 
11 b = rslvCramersRule(&B, cols, &BA, colsA, &BAi, &Bret);
12 try std.testing.expectEqual(true, b);
13 std.debug.print("Found Cramer's Rule Return Values:\n", .{});
14 
15 //Found Cramer's Rule Return Values:
16 //x: 2.0e+00 y: -1.0e+00
17 
18 prntXvec(&Bret);
19 try std.testing.expectEqual(true, isEquF32(expX1, Bret[0], true));
20 try std.testing.expectEqual(true, isEquF32(expX2, Bret[1], true));
</pre>