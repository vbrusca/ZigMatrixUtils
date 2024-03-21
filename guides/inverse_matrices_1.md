#Inverse of Matrices

To find the inverse of a matrix we need to reduce it to row eschelon form
while tracking the matrix permutations on an identity matrix. This can all be handled with one call to the <b>rdcXmtx</b> function.
<br>
<br>
Arguments 5 and 6 of the function call are a Boolean value that indicates if an identity matrix has been provided 

//"XMTX: ELA - Larson, Edwards: 2.3 Example 2 test"
</pre>
//Find the inverse of A and verify it is correct.
//A = 1  4
//   -1  3
//
//A^-1 = -3 -4
//        1  1
//
var A: [4]f32 = .{ 1, 4, -1, -3 };
var I: [4]f32 = .{ 1, 0, 0, 1 };
var B: [4]f32 = .{ 0, 0, 0, 0 };
var exp: [4]f32 = .{ -3, -4, 1, 1 };
var res: bool = false;

var sclr: f32 = 0.0;
res = rdcXmtx(&A, 2, false, &B, true, &I, 2, false, &sclr);
try std.testing.expectEqual(true, res);

prntXmtx(&B, 2);
prntNl();

prntXmtx(&I, 2);
prntNl();
try std.testing.expectEqual(true, equXmtx(&exp, &I));
</pre>