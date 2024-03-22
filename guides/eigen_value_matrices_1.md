# Eigen Values and Matrices

Calculating Eigen values for 2x2 matries is supported by using the EigVal2 structure. The fields of the structure are listed subsequently.

<pre>
EigVal2:
    Used in aiding the calculation of Eigen Values for a 2x2 matrix.

Fields:
    lamdExp: [5]f32 = .{0, 0, 0, 0, 0},
    polyExp: [3]f32 = .{0, 0, 0},
    eignVals: [2]f32 = .{0, 0},
    eignVec1: [2]f32 = .{0, 0},
    eignVec2: [2]f32 = .{0, 0},
</pre>

The main purpose of the <b>EigenVal2</b> structure is to hold pertinent information during the process of calculating the Eigen values of the given matrix.

<!-- "XMTX: ELA - Larson, Edwards: 3.4 Example 1, 3, 5 test" -->
<pre>
01 //Let A = |1,4 ,2,3|   x1 = |1, 1| x2 = |2, -1|
02 //Verify that l1 is an eigenvalue of A corresponding to x1 and that l2 = -1 is an eigenvalue of A corresponding to x2.
03 var mtx: [4]f32 = .{ 1, 4, 2, 3 };
04 var evs2: EigVal2 = .{};
05 var b: bool = false;
06 const expEv1: f32 = 5.0;
07 const expEv2: f32 = -1.0;
08 
09 b = fndEigVal2(&mtx, &evs2);
10 try std.testing.expectEqual(true, b);
11 std.debug.print("Found Eigen Values:\n", .{});
12 prntXvec(&evs2.eignVals);
13 
14 //Found Eigen Values:
15 //x: 5.0e+00 y: -1.0e+00
16 
17 try std.testing.expectEqual(true, isEquF32(expEv1, evs2.eignVals[0], true));
18 try std.testing.expectEqual(true, isEquF32(expEv2, evs2.eignVals[1], true));
</pre>