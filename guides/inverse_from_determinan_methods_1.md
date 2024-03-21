<!-- "XMTX: getInvFromDet2 test" -->
<pre>
var m1: [4]f32 = .{ 2, 1, 7, 4 };
var invM1: [4]f32 = .{ 4, -1, -7, 2 };
var res: [4]f32 = std.mem.zeroes([4]f32); //.{};
const detM1: f32 = detXmtx2(&m1);
const b: bool = getInvFromDet2(&m1, detM1, &res);
try std.testing.expectEqual(true, b);
try std.testing.expectEqual(true, equXmtx(&res, &invM1));
</pre>

<!-- "XMTX: getInvFromDet3 test" -->
<pre>
//A = 1 2 3
//    4 5 6
//    7 2 9

//detA = (1)(33)          - (2)(-6)          + (3)(-27)
//detA = (1)(det 5 6 2 9) - (2)(det 4 6 7 9) + (3)(det 4 5 7 2)
//detA = -36;

//A^-1 = -11/12  1/3  1/12
//       -1/6    1/3  -1/6
//       3/4     -1/3 1/12

var A: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 2, 9 };
var res: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
const detA: f32 = -36;
var invA: [9]f32 = .{ (-11.0 / 12.0), (1.0 / 3.0), (1.0 / 12.0), (-1.0 / 6.0), (1.0 / 3.0), (-1.0 / 6.0), (3.0 / 4.0), (-1.0 / 3.0), (1.0 / 12.0) };
const b: bool = getInvFromDet3(&A, detA, &res);
try std.testing.expectEqual(true, b);
try std.testing.expectEqual(true, equXmtx(&res, &invA));
</pre>

<!-- "XMTX: getInvFromDet4 test" -->
<pre>
var ret: [16]f32 = std.mem.zeroes([16]f32); //.{};
var mtx: [16]f32 = .{ 1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1 };
var expMtx: [16]f32 = .{ (1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0), (-1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0), (-1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0), (-1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0), (-1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0) };
const det: f32 = detXmtx4(&mtx);
const expDet: f32 = -16;
const b: bool = getInvFromDet4(&mtx, det, &ret);
const cols: usize = 4;

std.debug.print("getInvFromDet4 det: {} expDet: {}\n", .{ det, expDet });
prntXmtx(&ret, cols);
prntNl();

std.debug.print("getInvFromDet4 det:\n", .{});
prntXmtx(&expMtx, cols);
prntNl();

try std.testing.expectEqual(true, b);
try std.testing.expectEqual(true, isEquF32(expDet, det, true));
try std.testing.expectEqual(true, equXmtx(&expMtx, &ret));
</pre>
