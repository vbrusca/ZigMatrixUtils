# Encryption Using Matrices

<!-- //"XMTX: ELA - Larson, Edwards: 2.5 Example 5, 6 test" -->
<pre>
//Matrix Cryptography
//0 = _     6 = F       12 = L      18 = R      24 = X
//1 = A     7 = G       13 = M      19 = S      25 = Y
//2 = B     8 = H       14 = N      20 = T      26 = Z
//3 = C     9 = I       15 = O      21 = U
//4 = D     10 = J      16 = P      22 = V
//5 = E     11 = K      17 = Q      23 = W

//Clear text matrices
//[13  5  5] [20  0 13] [5  0 13] [15 14  4] [1 25  0]
//  M  E  E    T  _  M   E  _  M    O  N  D   A  Y  _

//A = 1  -2   2
//   -1   1   3
//    1  -1  -4
var A: [9]f32 = .{ 1, -2, 2, -1, 1, 3, 1, -1, -4 };
var tmp3: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
var idt3: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
var b: bool = false;
var w1: [3]f32 = .{ 13, 5, 5 };
var w2: [3]f32 = .{ 20, 0, 13 };
var w3: [3]f32 = .{ 5, 0, 13 };
var w4: [3]f32 = .{ 15, 14, 4 };
var w5: [3]f32 = .{ 1, 25, 0 };

var eW1: [3]f32 = .{ 0, 0, 0 };
var eW2: [3]f32 = .{ 0, 0, 0 };
var eW3: [3]f32 = .{ 0, 0, 0 };
var eW4: [3]f32 = .{ 0, 0, 0 };
var eW5: [3]f32 = .{ 0, 0, 0 };

var cW1: [3]f32 = .{ 0, 0, 0 };
var cW2: [3]f32 = .{ 0, 0, 0 };
var cW3: [3]f32 = .{ 0, 0, 0 };
var cW4: [3]f32 = .{ 0, 0, 0 };
var cW5: [3]f32 = .{ 0, 0, 0 };

//Get inverse matrix
var sclr: f32 = 0.0;
b = rdcXmtx(&A, 3, false, &tmp3, true, &idt3, 3, false, &sclr);
try std.testing.expectEqual(true, b);

prntNl();
std.debug.print("XMTX: A:\n", .{});
prntXmtx(&A, 3);

//XMTX: A:
//0: x: 1.0e+00 y: -2.0e+00 z: 2.0e+00
//1: x: -1.0e+00 y: 1.0e+00 z: 3.0e+00
//2: x: 1.0e+00 y: -1.0e+00 z: -4.0e+00

prntNl();
std.debug.print("XMTX: TMP3:\n", .{});
prntXmtx(&tmp3, 3);

//XMTX: TMP3:
//0: x: 1.0e+00 y: 0.0e+00 z: 0.0e+00
//1: x: -0.0e+00 y: 1.0e+00 z: 0.0e+00
//2: x: -0.0e+00 y: -0.0e+00 z: 1.0e+00 

prntNl();
std.debug.print("XMTX: IDT3:\n", .{});
prntXmtx(&idt3, 3);

//XMTX: IDT3:
//0: x: -1.0e+00 y: -1.0e+01 z: -8.0e+00
//1: x: -1.0e+00 y: -6.0e+00 z: -5.0e+00
//2: x: -0.0e+00 y: -1.0e+00 z: -1.0e+00

A = .{ 1, -2, 2, -1, 1, 3, 1, -1, -4 };

//Encode word 1
b = tmsXmtx(&w1, 3, &A, 3, &eW1, 3);
try std.testing.expectEqual(true, b);

prntNl();
std.debug.print("XMTX: W1:\n", .{});
prntXmtx(&w1, 3);

prntNl();
std.debug.print("XMTX: ENC_W1:\n", .{});
prntXmtx(&eW1, 3);

//Encode word 2
b = tmsXmtx(&w2, 3, &A, 3, &eW2, 3);
try std.testing.expectEqual(true, b);

prntNl();
std.debug.print("XMTX: W2:\n", .{});
prntXmtx(&w2, 3);

prntNl();
std.debug.print("XMTX: ENC_W2:\n", .{});
prntXmtx(&eW2, 3);

//Encode word 3
b = tmsXmtx(&w3, 3, &A, 3, &eW3, 3);
try std.testing.expectEqual(true, b);

prntNl();
std.debug.print("XMTX: W3:\n", .{});
prntXmtx(&w3, 3);

prntNl();
std.debug.print("XMTX: ENC_W3:\n", .{});
prntXmtx(&eW3, 3);

//Encode word 4
b = tmsXmtx(&w4, 3, &A, 3, &eW4, 3);
try std.testing.expectEqual(true, b);

prntNl();
std.debug.print("XMTX: W4:\n", .{});
prntXmtx(&w4, 3);

prntNl();
std.debug.print("XMTX: ENC_W4:\n", .{});
prntXmtx(&eW4, 3);

//Encode word 5
b = tmsXmtx(&w5, 3, &A, 3, &eW5, 3);
try std.testing.expectEqual(true, b);

prntNl();
std.debug.print("XMTX: W5:\n", .{});
prntXmtx(&w5, 3);

prntNl();
std.debug.print("XMTX: ENC_W5:\n", .{});
prntXmtx(&eW5, 3);

//Decode word 1
b = tmsXmtx(&eW1, 3, &idt3, 3, &cW1, 3);
try std.testing.expectEqual(true, b);

prntNl();
std.debug.print("XMTX: DEC_W1:\n", .{});
prntXmtx(&cW1, 3);

//Decode word 2
b = tmsXmtx(&eW2, 3, &idt3, 3, &cW2, 3);
try std.testing.expectEqual(true, b);

prntNl();
std.debug.print("XMTX: DEC_W2:\n", .{});
prntXmtx(&cW2, 3);

//Decode word 3
b = tmsXmtx(&eW3, 3, &idt3, 3, &cW3, 3);
try std.testing.expectEqual(true, b);

prntNl();
std.debug.print("XMTX: DEC_W3:\n", .{});
prntXmtx(&cW3, 3);

//Decode word 4
b = tmsXmtx(&eW4, 3, &idt3, 3, &cW4, 3);
try std.testing.expectEqual(true, b);

prntNl();
std.debug.print("XMTX: DEC_W4:\n", .{});
prntXmtx(&cW4, 3);

//Decode word 5
b = tmsXmtx(&eW5, 3, &idt3, 3, &cW5, 3);
try std.testing.expectEqual(true, b);

prntNl();
std.debug.print("XMTX: DEC_W5:\n", .{});
prntXmtx(&cW5, 3);

try std.testing.expectEqual(true, equXmtx(&w1, &cW1));
try std.testing.expectEqual(true, equXmtx(&w2, &cW2));
try std.testing.expectEqual(true, equXmtx(&w3, &cW3));
try std.testing.expectEqual(true, equXmtx(&w4, &cW4));
try std.testing.expectEqual(true, equXmtx(&w5, &cW5));
</pre>
<!--
XMTX: W1:
0: x: 1.3e+01 y: 5.0e+00 z: 5.0e+00

XMTX: ENC_W1:
0: x: 1.3e+01 y: -2.6e+01 z: 2.1e+01 

XMTX: W2:
0: x: 2.0e+01 y: 0.0e+00 z: 1.3e+01

XMTX: ENC_W2:
0: x: 3.3e+01 y: -5.3e+01 z: -1.2e+01 

XMTX: W3:
0: x: 5.0e+00 y: 0.0e+00 z: 1.3e+01

XMTX: ENC_W3:
0: x: 1.8e+01 y: -2.3e+01 z: -4.2e+01 

XMTX: W4:
0: x: 1.5e+01 y: 1.4e+01 z: 4.0e+00

XMTX: ENC_W4:
0: x: 5.0e+00 y: -2.0e+01 z: 5.6e+01

XMTX: W5:
0: x: 1.0e+00 y: 2.5e+01 z: 0.0e+00

XMTX: ENC_W5:
0: x: -2.4e+01 y: 2.3e+01 z: 7.7e+01

XMTX: DEC_W1:
0: x: 1.3e+01 y: 5.0e+00 z: 5.0e+00 

XMTX: DEC_W2:
0: x: 2.0e+01 y: 0.0e+00 z: 1.3e+01

XMTX: DEC_W3:
0: x: 5.0e+00 y: 0.0e+00 z: 1.3e+01

XMTX: DEC_W4:
0: x: 1.5e+01 y: 1.4e+01 z: 4.0e+00

XMTX: DEC_W5:
0: x: 1.0e+00 y: 2.5e+01 z: 0.0e+00
-->