# Zig Matrix Utils
An open source matrix utility library written in Zig.
The full Zig project is included in this repo.

## Zig Build Version
This project was built against Zig version "zig-windows-x86-0.12.0-dev.3284+153ba46a5".

## Developers
Victor Brusca<br>
Carlo Bruscani

## Documentation
I'll provide a viewable but not 100% usable link to the Zig Matrix Utils documentation here. For proper use clone or download the repo
and view the documentation site locally.<br>
[Zig Generated API Docs](https://htmlpreview.github.io/?https://github.com/vbrusca/ZigMatrixUtils/blob/main/docs/index.html)

## Source Material
This project was built using the following books as a basis.<br>
1. Elementary Linear Algebra by Larson, Edwards
2. Mathematics for 3D Game Programming and Computer Graphics 3rd Edition by Eric Lengyel

## Project Goals
1. Continuous development to complete the material covered in both books regarding linear algebra and matrix manipulations.<br>
2. To keep the project up to date with new versions of Zig as the language matures.
3. To complete and refine the code documentation.

## Running Unit Tests
You can run the unit tests from inside the project with the following command.<br>
<pre>
zig test ./src/XmtxUtils.zig
</pre>
There are over 180 test ran to verify functionality. Feel free to think of them as demonstrations of the associated functions.

## Project Build Commands
How to build an exe (not used).
<pre>
zig build-exe -femit-docs ./src/main.zig -O ReleaseSmall -fstrip -fsingle-threaded -femit-bin="zig-out/bin/main.exe"
</pre>

How to build a static library.
<pre>
zig build-lib -femit-docs ./src/XmtxUtils.zig -O ReleaseSmall -fstrip -fsingle-threaded -femit-bin="zig-out/lib/XmtxUtils.lib"
</pre>

How to build an object.
<pre>
zig build-obj -femit-docs ./src/XmtxUtils.zig -O ReleaseSmall -fstrip -fsingle-threaded -femit-bin="zig-out/lib/XmtxUtils.obj"
</pre>

How to build a dynamic library.
<pre>
zig build-lib -femit-docs ./src/XmtxUtils.zig -lc -dynamic -isystem -fstrip -fsingle-threaded -femit-bin="zig-out/lib/XmtxUtils.dll"
</pre>

## Rough Example of Usage
<pre>
test "XMTX: MF3D - Lengyel: Theorem 3.21 test" {
    prntNl();
    //Let F be an n X n matrix and let the entries of n X n matrix G be defined as
    //Gij = Cji(F) * (1 / detF)
    //where Cji(F) is the cofactor of (F^Tij) then G = F^-1
    //Gij = Cij(F^T) * (1 / detF)

    //Cij(H) = (-1)^(i + j)detH
    std.debug.print("Test 1:\n", .{});
    const alloc: std.mem.Allocator = std.testing.allocator;
    var F: [4]f32 = .{ 5, 6, 8, 9 };
    var invF: [4]f32 = .{ 0, 0, 0, 0 };
    var cols: usize = 2;

    std.debug.print("F Matrix:\n", .{});
    clnXmtx(&F);
    prntXmtx(&F, cols);

    const detF = try detXmtx(&F, cols, &alloc, 0);
    std.debug.print("detF = {}\n", .{detF});

    F = .{ 5, 6, 8, 9 };
    var idtF: [4]f32 = .{ 1, 0, 0, 1 };
    var sclr: f32 = 0.0;
    var b: bool = rdcXmtxInl(&F, cols, false, true, &idtF, 2, false, &sclr);
    try std.testing.expectEqual(true, b);

    std.debug.print("Calculated inverse F (should be identity matrix):\n", .{});
    clnXmtx(&F);
    prntXmtx(&F, cols);
    try std.testing.expectEqual(true, isIdtXmtx(&F, cols));

    std.debug.print("Calculated inverse F (should be inverse matrix):\n", .{});
    clnXmtx(&idtF);
    prntXmtx(&idtF, cols);

    F = .{ 5, 6, 8, 9 };
    b = getInvFromDet2(&F, detF, &invF);
    try std.testing.expectEqual(true, b);

    std.debug.print("Generated inverse from detF (should be inverse matrix):\n", .{});
    clnXmtx(&invF);
    prntXmtx(&invF, cols);
    try std.testing.expectEqual(true, equXmtx(&idtF, &invF));

    std.debug.print("Test 2:\n", .{});
    var F2: [9]f32 = .{ 2, 3, 8, 6, 0, -3, -1, 3, 2 };
    var invF2: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    cols = 3;

    std.debug.print("F2 Matrix:\n", .{});
    clnXmtx(&F2);
    prntXmtx(&F2, cols);

    const detF2 = try detXmtx(&F2, cols, &alloc, 0);
    std.debug.print("detF2 = {}\n", .{detF2});

    F2 = .{ 2, 3, 8, 6, 0, -3, -1, 3, 2 };
    var idtF2: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    sclr = 0.0;
    b = rdcXmtxInl(&F2, cols, false, true, &idtF2, 3, false, &sclr);
    try std.testing.expectEqual(true, b);

    std.debug.print("Calculated inverse F2 (should be identity matrix):\n", .{});
    clnXmtx(&F2);
    prntXmtx(&F2, cols);
    try std.testing.expectEqual(true, isIdtXmtx(&F2, cols));

    std.debug.print("Calculated inverse F2 (should be inverse matrix):\n", .{});
    clnXmtx(&idtF2);
    prntXmtx(&idtF2, cols);

    F2 = .{ 2, 3, 8, 6, 0, -3, -1, 3, 2 };
    b = getInvFromDet3(&F2, detF2, &invF2);
    try std.testing.expectEqual(true, b);

    std.debug.print("Generated inverse from detF2 (should be inverse matrix):\n", .{});
    clnXmtx(&invF2);
    prntXmtx(&invF2, cols);
    try std.testing.expectEqual(true, equXmtx(&idtF2, &invF2));
}
</pre>

#Introduction to the XmtxUtils Library

## Matrix Conversion to Reduced Row Eschelon Form

The library works mainly with arrays of float 32 to represent a vector or an array.
The library is primarily designed to work with smaller matrices 2x2, 3x3, 4x4.
Some operations can be performed on arbitrarily sized matrices or vectors when it makes sense.
<br>
<br>
In this example we'll convert a matrix to reduced row eschelon form.
First we start by declaring some matrices.

<pre>
01  var m1: [12]f32 = .{ 1, -2, 3, 9, -1, 3, 0, -4, 2, -5, 5, 17 };
02  var retM1: [12]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
03  var idtM1: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
04  const dim: usize = 3; //overrides the has augment column difference of 1 and controls the zero row check
05  const cols: usize = 4;
06  var b: bool = true; //holds the return value
07
08  const hasAug: bool = true; //toggles isAugmented flag for the reduction function
09
10  const hasIdt: bool = true; //indicates an matrix was provided to calculate and hold the inverse of m1.
11
12  const triag: bool = false; //A Boolean value indicating if the reduction operation should stop when the matrix is triangular.
13
14  var sclr: f32 = 0.0;
</pre>

On lines 1 we declare a 3x4 matrix containing out test data.
On line 2 we declare a en empty matrix used to hold the results and on line
3 we declare an empty matrix ment to eventually hold the identity matrix.
Next we define the dimension of the input matrix, and subsequently the column count.
Then we set up Boolean flags to handle some function arguments. Lastly, a scalar variable to track multiplication against the matrix. This comes in handy in some cases but isn't important here.

<pre>
15  b = rdcXmtx(&m1, cols, hasAug, &retM1, hasIdt, &idtM1, dim, triag, &sclr);
16
17  std.debug.print("Matrix M1:\n", .{});
18  prntXmtx(&m1, cols);
19  prntNl();
20  
21  //Initial matrix state
22  //Matrix M1:
23  //0: x: 1.0e+00 y: -2.0e+00 z: 3.0e+00 w: 9.0e+00
24  //1: x: -1.0e+00 y: 3.0e+00 z: 0.0e+00 w: -4.0e+00
25  //2: x: 2.0e+00 y: -5.0e+00 z: 5.0e+00 w: 1.7e+01  
26  
</pre>

Next on line 15 we call the rdcXmtx function and reduce the matrix to row eschelon form. The result of the operation is stored in b. We're going to look right into the matrices and ignore that value but this is a good spot for an assertion if your writing a unit test. On line 18 we print the matrix and the contents are listed on lines 23 - 25.

<pre>
27  std.debug.print("Matrix Ret:\n", .{});
28  prntXmtx(&retM1, cols);
29  prntNl();
30  
31  //Reduced Row Eschelon Form
32  //Matrix Ret:
33  //0: x: 1.0e+00 y: 0.0e+00 z: 0.0e+00 w: 1.0e+00
34  //1: x: 0.0e+00 y: 1.0e+00 z: 0.0e+00 w: -1.0e+00
35  //2: x: 0.0e+00 y: 0.0e+00 z: 1.0e+00 w: 2.0e+00
36  
</pre>

The next matrix print out shows the reduced form and solutions.

<pre>
37  std.debug.print("Matrix Inv:\n", .{});
38  prntXmtx(&idtM1, dim);
39  prntNl();
40  
41  //Inverse matrix
42  //Matrix Inv:
43  //0: x: 7.5e+00 y: -2.5e+00 z: -4.5e+00
44  //1: x: 2.5e+00 y: -5.0e-01 z: -1.5e+00
45  //2: x: -5.0e-01 y: 5.0e-01 z: 5.0e-01 
46
</pre>

The matrix inverse is calculated from the provided identity matrix.

<pre>
47  cpyLessXmtx(&retM1, &idtM1, cols, dim);
48  std.debug.print("Copy Matrix M1:\n", .{});
49  prntXmtx(&idtM1, dim);
50  clnXmtx(&idtM1);
51  prntNl();
52
53  //Partial copy of the resulting row eschelon form matrix
54  //into an empty matrix, should be the identity matrix
55  //Copy Matrix M1:
56  //0: x: 1.0e+00 y: 0.0e+00 z: 0.0e+00
57  //1: x: 0.0e+00 y: 1.0e+00 z: 0.0e+00
58  //2: x: 0.0e+00 y: 0.0e+00 z: 1.0e+00
</pre>

Lastly we copy of a part of the solution matrix, the 3x3 non-augmentes matrix
and show that it is the identity matrix.
