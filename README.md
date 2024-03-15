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
