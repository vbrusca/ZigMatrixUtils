# Zig Matrix Utils
An open source matrix utility library written in Zig.
The full Zig project is included in this repo.

## Zig Build Version
This project was built against Zig version "zig-windows-x86-0.12.0-dev.3284+153ba46a5".

## Developers
Victor Brusca<br>
Carlo Bruscani

## Documentation
I'll provide a viewable but not 100% usable link to the Zig Matrix Utils documentation here. Links don't work but you can browse the functions and their signature. For proper use clone or download the repo and view the documentation site locally.<br> [Zig Generated API Docs](https://htmlpreview.github.io/?https://github.com/vbrusca/ZigMatrixUtils/blob/main/docs/index.html)

## Source Material
This project was built using the following books as a basis.<br>
1. Elementary Linear Algebra by Larson, Edwards
2. Mathematics for 3D Game Programming and Computer Graphics 3rd Edition by Eric Lengyel

## Project Goals
1. LT: Continuous development to complete the material covered in both books regarding linear algebra and matrix manipulations.<br>
2. LT: To keep the project up to date with new versions of Zig as the language matures.
3. LT: To complete and refine the code documentation.
4. ST: Version 0.6: Add basic matrix scaling functions.
5. ST: Version 0.6: Add execution time tracking to the remaining theorem/problem unit tests.
6. ST: Version 0.6: Add support for fast floating point math and check performance.

## Running Unit Tests
You can run the full set of unit tests from inside the project with the following command.<br>
<pre>
zig test ./src/main.zig
</pre>
There are over 240 test ran to verify functionality. Feel free to think of them as demonstrations of the associated functions. You can also capture the test output with the following command on DOS terminals. You'll have to search around and find an equivalent command if you are on MacOS, Linux, or Unix for your respective shell.
<pre>
zig test ./src/XmtxUtils.zig > all_test_output.txt 2>&1
</pre>
Currently the library sets each module to use fast floating point math. This has already shown an impact in the performance of different functions in the function execution time list. If there is some instability in floating point math just comment out this line in the header of main.zig and XmtxUtils.zip.
<pre>
comptime {
    @setFloatMode(std.builtin.FloatMode.Optimized);
}
</pre>

## Project Build Commands
How to build an exe (NOT USED). There is no real main code for this librar currently.
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

## Guides

You can look into specific use cases for the library in the "Guides" section.
**Please note that the guides may be slightly out of sync with regard to the latest version of the library as it develops. Those gaps will be closed periodically over time and as a better way of tracking the guides associated with code changes evolves. The guides associate library functions with vector and matrix actions etc. You can use them as a loose example of how to use the library.

[Guides](https://github.com/vbrusca/ZigMatrixUtils/tree/main/guides)

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

## Execution Times

Execution times are calculated after the last unit test that uses them. You can find the execution time summary in the test output by searching for the term "Function Execution Times List". An example of the execution time summary is shown below.

<pre>
Function Execution Times List:
 absF32:	Count: 1.0e+00	Avg: 0.000ms 200.000ns
 absF32Ref:	Count: 1.0e+00	Avg: 0.000ms 100.000ns
 absF32Ret:	Count: 1.0e+00	Avg: 0.000ms 100.000ns
 absXmtx:	Count: 2.0e+00	Avg: 0.000ms 300.000ns
 addSclMulXmtxRows:	Count: 1.0e+00	Avg: 0.001ms 500.000ns
 addSclMulXmtxRowsInl:	Count: 1.0e+00	Avg: 0.000ms 200.000ns
 addSclXmtxRows:	Count: 1.0e+00	Avg: 0.000ms 300.000ns
 addSclXmtxRowsInl:	Count: 1.0e+00	Avg: 0.005ms 4900.000ns
 addXvec:	Count: 1.0e+00	Avg: 0.000ms 200.000ns
 adjXmtx3:	Count: 1.0e+00	Avg: 0.012ms 12200.000ns
 adjXmtx4:	Count: 1.0e+00	Avg: 0.013ms 12600.000ns
 aglBtwnXvec:	Count: 1.0e+00	Avg: 0.009ms 9400.000ns
 altXmtxRows:	Count: 1.0e+00	Avg: 0.001ms 500.000ns
 ...
</pre>

You can see the function name, the usage count, and the average execution time in ms and ns. To effectively turn off this memory allocation set <b>MAX_EXEC_TIMES</b> equal to 1.
