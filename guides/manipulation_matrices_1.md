# Manipulation of Matrices

There are a number of functions that help with matrix manipulations.
Some of them are listed as follows. This set of functions is particularly useful for Gaussian elimination.

<ul>
    <li>addSclMulXmtxRows</li>
    <li>addSclMulXmtxRowsInl</li>
    <li>addSclXmtxRows</li>
    <li>addSclXmtxRowsInl</li>
    <li>altXmtxRows</li>
    <li>altXmtxRowsInl</li>
</ul>

The <b>addSclMulXmtxRows</b> and <b>addSclMulXmtxRowsInl</b> functions add a scalar multiple of one matrix row to another. In the example below we add a scalar multiple of one matrix row to another and then check the results against a known matrix, exp. The function takes the source row, the destination row, a scalar multiplier, the number of columns in the matrix, the source matrix, and the destination matrix as aruments, line 12. The inline version of this function performs the matrix manipulation on the source matrix, with no return matrix specified.

<!-- //"XMTX: addSclMulXmtxRows test" -->
<pre>
01 var mtx: [9]f32 = .{ 1, 1, 1, 2, 2, 2, 3, 3, 3 };
02 var res: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
03 var exp: [9]f32 = .{ 1, 1, 1, 2, 2, 2, 7, 7, 7 };
04 
05 prntXmtx(&mtx, 3);
06 
07 //0: x: 1.0e+00 y: 1.0e+00 z: 1.0e+00 
08 //1: x: 2.0e+00 y: 2.0e+00 z: 2.0e+00 
09 //2: x: 3.0e+00 y: 3.0e+00 z: 3.0e+00 
10 
11 const start = try Instant.now();
12 addSclMulXmtxRows(1, 2, 2, 3, &mtx, &res);
13 const end = try Instant.now();
14 const elapsed1: f64 = @floatFromInt(end.since(start));
15 std.debug.print("\naddSclMulXmtxRows: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
16 addExecTime("addSclMulXmtxRows", elapsed1);
17 
18 prntXmtx(&res, 3);
19 
20 //addSclMulXmtxRows: Time elapsed is: 0.001ms, 500.000ns
21 //0: x: 1.0e+00 y: 1.0e+00 z: 1.0e+00 
22 //1: x: 2.0e+00 y: 2.0e+00 z: 2.0e+00 
23 //2: x: 7.0e+00 y: 7.0e+00 z: 7.0e+00 
24 
25 try std.testing.expectEqual(true, equXvec(&exp, &res));
</pre>

To add a scalar to a matrix row use either the <b>addSclXmtxRows</b> or the <b>addSclXmtxRowsInl</b> function. An example is as follows. Note the arguments to the function are the source row, the scalar amount, the number of columns in the matrix, a source matrix, and a destination matrix, line 11. The inline version of this function performs the matrix manipulation on the source matrix, with no return matrix specified.

<!-- //"XMTX: addSclXmtxRows test" -->
<pre>
01 var mtx: [9]f32 = .{ 1, 1, 1, 1, 1, 1, 1, 1, 1 };
02 var res: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
03 var exp: [9]f32 = .{ 4, 4, 4, 1, 1, 1, 1, 1, 1 };
04 prntXmtx(&mtx, 3);
05 
06 //0: x: 1.0e+00 y: 1.0e+00 z: 1.0e+00 
07 //1: x: 1.0e+00 y: 1.0e+00 z: 1.0e+00 
08 //2: x: 1.0e+00 y: 1.0e+00 z: 1.0e+00 
09 
10 const start = try Instant.now();
11 addSclXmtxRows(0, 3, 3, &mtx, &res);
12 const end = try Instant.now();
13 const elapsed1: f64 = @floatFromInt(end.since(start));
14 std.debug.print("\naddSclXmtxRows: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
15addExecTime("addSclXmtxRows", elapsed1);
16
17 //addSclXmtxRowsInl: Time elapsed is: 0.005ms, 4900.000ns
18 //0: x: 4.0e+00 y: 4.0e+00 z: 4.0e+00 
19 //1: x: 1.0e+00 y: 1.0e+00 z: 1.0e+00 
20 //2: x: 1.0e+00 y: 1.0e+00 z: 1.0e+00 
21 
22 prntXmtx(&res, 3);
23 try std.testing.expectEqual(true, equXvec(&exp, &res));
</pre>

To alternate two rows of a matrix use the <b>altXmtxRows</b> or the <b>altXmtxRowsInl</b> function. The function takes the source row, the destination row, the number of columns in the matrix, and a souce and destination matrix as arguments. The inline version of this function makes the changes to the source matrix and not a target output matrix.

<!-- //"XMTX: altXmtxRows test" -->
<pre>
01 var mtx: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
02 var exp: [9]f32 = .{ 7, 8, 9, 4, 5, 6, 1, 2, 3 };
03 var alt: [9]f32 = .{ 1, 1, 1, 1, 1, 1, 1, 1, 1 };
04 prntXmtx(&mtx, 3);
05 
06 //0: x: 1.0e+00 y: 2.0e+00 z: 3.0e+00 
07 //1: x: 4.0e+00 y: 5.0e+00 z: 6.0e+00 
08 //2: x: 7.0e+00 y: 8.0e+00 z: 9.0e+00 
09 
10 const start = try Instant.now();
11 altXmtxRows(0, 2, 3, &mtx, &alt);
12 const end = try Instant.now();
13 const elapsed1: f64 = @floatFromInt(end.since(start));
14 std.debug.print("\naltXmtxRows: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
15 addExecTime("altXmtxRows", elapsed1);
16 
17 prntXmtx(&alt, 3);
18 
19 //altXmtxRows: Time elapsed is: 0.001ms, 500.000ns
20 //0: x: 7.0e+00 y: 8.0e+00 z: 9.0e+00 
21 //1: x: 4.0e+00 y: 5.0e+00 z: 6.0e+00 
22 //2: x: 1.0e+00 y: 2.0e+00 z: 3.0e+00 
23 
24 try std.testing.expectEqual(true, equXvec(&exp, &alt));
</pre>