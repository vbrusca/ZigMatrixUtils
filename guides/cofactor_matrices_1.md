# Cofactor of Matrices

In order to calculate the cofactors of a given X by X marix you need to calculate the determinant of a number of (X - 1) by (X - 1) matrices.
<br>
<br>
The initialization of some demo data starts on line 1 with the <b>B</b> variable. This is the matrix we're going to calculate cofactors for. The answer to the problem is stored in the <b>expB</b> matrix and used for a final test.
The <b>cofSign</b> matrix holds the sign of the determinant at the given cofactor matrix location. The <b>cofB</b> matrix holds the resulting cofactor value and sign and the <b>retB</b> matrix is used to specify the inner matrix the cofactor is based on.
<br>
<br> 
The next few variables are used in the cofactor matrix creation loop. The <b>b</b> bool is used to track the result of certain function calls. The <b>row</b> and <b>col</b> variable holds the current row and col of the loop. The last two variables are the counts for rows and cols in this example.

<!-- "XMTX: ELA - Larson, Edwards: 3.1 Example 1, 2, 3, 4, 5 test" -->
<pre>
01 var B: [9]f32 = .{ 0, 2, 1, 3, -1, 2, 4, 0, 1 };
02 var expB: [9]f32 = .{ -1, 5, 4, -2, -4, 8, 5, 3, -6 };
03 var cofSign: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
04 var cofB: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
05 var retB: [4]f32 = .{ 0, 0, 0, 0 };
06 var b: bool = false;
07 var row: usize = 0;
08 var col: usize = 0;
09 const rows: usize = 3;
10 const cols: usize = 3;
</pre>

The next segment of code shows a nested for loop that is used to determine the cofactor and sign of each position in a cofactor matrix. First we store the sign of the current matrix position on line 17. Determine the cofactor matrix for the current <b>row</b> and <b>col</b>, line 19. The final value of the cofactor is set on line 28 using the deerminant of the cofactor matrix and the cofactor sign. Each iteration has some logging, line 29, and a test is done on line 32 to make sure the calculation works as expected.

<pre>
11 while (row < rows) : (row += 1) {
12     col = 0;
13     while (col < cols) : (col += 1) {
14         prntNl();
15         std.debug.print("Row: {} Col: {}\n", .{ row, col });
16 
17         cofSign[(row * cols) + col] = cofXmtxSign(row, col, true);
18 
19         b = cofXmtx(&B, cols, row, col, &retB, 2);
20         try std.testing.expectEqual(true, b);
21 
22         std.debug.print("XMTX: B:\n", .{});
23         prntXmtx(&B, cols);
24 
25         std.debug.print("XMTX: retB:\n", .{});
26         prntXmtx(&retB, 2);
27 
28         cofB[(row * cols) + col] = detXmtx2(&retB) * cofSign[(row * cols) + col];
29         std.debug.print("Sign: {} Val: {}\n", .{ cofSign[(row * cols) + col], cofB[(row * cols) + col] });
30     }
31 }
32 try std.testing.expectEqual(true, equXmtx(&cofB, &expB));
</pre>

The output of the above while loops is shown below.

<pre>
Row: 0 Col: 0
XMTX: retB:
0: x: -1.0e+00 y: 2.0e+00
1: x: 0.0e+00 y: 1.0e+00
Sign: 1.0e+00 Val: -1.0e+00

Row: 0 Col: 1
XMTX: retB:
0: x: 3.0e+00 y: 2.0e+00
1: x: 4.0e+00 y: 1.0e+00
Sign: -1.0e+00 Val: 5.0e+00

Row: 0 Col: 2
XMTX: retB:
0: x: 3.0e+00 y: -1.0e+00
1: x: 4.0e+00 y: 0.0e+00 
Sign: 1.0e+00 Val: 4.0e+00

Row: 1 Col: 0
XMTX: retB:
0: x: 2.0e+00 y: 1.0e+00
1: x: 0.0e+00 y: 1.0e+00 
Sign: -1.0e+00 Val: -2.0e+00

Row: 1 Col: 1
XMTX: retB:
0: x: 0.0e+00 y: 1.0e+00
1: x: 4.0e+00 y: 1.0e+00
Sign: 1.0e+00 Val: -4.0e+00

Row: 1 Col: 2
XMTX: retB:
0: x: 0.0e+00 y: 2.0e+00
1: x: 4.0e+00 y: 0.0e+00
Sign: -1.0e+00 Val: 8.0e+00

Row: 2 Col: 0
XMTX: retB:
0: x: 2.0e+00 y: 1.0e+00
1: x: -1.0e+00 y: 2.0e+00
Sign: 1.0e+00 Val: 5.0e+00

Row: 2 Col: 1
XMTX: retB:
0: x: 0.0e+00 y: 1.0e+00 
1: x: 3.0e+00 y: 2.0e+00
Sign: -1.0e+00 Val: 3.0e+00

Row: 2 Col: 2
XMTX: retB:
0: x: 0.0e+00 y: 2.0e+00
1: x: 3.0e+00 y: -1.0e+00
Sign: 1.0e+00 Val: -6.0e+00
</pre>

Note the values found here and the epexed values in the <b>expB</b> matrix.