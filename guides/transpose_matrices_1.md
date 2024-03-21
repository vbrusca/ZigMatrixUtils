#Transpose of Matrices

The transpose of a matrix can be calculated using the <b>trnXmtx</b> or <b>trnXmtxRect</b> functions, lines 31 - 32, 33 - 34, respectively.

//"XMTX: ELA - Larson, Edwards: 2.2 Example 8 test"
<pre>
01 var m1: [2]f32 = .{ 2, 8 };
02 var r1: [2]f32 = .{ 0, 0 };
03 var e1: [2]f32 = .{ 2, 8 };
04 //A = 2 => 2 8
05 //    8
06 //
07 
08 var m2: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
09 var r2: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
10 var e2: [9]f32 = .{ 1, 4, 7, 2, 5, 8, 3, 6, 9 };
11 //B = 1  2  3 => 1  4  7
12 //    4  5  6    2  5  8
13 //    7  8  9    3  6  9
14 //
15 
16 var m3: [9]f32 = .{ 1, 2, 0, 2, 1, 0, 0, 0, 1 };
17 var r3: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
18 var e3: [9]f32 = .{ 1, 2, 0, 2, 1, 0, 0, 0, 1 };
19 //C = 1  2  0 => 1  2  0
20 //    2  1  0    2  1  0
21 //    0  0  1    0  0  1
22 //
23 
24 var m4: [6]f32 = .{ 0, 1, 2, 4, 1, -1 };
25 var r4: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
26 var e4: [6]f32 = .{ 0, 2, 1, 1, 4, -1 };
27 //D = 0  1 => 0  2  1
28 //    2  4    1  4 -1
29 //    1 -1
30
31 trnXmtx(&m1, 1, &r1);
32 trnXmtx(&m2, 3, &r2);
33 trnXmtxRect(&m3, 3, &r3, 3);
34 trnXmtxRect(&m4, 2, &r4, 3);
35 
36 prntXmtx(&m4, 2);
37 prntNl();
38 
39 //0: x: 0.0e+00 y: 1.0e+00
40 //1: x: 2.0e+00 y: 4.0e+00
41 //2: x: 1.0e+00 y: -1.0e+00
42 
43 prntXmtx(&r4, 3);
44 prntNl();
45 
45 //0: x: 0.0e+00 y: 2.0e+00 z: 1.0e+00
46 //1: x: 1.0e+00 y: 4.0e+00 z: -1.0e+00
47 
48 const b1: bool = equXmtx(&r1, &e1);
49 const b2: bool = equXmtx(&r2, &e2);
50 const b3: bool = equXmtx(&r3, &e3);
51 const b4: bool = equXmtx(&r4, &e4);
</pre>

Note that one tranpose matrix function, <b>trnXmtx</b> assumes the matrices provided are square based on the input matrix. The second version of the function, <b>trnXmtxRect</b> allows you to specify the column count of the matrix arguments. 