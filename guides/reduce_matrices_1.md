# Reduction of Matrices

## Matrix Conversion to Reduced Row Eschelon Form

The library works mainly with arrays of float 32 to represent a vector or an array. The library is primarily designed to work with smaller matrices 2x2, 3x3, 4x4. Some operations can be performed on arbitrarily sized matrices or vectors when it makes sense.
<br>
<br>
In this example we'll convert a matrix to reduced row eschelon form.
First we start by declaring some matrices.

<!-- //"XMTX: ELA - Larson, Edwards: 1.2 Example 3 test" -->
<pre>
01  var m1: [12]f32 = .{ 1, -2, 3, 9, -1, 3, 0, -4, 2, -5, 5, 17 };
02  var retM1: [12]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
03  var idtM1: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };

04  const dim: usize = 3; //overrides the has augment column difference of 1 and controls the zero row check

05  const cols: usize = 4;
06  var b: bool = true; //toggles isAugmented flag for the reduction function
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

Next on line 15 we call the rdcXmtx function and reduce the matrix to row eschelon form. The result of the operation is stored in b. We're going to look right into the matrices and ignore that value but this is a good spot for an assertion if your writing a unit test. On line 18 we print the matrix and the contents are listed on lines 23 - 25. Note that the output is from the <b>printXmtx</b> function.

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
49  clnXmtx(&idtM1);
50  prntXmtx(&idtM1, dim);
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
and show that it is the identity matrix. Note the call on line 49 to the <b>clnXmtx</b> function used to clean up the matrix entry values rounding them to the nearest significance.