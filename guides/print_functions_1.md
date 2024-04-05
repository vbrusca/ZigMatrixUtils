## Vector, Matrix, and Other Print Functions

There are a number of functions available to assist with printing out vectors and matrices to the terminal. To print out information to the terminal the following methods are available. We'll review each entry in more detail subsquently but we'll list the functions here to start things off.

<ul>
    <li>prntNl</li>
    <li>prntNlStr</li>
    <li>prntNlStrArgs</li>
</ul>

The following examples demonstrate the usage of the simple print functions included in the library.

<!-- //"XMTX: prntNl test" -->
<pre>
01 prntNl();
</pre>

<!-- //"XMTX: prntNlStr test" -->
<pre>
01 prntNlStr("prntNlStr: Hello World");
</pre>

<!-- //"XMTX: prntNlStrArgs test" -->
<pre>
01 try prntNlStrArgs("prntNlStrArgs: {}", .{0});
</pre>

If you need to print a vector to the terminal the following functions can be used.

<ul>
    <li>prntXvec</li>
    <li>prntXvecNl</li>
</ul>

The following examples demonstrate printing vectors to the terminal.

<!-- //"XMTX: prntXvec test" -->
<pre>
01 var v1: [3]f32 = .{ 3, 3, 3 };
02 prntXvec(&v1);
</pre>

<!-- //"XMTX: prntXvecNl test" -->
<pre>
01 var v1: [3]f32 = .{ 3, 3, 3 };
02 prntXvecNl(&v1);
</pre>

If you need to print a matrix to the terminal the following functions can be used.

<ul>
    <li>prntXmtx</li>
    <li>prntXmtxNl</li>
</ul>

The following examples demonstrate printing a matrix to the terminal.

<!-- //"XMTX: prntXmtx test" -->
<pre>
01 var v1: [9]f32 = .{ 3, 3, 3, 0, 0, 0, 1, 1, 1 };
02 prntXmtx(&v1, 3);
</pre>

<!-- //"XMTX: prntXmtx test" -->
<pre>
01 var v1: [9]f32 = .{ 3, 3, 3, 0, 0, 0, 1, 1, 1 };
02 prntXmtx(&v1, 3);
</pre>

The remaining print functions are used less frequently but can come in handy from time to time. I'll list them here.

<ul>
    <li>prntPolyExp</li>
    <li>EigVal2.prnt</li>
    <li>EigVal2.prntNl</li>
    <li>EigVal3.prnt</li>
    <li>EigVal3.prntNl</li>
    <li>ExecTime.prnt</li>
    <li>ExecTime.prntNl</li>
</ul>

In all cases using the print functions is very straight forward. Note that functions that end with an "Nl" are used to print a new line before printing the desired output.