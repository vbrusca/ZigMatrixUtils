# Polynomial Order

In some cases during certain matrix manipulations, Eigen Values, we need to work with characteristic equations of a certain order. Those equations are represented as an array of floating point values that act as the values of the polynomial at that position.
<br>
<br>
To find the order of a given polynomial represented array of floating point values you can use the <b>getPolyOrder</b> function. The polynomial order is always one less than the number of full terms in the polynomial.

<!-- //"XMTX: getPolyOrder test" -->
<pre>
01 var pOdr1: [2]f32 = .{ 0, 0 };          //X, a
02 var pOdr2: [3]f32 = .{ 0, 0, 0 };       //X^2, X, a
03 var pOdr3: [4]f32 = .{ 0, 0, 0, 0 };    //X^3, X^2, X, a
04 const ans1 = getPolyOrder(&pOdr1);
05 const ans2 = getPolyOrder(&pOdr2);
06 const ans3 = getPolyOrder(&pOdr3);
07 const exp1: f32 = 1;
08 const exp2: f32 = 2;
09 const exp3: f32 = 3;
10 try std.testing.expectEqual(ans1, exp1);
11 try std.testing.expectEqual(ans2, exp2);
12 try std.testing.expectEqual(ans3, exp3);
</pre>