#Difference of Vectors 1

To find the difference between two vectors use the <b>diff1Xvec</b> function
which takes the difference of <b>m2 - m1</b> and stores the result in <b>m1</b>.

//"XMTX: ELA - Larson, Edwards: 2.1 Example 3 test"
<pre>
01 var m1: [9]f32 = .{ 3, 6, 12, -9, 0, -3, 6, 3, 6 };
02 var m2: [9]f32 = .{ 2, 0, 0, 1, -4, 3, -1, 3, 2 };
03 diff1Xvec(&m1, &m2);
04 var exp1: [9]f32 = .{ 1, 6, 12, -10, 4, -6, 7, 0, 4 };
05 clnXmtx(&m1);
06 const b: bool = equXmtx(&exp1, &m1);
</pre>

Notice that we clean the resulting matrix, <b>m1</b>, before using it as an argument in the <b>equXmtx</b> function.