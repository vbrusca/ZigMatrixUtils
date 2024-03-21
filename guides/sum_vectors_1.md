#Summation of Vectors

In this section we'll take a look at how to add multiple vectors.
In this case we're adding the values of vectors m1, m2, and m3 together and storing the result in the m3 matrix. 

//"XMTX: ELA - Larson, Edwards: 2.1 Example 2 test"
<pre>
01 var m1: [4]f32 = .{ -1, 2, 0, 1 };
02 var m2: [4]f32 = .{ 1, 3, -1, 2 };
03 var m3: [4]f32 = .{ 0, 0, 0, 0 };
04 sum2Xvec(&m3, &m1, &m2);
05 var exp1: [4]f32 = .{ 0, 5, -1, 3 };
06 const b: bool = equXmtx(&exp1, &m3);
</pre>

Notice that we check to see if the two vectors are equal by using the <b>equXmtx</b> function and thinking of the matrices <b>m3</b> and <b>exp1</b> as two single row vectors.