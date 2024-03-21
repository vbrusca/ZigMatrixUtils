# Encryption Using Matrices

You can use matrices to encrypt messages. Start by defining an encoding matrix for the English alphabet's capitalized letters, lines 2 - 7. The encoded, clear text matrices are described on lines 9 - 11. The description and variables for the encryption marix are on lines 13 - 18. The return value of matrices, b, is declared on line 19.
<br>
<br>
Variables to hold the letters to encrypt are declared on lines 20 - 24 while variables to hold their encrypted values are declared on lines 26 - 30. The last group of variables, lines 32 - 36, are for decrypting values and checking the results. The encryptiong is based on our ability to multiply by a matrix that has an inverse knowing that we can reverse the process by using the inverse and returning to the original matrix. 

<!-- //"XMTX: ELA - Larson, Edwards: 2.5 Example 5, 6 test" -->
<pre>
001 //Matrix Cryptography
002 //0 = _     6 = F       12 = L      18 = R      24 = X
003 //1 = A     7 = G       13 = M      19 = S      25 = Y
004 //2 = B     8 = H       14 = N      20 = T      26 = Z
005 //3 = C     9 = I       15 = O      21 = U
006 //4 = D     10 = J      16 = P      22 = V
007 //5 = E     11 = K      17 = Q      23 = W
008 
009 //Clear text matrices
010 //[13  5  5] [20  0 13] [5  0 13] [15 14  4] [1 25  0]
011 //  M  E  E    T  _  M   E  _  M    O  N  D   A  Y  _
012 
013 //A = 1  -2   2
014 //   -1   1   3
015 //    1  -1  -4
016 var A: [9]f32 = .{ 1, -2, 2, -1, 1, 3, 1, -1, -4 };
017 var tmp3: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
018 var idt3: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
019 var b: bool = false;
020 var w1: [3]f32 = .{ 13, 5, 5 };
021 var w2: [3]f32 = .{ 20, 0, 13 };
022 var w3: [3]f32 = .{ 5, 0, 13 };
023 var w4: [3]f32 = .{ 15, 14, 4 };
024 var w5: [3]f32 = .{ 1, 25, 0 };
025 
026 var eW1: [3]f32 = .{ 0, 0, 0 };
027 var eW2: [3]f32 = .{ 0, 0, 0 };
028 var eW3: [3]f32 = .{ 0, 0, 0 };
029 var eW4: [3]f32 = .{ 0, 0, 0 };
030 var eW5: [3]f32 = .{ 0, 0, 0 };
031 
032 var cW1: [3]f32 = .{ 0, 0, 0 };
033 var cW2: [3]f32 = .{ 0, 0, 0 };
034 var cW3: [3]f32 = .{ 0, 0, 0 };
035 var cW4: [3]f32 = .{ 0, 0, 0 };
036 var cW5: [3]f32 = .{ 0, 0, 0 };
037
038 //Get inverse matrix
039 var sclr: f32 = 0.0;
040 b = rdcXmtx(&A, 3, false, &tmp3, true, &idt3, 3, false, &sclr);
041 try std.testing.expectEqual(true, b);
042 
043 prntNl();
044 std.debug.print("XMTX: A:\n", .{});
045 prntXmtx(&A, 3);
046 
047 //XMTX: A:
048 //0: x: 1.0e+00 y: -2.0e+00 z: 2.0e+00
049 //1: x: -1.0e+00 y: 1.0e+00 z: 3.0e+00
050 //2: x: 1.0e+00 y: -1.0e+00 z: -4.0e+00
051 
052 prntNl();
053 std.debug.print("XMTX: TMP3:\n", .{});
054 prntXmtx(&tmp3, 3);
055 
056 //XMTX: TMP3:
057 //0: x: 1.0e+00 y: 0.0e+00 z: 0.0e+00
058 //1: x: -0.0e+00 y: 1.0e+00 z: 0.0e+00
059 //2: x: -0.0e+00 y: -0.0e+00 z: 1.0e+00 
060 
061 prntNl();
062 std.debug.print("XMTX: IDT3:\n", .{});
063 prntXmtx(&idt3, 3);
064 
065 //XMTX: IDT3:
066 //0: x: -1.0e+00 y: -1.0e+01 z: -8.0e+00
067 //1: x: -1.0e+00 y: -6.0e+00 z: -5.0e+00
068 //2: x: -0.0e+00 y: -1.0e+00 z: -1.0e+00
069 
070 A = .{ 1, -2, 2, -1, 1, 3, 1, -1, -4 };
071 
</pre>

At the end of the first section we have the matrix to multiply by, <b>A</b>, and its inverse, <b>idt3</b>.
In the next section we'll be encrypting the letters one word at a time for the 5 example words we have. We'll go over one example here. The first word is encrypted on line 73 by multiplying it by the matrix <b>A</b>, storing the results in the <b>eW1</b> matrix.

<pre>
072 //Encode word 1
073 b = tmsXmtx(&w1, 3, &A, 3, &eW1, 3);
074 try std.testing.expectEqual(true, b);
075 
076 prntNl();
077 std.debug.print("XMTX: W1:\n", .{});
078 prntXmtx(&w1, 3);
079 
080 //XMTX: W1:
081 //0: x: 1.3e+01 y: 5.0e+00 z: 5.0e+00
082 
083 prntNl();
084 std.debug.print("XMTX: ENC_W1:\n", .{});
085 prntXmtx(&eW1, 3);
086 
087 //XMTX: ENC_W1:
088 //0: x: 1.3e+01 y: -2.6e+01 z: 2.1e+01
089 
090 //Encode word 2
091 b = tmsXmtx(&w2, 3, &A, 3, &eW2, 3);
092 try std.testing.expectEqual(true, b);
093 
094 prntNl();
095 std.debug.print("XMTX: W2:\n", .{});
096 prntXmtx(&w2, 3);
097 
098 //XMTX: W2:
099 //0: x: 2.0e+01 y: 0.0e+00 z: 1.3e+01
100 
101 prntNl();
102 std.debug.print("XMTX: ENC_W2:\n", .{});
103 prntXmtx(&eW2, 3);
104 
105 //XMTX: ENC_W2:
106 //0: x: 3.3e+01 y: -5.3e+01 z: -1.2e+01 
107 
108 //Encode word 3
109 b = tmsXmtx(&w3, 3, &A, 3, &eW3, 3);
110 try std.testing.expectEqual(true, b);
111 
112 prntNl();
113 std.debug.print("XMTX: W3:\n", .{});
114 prntXmtx(&w3, 3);
115 
116 //XMTX: W3:
117 //0: x: 5.0e+00 y: 0.0e+00 z: 1.3e+01
118 
119 prntNl();
120 std.debug.print("XMTX: ENC_W3:\n", .{});
121 prntXmtx(&eW3, 3);
122 
123 //XMTX: ENC_W3:
124 //0: x: 1.8e+01 y: -2.3e+01 z: -4.2e+01
125 
126 //Encode word 4
127 b = tmsXmtx(&w4, 3, &A, 3, &eW4, 3);
128 try std.testing.expectEqual(true, b);
129 
130 prntNl();
131 std.debug.print("XMTX: W4:\n", .{});
132 prntXmtx(&w4, 3);
133 
134 //XMTX: W4:
135 //0: x: 1.5e+01 y: 1.4e+01 z: 4.0e+00
136 
137 prntNl();
138 std.debug.print("XMTX: ENC_W4:\n", .{});
139 prntXmtx(&eW4, 3);
140 
141 //XMTX: ENC_W4:
142 //0: x: 5.0e+00 y: -2.0e+01 z: 5.6e+01
143 
144 //Encode word 5
145 b = tmsXmtx(&w5, 3, &A, 3, &eW5, 3);
146 try std.testing.expectEqual(true, b);
147 
148 prntNl();
149 std.debug.print("XMTX: W5:\n", .{});
150 prntXmtx(&w5, 3);
151 
152 //XMTX: W5:
153 //0: x: 1.0e+00 y: 2.5e+01 z: 0.0e+00
154 
155 prntNl();
156 std.debug.print("XMTX: ENC_W5:\n", .{});
157 prntXmtx(&eW5, 3);
158 
159 //XMTX: ENC_W5:
160 //0: x: -2.4e+01 y: 2.3e+01 z: 7.7e+01
161 
</pre>

The reversal of the encryption process requires us to multiply the encoded matrix by the inverse of the encrypting matrix, <b>A</b>, that's the matrix <b>idt3</b>, with the results being stored in the <b>cW1 - cW5</b> matrices as shown below.

<pre>
162 //Decode word 1
163 b = tmsXmtx(&eW1, 3, &idt3, 3, &cW1, 3);
164 try std.testing.expectEqual(true, b);
165 
166 prntNl();
167 std.debug.print("XMTX: DEC_W1:\n", .{});
168 prntXmtx(&cW1, 3);
169 
170 //XMTX: DEC_W1:
171 //0: x: 1.3e+01 y: 5.0e+00 z: 5.0e+00
172 
173 //Decode word 2
174 b = tmsXmtx(&eW2, 3, &idt3, 3, &cW2, 3);
175 try std.testing.expectEqual(true, b);
176 
177 prntNl();
178 std.debug.print("XMTX: DEC_W2:\n", .{});
179 prntXmtx(&cW2, 3);
180 
181 //XMTX: DEC_W2:
182 //0: x: 2.0e+01 y: 0.0e+00 z: 1.3e+01
183 
184 //Decode word 3
185 b = tmsXmtx(&eW3, 3, &idt3, 3, &cW3, 3);
186 try std.testing.expectEqual(true, b);
187 
188 prntNl();
189 std.debug.print("XMTX: DEC_W3:\n", .{});
190 prntXmtx(&cW3, 3);
191 
192 //XMTX: DEC_W3:
193 //0: x: 5.0e+00 y: 0.0e+00 z: 1.3e+01
194 
195 //Decode word 4
196 b = tmsXmtx(&eW4, 3, &idt3, 3, &cW4, 3);
197 try std.testing.expectEqual(true, b);
198 
199 prntNl();
200 std.debug.print("XMTX: DEC_W4:\n", .{});
201 prntXmtx(&cW4, 3);
202 
203 //XMTX: DEC_W4:
204 //0: x: 1.5e+01 y: 1.4e+01 z: 4.0e+00
205 
206 //Decode word 5
207 b = tmsXmtx(&eW5, 3, &idt3, 3, &cW5, 3);
208 try std.testing.expectEqual(true, b);
209 
210 prntNl();
211 std.debug.print("XMTX: DEC_W5:\n", .{});
212 prntXmtx(&cW5, 3);
213 
214 //XMTX: DEC_W5:
215 //0: x: 1.0e+00 y: 2.5e+01 z: 0.0e+00
216 
217 try std.testing.expectEqual(true, equXmtx(&w1, &cW1));
218 try std.testing.expectEqual(true, equXmtx(&w2, &cW2));
219 try std.testing.expectEqual(true, equXmtx(&w3, &cW3));
220 try std.testing.expectEqual(true, equXmtx(&w4, &cW4));
221 try std.testing.expectEqual(true, equXmtx(&w5, &cW5));
</pre>

The process ends with checks to ensure the decrypted value matches the original clear text on lines 217 - 221.