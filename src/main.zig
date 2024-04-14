const std = @import("std");
const xmu = @import("./XmtxUtils.zig");

pub fn main() !void {
    // Prints to stderr (it's a shortcut based on `std.io.getStdErr()`)
    std.debug.print("All your {s} are belong to us.\n", .{"codebase"});

    // stdout is for the actual output of your application, for example if you
    // are implementing gzip, then only the compressed bytes should be sent to
    // stdout, not any debugging messages.
    const stdout_file = std.io.getStdOut().writer();
    var bw = std.io.bufferedWriter(stdout_file);
    const stdout = bw.writer();

    try stdout.print("Run `zig build test` to run the tests.\n", .{});

    var v1: [3]f32 = .{ 1, 2, 3 };
    var v2: [3]f32 = .{ 0, 0, 0 };
    var v3: [3]f32 = .{ 0, 0, 0 };
    xmu.cpyXvec(&v1, &v2);
    xmu.cpyXvec(&v1, &v3);
    _ = xmu.equXvec(&v1, &v2);
    _ = xmu.equXvec(&v1, &v3);

    try bw.flush(); // don't forget to flush!
}

test "simple test" {
    var v1: [3]f32 = .{ 1, 2, 3 };
    var v2: [3]f32 = .{ 0, 0, 0 };
    var v3: [3]f32 = .{ 1, 1, 1 };
    xmu.cpyXvec(&v1, &v2);
    xmu.cpyXvec(&v1, &v3);
    std.debug.print("simple test:\n", .{});
    try std.testing.expectEqual(true, xmu.equXvec(&v1, &v2));
    try std.testing.expectEqual(true, xmu.equXvec(&v1, &v3));
}

//
//
//
//
//
//
//
//
//
//
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//START EXTENDED PROBLEMS
//Elementary Linear Algebra - Larson, Edwards- 4th Edition
//Chapters 5.2 - 7.3
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

fn _diffInrPrdctR2(u: []f32, v: []f32) f32 {
    return ((u[0] * v[0]) + (2.0 * (u[1] * v[1])));
}

fn _notInrPrdctR3(u: []f32, v: []f32) f32 {
    return ((u[0] * v[0]) - (2.0 * (u[1] * v[1])) + (u[2] * v[2]));
}

fn _diffInrPrdctM22(u: []f32, v: []f32) f32 {
    //u_idx = |0, 1|   |u_11, u_12|
    //        |2, 3|   |u_21, u_22|
    return ((u[0] * v[0]) + (u[2] * v[2]) + (u[1] * v[1]) + (u[3] * v[3]));
}

fn _diffInrPrdctP2(p: []f32, q: []f32) f32 {
    if (p.len != q.len) {
        xmu.prntNlStr("_diffInrPrdctP2: Error, argument vectors p and q must be the same length.");
        return std.math.floatMin(f32);
    }

    var i: usize = 0;
    const l: usize = p.len;
    var ret: f32 = 0.0;
    while (i < l) : (i += 1) {
        ret += (p[i] * q[i]);
    }
    return ret;
}

fn _diffInrPrdctR3(u: []f32, v: []f32) f32 {
    return ((u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2]));
}

test "XMTX: ELA - Larson, Edwards: 5.2 Example 1 test" {
    //Example 1: pg 264

    var mtx3u: [3]f32 = .{ 1, 0, 0 };
    var mtx3v: [3]f32 = .{ 0, 1, 0 };
    var mtx3w: [3]f32 = .{ 0, 0, 1 };
    const mtx3c: f32 = 5.0;
    const alloc = std.testing.allocator;
    var b: bool = false;
    var exp: bool = true;

    b = try xmu.isInrPrdctSpc(&mtx3u, &mtx3v, &mtx3w, mtx3c, xmu.dotPrdXvec, &alloc);
    try std.testing.expectEqual(exp, b);
    xmu.prntNl();

    mtx3u = .{ 1, 0, 0 };
    mtx3v = .{ 1, 0, 0 };
    mtx3w = .{ 1, 0, 0 };
    b = false;
    exp = true;

    b = try xmu.isInrPrdctSpc(&mtx3u, &mtx3v, &mtx3w, mtx3c, xmu.dotPrdXvec, &alloc);
    try std.testing.expectEqual(exp, b);
    xmu.prntNl();

    var mtx4u: [4]f32 = .{ 1, 0, 0, 0 };
    var mtx4v: [4]f32 = .{ 0, 1, 0, 0 };
    var mtx4w: [4]f32 = .{ 0, 0, 1, 0 };
    const mtx4c: f32 = 7.0;

    b = try xmu.isInrPrdctSpc(&mtx4u, &mtx4v, &mtx4w, mtx4c, xmu.dotPrdXvec, &alloc);
    try std.testing.expectEqual(exp, b);
    xmu.prntNl();

    mtx4u = .{ 1, 0, 0, 0 };
    mtx4v = .{ 1, 0, 0, 0 };
    mtx4w = .{ 1, 0, 0, 0 };
    b = false;
    exp = true;

    b = try xmu.isInrPrdctSpc(&mtx4u, &mtx4v, &mtx4w, mtx4c, xmu.dotPrdXvec, &alloc);
    try std.testing.expectEqual(exp, b);
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.2 Example 2 test" {
    //Example 2: pg 265

    var mtx2u: [2]f32 = .{ 1, 0 };
    var mtx2v: [2]f32 = .{ 0, 1 };
    var mtx2w: [2]f32 = .{ 1, 1 };
    const mtx2c: f32 = 5.0;
    const alloc = std.testing.allocator;
    var b: bool = false;
    const exp: bool = true;

    b = try xmu.isInrPrdctSpc(&mtx2u, &mtx2v, &mtx2w, mtx2c, _diffInrPrdctR2, &alloc);
    try std.testing.expectEqual(exp, b);
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.2 Example 3 test" {
    //Example 3: pg 265

    var mtx3u: [3]f32 = .{ 1, 0, 0 };
    var mtx3v: [3]f32 = .{ 0, 1, 0 };
    var mtx3w: [3]f32 = .{ 0, 0, 1 };
    const mtx3c: f32 = 5.0;
    const alloc = std.testing.allocator;
    var b: bool = false;
    const exp: bool = false;

    b = try xmu.isInrPrdctSpc(&mtx3u, &mtx3v, &mtx3w, mtx3c, _notInrPrdctR3, &alloc);
    try std.testing.expectEqual(exp, b);
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.2 Example 4 test" {
    //Example 4: pg 266
    //An inner product on M 2x2 represented by V 4

    var mtx4u: [4]f32 = .{ 1, 0, 0, 0 };
    var mtx4v: [4]f32 = .{ 0, 1, 0, 0 };
    var mtx4w: [4]f32 = .{ 0, 0, 1, 0 };
    const mtx4c: f32 = 13.0;
    const alloc = std.testing.allocator;
    var b: bool = false;
    const exp: bool = true;

    b = try xmu.isInrPrdctSpc(&mtx4u, &mtx4v, &mtx4w, mtx4c, _diffInrPrdctM22, &alloc);
    try std.testing.expectEqual(exp, b);
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.2 Theorem 5.7 test" {
    //Theorem 5.7: pg 267
    var mtx3u: [3]f32 = .{ 1, 0, 0 };
    var mtx3v: [3]f32 = .{ 0, 1, 0 };
    var mtx3w: [3]f32 = .{ 0, 0, 1 };
    var mtx3t: [3]f32 = .{ 0, 0, 0 };
    const exp: f32 = 0.0;

    //Let u,v, w be vectors in an inner product space V, and let c be any real number.
    //1. <0, v> = <v, 0> = 0
    const v1 = xmu.dotPrdXvec3(&mtx3t, &mtx3v);
    const v2 = xmu.dotPrdXvec3(&mtx3v, &mtx3t);
    try std.testing.expectEqual(v1, v2);
    try std.testing.expectEqual(exp, v1);
    try std.testing.expectEqual(exp, v2);
    xmu.prntNl();

    //2. <u + v, w> = <u,w> + <v,w>
    xmu.clrXvec(&mtx3t);
    xmu.sum2Xvec(&mtx3t, &mtx3u, &mtx3v);
    const v3 = xmu.dotPrdXvec3(&mtx3t, &mtx3w);
    const v4 = xmu.dotPrdXvec3(&mtx3u, &mtx3w);
    const v5 = xmu.dotPrdXvec3(&mtx3v, &mtx3w);
    const v6 = (v4 + v5);
    try std.testing.expectEqual(v3, v6);
    xmu.prntNl();

    //3. <u,cv> = c<u,v>
    const c: f32 = 15.0;
    xmu.clrXvec(&mtx3t);
    xmu.cpyXvec(&mtx3v, &mtx3t);
    xmu.mulXvec(&mtx3t, c);
    const v7 = xmu.dotPrdXvec3(&mtx3u, &mtx3t);
    const v8 = xmu.dotPrdXvec3(&mtx3u, &mtx3v);
    const v9 = (c * v8);
    try std.testing.expectEqual(v7, v9);
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.2 Def of Norm, Distance, Angle test" {
    //Definition of Norm, Distance, and Angle: pg 267
    var mtx2u: [2]f32 = .{ 1, 0 };
    var mtx2v: [2]f32 = .{ 0, 1 };
    var mtx2t: [2]f32 = .{ 0, 0 };
    var exp: f32 = ((mtx2u[0] * mtx2u[0]) + (2.0 * (mtx2u[1] * mtx2u[1])));
    var b: f32 = -1.0;

    //1. The norm, length, magnitude of u is ||u|| = sqrt(<u,u>)
    // use prdct = ((u[0] * v[0]) + (2.0 * (u[1] * v[1])));
    b = xmu.magInrPrdctSpcXvec(&mtx2u, _diffInrPrdctR2);
    xmu.prntNlStrArgs("Found {} = {}", .{ exp, b });
    try std.testing.expectEqual(exp, b);
    xmu.prntNl();

    //2. The distance between u and v is d(u,v) = ||u - v||
    exp = 1.73205077e+00;
    xmu.diff2Xvec(&mtx2t, &mtx2u, &mtx2v);
    b = xmu.magInrPrdctSpcXvec(&mtx2t, _diffInrPrdctR2);
    xmu.prntNlStrArgs("Found {} = {}", .{ exp, b });
    try std.testing.expectEqual(exp, b);
    xmu.prntNl();

    //3. The angle between two vectors is cos theta = (<u,v> / (||u||||v||))
    exp = 1.57079625e+00;
    b = xmu.aglBtwnInrPrdctSpcXvec(&mtx2u, &mtx2v, _diffInrPrdctR2);
    xmu.prntNlStrArgs("Found {} = {}", .{ exp, b });
    try std.testing.expectEqual(exp, b);
    xmu.prntNl();

    exp = 1.57079625e+00;
    b = xmu.aglBtwnXvec(&mtx2u, &mtx2v);
    xmu.prntNlStrArgs("Found {} = {}", .{ exp, b });
    try std.testing.expectEqual(exp, b);
    xmu.prntNl();

    mtx2u = .{ 3, 4 };
    mtx2v = .{ 5, -12 };
    exp = 2.10330057e+00;
    b = xmu.aglBtwnXvec(&mtx2u, &mtx2v);
    xmu.prntNlStrArgs("Found {} = {}", .{ exp, b });
    try std.testing.expectEqual(exp, b);
    xmu.prntNl();

    //4. u and v are orthogonal <u,v> = 0
    mtx2u = .{ 1, 0 };
    mtx2v = .{ 0, 1 };
    exp = 0;
    b = _diffInrPrdctR2(&mtx2u, &mtx2v);
    xmu.prntNlStrArgs("Found {} = {}", .{ exp, b });
    try std.testing.expectEqual(exp, b);
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.2 Example 6 test" {
    //Example 4: pg 267
    //For polynomials p = A0 + A1*X + .. + An*Xn and q = B0 + B1*X + ... + Bn*Xn
    //<p,q> = A0*B0 + A1B1 + ... + AnBn
    //Let p(x) = 1 - 2*X^2, q(x) = 4 - 2*X + X^2, and r(x) = X + 2*X^2 be polynomials
    //in P2 and determine the following.

    //p(x) = |1,  0, -2|
    //q(x) = |4, -2,  1|
    //r(x) = |0,  1,  2|

    //a. <p,q>
    const p: [3]f32 = .{ 1, 0, -2 };
    const q: [3]f32 = .{ 4, -2, 1 };
    const prdct: *const fn (l: []f32, r: []f32) f32 = _diffInrPrdctP2;
    var v: f32 = prdct(@constCast(&p), @constCast(&q));
    var exp: f32 = 2.0;
    try std.testing.expectEqual(exp, v);
    xmu.prntNl();

    //b. <q,r>
    const r: [3]f32 = .{ 0, 1, 2 };
    v = prdct(@constCast(&q), @constCast(&r));
    exp = 0.0;
    try std.testing.expectEqual(exp, v);
    xmu.prntNl();

    //c. ||q||
    v = xmu.magInrPrdctSpcXvec(@constCast(&q), prdct);
    exp = std.math.sqrt(21.0);
    try std.testing.expectEqual(exp, v);
    xmu.prntNl();

    //d. d(p,q) = ||p - q||;
    var t: [3]f32 = .{ 0, 0, 0 };
    xmu.diff2Xvec(@constCast(&t), @constCast(&p), @constCast(&q));
    v = xmu.magInrPrdctSpcXvec(&t, prdct);
    exp = std.math.sqrt(22.0);
    try std.testing.expectEqual(exp, v);
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.2 Theorem 5.8 test" {
    //Theorem 5.8: pg 270
    //Let u and v be vectors in an inner product space V.
    //1. Cauchy-Schwarz Inequality: |<u,v>| <= ||u|| * ||v||
    const prdct: *const fn (l: []f32, r: []f32) f32 = _diffInrPrdctR3;
    var u: [3]f32 = .{ 1, 2, 3 };
    var v: [3]f32 = .{ 2, 0, 1 };
    var t: [3]f32 = .{ 0, 0, 0 };
    var v1: f32 = xmu.absF32(prdct(&u, &v));
    var v2: f32 = xmu.magInrPrdctSpcXvec(&u, prdct);
    var v3: f32 = xmu.magInrPrdctSpcXvec(&v, prdct);
    var v4: f32 = (v2 * v3);
    try std.testing.expectEqual(true, (v1 <= v4));
    xmu.prntNl();

    //2. Triangle Inequality: ||u + v|| <= ||u|| + ||v||
    xmu.sum2Xvec(&t, &u, &v);
    v1 = xmu.magInrPrdctSpcXvec(&t, prdct);
    v4 = (v2 + v3);
    try std.testing.expectEqual(true, (v1 <= v4));
    xmu.prntNl();

    //3. Pythagorean Theorem: u and v are orthogonal if and only if ||u + v||^2 = ||u||^2 + ||v||^2
    v1 = (v1 * v1);
    v4 = (v2 * v2) + (v3 * v3);
    try std.testing.expectEqual(true, (v1 != v4));
    xmu.prntNl();

    //4.
    u = .{ 1, 0, 0 };
    v = .{ 0, 1, 0 };
    t = .{ 0, 0, 0 };
    xmu.sum2Xvec(&t, &u, &v);
    v1 = xmu.magInrPrdctSpcXvec(&t, prdct);
    v2 = xmu.magInrPrdctSpcXvec(&u, prdct);
    v3 = xmu.magInrPrdctSpcXvec(&v, prdct);
    v1 = (v1 * v1);
    v4 = (v2 * v2) + (v3 * v3);

    var w: [2]f32 = .{ v1, v4 };
    xmu.clnXvec(&w);
    xmu.prntNlStrArgs("Found: v1: {} v4: {}", .{ v1, v4 });
    xmu.prntNlStrArgs("Found: w[0]: {} w[1]: {}", .{ w[0], w[1] });
    try std.testing.expectEqual(true, (w[0] == w[1]));
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.2 Example 9 test" {
    //Example 9: pg 271
    //Finding the orthogonal projection of u onto v
    //In R^2 given u=(4,2) and v=(3,4), find the orthogonal projection of u onto v.
    var u: [2]f32 = .{ 4, 2 };
    var v: [2]f32 = .{ 3, 4 };
    const projUontoV: []f32 = xmu.projXvec_VecP_Onto_VecQ_InrPrdctSpc(&u, &v, xmu.dotPrdXvec);
    var exp: [2]f32 = .{ (12.0 / 5.0), (16.0 / 5.0) };

    xmu.clnXvec(projUontoV);
    xmu.clnXvec(&exp);
    try std.testing.expectEqual(true, xmu.equXvec(projUontoV, &exp));
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.2 Example 10 test" {
    //Example 10: pg 272
    //Finding orthogonal projection in R^3
    //Use the Euclidean inner product in R^3 to find the orthogonal projection of u=(6, 2, 4)
    //onto v=(1, 2, 0)
    const u: [3]f32 = .{ 6, 2, 4 };
    const v: [3]f32 = .{ 1, 2, 0 };
    const exp: [3]f32 = .{ 2, 4, 0 };
    const projUontoV: []f32 = xmu.projXvec_VecP_Onto_VecQ_InrPrdctSpc(@constCast(&u), @constCast(&v), xmu.dotPrdXvec);

    xmu.prntNlStr("Expected:");
    xmu.prntXvec(@constCast(&exp));
    xmu.prntNl();

    xmu.prntNlStr("Projected:");
    xmu.prntXvec(projUontoV);
    xmu.prntNl();

    try std.testing.expectEqual(true, xmu.equXvec(projUontoV, @constCast(&exp)));
    xmu.prntNl();
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//STOP ENTENDED PROBLEMS
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
