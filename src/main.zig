//! The main execution point for running unit tests. This isn't used by the library
//! except for unit testing.
//!
//! Developers: Carlo Bruscani, Victor Brusca
//! 04/17/2024
//!

const std = @import("std");
const xmu = @import("./XmtxUtils.zig");

comptime {
    @setFloatMode(std.builtin.FloatMode.optimized);
}

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
    //version 1
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

    //version 2
    var lu2 = [3]f32{ 6, 2, 4 };
    var lv2 = [3]f32{ 1, 2, 0 };
    var lprojUontoV: []f32 = xmu.projXvec_VecP_Onto_VecQ_InrPrdctSpc(&lu2, &lv2, xmu.dotPrdXvec);
    try std.testing.expectEqual(true, xmu.equXvec(lprojUontoV, @constCast(&exp)));
    xmu.prntNl();

    lprojUontoV = xmu.projXvec_VecP_Onto_VecQ_InrPrdctSpc(lu2[0..3], lv2[0..3], xmu.dotPrdXvec);
    try std.testing.expectEqual(true, xmu.equXvec(lprojUontoV, @constCast(&exp)));
    xmu.prntNl();

    const llu2: []f32 = lu2[0..3];
    const llv2: []f32 = lv2[0..3];
    lprojUontoV = xmu.projXvec_VecP_Onto_VecQ_InrPrdctSpc(llu2, llv2, xmu.dotPrdXvec);
    try std.testing.expectEqual(true, xmu.equXvec(lprojUontoV, @constCast(&exp)));
    xmu.prntNl();

    var lllu2: []f32 = lu2[0..3];
    var lllv2: []f32 = lv2[0..3];
    lprojUontoV = xmu.projXvec_VecP_Onto_VecQ_InrPrdctSpc(lllu2, lllv2, xmu.dotPrdXvec);
    try std.testing.expectEqual(true, xmu.equXvec(lprojUontoV, @constCast(&exp)));
    xmu.prntNl();

    lllu2 = lu2[0..3];
    lllv2 = lv2[0..3];
}

test "XMTX: ELA - Larson, Edwards: 5.2 Theorem 5.9 test" {
    //Theorem 5.9: pg 273
    //Let u and v be two vectors in an inner product space V, such that v != 0 vector. Then...
    //d(u, projUontoV) < d(u, c*v), iff c != <u,v>/<v,v>

    var u = [3]f32{ 6, 2, 4 };
    var v = [3]f32{ 1, 2, 0 };
    var cv = [3]f32{ 1, 2, 0 };
    const n = xmu.dotPrdXvec(&u, &v);
    const d = xmu.dotPrdXvec(&v, &v);
    const c = n / d;

    xmu.prntNlStrArgs("Found c: {}", .{c});
    const alloc = std.testing.allocator;
    const projUontoV = xmu.projXvec_VecP_Onto_VecQ_InrPrdctSpc(&u, &v, xmu.dotPrdXvec);
    xmu.mulXvec(&cv, c);

    xmu.prntNlStr("Found cv:");
    xmu.prntXmtxNl(&cv, 3);
    const v1 = try xmu.dstInrPrdctSpcXvec(&u, @constCast(projUontoV), xmu.dotPrdXvec, &alloc);
    const v2 = try xmu.dstInrPrdctSpcXvec(&u, &cv, xmu.dotPrdXvec, &alloc);

    xmu.prntNlStrArgs("Found v1: {}", .{v1});
    xmu.prntNlStrArgs("Found v2: {}", .{v2});
    try std.testing.expectEqual(v1, v2);
    xmu.prntNl();
}

fn _diffInrPrdPrdctProblem3(u: []f32, v: []f32) f32 {
    return ((3 * u[0] * v[0]) + (u[1] * v[1]));
}

fn _diffInrPrdPrdctProblem7(u: []f32, v: []f32) f32 {
    return ((u[0] * v[0]) + (2.0 * u[1] * v[1]) + (u[2] * v[2]));
}

test "XMTX: ELA - Larson, Edwards: 5.2 Problem 1, 3, 5, 7 test" {
    //Chapter 5: Section 5.2: Problem 1, 3, 5, 7: pg 274
    //Find (a) <u,v>, (b) ||u||, and (c) d(u, v) for the given inner product defined in R^n.
    const alloc = std.testing.allocator;

    //1: u=(3, 4), v=(5, -12), <u,v>=(dotPrdXvec(u, v))
    //Find: (a) <u,v>, (b) ||u||, (c) d(u, v)
    //Exp: (a) -33, (b) 5, (c) 2 * sqrt(65)
    //(a)
    var u: [2]f32 = .{ 3, 4 };
    var v: [2]f32 = .{ 5, -12 };
    var val: f32 = xmu.inrPrdct(&u, &v, xmu.dotPrdXvec);
    var exp: f32 = -33.0;
    xmu.prntNlStrArgs("1a Found val: {}", .{val});
    xmu.prntNlStrArgs("1a Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //(b)
    val = xmu.magInrPrdctSpcXvec(&u, xmu.dotPrdXvec);
    exp = 5.0;
    xmu.prntNlStrArgs("1b Found val: {}", .{val});
    xmu.prntNlStrArgs("1b Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //(c)
    val = try xmu.dstInrPrdctSpcXvec(&u, &v, xmu.dotPrdXvec, &alloc);
    exp = (2.0 * std.math.sqrt(65.0));
    xmu.prntNlStrArgs("1c Found val: {}", .{val});
    xmu.prntNlStrArgs("1c Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //3: u=(-4, 3), v=(0,5), <u,v>=((3*u1*v1) + (u2*v2))
    //Find: (a) <u,v>, (b) ||u||, and (c) d(u, v)
    //Exp: (a)-15, (b) sqrt(57), (c) 2*sqrt(13)
    //(a)
    u = [2]f32{ -4, 3 };
    v = [2]f32{ 0, 5 };
    val = xmu.inrPrdct(&u, &v, _diffInrPrdPrdctProblem3);
    exp = 15.0;
    xmu.prntNlStrArgs("3a Found val: {}", .{val});
    xmu.prntNlStrArgs("3a Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //(b)
    val = xmu.magInrPrdctSpcXvec(&u, _diffInrPrdPrdctProblem3);
    exp = std.math.sqrt(57.0);
    xmu.prntNlStrArgs("3b Found val: {}", .{val});
    xmu.prntNlStrArgs("3b Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //(c)
    val = try xmu.dstInrPrdctSpcXvec(&u, &v, _diffInrPrdPrdctProblem3, &alloc);
    exp = (2.0 * std.math.sqrt(13.0));
    xmu.prntNlStrArgs("3c Found val: {}", .{val});
    xmu.prntNlStrArgs("3c Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //5. u=(0,9,4), v=(9,-2,-4), <u,v>=(dotPrdXvec(u,v))
    //Find: (a) <u,v>, (b) ||u||, (c) d(u, v)
    //Exp: (a) -34, (b) sqrt(97), (c) sqrt(266)
    //(a)
    var lu3: [3]f32 = .{ 0, 9, 4 };
    var lv3: [3]f32 = .{ 9, -2, -4 };
    val = xmu.inrPrdct(&lu3, &lv3, xmu.dotPrdXvec);
    exp = -34.0;
    xmu.prntNlStrArgs("5a Found val: {}", .{val});
    xmu.prntNlStrArgs("5a Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //(b)
    val = xmu.magInrPrdctSpcXvec(&lu3, xmu.dotPrdXvec);
    exp = std.math.sqrt(97.0);
    xmu.prntNlStrArgs("5b Found val: {}", .{val});
    xmu.prntNlStrArgs("5b Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //(c)
    val = try xmu.dstInrPrdctSpcXvec(&lu3, &lv3, xmu.dotPrdXvec, &alloc);
    exp = std.math.sqrt(266.0);
    xmu.prntNlStrArgs("5c Found val: {}", .{val});
    xmu.prntNlStrArgs("5c Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //7: u=(1,1,1), v=(2,5,2), <u,v>=((u1*v1) + (2*u2*v2) + (u3*v3))
    //Find: (a) <u,v>, (b) ||u||, and (c) d(u, v)
    //Exp: (a)14, (b) 2, (c) sqrt(34)
    //(a)
    lu3 = [3]f32{ 1, 1, 1 };
    lv3 = [3]f32{ 2, 5, 2 };
    val = xmu.inrPrdct(&lu3, &lv3, _diffInrPrdPrdctProblem7);
    exp = 14.0;
    xmu.prntNlStrArgs("7a Found val: {}", .{val});
    xmu.prntNlStrArgs("7a Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //(b)
    val = xmu.magInrPrdctSpcXvec(&lu3, _diffInrPrdPrdctProblem7);
    exp = 2.0;
    xmu.prntNlStrArgs("7b Found val: {}", .{val});
    xmu.prntNlStrArgs("7b Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //(c)
    val = try xmu.dstInrPrdctSpcXvec(&lu3, &lv3, _diffInrPrdPrdctProblem7, &alloc);
    exp = std.math.sqrt(34.0);
    xmu.prntNlStrArgs("7c Found val: {}", .{val});
    xmu.prntNlStrArgs("7c Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();
}

fn _diffInrPrdPrdctProblem13(u: []f32, v: []f32) f32 {
    return ((2.0 * u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2]) + (2.0 * u[3] * v[3]));
}

test "XMTX: ELA - Larson, Edwards: 5.2 Problem 13 test" {
    //Chapter 5: Section 5.2: Problem 13: pg 274
    //<a,b>=((2.0*A11*B11) + (A12*B12) + (2.0*A22*B22))

    //13: A=|-1  3|     B=|0 -2|
    //      | 4 -2|       |1  1|
    //Find: (a)<A,B>, (b)||A||, (c)d(A,B)
    //Exp: (a) -6, (b) sqrt(35.0), (c) 3.0*sqrt(6.0)
    //(a)
    var la = [4]f32{ -1, 3, 4, -2 };
    var lb = [4]f32{ 0, -2, 1, 1 };
    var val: f32 = xmu.inrPrdct(&la, &lb, _diffInrPrdPrdctProblem13);
    var exp: f32 = -6.0;
    var alloc = std.testing.allocator;
    xmu.prntNlStrArgs("13a Found val: {}", .{val});
    xmu.prntNlStrArgs("13a Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //(b)
    val = xmu.magInrPrdctSpcXvec(&la, _diffInrPrdPrdctProblem13);
    exp = std.math.sqrt(35.0);
    xmu.prntNlStrArgs("13b Found val: {}", .{val});
    xmu.prntNlStrArgs("13b Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //(c)
    val = try xmu.dstInrPrdctSpcXvec(&la, &lb, _diffInrPrdPrdctProblem13, &alloc);
    exp = (3.0 * std.math.sqrt(6.0));
    xmu.prntNlStrArgs("13c Found val: {}", .{val});
    xmu.prntNlStrArgs("13c Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.2 Problem 15 test" {
    //Chapter 5: Section 5.2: Problem 15: pg 274
    //<p,q> = a0b0 + a1b1 + a2b2

    //15: p(x)=(1 - x + (3.0 * x^2)), q(x)=(x - x^2)
    //Find: (a)<p,q>, (b)||p||, (c)d(p,q)
    //Exp: (a) -4, (b) sqrt(11), (c) sqrt(21)
    //(a)
    var p = [3]f32{ 1, -1, 3 };
    var q = [3]f32{ 0, 1, -1 };
    var val = xmu.inrPrdct(&p, &q, xmu.dotPrdXvec);
    var exp: f32 = -4;
    var alloc = std.testing.allocator;

    xmu.prntNlStrArgs("15a Found val: {}", .{val});
    xmu.prntNlStrArgs("15a Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //(b)
    val = xmu.magInrPrdctSpcXvec(&p, xmu.dotPrdXvec);
    exp = std.math.sqrt(11.0);
    xmu.prntNlStrArgs("15b Found val: {}", .{val});
    xmu.prntNlStrArgs("15b Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //(c)
    val = try xmu.dstInrPrdctSpcXvec(&p, &q, xmu.dotPrdXvec, &alloc);
    exp = std.math.sqrt(21.0);
    xmu.prntNlStrArgs("15c Found val: {}", .{val});
    xmu.prntNlStrArgs("15c Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();
}

fn _diffInrPrdPrdctProblem15(u: []f32, v: []f32) f32 {
    return ((u[0] * v[0]) - (2.0 * u[1] * v[1]) + (u[2] * v[2]));
}

test "XMTX: ELA - Larson, Edwards: 5.2 Problem 17 test" {
    //Chapter 5: Section 5.2: Problem 17: pg 274
    //Prove that the given function is an inner product.

    //17: <u,v>=(u1 * v1) - (2.0 * u2 * v2) + (u3 * v3)
    //Find: if <u,v> is an inner product space
    //Exp: false
    const u = [3]f32{ 1, 0, 0 };
    const v = [3]f32{ 0, 1, 0 };
    const w = [3]f32{ 0, 0, 1 };
    const c: f32 = 5.0;
    const exp: bool = false;
    const alloc = std.testing.allocator;
    const val: bool = try xmu.isInrPrdctSpc(@constCast(&u), @constCast(&v), @constCast(&w), c, _diffInrPrdPrdctProblem15, &alloc);
    xmu.prntNlStrArgs("17 Found val: {}", .{val});
    xmu.prntNlStrArgs("17 Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.2 Problem 19 test" {
    //Chapter 5: Section 5.2: Problem 19: pg 274
    //Prove that the given function is an inner product.

    //19: <a,b>=((2.0*A11*B11) + (A12*B12) + (2.0*A22*B22))
    //From problem 13.
    //13: A=|-1  3|     B=|0 -2|
    //      | 4 -2|       |1  1|
    var u = [4]f32{ -1, 3, 4, -2 };
    var v = [4]f32{ 0, -2, 1, 1 };
    var w = [4]f32{ 0, 0, 0, 0 };
    const c: f32 = 7.0;
    const exp: bool = true;
    const alloc = std.testing.allocator;
    const val: bool = try xmu.isInrPrdctSpc(@constCast(&u), @constCast(&v), @constCast(&w), c, _diffInrPrdPrdctProblem13, &alloc);
    xmu.prntNlStrArgs("19 Found val: {}", .{val});
    xmu.prntNlStrArgs("19 Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.2 Problem 31 test" {
    //Chapter 5: Section 5.2: Problem 31: pg 274
    //31:
    //(a) Verify the Cauchy-Schwarz Inequality.
    //(b) Verify the Triangle Inquality.

    //(a)
    var u: [2]f32 = .{ 5, 12 };
    var v: [2]f32 = .{ 3, 4 };
    const prdct: *const fn (l: []f32, r: []f32) f32 = xmu.dotPrdXvec;
    const alloc = std.testing.allocator;
    var exp: bool = true;
    var val: bool = xmu.tstCauchySchwarzIneq(&u, &v, prdct);
    xmu.prntNlStrArgs("31a Found val: {}", .{val});
    xmu.prntNlStrArgs("31a Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //(b)
    exp = true;
    val = try xmu.tstTriangleIneq(&u, &v, prdct, &alloc);
    xmu.prntNlStrArgs("31b Found val: {}", .{val});
    xmu.prntNlStrArgs("31b Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.2 Problem 33 test" {
    //Chapter 5: Section 5.2: Problem 33: pg 274
    //33:
    //(a) Verify the Cauchy-Schwarz Inequality.
    //(b) Verify the Triangle Inquality.
    //<p,q> = a0*b0 + a1*b1 + a2*b2 = dotPrdctXvec
    //p(x) = 2x, q(x) = 3.0 * x^2 + 1.0

    //(a)
    var u: [3]f32 = .{ 0, 2, 0 };
    var v: [3]f32 = .{ 1, 0, 3 };
    const prdct: *const fn (l: []f32, r: []f32) f32 = xmu.dotPrdXvec;
    const alloc = std.testing.allocator;
    var exp: bool = true;
    var val: bool = xmu.tstCauchySchwarzIneq(&u, &v, prdct);
    xmu.prntNlStrArgs("33a Found val: {}", .{val});
    xmu.prntNlStrArgs("33a Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //(b)
    exp = true;
    val = try xmu.tstTriangleIneq(&u, &v, prdct, &alloc);
    xmu.prntNlStrArgs("33b Found val: {}", .{val});
    xmu.prntNlStrArgs("33b Found exp: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.2 Problem 41 test" {
    //Chapter 5: Section 5.2: Problem 41: pg 275
    //Find: (a) projUontoV, (b) projVontoU

    //(a)
    var u: [2]f32 = .{ 1, 2 };
    var v: [2]f32 = .{ 2, 1 };
    var val: [2]f32 = xmu.projXvec_VecP_Onto_VecQ(@constCast(&u), @constCast(&v))[0..2].*;
    var exp: [2]f32 = .{ (8.0 / 5.0), (4.0 / 5.0) };
    xmu.prntNlStrArgs("41a Found val: {any}", .{val});
    xmu.prntNlStrArgs("41a Found exp: {any}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();

    //(b)
    u = [2]f32{ 1, 2 };
    v = [2]f32{ 2, 1 };
    val = xmu.projXvec_VecP_Onto_VecQ(@constCast(&v), @constCast(&u))[0..2].*;
    exp = .{ (4.0 / 5.0), (8.0 / 5.0) };
    xmu.prntNlStrArgs("41b Found val: {any}", .{val});
    xmu.prntNlStrArgs("41b Found exp: {any}", .{exp});
    try std.testing.expectEqual(val, exp);
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.3 Example 1 test" {
    //Chapter 5: Section 5.3: Example 1: pg 276
    const cols: usize = 3.0;
    const alloc = std.testing.allocator;
    var S: [9]f32 = .{ (1.0 / std.math.sqrt(2.0)), (1.0 / std.math.sqrt(2.0)), (0.0), (std.math.sqrt(2.0) / -6.0), (std.math.sqrt(2.0) / 6.0), ((2.0 * std.math.sqrt(2.0)) / 3.0), (2.0 / 3.0), (-2.0 / 3.0), (1.0 / 3.0) };
    var retS: [9]f32 = .{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    var idtS: [9]f32 = .{ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
    var trnS: [9]f32 = .{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    const b1: bool = try xmu.isOrthXmtx(&S, cols, &retS, &idtS, &trnS);
    const b2: bool = try xmu.isOrthogonalXmtx(&S, cols, &alloc);
    const b3: bool = try xmu.isOrthonormalXmtx(&S, cols, &alloc);
    xmu.prntNlStrArgs("Found isOrthXmtx: {}, isOrthogonalXmtx: {}, isOrthonormalXmtx: {}", .{ b1, b2, b3 });
    try std.testing.expectEqual(b1, b2);
    try std.testing.expectEqual(b1, b3);
    try std.testing.expectEqual(b2, b3);
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.3 Example 6 test" {
    //Chapter 5: Section 5.3: Example 6: page 282 - 283
    var basis: [4]f32 = .{ 1, 1, 0, 1 };
    const cols: usize = 2;
    const alloc = std.testing.allocator;
    const prdct: *const fn (l: []f32, r: []f32) f32 = xmu.dotPrdXvec;
    var res: [4]f32 = .{ 0, 0, 0, 0 };
    const exp: [4]f32 = .{ std.math.sqrt(2.0) / 2.0, std.math.sqrt(2.0) / 2.0, -std.math.sqrt(2.0) / 2.0, std.math.sqrt(2.0) / 2.0 };

    try xmu.gramSchmidtOthonormal(&basis, &res, cols, &alloc, prdct);

    xmu.prntNlStr("Basis:");
    xmu.prntXmtxNl(&basis, cols);
    xmu.prntNlStr("Res:");
    xmu.prntXmtxNl(&res, cols);
    xmu.prntNlStr("Exp:");
    xmu.prntXmtxNl(@constCast(&exp), cols);

    try std.testing.expectEqual(true, xmu.equXvecWrkr(&res, @constCast(&exp), false));
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.3 Example 7 test" {
    //Chapter 5: Section 5.3: Example 7: page 283 - 284
    var basis: [9]f32 = .{ 1, 1, 0, 1, 2, 0, 0, 1, 2 };
    const cols: usize = 3;
    const alloc = std.testing.allocator;
    const prdct: *const fn (l: []f32, r: []f32) f32 = xmu.dotPrdXvec;
    var res: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const exp: [9]f32 = .{ std.math.sqrt(2.0) / 2.0, std.math.sqrt(2.0) / 2.0, 0.0, -std.math.sqrt(2.0) / 2.0, std.math.sqrt(2.0) / 2.0, 0.0, 0.0, 0.0, 1.0 };

    try xmu.gramSchmidtOthonormal(&basis, &res, cols, &alloc, prdct);

    xmu.prntNlStr("Basis:");
    xmu.prntXmtxNl(&basis, cols);
    xmu.prntNlStr("Res:");
    xmu.prntXmtxNl(&res, cols);
    xmu.prntNlStr("Exp:");
    xmu.prntXmtxNl(@constCast(&exp), cols);

    try std.testing.expectEqual(true, xmu.equXvecWrkr(&res, @constCast(&exp), false));
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.3 Example 8 test" {
    //Chapter 5: Section 5.3: Example 8: page 284
    var basis: [6]f32 = .{ 0, 1, 0, 1, 1, 1 };
    const cols: usize = 3;
    const alloc = std.testing.allocator;
    const prdct: *const fn (l: []f32, r: []f32) f32 = xmu.dotPrdXvec;
    var res: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    const exp: [6]f32 = .{ 0, 1, 0, std.math.sqrt(2.0) / 2.0, 0, std.math.sqrt(2.0) / 2.0 };

    try xmu.gramSchmidtOthonormal(&basis, &res, cols, &alloc, prdct);

    xmu.prntNlStr("Basis:");
    xmu.prntXmtxNl(&basis, cols);
    xmu.prntNlStr("Res:");
    xmu.prntXmtxNl(&res, cols);
    xmu.prntNlStr("Exp:");
    xmu.prntXmtxNl(@constCast(&exp), cols);

    try std.testing.expectEqual(true, xmu.equXvecWrkr(&res, @constCast(&exp), false));
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.3 Example 10 test" {
    //Chapter 5: Section 5.3: Example 10: page 285 - 286
    var basis: [8]f32 = .{ -2, 2, 1, 0, 1, -8, 0, 1 };
    const cols: usize = 4;
    const alloc = std.testing.allocator;
    const prdct: *const fn (l: []f32, r: []f32) f32 = xmu.dotPrdXvec;
    var res: [8]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0 };
    const exp: [8]f32 = .{ -2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0, 0, -3.0 / std.math.sqrt(30.0), -4.0 / std.math.sqrt(30.0), 2.0 / std.math.sqrt(30.0), 1.0 / std.math.sqrt(30.0) };

    try xmu.gramSchmidtOthonormal(&basis, &res, cols, &alloc, prdct);

    xmu.prntNlStr("Basis:");
    xmu.prntXmtxNl(&basis, cols);
    xmu.prntNlStr("Res:");
    xmu.prntXmtxNl(&res, cols);
    xmu.prntNlStr("Exp:");
    xmu.prntXmtxNl(@constCast(&exp), cols);

    try std.testing.expectEqual(true, xmu.equXvecWrkr(&res, @constCast(&exp), false));
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.3 Problem 1, 3, 5, 7, 9, 11 test" {
    //Chapter 5: Section 5.3: Problem 1, 3, 5, 7, 9, 11: pg 286
    //For problems 1 - 11 determine if the et of vectors in R^n is orthogonal, orthonormal, or neither.
    //(1) neither
    var mtx1: [4]f32 = .{ -4, 6, 5, 0 };
    var cols1: usize = 2;
    const alloc = std.testing.allocator;
    var b1: bool = false;
    var b2: bool = false;
    var exp1: bool = false;
    var exp2: bool = false;

    b1 = try xmu.isOrthogonalXmtx(@constCast(&mtx1), cols1, &alloc);
    b2 = try xmu.isOrthonormalXmtx(@constCast(&mtx1), cols1, &alloc);
    try std.testing.expectEqual(exp1, b1);
    try std.testing.expectEqual(exp2, b2);
    xmu.prntNl();

    //(3) orthonormal
    mtx1 = [4]f32{ 3.0 / 5.0, 4.0 / 5.0, -4.0 / 5.0, 3.0 / 5.0 };
    cols1 = 2;
    exp1 = true;
    exp2 = true;
    b1 = try xmu.isOrthogonalXmtx(@constCast(&mtx1), cols1, &alloc);
    b2 = try xmu.isOrthonormalXmtx(@constCast(&mtx1), cols1, &alloc);
    try std.testing.expectEqual(exp1, b1);
    try std.testing.expectEqual(exp2, b2);
    xmu.prntNl();

    //(5) orthogonal
    var mtx2: [9]f32 = .{ 4, -1, 1, -1, 0, 4, -4, -17, -1 };
    cols1 = 3;
    exp1 = true;
    exp2 = false;
    b1 = try xmu.isOrthogonalXmtx(@constCast(&mtx2), cols1, &alloc);
    b2 = try xmu.isOrthonormalXmtx(@constCast(&mtx2), cols1, &alloc);
    try std.testing.expectEqual(exp1, b1);
    try std.testing.expectEqual(exp2, b2);
    xmu.prntNl();

    //(7) neither
    mtx2 = .{ std.math.sqrt(2.0) / 3.0, 0, std.math.sqrt(2.0) / -6.0, 0, (2.0 * std.math.sqrt(5.0)) / 5.0, std.math.sqrt(5.0) / -5.0, std.math.sqrt(5.0) / 5.0, 0, 1.0 / 2.0 };
    cols1 = 3;
    exp1 = false;
    exp2 = false;
    b1 = try xmu.isOrthogonalXmtx(@constCast(&mtx2), cols1, &alloc);
    b2 = try xmu.isOrthonormalXmtx(@constCast(&mtx2), cols1, &alloc);
    try std.testing.expectEqual(exp1, b1);
    try std.testing.expectEqual(exp2, b2);
    xmu.prntNl();

    //(9) orthonormal
    var mtx3: [12]f32 = .{ std.math.sqrt(2.0) / 2.0, 0, 0, std.math.sqrt(2.0) / 2.0, 0, std.math.sqrt(2.0) / 2.0, std.math.sqrt(2.0) / 2.0, 0, -0.5, 0.5, -0.5, 0.5 };
    cols1 = 4;
    exp1 = true;
    exp2 = true;
    b1 = try xmu.isOrthogonalXmtx(@constCast(&mtx3), cols1, &alloc);
    b2 = try xmu.isOrthonormalXmtx(@constCast(&mtx3), cols1, &alloc);
    try std.testing.expectEqual(exp1, b1);
    try std.testing.expectEqual(exp2, b2);
    xmu.prntNl();

    //(11) orthonormal
    var mtx4: [16]f32 = .{ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };
    cols1 = 4;
    exp1 = true;
    exp2 = true;
    b1 = try xmu.isOrthogonalInrPrdctSpcXmtx(@constCast(&mtx4), cols1, xmu.dotPrdXvec, &alloc);
    b2 = try xmu.isOrthonormalInrPrdctSpcXmtx(@constCast(&mtx4), cols1, xmu.dotPrdXvec, &alloc);
    try std.testing.expectEqual(exp1, b1);
    try std.testing.expectEqual(exp2, b2);
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.3 Problem 13, 15, 17 test" {
    //Chapter 5: Section 5.3: Problem 13, 15, 17: pg 287
    //Find the coordinates of x relative to the orthonormal basis in R^n.

    //(13)
    //B = {
    //        ( ((-2.0 / 13.0) * std.math.sqrt(13)), ((3.0 / 13.0) * std.math.sqrt(13)) ),
    //        ( ((3.0 / 13.0) * std.math.sqrt(13)), ((2.0 / 13.0) * std.math.sqrt(13)) )
    //    }
    //x = (1,2)
    const alloc = std.testing.allocator;
    var B1: [4]f32 = .{ ((-2.0 * std.math.sqrt(13.0)) / 13.0), ((3.0 * std.math.sqrt(13.0)) / 13.0), ((3.0 * std.math.sqrt(13.0)) / 13.0), ((2.0 * std.math.sqrt(13.0)) / 13.0) };
    var x1: [2]f32 = .{ 1.0, 2.0 };
    var cols: usize = 2;
    var res1: [2]f32 = .{ 0.0, 0.0 };
    var exp1: [2]f32 = .{ ((4.0 * std.math.sqrt(13.0)) / 13.0), ((7.0 * std.math.sqrt(13.0)) / 13.0) };

    try xmu.coordRelOrthBasisInrPrdctSpc(&B1, cols, &x1, &res1, &alloc, xmu.dotPrdXvec);

    xmu.prntNlStr("Problem 13:");
    xmu.prntNlStr("Basis B1:");
    xmu.prntXmtxNl(&B1, cols);
    xmu.prntNlStr("Vector x1:");
    xmu.prntXvecNl(&x1);
    xmu.prntNlStrArgs("Cols: {}", .{cols});
    xmu.prntNlStr("Res1:");
    xmu.prntXvecNl(&res1);
    xmu.prntNlStr("Exp1:");
    xmu.prntXvecNl(&exp1);

    try std.testing.expectEqual(true, xmu.equXvecWrkr(&exp1, &res1, false));
    xmu.prntNl();

    //(15)
    var B2: [9]f32 = .{ ((1.0 / 10.0) * std.math.sqrt(10)), 0, ((3.0 / 10.0) * std.math.sqrt(10)), 0, 1, 0, ((-3.0 / 10.0) * std.math.sqrt(10)), 0, ((1.0 / 10.0) * std.math.sqrt(10)) };
    var x2: [3]f32 = .{ 2, -2, 1 };
    cols = 3;
    var res2: [3]f32 = .{ 0.0, 0.0, 0.0 };
    var exp2: [3]f32 = .{ ((1.0 / 2.0) * std.math.sqrt(10)), -2.0, ((-1.0 / 2.0) * std.math.sqrt(10)) };

    try xmu.coordRelOrthBasisInrPrdctSpc(&B2, cols, &x2, &res2, &alloc, xmu.dotPrdXvec);

    xmu.prntNlStr("Problem 15:");
    xmu.prntNlStr("Basis B2:");
    xmu.prntXmtxNl(&B2, cols);
    xmu.prntNlStr("Vector x2:");
    xmu.prntXvecNl(&x2);
    xmu.prntNlStrArgs("Cols: {}", .{cols});
    xmu.prntNlStr("Res2:");
    xmu.prntXvecNl(&res2);
    xmu.prntNlStr("Exp2:");
    xmu.prntXvecNl(&exp2);

    try std.testing.expectEqual(true, xmu.equXvecWrkr(&exp2, &res2, false));
    xmu.prntNl();

    //(17)
    var B3: [9]f32 = .{ (3.0 / 5.0), (4.0 / 5.0), 0, (-4.0 / 5.0), (3.0 / 5.0), 0, 0, 0, 1 };
    var x3: [3]f32 = .{ 5, 10, 15 };
    cols = 3;
    var res3: [3]f32 = .{ 0.0, 0.0, 0.0 };
    var exp3: [3]f32 = .{ 11.0, 2.0, 15.0 };

    try xmu.coordRelOrthBasisInrPrdctSpc(&B3, cols, &x3, &res3, &alloc, xmu.dotPrdXvec);

    xmu.prntNlStr("Problem 17:");
    xmu.prntNlStr("Basis B3:");
    xmu.prntXmtxNl(&B3, cols);
    xmu.prntNlStr("Vector x3:");
    xmu.prntXvecNl(&x3);
    xmu.prntNlStrArgs("Cols: {}", .{cols});
    xmu.prntNlStr("Res3:");
    xmu.prntXvecNl(&res3);
    xmu.prntNlStr("Exp3:");
    xmu.prntXvecNl(&exp3);

    try std.testing.expectEqual(true, xmu.equXvecWrkr(&exp3, &res3, false));
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.3 Problem 19, 21, 23, 25 test" {
    //Chapter 5: Section 5.3: Problem 19, 21, 23, 25: pg 287
    //Use Gram Schmidt orthonormalization process to transform the given basis for R^n into an orthonormal
    //basis. Use the Euclidean inner product for R^n.

    //(19)
    //B = {(3,4), (1,0)}
    //exp = {
    //      (3/5, 4/5),
    //      (4/5, -3/5)
    //}
    const alloc = std.testing.allocator;
    var B1: [4]f32 = .{ 3, 4, 1, 0 };
    var cols: usize = 2;
    var res1: [4]f32 = .{ 0, 0, 0, 0 };
    var exp1: [4]f32 = .{ (3.0 / 5.0), (4.0 / 5.0), (4.0 / 5.0), (-3.0 / 5.0) };

    try xmu.gramSchmidtOthonormal(&B1, &res1, cols, &alloc, xmu.dotPrdXvec);

    xmu.prntNlStr("Problem 19:");
    xmu.prntNlStr("Basis B1:");
    xmu.prntXmtxNl(&B1, cols);
    xmu.prntNlStrArgs("Cols: {}", .{cols});
    xmu.prntNlStr("Res1:");
    xmu.prntXvecNl(&res1);
    xmu.prntNlStr("Exp1:");
    xmu.prntXvecNl(&exp1);

    try std.testing.expectEqual(true, xmu.equXvecWrkr(&exp1, &res1, false));
    xmu.prntNl();

    //(21)
    //B = {(1,-1), (1,1)}
    //exp = {
    //      (sqrt(2)/2, -sqrt(2)/2),
    //      (sqrt(2)/2, sqrt(2)/2)
    //}
    var B2: [4]f32 = .{ 1, -1, 1, 1 };
    cols = 2;
    var res2: [4]f32 = .{ 0, 0, 0, 0 };
    var exp2: [4]f32 = .{ ((std.math.sqrt(2.0)) / 2.0), ((-1.0 * std.math.sqrt(2.0)) / 2.0), ((std.math.sqrt(2.0)) / 2.0), ((std.math.sqrt(2.0)) / 2.0) };

    try xmu.gramSchmidtOthonormal(&B2, &res2, cols, &alloc, xmu.dotPrdXvec);

    xmu.prntNlStr("Problem 21:");
    xmu.prntNlStr("Basis B2:");
    xmu.prntXmtxNl(&B2, cols);
    xmu.prntNlStrArgs("Cols: {}", .{cols});
    xmu.prntNlStr("Res2:");
    xmu.prntXvecNl(&res2);
    xmu.prntNlStr("Exp2:");
    xmu.prntXvecNl(&exp2);

    try std.testing.expectEqual(true, xmu.equXvecWrkr(&exp2, &res2, false));
    xmu.prntNl();

    //(23)
    //B = {(4,-3,0), (1,2,0), (0,0,4)}
    //exp = {
    //      (4/5, -3/5, 0),
    //      (3/5, 4/5, 0),
    //      (0, 0, 1)
    //}
    var B3: [9]f32 = .{ 4, -3, 0, 1, 2, 0, 0, 0, 4 };
    cols = 3;
    var res3: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var exp3: [9]f32 = .{ 4.0 / 5.0, -3.0 / 5.0, 0, 3.0 / 5.0, 4.0 / 5.0, 0, 0, 0, 1.0 };

    try xmu.gramSchmidtOthonormal(&B3, &res3, cols, &alloc, xmu.dotPrdXvec);

    xmu.prntNlStr("Problem 23:");
    xmu.prntNlStr("Basis B3:");
    xmu.prntXmtxNl(&B3, cols);
    xmu.prntNlStrArgs("Cols: {}", .{cols});
    xmu.prntNlStr("Res3:");
    xmu.prntXvecNl(&res3);
    xmu.prntNlStr("Exp3:");
    xmu.prntXvecNl(&exp3);

    try std.testing.expectEqual(true, xmu.equXvecWrkr(&exp3, &res3, false));
    xmu.prntNl();

    //(25)
    //B = {(0,1,1), (1,1,0), (1,0,1)}
    //exp = {
    //      (0, sqrt(2)/2, sqrt(2)/2),
    //      (sqrt(6)/3, sqrt(6)/6, -sqrt(6)/6),
    //      (sqrt(3)/3, -sqrt(3)/3, sqrt(3)/3)
    //}
    var B4: [9]f32 = .{ 0, 1, 1, 1, 1, 0, 1, 0, 1 };
    cols = 3;
    var res4: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var exp4: [9]f32 = .{ 0, ((std.math.sqrt(2.0)) / 2.0), ((std.math.sqrt(2.0)) / 2.0), ((std.math.sqrt(6.0)) / 3.0), ((std.math.sqrt(6.0)) / 6.0), ((-1.0 * std.math.sqrt(6.0)) / 6.0), ((std.math.sqrt(3.0)) / 3.0), ((-1.0 * std.math.sqrt(3.0)) / 3.0), ((std.math.sqrt(3.0)) / 3.0) };

    try xmu.gramSchmidtOthonormal(&B4, &res4, cols, &alloc, xmu.dotPrdXvec);

    xmu.prntNlStr("Problem 25:");
    xmu.prntNlStr("Basis B4:");
    xmu.prntXmtxNl(&B4, cols);
    xmu.prntNlStrArgs("Cols: {}", .{cols});
    xmu.prntNlStr("Res4:");
    xmu.prntXvecNl(&res4);
    xmu.prntNlStr("Exp4:");
    xmu.prntXvecNl(&exp4);

    try std.testing.expectEqual(true, xmu.equXvecWrkr(&exp4, &res4, false));
    xmu.prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.3 Problem 27, 29, 31 test" {
    //Chapter 5: Section 5.3: Problem 27, 29, 31: pg 287
    //Use Gram Schmidt orthonormalization process to transform the given basis for a subspace of R^n
    //into an orthonormal basis for the subspace.

    //(27)
    //B = {(-8, 3, 5)}
    //exp = {
    //      (-4.0 * sqrt(2.0)) / 7.0,
    //      ( 3.0 * sqrt(2.0)) / 14.0,
    //      ( 5.0 * sqrt(2.0)) / 14.0
    //}
    const alloc = std.testing.allocator;
    var cols: usize = 3;
    var B1: [3]f32 = .{ -8, 3, 5 };
    var res1: [3]f32 = .{ 0, 0, 0 };
    var exp1: [3]f32 = .{ (-4.0 * std.math.sqrt(2.0)) / 7.0, (3.0 * std.math.sqrt(2.0)) / 14.0, (5.0 * std.math.sqrt(2.0)) / 14.0 };

    try xmu.gramSchmidtOthonormal(&B1, &res1, cols, &alloc, xmu.dotPrdXvec);

    xmu.prntNlStr("Problem 27:");
    xmu.prntNlStr("Basis B1:");
    xmu.prntXmtxNl(&B1, cols);
    xmu.prntNlStrArgs("Cols: {}", .{cols});
    xmu.prntNlStr("Res1:");
    xmu.prntXvecNl(&res1);
    xmu.prntNlStr("Exp1:");
    xmu.prntXvecNl(&exp1);

    try std.testing.expectEqual(true, xmu.equXvecWrkr(&exp1, &res1, false));
    xmu.prntNl();

    //(29)
    //B = {(3, 4, 0), (1, 0, 0)}
    //exp = {
    //      (3.0/5.0, 4.0/5.0, 0.0),
    //      (4.0/5.0, -3.0/5.0, 0.0)
    //}
    cols = 3;
    var B2: [6]f32 = .{ 3, 4, 0, 1, 0, 0 };
    var res2: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    var exp2: [6]f32 = .{ 3.0 / 5.0, 4.0 / 5.0, 0.0, 4.0 / 5.0, -3.0 / 5.0, 0.0 };

    try xmu.gramSchmidtOthonormal(&B2, &res2, cols, &alloc, xmu.dotPrdXvec);

    xmu.prntNlStr("Problem 29:");
    xmu.prntNlStr("Basis B2:");
    xmu.prntXmtxNl(&B2, cols);
    xmu.prntNlStrArgs("Cols: {}", .{cols});
    xmu.prntNlStr("Res2:");
    xmu.prntXvecNl(&res2);
    xmu.prntNlStr("Exp2:");
    xmu.prntXvecNl(&exp2);

    try std.testing.expectEqual(true, xmu.equXvecWrkr(&exp2, &res2, false));
    xmu.prntNl();

    //(31)
    //B = {(1, 2, -1, 0), (2, 2, 0, 1)}
    //exp = {
    //      std.math.sqrt(6.0)/6.0,
    //      std.math.sqrt(6.0)/3.0,
    //      (-1.0 * std.math.sqrt(6.0))/6.0,
    //      0,
    //      std.math.sqrt(3.0)/3.0,
    //      0,
    //      std.math.sqrt(3.0)/3.0,
    //      std.math.sqrt(3.0)/3.0
    //}

    cols = 4;
    var B3: [8]f32 = .{ 1, 2, -1, 0, 2, 2, 0, 1 };
    var res3: [8]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0 };
    var exp3: [8]f32 = .{ std.math.sqrt(6.0) / 6.0, std.math.sqrt(6.0) / 3.0, (-1.0 * std.math.sqrt(6.0)) / 6.0, 0, std.math.sqrt(3.0) / 3.0, 0, std.math.sqrt(3.0) / 3.0, std.math.sqrt(3.0) / 3.0 };

    try xmu.gramSchmidtOthonormal(&B3, &res3, cols, &alloc, xmu.dotPrdXvec);

    xmu.prntNlStr("Problem 31:");
    xmu.prntNlStr("Basis B3:");
    xmu.prntXmtxNl(&B3, cols);
    xmu.prntNlStrArgs("Cols: {}", .{cols});
    xmu.prntNlStr("Res3:");
    xmu.prntXvecNl(&res3);
    xmu.prntNlStr("Exp3:");
    xmu.prntXvecNl(&exp3);

    try std.testing.expectEqual(true, xmu.equXvecWrkr(&exp3, &res3, false));
    xmu.prntNl();
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//STOP ENTENDED PROBLEMS
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
