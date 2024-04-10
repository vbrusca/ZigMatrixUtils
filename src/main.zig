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

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//STOP ENTENDED PROBLEMS
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
