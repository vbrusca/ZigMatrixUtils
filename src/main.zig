const std = @import("std");
const xUtils = @import("./XmtxUtils.zig");

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
    xUtils.cpyXvec(&v1, &v2);
    xUtils.cpyXvec(&v1, &v3);
    _ = xUtils.equXvec(&v1, &v2);
    _ = xUtils.equXvec(&v1, &v3);

    try bw.flush(); // don't forget to flush!
}

test "simple test" {
    var v1: [3]f32 = .{ 1, 2, 3 };
    var v2: [3]f32 = .{ 0, 0, 0 };
    var v3: [3]f32 = .{ 1, 1, 1 };
    xUtils.cpyXvec(&v1, &v2);
    xUtils.cpyXvec(&v1, &v3);
    std.debug.print("simple test:\n", .{});
    try std.testing.expectEqual(true, xUtils.equXvec(&v1, &v2));
    try std.testing.expectEqual(true, xUtils.equXvec(&v1, &v3));
}
