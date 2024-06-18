//! The XmtxUtils.zig file contains a library of linear algebra functions
//! used to solve different problems in linear algebra.
//!
//! The work here is based on the following texts.
//!
//!
//! Elementary Linear Algebra by Larson, Edwards
//!
//! Mathematics for 3D Game Programming and Computer Graphics 3rd Edition by Eric Lengyel<br/>
//!
//!
//! Developers: Carlo Bruscani, Victor Brusca
//! 04/17/2024
//!
//! How to run tests: zig test .\src\XmtxUtils.zig
//!
//! How to build an exe (not used): zig build-exe -femit-docs ./src/main.zig -O ReleaseSmall -fstrip -fsingle-threaded -femit-bin="zig-out/bin/main.exe"
//!
//! How to build a static library: zig build-lib -femit-docs ./src/XmtxUtils.zig -O ReleaseSmall -fstrip -fsingle-threaded -femit-bin="zig-out/lib/XmtxUtils.lib"
//!
//! How to build an object: zig build-obj -femit-docs ./src/XmtxUtils.zig -O ReleaseSmall -fstrip -fsingle-threaded -femit-bin="zig-out/lib/XmtxUtils.obj"
//!
//! How to build a dynamic library: zig build-lib -femit-docs ./src/XmtxUtils.zig -lc -dynamic -isystem -fstrip -fsingle-threaded -femit-bin="zig-out/lib/XmtxUtils.dll"
//!

const std = @import("std");

comptime {
    @setFloatMode(std.builtin.FloatMode.optimized);
}

const time = std.time;
const Instant = time.Instant;

///The percision to use when determining if a number is zero or when determining if two f32 values are equal.
pub const ZERO_F32: f32 = 0.001;

///A Boolean indicating if the library should use exact comparison of f32 values, or precision comparison using the value of ZERO_F32.
pub const COMPARE_MODE_EXACT: bool = false;

///Error to use when invalid array lengths are encountered.
pub const Error = error{ InvalidLengths, OperationFailed, DivideByZero, NonInvertibleMatrix, AugmenedMatrixRequired, UnsupportedMarixDimension };

///Identity for a 1x1 matrix.
pub const iM1: [1]f32 = .{1};

///Identity for a 2x2 matrix.
pub const iM2: [4]f32 = .{ 1, 0, 0, 1 };

///Identity for a 3x3 matrix.
pub const iM3: [3]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };

///Identity for a 4x4 matrix.
pub const iM4: [4]f32 = .{ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };

///i^ in 1x1 space
pub const iH1: [1]f32 = .{1};

///i^ in 2x2 space
pub const iH2: [2]f32 = .{ 1, 0 };

///j^ in 2x2 space
pub const jH2: [2]f32 = .{ 0, 1 };

///i^ in 3x3 space
pub const iH3: [3]f32 = .{ 1, 0, 0 };

///j^ in 3x3 space
pub const jH3: [3]f32 = .{ 0, 1, 0 };

///k^ in 3x3 space
pub const kH3: [3]f32 = .{ 0, 0, 1 };

///i^ in 4x4 space
pub const iH4: [4]f32 = .{ 1, 0, 0, 0 };

///j^ in 4x4 space
pub const jH4: [4]f32 = .{ 0, 1, 0, 0 };

///k^ in 4x4 space
pub const kH4: [4]f32 = .{ 0, 0, 1, 0 };

///l^ in 4x4 space
pub const lH4: [4]f32 = .{ 0, 0, 0, 1 };

///An enumeration that describes the solution types, discriminant, of a second order polynomial.
pub const POLY2_SOL_TYPE = enum { NONE, ONE_REAL, TWO_REAL };

///An enumeration that describes the solution types, discriminant, of a third order polynomial.
pub const POLY3_SOL_TYPE = enum { NONE, ONE_REPEATED_REAL, THREE_DISTINCT_REALS, ONE_REAL_TWO_IMAGINARY };

///An enumeration of matrix operations that are used by the idnfXmtx and procXmtx functions.
pub const MTX_OPS = enum { MTX_MUL, MTX_DIV, MTX_ADD, MTX_SUB, MTX_PRNT, MTX_NRM, MTX_ABS, MTX_IS_LIN_INDP, MTX_IS_INVERTIBLE, MTX_IS_ZERO };

///An enumeration that is used to describe the type of basis of a 3D set of vectors.
pub const BASIS_HAND = enum { RIGHT, LEFT, ERROR_ZERO, ERROR_INVALID_MATRIX };

///Used in aiding the calculation of Eigen Values for a 2x2 matrix.
pub const EigVal2 = struct {
    lamdExp: [5]f32 = .{ 0, 0, 0, 0, 0 }, //a, -&, d, -&, (-cd)
    polyExp: [3]f32 = .{ 0, 0, 0 }, //a&^2 + b& + c = 0
    eignVals: [2]f32 = .{ 0, 0 },
    eignVec1: [2]f32 = .{ 0, 0 },
    eignVec2: [2]f32 = .{ 0, 0 },

    ///Supports printing out the Eigen value 2 structure in a readable form.
    pub fn prnt(self: EigVal2) void {
        std.debug.print("\n(a + -&)(d + -&) + (-bc)", .{});
        std.debug.print("\n({} + -&)({} + -&) + ({})", .{ self.lamdExp[0], self.lamdExp[2], self.lamdExp[4] });
        std.debug.print("\na&^2 + b& + c = 0", .{});
        std.debug.print("\n{}&^2 + {}& + {} = 0", .{ self.polyExp[0], self.polyExp[1], self.polyExp[2] });
        std.debug.print("\nEigen Values {}, {} = 0", .{ self.eignVals[0], self.eignVals[1] });
    }

    ///Supports printing out the Eigen value 2 structure in a readable form, following a new line.
    pub fn prntNl(self: EigVal2) void {
        std.debug.print("\n");
        prnt(self);
    }
};

///Used in aiding the calculation of Eigen Values for a 3x3 matrix.
pub const EigVal3 = struct {
    // (a - &)*|(e - &), f, h, (i - &)| + (-b)*|d, f, g, (i - &)| + (c)*|d, (e - &), g, h|
    // 0    1    2   3   4  5  6    7     8     9 10 11  12   13    14   15 16   17  18 19
    lamdExp: [20]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    polyExp: [4]f32 = .{ 0, 0, 0, 0 }, //a&^3 + b&^2 + c& + d = 0
    eignVals: [3]f32 = .{ 0, 0, 0 },
    eignVec1: [3]f32 = .{ 0, 0, 0 },
    eignVec2: [3]f32 = .{ 0, 0, 0 },
    eignVec3: [3]f32 = .{ 0, 0, 0 },

    ///Supports printing out the Eigen value 3 structure in a readable form.
    pub fn prnt(self: EigVal3) void {
        std.debug.print("\n(a - &)*|(e - &), f, h, (i - &)| + (-b)*|d, f, g, (i - &)| + (c)*|d, (e - &), g, h|", .{});
        std.debug.print("\n({} - &)*|({} - &), {}, {}, ({} - &)| + (-{})*|{}, {}, {}, ({} - &)| + ({})*|{}, ({} - &), {}, {}|", .{ self.lamdExp[0], self.lamdExp[2], self.lamdExp[4], self.lamdExp[5], self.lamdExp[6], self.lamdExp[8], self.lamdExp[9], self.lamdExp[10], self.lamdExp[11], self.lamdExp[12], self.lamdExp[14], self.lamdExp[15], self.lamdExp[16], self.lamdExp[18], self.lamdExp[19] });
        std.debug.print("\na&^3 + b&^2 + c& + d = 0", .{});
        std.debug.print("\n{}&^3 + {}&^2 + {}& + {} = 0", .{ self.polyExp[0], self.polyExp[1], self.polyExp[2], self.polyExp[3] });
        std.debug.print("\nEigen Values {}, {}, {} = 0", .{ self.eignVals[0], self.eignVals[1], self.eignVals[2] });
    }

    ///Supports printing out the Eigen value 3 structure in a readable form, following a new line.
    pub fn prntNl(self: EigVal3) void {
        std.debug.print("\n");
        prnt(self);
    }
};

///Used to describe a function's average execution time.
pub const ExecTime = struct {
    fncName: []const u8,
    fncTime: f64,
    fncCnt: f64,
    fncAvg: f64,

    ///Supports printing out the ExecTime in a readable form.
    pub fn prnt(self: ExecTime) void {
        std.debug.print("{s}:\tCount: {}\tAvg: {d:.3}ms {d:.3}ns", .{ self.fncName, self.fncCnt, self.fncAvg / time.ns_per_ms, self.fncAvg });
    }

    ///Supports printing out the ExecTime in a readable form, following a new line.
    pub fn prntNl(self: ExecTime) void {
        std.debug.print("\n");
        prnt(self);
    }
};

///The current count of recorded function execution times.
pub var execTimesCnt: usize = 0;

///The maximum number of function execution times that can be recorded.
pub const MAX_EXEC_TIMES: usize = 1024;

///Used to hold a series of function execution times.
pub var execTimes: [MAX_EXEC_TIMES]ExecTime = undefined;

///A Boolean that controls verbose logging.
pub var VERBOSE: bool = false;

///A function to add a new ExecTime struct into the execTimes array.
///
///  key = A string representing the function name.
///
///  value = An f64 value representing the execution time.
///
pub fn addExecTime(key: []const u8, value: f64) void {
    if (execTimesCnt < MAX_EXEC_TIMES) {
        execTimes[execTimesCnt] = .{
            .fncName = key,
            .fncTime = value,
            .fncCnt = 1.0,
            .fncAvg = value,
        };
        execTimesCnt += 1;
    }

    if (VERBOSE) {
        prntNlStrArgs("Added a new entry, maxExecTime entry count is now {}.", .{execTimesCnt});
    }
}

///Returns the absolute value of the provided argument.
///
///  v = The f32 value to calculate the absolute value of.
///
///  returns = The absolute value of the argument v.
///
pub fn absF32(v: f32) f32 {
    if (v < 0.0) {
        return (v * -1.0);
    }
    return v;
}

test "XMTX: absF32 test" {
    const v: f32 = 7;

    const start = try Instant.now();
    const res: f32 = absF32(-7);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nabsF32: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("absF32", elapsed1);

    try std.testing.expectEqual(v, res);
    prntNl();
}

///Calculates and stores the absolute value of the given f32 argument v in v.
///
///  v = A reference to the f32 vaiable to calculate and store the absolute value of.
///
pub fn absF32Ref(v: *f32) void {
    if (v.* < 0.0) {
        v.* *= -1.0;
    }
}

test "XMTX: absF32Ref test" {
    const exp: f32 = 7;
    var v: f32 = -7;
    std.debug.print("\nBefore Exp = {}, V = {}", .{ exp, v });

    const start = try Instant.now();
    absF32Ref(&v);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nabsF32Ref: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("absF32Ref", elapsed1);

    std.debug.print("\nAfter Exp = {}, V = {}", .{ exp, v });
    try std.testing.expectEqual(exp, v);
    prntNl();
}

///Calculates and stores the absolute value of the given f32 argument v into the function argument ret.
///
///  v = A constant reference to the f32 value to calculate the absolute value of.
///
///  ret = A reference to the f32 variable that will store the resulting absolute value.
///
pub fn absF32Ret(v: *const f32, ret: *f32) void {
    ret.* = absF32(v.*);
}

test "XMTX: absF32Ret test" {
    const exp: f32 = 7;
    var v: f32 = -7;
    var ret: f32 = 0;
    std.debug.print("\nBefore Exp = {}, V = {}, Ret = {}", .{ exp, v, ret });

    const start = try Instant.now();
    absF32Ret(&v, &ret);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nabsF32Ret: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("absF32Ret", elapsed1);

    std.debug.print("\nAfter Exp = {}, V = {}, Ret = {}", .{ exp, v, ret });
    try std.testing.expectEqual(false, (exp == v));
    try std.testing.expectEqual(true, (exp == ret));
    prntNl();
}

///Returns the order of the polynomial expression represented as an f32 array.
///
///  polyExp = The polynomial expression to process, represented as an f32 array.
///
///  returns = The order of the polynomial expression.
///
///  NOTE: Description of how polynomials are represented as arrays.
///
///  polyExp of len 2 = a& + b = 0                 ORDER 1
///
///                     [0]  [1]
///
///  polyExp of len 3 = a&^2 + b& + c = 0          ORDER 2
///
///                     [0]    [1]  [2]
///
///  polyExp of len 4 = a&^3 + b&^2 + c& + d = 0   ORDER 3
///
///                     [0]    [1]    [2]  [3]
///
pub fn getPolyOrder(polyExp: []f32) f32 {
    const l: usize = polyExp.len;
    if ((l - 1) > 0) {
        return @floatFromInt(l - 1);
    }
    return 0.0;
}

test "XMTX: getPolyOrder test" {
    var pOdr1: [2]f32 = .{ 0, 0 };
    var pOdr2: [3]f32 = .{ 0, 0, 0 };
    var pOdr3: [4]f32 = .{ 0, 0, 0, 0 };

    var start = try Instant.now();
    const ans1 = getPolyOrder(&pOdr1);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ngetPolyOrder #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getPolyOrder", elapsed1);

    start = try Instant.now();
    const ans2 = getPolyOrder(&pOdr2);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ngetPolyOrder #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getPolyOrder", elapsed1);

    start = try Instant.now();
    const ans3 = getPolyOrder(&pOdr3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ngetPolyOrder #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getPolyOrder", elapsed1);

    const exp1: f32 = 1;
    const exp2: f32 = 2;
    const exp3: f32 = 3;
    try std.testing.expectEqual(ans1, exp1);
    try std.testing.expectEqual(ans2, exp2);
    try std.testing.expectEqual(ans3, exp3);
    prntNl();
}

///Sets the order of the polynomial expression, polyExp - represented as an array of f32 values, to the return argument reference ret.
///
///  polyExp = The constant reference to the polynomial expression to process, represented as an f32 array.
///
///  ret = A reference to the f32 variable that will hold the value of the polynomial order.
///
///  Description of how polynomials are represented as arrays.
///
///  polyExp of len 2 = a& + b = 0                 ORDER 1
///
///                     [0]  [1]
///
///  polyExp of len 3 = a&^2 + b& + c = 0          ORDER 2
///
///                     [0]    [1]  [2]
///
///  polyExp of len 4 = a&^3 + b&^2 + c& + d = 0   ORDER 3
///
///                     [0]    [1]    [2]  [3]
///
pub fn getPolyOrderRet(polyExp: []f32, ret: *f32) void {
    ret.* = getPolyOrder(polyExp);
}

test "XMTX: getPolyOrderRet test" {
    var pOdr1: [2]f32 = .{ 0, 0 };
    var pOdr2: [3]f32 = .{ 0, 0, 0 };
    var pOdr3: [4]f32 = .{ 0, 0, 0, 0 };

    var ans1: f32 = -1;
    var start = try Instant.now();
    getPolyOrderRet(&pOdr1, &ans1);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ngetPolyOrderRet #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getPolyOrderRet", elapsed1);

    var ans2: f32 = -1;
    start = try Instant.now();
    getPolyOrderRet(&pOdr2, &ans2);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ngetPolyOrderRet #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getPolyOrderRet", elapsed1);

    var ans3: f32 = -1;
    start = try Instant.now();
    getPolyOrderRet(&pOdr3, &ans3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ngetPolyOrderRet #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getPolyOrderRet", elapsed1);

    const exp1: f32 = 1;
    const exp2: f32 = 2;
    const exp3: f32 = 3;
    try std.testing.expectEqual(ans1, exp1);
    try std.testing.expectEqual(ans2, exp2);
    try std.testing.expectEqual(ans3, exp3);
    prntNl();
}

///Returns the first Z factors of the number x where Z is equal to the length of the array argument ret.
///
///  x   = The f32 number, should hold an integer value, that factors are searched for.
///
///  ret = The array to store the factors, each factor is stored in its
///        associated array index, i.e. index = 0 would have a value of 1
///        if 1 is a factor of x.
///
///  returns = The number of factors found and stored in the ret argument.
///
pub fn fndFactorsOf(x: f32, ret: []f32) f32 {
    return fndFactorsOfRef(&x, ret);
}

test "XMTX: fndFactorsOf test" {
    const a: f32 = 6;
    var ret: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    const exp: [6]f32 = .{ 1, 2, 3, 0, 0, 6 };

    const start = try Instant.now();
    const fnd: f32 = fndFactorsOf(a, &ret);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nfndFactorsOf: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("fndFactorsOf", elapsed1);

    std.debug.print("\nGiven X = {} we have found {} factors.", .{ a, fnd });
    const l: usize = ret.len;
    var i: usize = 0;
    var cnt: usize = 0;
    while (i < l) : (i += 1) {
        if (ret[i] != 0) {
            cnt += 1;
            std.debug.print("\n{} is a factor of {}", .{ ret[i], a });
            try std.testing.expectEqual(exp[i], ret[i]);
        }
    }
    try std.testing.expectEqual(cnt, @as(usize, @intFromFloat(fnd)));
    prntNl();
}

///Returns the first Z factors of the number x where Z is equal to the length of the array argument ret.
///
///  x   = The f32 number, should hold an integer value, that factors are searched for.
///
///  ret = The array to store the factors, each factor is stored in its
///        associated array index, i.e. index = 0 would have a value of 1
///        if 1 is a factor of x.
///
///  returns = The number of factors found and stored in the ret argument.
///
pub fn fndFactorsOfRef(x: *const f32, ret: []f32) f32 {
    var i: usize = 0;
    const l: usize = ret.len;
    var fnd: f32 = 0;
    const lx: f32 = absF32(x.*);
    const t: usize = @intFromFloat(lx);

    while (i < t) : (i += 1) {
        if (i < l) {
            if (VERBOSE) {
                prntNlStrArgs("t = {}, (i + 1) = {}, (t % (i + 1)) = {}", .{ t, (i + 1), (t % (i + 1)) });
            }

            if ((t % (i + 1)) == 0) {
                ret[i] = @floatFromInt(i + 1);
                fnd += 1;
            } else {
                ret[i] = 0;
            }
        } else {
            break;
        }
    }
    return fnd;
}

test "XMTX: fndFactorsOfRef test" {
    const a: f32 = 6;
    var ret: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    const exp: [6]f32 = .{ 1, 2, 3, 0, 0, 6 };

    const start = try Instant.now();
    const fnd: f32 = fndFactorsOfRef(&a, &ret);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nfndFactorsOfRef: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("fndFactorsOfRef", elapsed1);

    std.debug.print("\nGiven X = {} we have found {} factors.", .{ a, fnd });
    const l: usize = ret.len;
    var i: usize = 0;
    var cnt: usize = 0;
    while (i < l) : (i += 1) {
        if (ret[i] != 0) {
            cnt += 1;
            std.debug.print("\n{} is a factor of {}", .{ ret[i], a });
            try std.testing.expectEqual(exp[i], ret[i]);
        }
    }
    try std.testing.expectEqual(cnt, @as(usize, @intFromFloat(fnd)));
    prntNl();
}

///Returns the first Z factors of the number x where Z = ret.len.
///
///  x   = The f32 number, should hold an integer value, that factors are searched for.
///
///  ret = The array to store the factors, each factor is stored in its
///        associated array index, i.e. index = 0 would have a value of 1
///        if 1 is a factor of x.
///
///  cnt = A reference to the variable that holds the number of factors found and stored in the ret argument.
///
pub fn fndFactorsOfRet(x: *const f32, ret: []f32, cnt: *f32) void {
    cnt.* = fndFactorsOfRef(x, ret);
}

test "XMTX: fndFactorsOfRet test" {
    const a: f32 = 6;
    var ret: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    const exp: [6]f32 = .{ 1, 2, 3, 0, 0, 6 };
    var fnd: f32 = -1;

    const start = try Instant.now();
    fndFactorsOfRet(&a, &ret, &fnd);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nfndFactorsOfRet: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("fndFactorsOfRet", elapsed1);

    std.debug.print("\nGiven X = {} we have found {} factors.", .{ a, fnd });
    const l: usize = ret.len;
    var i: usize = 0;
    var cnt: usize = 0;
    while (i < l) : (i += 1) {
        if (ret[i] != 0) {
            cnt += 1;
            std.debug.print("\n{} is a factor of {}", .{ ret[i], a });
            try std.testing.expectEqual(exp[i], ret[i]);
        }
    }
    try std.testing.expectEqual(cnt, @as(usize, @intFromFloat(fnd)));
    prntNl();
}

///Stores the order 1 polynomial and returns the remainder after a polynomial of order 1 divides a polynomial of order 2.
///
///  poly1 = A pointer to an order 1 polynomial.
///
///  poly2 = A pointer to an order 2 polynomial.
///
///  retPoly1 = A pointer to a return polynomial of order 1.
///
///  returns = The remainder after the synthetic division.
///
pub fn synthDivPoly1IntoPoly2(poly1: *const [2]f32, poly2: *const [3]f32, retPoly1: *[2]f32) f32 {
    //  divPoly1    | poly2
    //              | wrkPoly2
    //              |__________________
    //              | ansPoly2

    var wrkPoly2: [3]f32 = .{ 0, 0, 0 }; //handle work under the quotient poly
    var ansPoly2: [3]f32 = .{ 0, 0, 0 }; //holds the calculated answer
    var ansPoly2Idx: usize = 0; //the current index of the answer polynomial
    var divPoly1: f32 = 0;

    //Step 1: find answer value one
    divPoly1 = (-1.0 * poly1[1]); //constant value of poly1
    ansPoly2[0] = poly2[0]; //initial copy down
    ansPoly2Idx = 1; //starts on index 1 because of the initial copy down

    //Step 2: find answer value two
    wrkPoly2[1] = (divPoly1 * ansPoly2[0]);
    ansPoly2[1] = (poly2[1] + wrkPoly2[1]);

    //Step 3: find answer value three
    wrkPoly2[2] = (divPoly1 * ansPoly2[1]);
    ansPoly2[2] = (poly2[2] + wrkPoly2[2]);

    retPoly1[0] = ansPoly2[0];
    retPoly1[1] = ansPoly2[1];
    return ansPoly2[2];
}

test "XMTX: synthDivPoly1IntoPoly2 test" {
    var poly2Exp: [3]f32 = .{ 1, 5, 6 };
    var poly1Exp: [2]f32 = .{ 1, -1 };
    var retPoly1Exp: [2]f32 = .{ 0, 0 };
    var exp: f32 = 12.0;

    var start = try Instant.now();
    var rmnd: f32 = synthDivPoly1IntoPoly2(&poly1Exp, &poly2Exp, &retPoly1Exp);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nsynthDivPoly1IntoPoly2 #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("synthDivPoly1IntoPoly2", elapsed1);

    std.debug.print("\nAAA Found remainder: {}", .{rmnd});
    prntXvecNl(&retPoly1Exp);
    prntNl();
    try std.testing.expectEqual(exp, rmnd);

    poly2Exp = .{ 1, 1, -2 };
    poly1Exp = .{ 1, 2 };
    retPoly1Exp = .{ 0, 0 };
    exp = 0.0;

    start = try Instant.now();
    rmnd = synthDivPoly1IntoPoly2(&poly1Exp, &poly2Exp, &retPoly1Exp);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nsynthDivPoly1IntoPoly2 #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("synthDivPoly1IntoPoly2", elapsed1);

    std.debug.print("\nBBB Found remainder: {}", .{rmnd});
    prntXvecNl(&retPoly1Exp);
    prntNl();
    try std.testing.expectEqual(exp, rmnd);
    prntNl();
}

///Stores the order 1 polynomial and remainder after a polynomial of order 1 divides a polynomial of order 2 in the retPoly1 function argument.
///
///  poly1 = A pointer to an order 1 polynomial.
///
///  poly2 = A pointer to an order 2 polynomial.
///
///  retPoly1 = A pointer to a return polynomial of order 1.
///
///  ret = Stores the remainder after the synthetic division.
///
pub fn synthDivPoly1IntoPoly2Ret(poly1: *const [2]f32, poly2: *const [3]f32, retPoly1: *[2]f32, ret: *f32) void {
    ret.* = synthDivPoly1IntoPoly2(poly1, poly2, retPoly1);
}

test "XMTX: synthDivPoly1IntoPoly2Ret test" {
    var poly2Exp: [3]f32 = .{ 1, 5, 6 };
    var poly1Exp: [2]f32 = .{ 1, -1 };
    var retPoly1Exp: [2]f32 = .{ 0, 0 };
    var exp: f32 = 12.0;
    var rmnd: f32 = -1;

    var start = try Instant.now();
    synthDivPoly1IntoPoly2Ret(&poly1Exp, &poly2Exp, &retPoly1Exp, &rmnd);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nsynthDivPoly1IntoPoly2Ret #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("synthDivPoly1IntoPoly2Ret", elapsed1);

    std.debug.print("\nAAA Found remainder: {}", .{rmnd});
    prntXvecNl(&retPoly1Exp);
    prntNl();
    try std.testing.expectEqual(exp, rmnd);

    poly2Exp = .{ 1, 1, -2 };
    poly1Exp = .{ 1, 2 };
    retPoly1Exp = .{ 0, 0 };
    exp = 0.0;
    rmnd = -1;

    start = try Instant.now();
    synthDivPoly1IntoPoly2Ret(&poly1Exp, &poly2Exp, &retPoly1Exp, &rmnd);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nsynthDivPoly1IntoPoly2Ret #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("synthDivPoly1IntoPoly2Ret", elapsed1);

    std.debug.print("\nBBB Found remainder: {}", .{rmnd});
    prntXvecNl(&retPoly1Exp);
    prntNl();
    try std.testing.expectEqual(exp, rmnd);
    prntNl();
}

///Stores the order 2 polynomial and returns the remainder after a polynomial of order 1 divides a polynomial of order 3.
///
///  poly1 = A pointer to an order 1 polynomial.
///
///  poly3 = A pointer to an order 3 polynomial.
///
///  retPoly2 = A pointer to a return polynomial of order 2.
///
///  returns = The remainder after the synthetic division.
///
pub fn synthDivPoly1IntoPoly3(poly1: *const [2]f32, poly3: *const [4]f32, retPoly2: *[3]f32) f32 {
    //  divPoly1    | poly3
    //              | wrkPoly3
    //              |__________________
    //              | ansPoly3

    var wrkPoly3: [4]f32 = .{ 0, 0, 0, 0 }; //handle work under the quotient poly
    var ansPoly3: [4]f32 = .{ 0, 0, 0, 0 }; //holds the calculated answer
    var ansPoly3Idx: usize = 0; //the current index of the answer polynomial
    var divPoly1: f32 = 0;

    //Step 1: find answer value one
    divPoly1 = (-1.0 * poly1[1]); //constant value of poly1
    ansPoly3[0] = poly3[0]; //initial copy down
    ansPoly3Idx = 1; //starts on index 1 because of the initial copy down

    //Step 2: find answer value two
    wrkPoly3[1] = (divPoly1 * ansPoly3[0]);
    ansPoly3[1] = (poly3[1] + wrkPoly3[1]);

    //Step 3: find answer value three
    wrkPoly3[2] = (divPoly1 * ansPoly3[1]);
    ansPoly3[2] = (poly3[2] + wrkPoly3[2]);

    //Step 4: find answer value four
    wrkPoly3[3] = (divPoly1 * ansPoly3[2]);
    ansPoly3[3] = (poly3[3] + wrkPoly3[3]);

    retPoly2[0] = ansPoly3[0];
    retPoly2[1] = ansPoly3[1];
    retPoly2[2] = ansPoly3[2];
    return ansPoly3[3];
}

test "XMTX: synthDivPoly1IntoPoly3 test" {
    var poly3Exp: [4]f32 = .{ 2, -3, 4, -1 };
    var poly1Exp: [2]f32 = .{ 1, 1 };
    var retPoly2Exp: [3]f32 = .{ 0, 0, 0 };
    const exp: f32 = -10.0;

    const start = try Instant.now();
    const rmnd: f32 = synthDivPoly1IntoPoly3(&poly1Exp, &poly3Exp, &retPoly2Exp);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nsynthDivPoly1IntoPoly3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("synthDivPoly1IntoPoly3", elapsed1);

    prntNl();
    std.debug.print("\nAAA Found remainder: {}", .{rmnd});
    std.debug.print("\nsynthDivPoly1IntoPoly3 Answers: {any}", .{retPoly2Exp});

    const pexp: [3]f32 = .{ 2.0, -5.0, 9.0 };
    prntXvecNl(&retPoly2Exp);
    prntNl();
    try std.testing.expectEqual(exp, rmnd);
    try std.testing.expectEqual(pexp[0], retPoly2Exp[0]);
    try std.testing.expectEqual(pexp[1], retPoly2Exp[1]);
    try std.testing.expectEqual(pexp[2], retPoly2Exp[2]);

    //TODO: add more tests
    //Possible second test
    //(x + 2) div 3, 1, 1, -5 result => 3, -5, 11, -27

    prntNl();
}

///Stores the order 2 polynomial and remainder after a polynomial of order 1 divides a polynomial of order 3 in the retPoly2 and ret function arguments.
///
///  poly1 = A pointer to an order 1 polynomial.
///
///  poly3 = A pointer to an order 3 polynomial.
///
///  retPoly2 = A pointer to a return polynomial of order 2.
///
///  ret = Stores the remainder after the synthetic division.
///
pub fn synthDivPoly1IntoPoly3Ret(poly1: *const [2]f32, poly3: *const [4]f32, retPoly2: *[3]f32, ret: *f32) void {
    ret.* = synthDivPoly1IntoPoly3(poly1, poly3, retPoly2);
}

test "XMTX: synthDivPoly1IntoPoly3Ret test" {
    var poly3Exp: [4]f32 = .{ 2, -3, 4, -1 };
    var poly1Exp: [2]f32 = .{ 1, 1 };
    var retPoly2Exp: [3]f32 = .{ 0, 0, 0 };
    const exp: f32 = -10.0;
    var rmnd: f32 = -1;

    const start = try Instant.now();
    synthDivPoly1IntoPoly3Ret(&poly1Exp, &poly3Exp, &retPoly2Exp, &rmnd);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nsynthDivPoly1IntoPoly3Ret: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("synthDivPoly1IntoPoly3Ret", elapsed1);

    prntNl();
    std.debug.print("\nAAA Found remainder: {}", .{rmnd});
    std.debug.print("\nsynthDivPoly1IntoPoly3 Answers: {any}", .{retPoly2Exp});

    const pexp: [3]f32 = .{ 2.0, -5.0, 9.0 };
    prntXvecNl(&retPoly2Exp);
    prntNl();
    try std.testing.expectEqual(exp, rmnd);
    try std.testing.expectEqual(pexp[0], retPoly2Exp[0]);
    try std.testing.expectEqual(pexp[1], retPoly2Exp[1]);
    try std.testing.expectEqual(pexp[2], retPoly2Exp[2]);

    //TODO: add more tests
    //Possible second test
    //(x + 2) div 3, 1, 1, -5 result => 3, -5, 11, -27

    prntNl();
}

///A function used to print out polynomials.
///
///  polyExp = The polynomial expression to print out.
///
pub fn prntPolyExp(polyExp: []f32) void {
    const l: usize = polyExp.len;
    if (l == 4) {
        prntNlStrArgs("{}x^3 + {}x^2 + {}x + {}", .{ polyExp[0], polyExp[1], polyExp[2], polyExp[3] });
    } else if (l == 3) {
        prntNlStrArgs("{}x^2 + {}x + {}", .{ polyExp[0], polyExp[1], polyExp[2] });
    } else if (l == 2) {
        prntNlStrArgs("{}x + {}", .{ polyExp[0], polyExp[1] });
    } else if (l == 1) {
        prntNlStrArgs("{}", .{polyExp[0]});
    } else {
        prntNlStrArgs("{any}", .{polyExp});
    }
}

test "XMTX: prntPolyExp test" {
    var p4: [4]f32 = .{ 4, 3, 2, 1 };
    var p3: [3]f32 = .{ 3, 2, 1 };
    var p2: [2]f32 = .{ 2, 1 };
    var p1: [1]f32 = .{1};
    prntPolyExp(&p4);
    prntPolyExp(&p3);
    prntPolyExp(&p2);
    prntPolyExp(&p1);
    prntNl();
}

///Returns the value of the order 3 polynomial resolved with the provided value of x.
///
///  x = The value to use for x.
///
///  polyExp3 = The 3rd order polynomial expression to resolve.
///
///  returns = The value of the polynomial resolved for the given value of x.
///
pub fn rslvPoly3(x: f32, polyExp3: *const [4]f32) f32 {
    return rslvPoly3Ref(&x, polyExp3);
}

test "XMTX: rslvPoly3 test" {
    //2x^3 + 3x^2 – 11x – 6 = 0
    //x = 2, -1/2, and -3
    var p1: [4]f32 = .{ 2, 3, -11, -6 };
    var exp: f32 = 0.0;
    var arg: f32 = 2.0;
    var ans: f32 = -1.0;

    var start = try Instant.now();
    ans = rslvPoly3(arg, &p1);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3 #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = -0.5;
    ans = -1.0;

    start = try Instant.now();
    ans = rslvPoly3(arg, &p1);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3 #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = -3.0;
    ans = -1.0;

    start = try Instant.now();
    ans = rslvPoly3(arg, &p1);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3 #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3", elapsed1);

    try std.testing.expectEqual(exp, ans);
    //x^3 – 2x^2 – x + 2 = 0
    //x = 1, -1, and 2
    var p2: [4]f32 = .{ 1, -2, -1, 2 };
    exp = 0.0;
    arg = 1.0;
    ans = -1.0;

    start = try Instant.now();
    ans = rslvPoly3(arg, &p2);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3 #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = -1.0;
    ans = -1.0;

    start = try Instant.now();
    ans = rslvPoly3(arg, &p2);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3 #4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = 2.0;
    ans = -1.0;

    start = try Instant.now();
    ans = rslvPoly3(arg, &p2);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3 #4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3", elapsed1);

    try std.testing.expectEqual(exp, ans);
    //x^3 − 6x^2 + 11x – 6 = 0
    //x = 1, x = 2 and x = 3
    var p3: [4]f32 = .{ 1, -6, 11, -6 };
    exp = 0.0;
    arg = 1.0;
    ans = -1.0;

    start = try Instant.now();
    ans = rslvPoly3(arg, &p3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3 #5: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = 2.0;
    ans = -1.0;

    start = try Instant.now();
    ans = rslvPoly3(arg, &p3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3 #6: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = 3.0;
    ans = -1.0;

    start = try Instant.now();
    ans = rslvPoly3(arg, &p3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3 #7: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3", elapsed1);

    try std.testing.expectEqual(exp, ans);
    prntNl();
}

///Stores the value of the order 3 polynomial resolved with the provided value of x.
///
///  x = A reference to the value to use for x.
///
///  polyExp3 = The 3rd order polynomial expression to resolve.
///
///  returns = The value of the polynomial resolved for the given value of x.
///
pub fn rslvPoly3Ref(x: *const f32, polyExp3: *const [4]f32) f32 {
    return (x.* * x.* * x.* * polyExp3[0]) + (x.* * x.* * polyExp3[1]) + (x.* * polyExp3[2]) + polyExp3[3];
}

test "XMTX: rslvPoly3Ref test" {
    //2x^3 + 3x^2 – 11x – 6 = 0
    //x = 2, -1/2, and -3
    var p1: [4]f32 = .{ 2, 3, -11, -6 };
    var exp: f32 = 0.0;
    var arg: f32 = 2.0;
    var ans: f32 = -1.0;

    var start = try Instant.now();
    ans = rslvPoly3Ref(&arg, &p1);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ref #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ref", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = -0.5;
    ans = -1.0;

    start = try Instant.now();
    ans = rslvPoly3Ref(&arg, &p1);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ref #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ref", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = -3.0;
    ans = -1.0;

    start = try Instant.now();
    ans = rslvPoly3Ref(&arg, &p1);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ref #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ref", elapsed1);

    try std.testing.expectEqual(exp, ans);
    //x^3 – 2x^2 – x + 2 = 0
    //x = 1, -1, and 2
    var p2: [4]f32 = .{ 1, -2, -1, 2 };
    exp = 0.0;
    arg = 1.0;
    ans = -1.0;

    start = try Instant.now();
    ans = rslvPoly3Ref(&arg, &p2);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ref #4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ref", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = -1.0;
    ans = -1.0;

    start = try Instant.now();
    ans = rslvPoly3Ref(&arg, &p2);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ref #5: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ref", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = 2.0;
    ans = -1.0;

    start = try Instant.now();
    ans = rslvPoly3Ref(&arg, &p2);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ref #6: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ref", elapsed1);

    try std.testing.expectEqual(exp, ans);
    //x^3 − 6x^2 + 11x – 6 = 0
    //x = 1, x = 2 and x = 3
    var p3: [4]f32 = .{ 1, -6, 11, -6 };
    exp = 0.0;
    arg = 1.0;
    ans = -1.0;

    start = try Instant.now();
    ans = rslvPoly3Ref(&arg, &p3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ref #7: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ref", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = 2.0;
    ans = -1.0;

    start = try Instant.now();
    ans = rslvPoly3Ref(&arg, &p3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ref #8: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ref", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = 3.0;
    ans = -1.0;

    start = try Instant.now();
    ans = rslvPoly3Ref(&arg, &p3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ref #9: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ref", elapsed1);

    try std.testing.expectEqual(exp, ans);
    prntNl();
}

///Stores the value of the order 3 polynomial resolved with the provided value of x.
///
///  x = A reference to the value to use for x.
///
///  polyExp3 = The 3rd order polynomial expression to resolve.
///
///  ret = Stores the resolved value of the polynomial resolved for the given value of x.
///
pub fn rslvPoly3Ret(x: *const f32, polyExp3: *const [4]f32, ret: *f32) void {
    ret.* = rslvPoly3Ref(x, polyExp3);
}

test "XMTX: rslvPoly3Ret test" {
    //2x^3 + 3x^2 – 11x – 6 = 0
    //x = 2, -1/2, and -3
    var p1: [4]f32 = .{ 2, 3, -11, -6 };
    var exp: f32 = 0.0;
    var arg: f32 = 2.0;
    var ans: f32 = -1;

    var start = try Instant.now();
    rslvPoly3Ret(&arg, &p1, &ans);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ret #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ret", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = -0.5;
    ans = -1.0;

    start = try Instant.now();
    rslvPoly3Ret(&arg, &p1, &ans);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ret #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ret", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = -3.0;
    ans = -1.0;

    start = try Instant.now();
    rslvPoly3Ret(&arg, &p1, &ans);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ret #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ret", elapsed1);

    try std.testing.expectEqual(exp, ans);
    //x^3 – 2x^2 – x + 2 = 0
    //x = 1, -1, and 2
    var p2: [4]f32 = .{ 1, -2, -1, 2 };
    exp = 0.0;
    arg = 1.0;
    ans = -1.0;

    start = try Instant.now();
    rslvPoly3Ret(&arg, &p2, &ans);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ret #4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ret", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = -1.0;
    ans = -1.0;

    start = try Instant.now();
    rslvPoly3Ret(&arg, &p2, &ans);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ret #5: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ret", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = 2.0;
    ans = -1.0;

    start = try Instant.now();
    rslvPoly3Ret(&arg, &p2, &ans);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ret #6: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ret", elapsed1);

    try std.testing.expectEqual(exp, ans);
    //x^3 − 6x^2 + 11x – 6 = 0
    //x = 1, x = 2 and x = 3
    var p3: [4]f32 = .{ 1, -6, 11, -6 };
    exp = 0.0;
    arg = 1.0;
    ans = -1.0;

    start = try Instant.now();
    rslvPoly3Ret(&arg, &p3, &ans);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ret #7: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ret", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = 2.0;
    ans = -1.0;

    start = try Instant.now();
    rslvPoly3Ret(&arg, &p3, &ans);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ret #8: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ret", elapsed1);

    try std.testing.expectEqual(exp, ans);
    exp = 0.0;
    arg = 3.0;
    ans = -1.0;

    start = try Instant.now();
    rslvPoly3Ret(&arg, &p3, &ans);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly3Ret #9: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly3Ret", elapsed1);

    try std.testing.expectEqual(exp, ans);
    prntNl();
}

///Returns the value of the order 2 polynomial resolved with the provided value of x.
///
///  x = The value to use for x.
///
///  polyExp2 = The 2nd order polynomial expression to resolve.
///
///  returns = The value of the polynomial resolved for the given value of x.
///
pub fn rslvPoly2(x: f32, polyExp2: *const [3]f32) f32 {
    return rslvPoly2Ref(&x, polyExp2);
}

test "XMTX: rslvPoly2 test" {
    //x 2 – 6 x - 16 = 0
    //( x – 8)( x + 2) = 0
    var p1: [3]f32 = .{ 1, -6, -16 };
    const exp: f32 = 0.0;
    var arg: f32 = 8.0;
    var ans: f32 = -1.0;

    var start = try Instant.now();
    ans = rslvPoly2(arg, &p1);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2 #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2", elapsed1);

    try std.testing.expectEqual(exp, ans);
    arg = -2.0;

    start = try Instant.now();
    ans = rslvPoly2(arg, &p1);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2 #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2", elapsed1);

    try std.testing.expectEqual(exp, ans);
    //x 2 + 6 x + 5 = 0
    //( x + 5)( x + 1) = 0
    var p2: [3]f32 = .{ 1, 6, 5 };
    arg = -5.0;

    start = try Instant.now();
    ans = rslvPoly2(arg, &p2);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2 #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2", elapsed1);

    try std.testing.expectEqual(exp, ans);
    arg = -1.0;

    start = try Instant.now();
    ans = rslvPoly2(arg, &p2);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2 #4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2", elapsed1);

    try std.testing.expectEqual(exp, ans);
    //x 2 – 16 = 0
    //( x + 4)( x - 4) = 0
    var p3: [3]f32 = .{ 1, 0, -16 };
    arg = -4.0;

    start = try Instant.now();
    ans = rslvPoly2(arg, &p3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2 #5: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2", elapsed1);

    try std.testing.expectEqual(exp, ans);
    arg = 4.0;

    start = try Instant.now();
    ans = rslvPoly2(arg, &p3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2 #6: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2", elapsed1);

    try std.testing.expectEqual(exp, ans);
    prntNl();
}

///Returns the value of the order 2 polynomial resolved with the provided value of x.
///
///  x = A reference to the value to use for x.
///
///  polyExp2 = The 2nd order polynomial expression to resolve.
///
///  returns = The value of the polynomial resolved for the given value of x.
///
pub fn rslvPoly2Ref(x: *const f32, polyExp2: *const [3]f32) f32 {
    return (x.* * x.* * polyExp2[0]) + (x.* * polyExp2[1]) + polyExp2[2];
}

test "XMTX: rslvPoly2Ref test" {
    //x 2 – 6 x - 16 = 0
    //( x – 8)( x + 2) = 0
    var p1: [3]f32 = .{ 1, -6, -16 };
    const exp: f32 = 0.0;
    var arg: f32 = 8.0;
    var ans: f32 = -1.0;

    var start = try Instant.now();
    ans = rslvPoly2Ref(&arg, &p1);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2Ref #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2Ref", elapsed1);

    try std.testing.expectEqual(exp, ans);
    arg = -2.0;

    start = try Instant.now();
    ans = rslvPoly2Ref(&arg, &p1);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2Ref #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2Ref", elapsed1);

    try std.testing.expectEqual(exp, ans);
    //x 2 + 6 x + 5 = 0
    //( x + 5)( x + 1) = 0
    var p2: [3]f32 = .{ 1, 6, 5 };
    arg = -5.0;

    start = try Instant.now();
    ans = rslvPoly2Ref(&arg, &p2);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2Ref #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2Ref", elapsed1);

    try std.testing.expectEqual(exp, ans);
    arg = -1.0;

    start = try Instant.now();
    ans = rslvPoly2Ref(&arg, &p2);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2Ref #4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2Ref", elapsed1);

    try std.testing.expectEqual(exp, ans);
    //x 2 – 16 = 0
    //( x + 4)( x - 4) = 0
    var p3: [3]f32 = .{ 1, 0, -16 };
    arg = -4.0;

    start = try Instant.now();
    ans = rslvPoly2Ref(&arg, &p3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2Ref #5: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2Ref", elapsed1);

    try std.testing.expectEqual(exp, ans);
    arg = 4.0;

    start = try Instant.now();
    ans = rslvPoly2Ref(&arg, &p3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2Ref #6: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2Ref", elapsed1);

    try std.testing.expectEqual(exp, ans);
    prntNl();
}

///Stores the value of the order 2 polynomial resolved with the provided value of x.
///
///  x = A reference to the value to use for x.
///
///  polyExp2 = The 2nd order polynomial expression to resolve.
///
///  ret = Stores the value of the polynomial resolved for the given value of x.
///
pub fn rslvPoly2Ret(x: *const f32, polyExp2: *const [3]f32, ret: *f32) void {
    ret.* = rslvPoly2Ref(x, polyExp2);
}

test "XMTX: rslvPoly2Ret test" {
    //x 2 – 6 x - 16 = 0
    //( x – 8)( x + 2) = 0
    var p1: [3]f32 = .{ 1, -6, -16 };
    const exp: f32 = 0.0;
    var arg: f32 = 8.0;
    var ans: f32 = -1;

    var start = try Instant.now();
    rslvPoly2Ret(&arg, &p1, &ans);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2Ret #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2Ret", elapsed1);

    try std.testing.expectEqual(exp, ans);
    arg = -2.0;
    ans = -1;

    start = try Instant.now();
    rslvPoly2Ret(&arg, &p1, &ans);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2Ret #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2Ret", elapsed1);

    try std.testing.expectEqual(exp, ans);
    //x 2 + 6 x + 5 = 0
    //( x + 5)( x + 1) = 0
    var p2: [3]f32 = .{ 1, 6, 5 };
    arg = -5.0;
    ans = -1;

    start = try Instant.now();
    rslvPoly2Ret(&arg, &p2, &ans);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2Ret #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2Ret", elapsed1);

    try std.testing.expectEqual(exp, ans);
    arg = -1.0;
    ans = -1;

    start = try Instant.now();
    rslvPoly2Ret(&arg, &p2, &ans);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2Ret #4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2Ret", elapsed1);

    try std.testing.expectEqual(exp, ans);
    //x 2 – 16 = 0
    //( x + 4)( x - 4) = 0
    var p3: [3]f32 = .{ 1, 0, -16 };
    arg = -4.0;
    ans = -1;

    start = try Instant.now();
    rslvPoly2Ret(&arg, &p3, &ans);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2Ret #5: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2Ret", elapsed1);

    try std.testing.expectEqual(exp, ans);
    arg = 4.0;
    ans = -1;

    start = try Instant.now();
    rslvPoly2Ret(&arg, &p3, &ans);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvPoly2Ret #6: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvPoly2Ret", elapsed1);

    try std.testing.expectEqual(exp, ans);
    prntNl();
}

///Finds the roots of the order 3 polynomial provided.
///
///  polyExp = An order 3 polynomial provided as the basis of this function.
///
///  factors = An array for storing the factors of the d term, used to find roots to the polynomial before synthetic division.
///
///  returns = Up to three values that are roots of the polynomial expression, using f32_nan for imaginary numbers which aren't currently supported.
///
pub fn rtsPoly3(polyExpIn: *const [4]f32, factors: []f32) [3]f32 {
    //p[0] = a*x^3, p[1] = b*x^2, p[2] = c*x, p[3] = d
    //check for two basic types of 3rd order polynomials one with a d and one without
    //prep local copy of the polynomial argument because we will alter the values if d is present
    var polyExpPrp: [4]f32 = .{ polyExpIn.*[0], polyExpIn.*[1], polyExpIn.*[2], polyExpIn.*[3] };
    var polyExp: *[4]f32 = &polyExpPrp;

    var ret: [3]f32 = .{ std.math.nan(f32), std.math.nan(f32), std.math.nan(f32) };
    const hasD: bool = (polyExp[3] != 0.0);
    const origD: f32 = polyExp[3];

    if (VERBOSE) {
        prntNlStrArgs("AAA Found {}x^3 + {}x^2 + {}x + {}", .{ polyExp[0], polyExp[1], polyExp[2], polyExp[3] });
    }

    if (polyExp[0] != 1.0) {
        if (polyExp[0] != 0.0) {
            polyExp[3] = (polyExp[3] / polyExp[0]);
            polyExp[2] = (polyExp[2] / polyExp[0]);
            polyExp[1] = (polyExp[1] / polyExp[0]);
            polyExp[0] = 1.0;
        } else {
            prntNlStr("!! Warning expected polynomial to have a degree 3 coefficient !!");
            ret = .{ std.math.nan(f32), std.math.nan(f32), std.math.nan(f32) };
            return ret;
        }
    }

    if (VERBOSE) {
        prntNlStrArgs("BBB Found {}x^3 + {}x^2 + {}x + {}", .{ polyExp[0], polyExp[1], polyExp[2], polyExp[3] });
    }

    if (!hasD) {
        //approach #1, if d = 0
        //solution #1 = 0
        //divide by x and solve the resulting quadratic equation
        var poly2: [3]f32 = .{ polyExp[0], polyExp[1], polyExp[2] };
        const poly2ret: [2]f32 = rtsPoly2(&poly2);
        ret[0] = 0.0;
        ret[1] = poly2ret[0];
        ret[2] = poly2ret[1];
        //std.debug.print("Found solutions in approach #1 {}, {}, {}\n", .{ret[0], ret[1], ret[2]});
        return ret;
    } else {
        //approach #2
        const lcnt: usize = @intFromFloat(absF32(origD)); //@as(usize, @intFromFloat(absF32(origD)));
        const factCnt: f32 = fndFactorsOf(origD, factors);
        _ = factCnt;

        var i: usize = 0;
        var j: usize = 0;
        var ans1: f32 = 0;
        var ans2: f32 = 0;
        //std.debug.print("Found {}x^3 + {}x^2 + {}x + {}\n", .{polyExp[0], polyExp[1], polyExp[2], polyExp[3]});
        //std.debug.print("Found {} factors of {}, {any}\n", .{factCnt, polyExp[3], factors});
        //std.debug.print("=====lcnt {}\n", .{lcnt});

        while (i < lcnt) : (i += 1) {
            if (i < factors.len and factors[i] != 0) {
                ans1 = rslvPoly3(factors[i], polyExp);
                ans2 = rslvPoly3((-1.0 * factors[i]), polyExp);
                //std.debug.print("factors[i] {} ans1 {} ans2 {}\n", .{factors[i], ans1, ans2});

                if (ans1 == 0) {
                    //found root
                    ret[j] = factors[i];
                    //std.debug.print("Found root {} at index {}\n", .{factors[i], j});
                    j += 1;
                }

                if (ans2 == 0) {
                    //found root
                    ret[j] = (-1.0 * factors[i]);
                    //std.debug.print("Found root {} at index {}\n", .{(-1.0 * factors[i]), j});
                    j += 1;
                }
            }
        }

        //If a root has been found use it and synthetic division
        //to find the remaining roots
        if (j > 0) {
            var poly1: [2]f32 = .{ 1, -1.0 * ret[0] };
            var retPoly2: [3]f32 = .{ 0, 0, 0 };
            const rmdr: f32 = synthDivPoly1IntoPoly3(&poly1, polyExp, &retPoly2);
            //std.debug.print("AAA Synth Div Return: {any}, Remainder: {}\n", .{retPoly2, rmdr});

            if (rmdr == 0) {
                const retPoly1: [2]f32 = rtsPoly2(&retPoly2);
                //std.debug.print("BBB Roots of Synth Div Return: {any}\n", .{retPoly1});
                ret[1] = retPoly1[0];
                ret[2] = retPoly1[1];
                return ret;
            } else {
                prntNlStr("!! Warning expected remainder from synthetic division to be 0.0 !!");
                ret = .{ std.math.nan(f32), std.math.nan(f32), std.math.nan(f32) };
                return ret;
            }
        } else {
            prntNlStr("!! Warning expected to find one root using factors of d !!");
            ret = .{ std.math.nan(f32), std.math.nan(f32), std.math.nan(f32) };
            return ret;
        }
    }
}

test "XMTX: rtsPoly3 test" {
    //2x^3 + 3x^2 – 11x – 6 = 0
    //x = 2, -1/2, and -3
    var p1: [4]f32 = .{ 2, 3, -11, -6 };
    var factors: [10]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    const start = try Instant.now();
    const ret: [3]f32 = rtsPoly3(&p1, &factors);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nrtsPoly3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rtsPoly3", elapsed1);

    prntPolyExp(&p1);
    std.debug.print("\nrtsPoly3 roots found: {any}", .{ret});

    const exp1: f32 = 2;
    const exp2: f32 = -0.5;
    const exp3: f32 = -3.0;
    try std.testing.expectEqual(exp1, ret[0]);
    try std.testing.expectEqual(exp2, ret[1]);
    try std.testing.expectEqual(exp3, ret[2]);
    prntNl();
}

///Finds the roots of the order 3 polynomial provided.
///
///  polyExp = An order 3 polynomial provided as the basis of this function.
///
///  factors = An array for storing the factors of the d term, used to find roots to the polynomial before synthetic division.
///
///  ret = Stores up to three values that are roots of the polynomial expression, using f32_nan for imaginary numbers which aren't currently supported.
///
pub fn rtsPoly3Ret(polyExp: *const [4]f32, factors: []f32, ret: *[3]f32) void {
    ret.* = rtsPoly3(polyExp, factors);
}

test "XMTX: rtsPoly3Ret test" {
    //2x^3 + 3x^2 – 11x – 6 = 0
    //x = 2, -1/2, and -3
    var p1: [4]f32 = .{ 2, 3, -11, -6 };
    var factors: [10]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var ret: [3]f32 = .{ 0, 0, 0 };

    const start = try Instant.now();
    rtsPoly3Ret(&p1, &factors, &ret);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nrtsPoly3Ret: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rtsPoly3Ret", elapsed1);

    prntPolyExp(&p1);
    std.debug.print("\nrtsPoly3 roots found: {any}", .{ret});

    const exp1: f32 = 2;
    const exp2: f32 = -0.5;
    const exp3: f32 = -3.0;
    try std.testing.expectEqual(exp1, ret[0]);
    try std.testing.expectEqual(exp2, ret[1]);
    try std.testing.expectEqual(exp3, ret[2]);

    prntNl();
}

///Finds the roots of the order 2 polynomial provided.
///
///  polyExp = An order 2 polynomial provided as the basis of this function.
///
///  returns = Up to two values that are roots of the polynomial expression.
///
pub fn rtsPoly2(polyExp: *const [3]f32) [2]f32 {
    //eigVal1 = ( -b + sqrt(b^2 - (4ac)) ) / (2a)
    var ret: [2]f32 = .{ std.math.nan(f32), std.math.nan(f32) };

    if (polyExp.len != 3) {
        prntNlStr("!! Warning expected a polynomial expression of order to to have a length of 3 !!");
        ret = .{ std.math.nan(f32), std.math.nan(f32) };
        return ret;
    }

    const descr: f32 = std.math.pow(f32, polyExp.*[1], 2) - (4 * polyExp[0] * polyExp[2]);
    const sqrt4ac: f32 = @sqrt(descr);
    const negB: f32 = (-1 * polyExp[1]);
    const twoA: f32 = (2 * polyExp[0]);

    //std.debug.print("\nDiscriminant: {}", .{descr});
    //std.debug.print("\n{} / {}", .{(negB + sqrt4ac), twoA});
    ret[0] = (negB + sqrt4ac) / twoA;

    //std.debug.print("\n{} / {}", .{(negB - sqrt4ac), twoA});
    ret[1] = (negB - sqrt4ac) / twoA;

    //std.debug.print("\n{}, {}", .{ret[0], ret[1]});
    return ret;
}

test "XMTX: rtsPoly2 test" {
    //x2 - 7x + 10 = 0 are x = 2 and x = 5
    var poly2Exp: [3]f32 = .{ 1, -7, 10 };

    const start = try Instant.now();
    var roots: [2]f32 = rtsPoly2(&poly2Exp);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nrtsPoly2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rtsPoly2", elapsed1);

    prntNl();
    prntXvecNl(&roots);

    const exp1: f32 = 5;
    const exp2: f32 = 2;
    try std.testing.expectEqual(exp1, roots[0]);
    try std.testing.expectEqual(exp2, roots[1]);
    prntNl();
}

///Finds the roots of the order 2 polynomial provided.
///
///  polyExp = An order 2 polynomial provided as the basis of this function.
///
///  ret = Stores up to two values that are roots of the polynomial expression.
///
pub fn rtsPoly2Ret(polyExp: *const [3]f32, ret: *[2]f32) void {
    ret.* = rtsPoly2(polyExp);
}

test "XMTX: rtsPoly2Ret test" {
    //x2 - 7x + 10 = 0 are x = 2 and x = 5
    var poly2Exp: [3]f32 = .{ 1, -7, 10 };
    var roots: [2]f32 = .{ 0, 0 };

    const start = try Instant.now();
    rtsPoly2Ret(&poly2Exp, &roots);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nrtsPoly2Ret: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rtsPoly2Ret", elapsed1);

    prntXvecNl(&roots);
    prntNl();

    const exp1: f32 = 5;
    const exp2: f32 = 2;
    try std.testing.expectEqual(exp1, roots[0]);
    try std.testing.expectEqual(exp2, roots[1]);

    prntNl();
}

///Returns an enumeration value describing the descriminant of the cubic polynomial.
///
///  polyExp = The cubic polynomial to process.
///
///  returns = The type of solution to this cubic polynomial.
///
pub fn dscrPoly3(polyExp: *const [4]f32) POLY3_SOL_TYPE {
    //d = 18abcd – 4b³d + b²c² – 4ac³ – 27a²d²
    //d > 0 => then the equation has three distinct real roots
    //d = 0 => hen the equation has a repeated root and all its roots are real
    //d < 0 => then the equation has one real root and two non-real complex conjugate roots.
    const descr: f32 = (18.0 * polyExp[0] * polyExp[1] * polyExp[2] * polyExp[3]) - (4.0 * polyExp[1] * polyExp[1] * polyExp[1] * polyExp[3]) + (polyExp[1] * polyExp[1] * polyExp[2] * polyExp[2]) - (4.0 * polyExp[0] * polyExp[2] * polyExp[2] * polyExp[2]) - (27.0 * polyExp[0] * polyExp[0] * polyExp[3] * polyExp[3]);
    if (isEquF32(descr, 0.0, COMPARE_MODE_EXACT)) {
        return POLY3_SOL_TYPE.ONE_REPEATED_REAL;
    } else if (descr > 0) {
        return POLY3_SOL_TYPE.THREE_DISTINCT_REALS;
    } else {
        return POLY3_SOL_TYPE.ONE_REAL_TWO_IMAGINARY;
    }
}

test "XMTX: dscrPoly3 test" {
    //2x^3 + 3x^2 – 11x – 6 = 0
    //x = 2, -1/2, and -3
    var p1: [4]f32 = .{ 2, 3, -11, -6 };

    const start = try Instant.now();
    const solt: POLY3_SOL_TYPE = dscrPoly3(&p1);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ndscrPoly3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("dscrPoly3", elapsed1);

    try std.testing.expectEqual(POLY3_SOL_TYPE.THREE_DISTINCT_REALS, solt);
    prntNl();
}

///Returns an enumeration value indicating the nature of the solutions for the given second order polynomial expression.
///
///  polyExp = A second order polynomial expression.
///
///  returns = A POLY2_SOL_TYPE enumeration entry that indicates the nature of the polynomial's solution.
///
pub fn dscrPoly2(poly2Exp: *[3]f32) POLY2_SOL_TYPE {
    const descr: f32 = std.math.pow(f32, poly2Exp[1], 2) - (4 * poly2Exp[0] * poly2Exp[2]);
    if (isEquF32(descr, 0.0, COMPARE_MODE_EXACT)) {
        return POLY2_SOL_TYPE.ONE_REAL;
    } else if (descr < 0.0) {
        return POLY2_SOL_TYPE.NONE;
    } else {
        return POLY2_SOL_TYPE.TWO_REAL;
    }
}

test "XMTX: dscrPoly2 test" {
    //x2 - 7x + 10 = 0 are x = 2 and x = 5
    var poly2Exp: [3]f32 = .{ 1, -7, 10 };
    var roots: [2]f32 = rtsPoly2(&poly2Exp);

    const start = try Instant.now();
    const solt: POLY2_SOL_TYPE = dscrPoly2(&poly2Exp);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ndscrPoly2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("dscrPoly2", elapsed1);

    try std.testing.expectEqual(POLY2_SOL_TYPE.TWO_REAL, solt);
    prntXvecNl(&roots);
    prntNl();

    const exp1: f32 = 5;
    const exp2: f32 = 2;
    try std.testing.expectEqual(exp1, roots[0]);
    try std.testing.expectEqual(exp2, roots[1]);

    prntNl();
}

///Find the eigen values for a 3x3 dimensional matrix.
///
///  mtx = A 3x3 matrix to find eigen values for.
///
///  evs = An EigVal3 pointer used for calculating the eigen values of the given matrix.
///
///  factors = An array for storing the factors of the d term, used to find roots to the polynomial before synthetic division.
///
///  returns = A boolean value indicating the success of the operation.
///
pub fn fndEigVal3(mtx: *const [9]f32, evs: *EigVal3, factors: []f32) bool {
    //Given 3x3 matrix A
    //A = a b c     0 1 2
    //    d e f     3 4 5
    //    g h i     6 7 8

    //detA = (a)*|e,f ,h,i| + (-b)*|d,f ,g,i| + (c)*|d,e ,g,h|

    //(A - &*i3) = (a - &) b       c
    //             d       (e - &) f
    //             g       h       (i - &)
    //eigDetA = det(A - &*I3) = (a - &)*|(e - &),f ,h,(i - &)| + (-b)*|d,f ,g,(i - &)| + (c)*|d,(e - &) ,g,h|
    //eigDetA = (-1)&^3 + (a + e + i)&^2 + (cg - ea - ia - ei + bd)& + (eig - bdi + bfg + cdh - ceg)

    //          (a - &)*|(e - &), f,    h, (i - &)| + (-b)*|d,      f, g, (i - &)| + (c)    *|d, (e - &), g, h|
    //lamdExp = 0    1    2   3   4     5  6    7     8     9       10 11  12  13    14       15 16   17  18 19

    //          (-1)&^3 + (a + e + i)&^2 + (cg - ea - ia - ei + bd)& + (eig - bdi + bfg + cdh - ceg)
    //polyExp = 0         1                2                           3
    //position  T         Z                X                           Y

    //if (mtx.len != 9) {
    //    std.debug.print("\n!! Warning fndEigVal3 requires 3x3 matrices !!", .{});
    //    return false;
    //}

    const lamda: f32 = std.math.nan(f32);
    evs.*.lamdExp[0] = mtx[0];
    evs.*.lamdExp[1] = -lamda;
    evs.*.lamdExp[2] = mtx[4];
    evs.*.lamdExp[3] = -lamda;
    evs.*.lamdExp[4] = mtx[5]; //f

    evs.*.lamdExp[5] = mtx[6];
    evs.*.lamdExp[6] = mtx[7];
    evs.*.lamdExp[7] = -lamda;
    evs.*.lamdExp[8] = -1.0 * mtx[1];
    evs.*.lamdExp[9] = mtx[3];

    evs.*.lamdExp[10] = mtx[5]; //f
    evs.*.lamdExp[11] = mtx[6]; //g
    evs.*.lamdExp[12] = mtx[8];
    evs.*.lamdExp[13] = -lamda;
    evs.*.lamdExp[14] = mtx[2];

    evs.*.lamdExp[15] = mtx[3]; //d
    evs.*.lamdExp[16] = mtx[4];
    evs.*.lamdExp[17] = -lamda;
    evs.*.lamdExp[18] = mtx[6];
    evs.*.lamdExp[19] = mtx[7];

    //Given 3x3 matrix A
    //A = a b c     0 1 2
    //    d e f     3 4 5
    //    g h i     6 7 8
    //          (-1)&^3 + (a + e + i)&^2 + (cg - ea - ia - ei + bd)& + (eig - bdi + bfg + cdh - ceg)
    //polyExp = 0         1                2                           3
    //position  T         Z                X                           Y
    evs.*.polyExp[0] = -1.0; //see position T

    //(a + e + i) -> (0 + 4 + 8)
    evs.*.polyExp[1] = mtx[0] + mtx[4] + mtx[8]; //see position Z

    //(cg - ea - ia - ei + bd) -> (2*6 - 4*0 - 8*0 - 4*8 + 1*3)
    evs.*.polyExp[2] = (mtx[2] * mtx[6]) - (mtx[4] * mtx[0]) - (mtx[8] * mtx[0]) - (mtx[4] * mtx[8]) + (mtx[1] * mtx[3]); //see position X

    //(eig - bdi + bfg + cdh - ceg) -> (4*8*6 - 1*3*8 + 1*5*6 + 2*3*7 - 2*4*6)
    evs.*.polyExp[3] = (mtx[4] * mtx[8] * mtx[6]) - (mtx[1] * mtx[3] * mtx[8]) + (mtx[1] * mtx[5] * mtx[6]) + (mtx[2] * mtx[3] * mtx[7]) - (mtx[2] * mtx[4] * mtx[6]); //see position Y
    evs.*.eignVals = rtsPoly3(&evs.*.polyExp, factors);
    return true;
}

test "XMTX: fndEigVal3 test" {
    //A = −2, −2,  4
    //    −4,  1,  2
    //     2,  2,  5
    var A: [9]f32 = .{ -2, -2, 4, -4, 1, 2, 2, 2, 5 };
    var factors: [128]f32 = std.mem.zeroes([128]f32); //.{0};
    var j: usize = 0;

    while (j < 128) : (j += 1) {
        factors[j] = 0;
    }

    var evs: EigVal3 = .{};
    var b: bool = false;

    const start = try Instant.now();
    b = fndEigVal3(&A, &evs, &factors);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nfndEigVal3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("fndEigVal3", elapsed1);

    try std.testing.expectEqual(true, b);
    std.debug.print("\nAAA Eigen Values Found:", .{});
    evs.prnt();

    prntNl();
}

///Find the eigen values for a 2x2 dimensional matrix.
///
///  mtx = A 2x2 matrix to find eigen values for.
///
///  evs = An EigVal2 pointer used for calculating the eigen values of the given matrix.
///
///  returns = A boolean value indicating the success of the operation.
///
pub fn fndEigVal2(mtx: *const [4]f32, evs: *EigVal2) bool {
    //Given 2x2 matrix A
    //A = a b
    //    c d

    //detA = ad - bc
    //eigDetA = (a - &)(d - &) + (-bc)
    //lamdExp    0   1  2   3    4
    //          ad + a(-&) + d(-&) + &^2 + (-bc)

    //polyExp = 2    1       1       0     2
    //position  Y    X       X       Z     Y

    //if (mtx.len != 4) {
    //    std.debug.print("\n!! Warning fndEigVal2 requires 2x2 matrices !!", .{});
    //    return false;
    //}

    const lamda: f32 = std.math.nan(f32);
    evs.*.lamdExp[0] = mtx[0];
    evs.*.lamdExp[1] = -lamda;
    evs.*.lamdExp[2] = mtx[3];
    evs.*.lamdExp[3] = -lamda;
    evs.*.lamdExp[4] = (-1 * mtx[1] * mtx[2]);

    evs.*.polyExp[0] = 1; //see position Z
    evs.*.polyExp[1] = (-1 * mtx[0]) + (-1 * mtx[3]); //see position X
    evs.*.polyExp[2] = (mtx[0] * mtx[3]) + evs.*.lamdExp[4]; //see position Y

    //eigVal1 = ( -b + sqrt(b^2 - (4ac)) ) / (2a)
    //eigVal2 = ( -b - sqrt(b^2 - (4ac)) ) / (2a)
    evs.*.eignVals = rtsPoly2(&evs.*.polyExp);
    return true;
}

test "XMTX: fndEigVal2 test" {
    var m1: [4]f32 = .{ 5, 6, 8, 9 };
    var evs: EigVal2 = .{};
    var b: bool = false;

    var start = try Instant.now();
    b = fndEigVal2(&m1, &evs);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nfndEigVal2 #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("fndEigVal2", elapsed1);

    try std.testing.expectEqual(true, b);
    std.debug.print("\nAAA Eigen Values Found:", .{});
    evs.prnt();
    try std.testing.expectEqual(@as(f32, 1.42111024e+01), evs.eignVals[0]);
    try std.testing.expectEqual(@as(f32, -2.11102485e-01), evs.eignVals[1]);

    m1 = .{ 1, 4, 2, 3 };
    evs = .{};
    b = false;

    start = try Instant.now();
    b = fndEigVal2(&m1, &evs);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nfndEigVal2 #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("fndEigVal2", elapsed1);

    try std.testing.expectEqual(true, b);
    std.debug.print("\nBBB Eigen Values Found:", .{});
    evs.prnt();
    try std.testing.expectEqual(@as(f32, 5.0), evs.eignVals[0]);
    try std.testing.expectEqual(@as(f32, -1.0), evs.eignVals[1]);

    prntNl();
}

///Cleans the given matrix by rounding f32 values to the nearest significance, ZERO_F32, resulting in clean 0.0, 1.0, etc, values.
///
///  mtx = The matrix to process.
///
///  returns = A clean matrix with values truncated.
///
pub fn clnXmtx(mtx: []f32) void {
    clnXvec(mtx);

    //const l: usize = mtx.len;
    //var i: usize = 0;
    //while (i < l) {
    //    if (absF32(@floor(mtx[i]) - mtx[i]) < ZERO_F32) {
    //        mtx[i] = @floor(mtx[i]);
    //    } else if (absF32(@ceil(mtx[i]) - mtx[i]) < ZERO_F32) {
    //        mtx[i] = @ceil(mtx[i]);
    //    }
    //
    //    if (mtx[i] == -0.0) {
    //        mtx[i] = 0.0;
    //    }
    //
    //    i += 1;
    //}
}

test "XMTX: clnXmtx test" {
    var m1: [3]f32 = .{ 1.00000001, 2.00000001, 3.00000001 };
    var m2: [3]f32 = .{ 1, 2, 3 };

    const start = try Instant.now();
    clnXmtx(&m1);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nclnXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("clnXmtx", elapsed1);

    try std.testing.expectEqual(true, equXvec(&m1, &m2));
    prntNl();
}

///Cleans the given vector by rounding float values to the nearest significance, ZERO_F32, resulting in clean 0.0, 1.0, etc, values.
///
///  vec = The vector to process.
///
///  returns = A clean vector with values truncated.
///
pub fn clnXvec(vec: []f32) void {
    const l: usize = vec.len;
    var i: usize = 0;
    while (i < l) {
        if (absF32(@floor(vec[i]) - vec[i]) < ZERO_F32) {
            vec[i] = @floor(vec[i]);
        } else if (absF32(@ceil(vec[i]) - vec[i]) < ZERO_F32) {
            vec[i] = @ceil(vec[i]);
        }

        if (vec[i] == -0.0) {
            vec[i] = 0.0;
        }

        i += 1;
    }
}

test "XMTX: clnXvec test" {
    var m1: [3]f32 = .{ 1.00000001, 2.00000001, 3.00000001 };
    var m2: [3]f32 = .{ 1, 2, 3 };

    const start = try Instant.now();
    clnXvec(&m1);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nclnXvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("clnXvec", elapsed1);

    try std.testing.expectEqual(true, equXvec(&m1, &m2));
    prntNl();
}

///Returns a Boolean value indicating if the two f32 numbers are equal using optional exact comparison or precision delta comparison.
///
///  l = The left-hand number in the comparison.
///
///  r = The right-hand number in the comparison.
///
///  comparisonModeExact = A Boolean value indicating if the comparison mode should be exact or within a specified precision delta.
///
///  returns = A Boolean value indicating if the two numbers provided are equal under the specified comparison mode.
///
pub fn isEquF32(l: f32, r: f32, compareModeExact: bool) bool {
    return isEquF32Ref(&l, &r, &compareModeExact);
    //if (compareModeExact) {
    //    if (l == r) {
    //        return true;
    //    } else {
    //        return false;
    //    }
    //} else {
    //    if (absF32(l - r) < ZERO_F32) {
    //        return true;
    //    } else {
    //        return false;
    //    }
    //}
}

test "XMTX: isEquF32 test" {
    var b: bool = false;

    var start = try Instant.now();
    b = isEquF32(1.0, 1.00001, false);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisEquF32 #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isEquF32", elapsed1);

    try std.testing.expectEqual(true, b);

    start = try Instant.now();
    b = isEquF32(1.0, 1.0, true);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisEquF32 #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isEquF32", elapsed1);

    try std.testing.expectEqual(true, b);

    prntNl();
}

///Returns a Boolean value indicating if the two f32 numbers are equal using optional exact comparison or precision delta comparison.
///
///  l = A reference to the left-hand number in the comparison.
///
///  r = A reference to the right-hand number in the comparison.
///
///  comparisonModeExact = A reference to a Boolean value indicating if the comparison mode should be exact or within a specified precision delta.
///
///  returns = A Boolean value indicating if the two numbers provided are equal under the specified comparison mode.
///
pub fn isEquF32Ref(l: *const f32, r: *const f32, compareModeExact: *const bool) bool {
    if (compareModeExact.*) {
        if (l.* == r.*) {
            return true;
        } else {
            return false;
        }
    } else {
        if (absF32(l.* - r.*) < ZERO_F32) {
            return true;
        } else {
            return false;
        }
    }
}

test "XMTX: isEquF32Ref test" {
    var f1: f32 = 1.000;
    var f2: f32 = 1.0001;
    var f3: f32 = 1.100;
    var b: bool = false;
    var cm: bool = false;

    var start = try Instant.now();
    b = isEquF32Ref(&f1, &f2, &cm);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisEquF32Ref #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isEquF32Ref", elapsed1);

    try std.testing.expectEqual(true, b);
    cm = true;

    start = try Instant.now();
    b = isEquF32Ref(&f1, &f2, &cm);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisEquF32Ref #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isEquF32Ref", elapsed1);

    try std.testing.expectEqual(false, b);
    cm = false;

    start = try Instant.now();
    b = isEquF32Ref(&f1, &f3, &cm);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisEquF32Ref #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isEquF32Ref", elapsed1);

    try std.testing.expectEqual(false, b);
    cm = true;

    start = try Instant.now();
    b = isEquF32Ref(&f1, &f3, &cm);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisEquF32Ref #4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isEquF32Ref", elapsed1);

    try std.testing.expectEqual(false, b);
    prntNl();
}

///Returns a Boolean value indicating if the two f32 numbers are equal using optional exact comparison or precision delta comparison.
///
///  l = A reference to the left-hand number in the comparison.
///
///  r = A reference to the right-hand number in the comparison.
///
///  comparisonModeExact = A reference to a Boolean value indicating if the comparison mode should be exact or within a specified precision delta.
///
///  ret = Stores a Boolean value indicating if the two numbers provided are equal under the specified comparison mode.
///
pub fn isEquF32Ret(l: *const f32, r: *const f32, compareModeExact: *const bool, ret: *bool) void {
    ret.* = isEquF32Ref(l, r, compareModeExact);
}

test "XMTX: isEquF32Ret test" {
    var f1: f32 = 1.000;
    var f2: f32 = 1.0001;
    var f3: f32 = 1.100;
    var b: bool = false;
    var cm: bool = false;

    var start = try Instant.now();
    isEquF32Ret(&f1, &f2, &cm, &b);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisEquF32Ret #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isEquF32Ret", elapsed1);

    try std.testing.expectEqual(true, b);
    cm = true;

    start = try Instant.now();
    isEquF32Ret(&f1, &f2, &cm, &b);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisEquF32Ret #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isEquF32Ret", elapsed1);

    try std.testing.expectEqual(false, b);
    cm = false;

    start = try Instant.now();
    isEquF32Ret(&f1, &f3, &cm, &b);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisEquF32Ret #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isEquF32Ret", elapsed1);

    try std.testing.expectEqual(false, b);
    cm = true;

    start = try Instant.now();
    isEquF32Ret(&f1, &f3, &cm, &b);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisEquF32Ret #4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isEquF32Ret", elapsed1);

    try std.testing.expectEqual(false, b);
    prntNl();
}

///Copies the mtx matrix argument into the ret matrix argument.
///
///  mtx = The matrix to copy into the ret matrix argument.
///
///  ret = The return matrix, the copy to be made.
///
pub fn cpyXmtx(mtx: []f32, ret: []f32) void {
    cpyXvec(mtx, ret);
}

test "XMTX: cpyXmtx test" {
    var m1: [3]f32 = .{ 1, 2, 3 };
    var m2: [3]f32 = .{ 0, 0, 0 };
    var m3: [3]f32 = .{ 1, 1, 1 };

    var start = try Instant.now();
    cpyXmtx(&m1, &m2);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncpyXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cpyXmtx", elapsed1);

    start = try Instant.now();
    cpyXmtx(&m1, &m3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ncpyXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cpyXmtx", elapsed1);

    try std.testing.expectEqual(true, equXvec(&m1, &m2));
    try std.testing.expectEqual(true, equXvec(&m1, &m3));

    prntNl();
}

///Copies the mtx matrix into the ret matrix argument if the column is less than cpyCols.
///
///  mtx = The matrix to copy into the ret matrix argument, less columns greater than or equal to cpyCols
///
///  ret = The return matrix to store the copy in.
///
///  mtxCols = The number of columns in the mtx matrix.
///
///  cpyCols = The number of columns to allow to be copied into the ret matrix, 0 .. (cpyCols - 1).
///
pub fn cpyLessXmtx(mtx: []f32, ret: []f32, mtxCols: usize, cpyCols: usize) void {
    const mtxRows: usize = mtx.len / mtxCols;
    var r: usize = 0;
    var c: usize = 0;
    while (r < mtxRows) : (r += 1) {
        c = 0;
        while (c < mtxCols) : (c += 1) {
            if (c < cpyCols) {
                ret[((r * cpyCols) + c)] = mtx[((r * mtxCols) + c)];
            }
        }
    }
}

test "XMTX: cpyLessXmtx test" {
    var m1: [4]f32 = .{ 1, 2, 3, 4 };
    var m2: [3]f32 = .{ 0, 0, 0 };
    var m3: [3]f32 = .{ 1, 2, 3 };

    const start = try Instant.now();
    cpyLessXmtx(&m1, &m2, 4, 3);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncpyLessXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cpyLessXmtx", elapsed1);

    try std.testing.expectEqual(true, equXvec(&m2, &m3));
    prntNl();
}

///Copies data from the matrix mtx at startCol, startRow to endCol, endRow to ret 0,0 to (endCol - startCol), (endRow - strtRow).
///
///  mtx = The data to copy into the ret matrix.
///
///  ret = The matrix that holds the return data.
///
///  strtCol = The starting column in the mtx matrix.
///
///  endCol = The ending column in the mtx matrix.
///
///  strtRow = The starting row in the mtx matrix.
///
///  endRow = The ending row in the mtx matrix.
///
///  mtxCols = The number of columns in the mtx matrix.
///
///  cpyCols = The number of columns in the ret matrix.
///
pub fn cpyLessColRowXmtx(mtx: []f32, ret: []f32, strtCol: usize, endCol: usize, strtRow: usize, endRow: usize, mtxCols: usize, cpyCols: usize) void {
    var r: usize = strtRow;
    var dstR: usize = 0;
    var c: usize = 0;
    var dstC: usize = 0;
    while (r < endRow) : (r += 1) {
        if (r >= strtRow) {
            c = strtCol;
            dstC = 0;
            while (c < endCol) : (c += 1) {
                if (c >= strtCol) {
                    ret[((dstR * cpyCols) + dstC)] = mtx[((r * mtxCols) + c)];
                    dstC += 1;
                }
            }
            dstR += 1;
        }
    }
}

test "XMTX: cpyLessColRowXmtx test" {
    var m1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var ret: [3]f32 = .{ 0, 0, 0 };
    var exp1: [3]f32 = .{ 1, 2, 3 };
    var exp2: [3]f32 = .{ 7, 8, 9 };

    var start = try Instant.now();
    cpyLessColRowXmtx(&m1, &ret, 0, 3, 0, 1, 3, 3);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncpyLessColRowXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cpyLessColRowXmtx", elapsed1);

    try std.testing.expectEqual(true, equXvec(&ret, &exp1));

    start = try Instant.now();
    cpyLessColRowXmtx(&m1, &ret, 0, 3, 2, 3, 3, 3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ncpyLessColRowXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cpyLessColRowXmtx", elapsed1);

    try std.testing.expectEqual(true, equXvec(&ret, &exp2));

    //using an array with slice argument
    start = try Instant.now();
    cpyLessColRowXmtx3(m1[0..9], &ret, 0, 3, 2, 3, 3, 3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ncpyLessColRowXmtx #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cpyLessColRowXmtx", elapsed1);

    try std.testing.expectEqual(true, equXvec(&ret, &exp2));
    prntNl();
}

///Copies data from the matrix mtx at strtSrcCol, strtSrcRow to endSrcCol, endSrcRow to ret strtRetCol, strtRetRow copying this many rows and columns into the region, (endCol - startCol), (endRow - strtRow).
///
///  mtxSrc = The src matrix data to copy into the ret matrix.
///
///  mtxSrcCols = The number of columns in the mtxSrc argument.
/// 
///  mtxRet = The matrix that holds the return data.
///
///  mtxRetCols = The number of columns in the mtxRet argument.
///
///  strtSrcCol = The starting column in the mtxSrc matrix.
///
///  endSrcCol = The end column in the mtxSrc matrix.
///
///  strtSrcRow = The starting row in the mtxSrc matrix.
///
///  endSrcRow = The ending row in the mtxSrc matrix.
///
///  strtRetCol = The starting column to copy the region into the retMtx matrix.
///
///  strtRetRow = The starting row to copy the region into the retMtx matrix.
/// 
pub fn cpyXmtxSqr(mtxSrc: []f32, mtxSrcCols: usize, mtxRet: []f32, mtxRetCols: usize, strtSrcCol: usize, endSrcCol: usize, strtSrcRow: usize, endSrcRow: usize, strtRetCol: usize, strtRetRow: usize) void {
    var r: usize = strtSrcRow;
    var dstR: usize = strtRetRow;
    var c: usize = 0;
    var dstC: usize = 0;
    while (r < endSrcRow) : (r += 1) {
        if (r >= strtSrcRow) {
            c = strtSrcCol;
            dstC = strtRetCol;
            while (c < endSrcCol) : (c += 1) {
                if (c >= strtSrcCol) {
                    mtxRet[((dstR * mtxRetCols) + dstC)] = mtxSrc[((r * mtxSrcCols) + c)];
                    dstC += 1;
                }
            }
            dstR += 1;
        }
    }
}

test "XMTX: cpyXmtxSqr test" {    
    // mtxSrc = | 1  2  3 |
    //          | 4  5  6 |
    //          | 7  8  9 |
    var mtxSrc: [9]f32 = .{1, 2, 3, 4, 5, 6, 7, 8, 9};
    const mtxSrcCols: usize = 3;
    var mtxRet: [4]f32 = .{1, 2, 3, 4};
    const mtxRetCols: usize = 2;
    const exp: [4]f32 = .{1, 2, 4, 5};

    const strtSrcCol: usize = 0;
    const strtSrcRow: usize = 0;
    const endSrcCol: usize = 2;
    const endSrcRow: usize = 2;
    const strtRetRow: usize = 0;
    const strtRetCol: usize = 0;

    cpyXmtxSqr(&mtxSrc, mtxSrcCols, &mtxRet, mtxRetCols, strtSrcCol, endSrcCol, strtSrcRow, endSrcRow, strtRetCol, strtRetRow);

    prntNlStr("MtxSrc:");
    prntXmtx(&mtxSrc, mtxSrcCols);
    prntNlStr("MtxRet:");
    prntXmtx(&mtxRet, mtxRetCols);
    prntNlStr("Exp:");
    prntXmtx(@constCast(&exp), mtxRetCols);

    try std.testing.expectEqual(true, equXvec(&mtxRet, @constCast(&exp)));
    prntNl();
}

///Copies data from the matrix mtx at startCol, startRow to endCol, endRow to ret 0,0 to (endCol - startCol), (endRow - strtRow).
///
///  mtx = The data to copy into the ret matrix.
///
///  ret = The matrix that holds the return data.
///
///  strtCol = The starting column in the mtx matrix.
///
///  endCol = The ending column in the mtx matrix.
///
///  strtRow = The starting row in the mtx matrix.
///
///  endRow = The ending row in the mtx matrix.
///
///  mtxCols = The number of colmns in the mtx matrix.
///
///  cpyCols = The number of columns in the ret matrix.
///
pub fn cpyLessColRowXmtx3(mtx: *const [9]f32, ret: []f32, strtCol: usize, endCol: usize, strtRow: usize, endRow: usize, mtxCols: usize, cpyCols: usize) void {
    const mtxRows: usize = mtx.len / mtxCols;
    _ = mtxRows;
    var r: usize = strtRow;
    var dstR: usize = 0;
    var c: usize = 0;
    var dstC: usize = 0;
    while (r < endRow) : (r += 1) {
        if (r >= strtRow) {
            c = strtCol;
            dstC = 0;
            while (c < endCol) : (c += 1) {
                if (c >= strtCol) {
                    ret[((dstR * cpyCols) + dstC)] = mtx[((r * mtxCols) + c)];
                    dstC += 1;
                }
            }
            dstR += 1;
        }
    }
}

test "XMTX: cpyLessColRowXmtx3 test" {
    var m1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var ret: [3]f32 = .{ 0, 0, 0 };
    var exp1: [3]f32 = .{ 1, 2, 3 };
    var exp2: [3]f32 = .{ 7, 8, 9 };

    var start = try Instant.now();
    cpyLessColRowXmtx3(&m1, &ret, 0, 3, 0, 1, 3, 3);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncpyLessColRowXmtx3 #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cpyLessColRowXmtx3", elapsed1);

    try std.testing.expectEqual(true, equXvec(&ret, &exp1));

    start = try Instant.now();
    cpyLessColRowXmtx3(&m1, &ret, 0, 3, 2, 3, 3, 3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ncpyLessColRowXmtx3 #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cpyLessColRowXmtx3", elapsed1);

    try std.testing.expectEqual(true, equXvec(&ret, &exp2));

    //using an array with slice argument
    start = try Instant.now();
    cpyLessColRowXmtx(m1[0..9], &ret, 0, 3, 2, 3, 3, 3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ncpyLessColRowXmtx3 #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cpyLessColRowXmtx3", elapsed1);

    try std.testing.expectEqual(true, equXvec(&ret, &exp2));
    prntNl();
}

///Copies the vec vector into the ret vector argument.
///
///  vec = The vector to copy into the ret vector.
///
///  ret = The vector to hold the new copy.
///
pub fn cpyXvec(vec: []f32, ret: []f32) void {
    const l: usize = vec.len;
    var i: usize = 0;
    while (i < l) {
        ret[i] = vec[i];
        i += 1;
    }
}

test "XMTX: cpyXvec test" {
    var v1: [3]f32 = .{ 1, 2, 3 };
    var v2: [3]f32 = .{ 0, 0, 0 };
    var v3: [3]f32 = .{ 1, 1, 1 };

    var start = try Instant.now();
    cpyXvec(&v1, &v2);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncpyXvec #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cpyXvec", elapsed1);

    try std.testing.expectEqual(true, equXvec(&v1, &v2));

    start = try Instant.now();
    cpyXvec(&v1, &v3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ncpyXvec #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cpyXvec", elapsed1);

    try std.testing.expectEqual(true, equXvec(&v1, &v3));
    prntNl();
}

///Returns a newly allocated matrix of the specified length.
///
///  length = The length of the newly created f32 matrix.
///
///  alloc = The allocator to use when creating the new matrix.
///
///  returns = The newly created matrix or an error code.
///
pub fn crtXmtxEz(length: usize, alloc: *const std.mem.Allocator) ![]f32 {
    return try alloc.*.alloc(f32, length);
}

test "XMTX: crtXmtxEz test" {
    const alloc: std.mem.Allocator = std.testing.allocator;
    const l: usize = 9;
    const cols: usize = 3;

    var start = try Instant.now();
    const m1: []f32 = try crtXmtxEz(l, &alloc);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncrtXmtxEz #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("crtXmtxEz", elapsed1);

    start = try Instant.now();
    const m2: []f32 = try crtXmtxEz(l, &alloc);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ncrtXmtxEz #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("crtXmtxEz", elapsed1);

    try std.testing.expectEqual(true, isSqrXmtx(m1, cols));
    try std.testing.expectEqual(false, isIdtXmtx(m2, cols));
    alloc.free(m1);
    alloc.free(m2);

    prntNl();
}

///Returns a newly created matrix of the specified length, using the provided allocator, with the provided number of columns. Optionaly filled with the identity matrix.
///
///  length = The length of the new matrix to create.
///
///  cols = The number of columns in the new matrix.
///
///  alloc = The allocator used to create the new matrix.
///
///  identity = A Boolean indicating if the resulting new matrix should be filled with the identity matrix for that order.
///
///  returns = The newly created matrix or an error code.
///
pub fn crtXmtx(length: usize, cols: usize, alloc: *const std.mem.Allocator, identity: bool) ![]f32 {
    var ret: []f32 = try alloc.*.alloc(f32, length);
    var i: usize = 0;
    if (identity) {
        while (i < length) {
            const r = i % cols;
            const c = i / cols;
            if (r == c) {
                ret[i] = 1;
            } else {
                ret[i] = 0;
            }
            i += 1;
        }
    }
    return ret;
}

test "XMTX: crtXmtx test" {
    const alloc: std.mem.Allocator = std.testing.allocator;
    const l: usize = 9;
    const cols: usize = 3;

    var start = try Instant.now();
    const m1: []f32 = try crtXmtx(l, cols, &alloc, false);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncrtXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("crtXmtx", elapsed1);

    start = try Instant.now();
    const m2: []f32 = try crtXmtx(l, cols, &alloc, true);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ncrtXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("crtXmtx", elapsed1);

    try std.testing.expectEqual(true, isSqrXmtx(m1, cols));
    try std.testing.expectEqual(true, isIdtXmtx(m2, cols));
    alloc.free(m1);
    alloc.free(m2);

    prntNl();
}

///Returns a newly created vector of the specified length using the provided allocator.
///
///  length = The length of the new vector to create.
///
///  alloc = The allocator used to create the new vector.
///
///  returns = The newly created vector or an error code.
///
pub fn crtXvec(length: usize, alloc: *const std.mem.Allocator) ![]f32 {
    return try alloc.*.alloc(f32, length);
}

test "XMTX: crtXvec test" {
    const alloc: std.mem.Allocator = std.testing.allocator;

    var start = try Instant.now();
    const v1: []f32 = try crtXmtx(9, 3, &alloc, true);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncrtXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("crtXmtx", elapsed1);

    start = try Instant.now();
    const v2: []f32 = try crtXvec(9, &alloc);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ncrtXvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("crtXvec", elapsed1);

    try std.testing.expectEqual(true, isIdtXmtx(v1, 3));
    try std.testing.expectEqual(false, isIdtXmtx(v2, 3));
    alloc.free(v1);
    alloc.free(v2);

    prntNl();
}

///Clears the values of the given matrix.
///
///  mtx = The matrix to reset to zero for all values.
///
pub fn clrXmtx(mtx: []f32) void {
    clrXvec(mtx);
    //const l = mtx.len;
    //var i: usize = 0;
    //while (i < l) : (i += 1) {
    //    mtx[i] = 0;
    //}
}

test "XMTX: clrXmtx test" {
    var mtx: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var exp: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    prntXmtxNl(&mtx, 3);
    prntXmtxNl(&mtx, 3);

    const start = try Instant.now();
    clrXmtx(&mtx);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nclrXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("clrXmtx", elapsed1);

    try std.testing.expectEqual(true, equXvec(&mtx, &exp));
    prntNl();
}

///Clears the values of the given vector.
///
///  vec = The vector to reset to zero for all values.
///
pub fn clrXvec(vec: []f32) void {
    const l = vec.len;
    var i: usize = 0;
    while (i < l) : (i += 1) {
        vec[i] = 0;
    }
}

test "XMTX: clrXvec test" {
    var mtx: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var exp: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    prntXmtxNl(&mtx, 3);
    prntXmtxNl(&mtx, 3);

    const start = try Instant.now();
    clrXvec(&mtx);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nclrXvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("clrXvec", elapsed1);

    try std.testing.expectEqual(true, equXvec(&mtx, &exp));
    prntNl();
}

///Copies the provided matrix to a newly created matrix and returns it.
///
///  mtx = The matrix to copy into a new matrix.
///
///  alloc = The allocator to use to create the new matrix.
///
///  returns = The newly created matrix of an error code.
///
pub fn cpyXmtxNew(mtx: []f32, alloc: *const std.mem.Allocator) ![]f32 {
    return cpyXvecNew(mtx, alloc);
}

test "XMTX: cpyXmtxNew test" {
    const alloc: std.mem.Allocator = std.testing.allocator;
    var m1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    const start = try Instant.now();
    const m2: []f32 = try cpyXmtxNew(&m1, &alloc);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncpyXmtxNew: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cpyXmtxNew", elapsed1);

    try std.testing.expectEqual(true, equXvec(&m1, m2));
    alloc.free(m2);
    prntNl();
}

///Copies the provided vector to a newly created vector and returns it.
///
///  mtx = The vector to copy into a new vector.
///
///  alloc = The allocator to use to create the new vector.
///
///  returns = The newly created vector of an error code.
///
pub fn cpyXvecNew(vec: []f32, alloc: *const std.mem.Allocator) ![]f32 {
    const l: usize = vec.len;
    var i: usize = 0;
    var ret: []f32 = try alloc.*.alloc(f32, l);
    while (i < l) {
        ret[i] = vec[i];
        i += 1;
    }
    return ret;
}

test "XMTX: cpyXvecNew test" {
    const alloc: std.mem.Allocator = std.testing.allocator;
    var v1: [3]f32 = .{ 1, 2, 3 };

    const start = try Instant.now();
    const v2: []f32 = try cpyXvecNew(&v1, &alloc);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncpyXvecNew: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cpyXvecNew", elapsed1);

    try std.testing.expectEqual(true, equXvec(&v1, v2));
    alloc.free(v2);

    prntNl();
}

///Returns the magnitude of the provided vector.
///
///  vec = The vector to calculate the magnitude of.
///
///  returns = The magnitude of the given vector.
///
pub fn magXvec(vec: []f32) f32 {
    const l: usize = vec.len;
    var i: usize = 0;
    var val: f32 = 0;
    while (i < l) {
        val += (vec[i] * vec[i]);
        i += 1;
    }
    return @sqrt(val);
}

test "XMTX: magXvec test" {
    var v1: [3]f32 = .{ 1, 2, 3 };
    const exp: f32 = 3.74165749e+00;

    const start = try Instant.now();
    const val: f32 = magXvec(&v1);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nmagXvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("magXvec", elapsed1);

    try std.testing.expectEqual(exp, val);
    prntNl();
}

///Returns a Boolean value indicating if the matrix has an inverse. Does not perform a deep check. Does not check for reduction or the determinant of the matrix mtx.
///
///  mtx = The matrix to check for an inverse.
///
///  cols = The number of columns in the specified matrix.
///
///  trnMtx = The transpose of the mtx matrix.
///
///  returns = A Boolean value indicating if the mtx matrix has an inverse.
///
pub fn hasInvXmtx(mtx: []f32, cols: usize, trnMtx: []f32) bool {
    if (isZeroXmtx(mtx, cols)) {
        prntNlStr("hasInvXmtx: Exit 1");
        return false;
    }

    if (!idnfXmtx(mtx, cols, MTX_OPS.MTX_IS_LIN_INDP)) {
        prntNlStr("hasInvXmtx: Exit 2");
        return false;
    }

    if (!isSqrXmtx(mtx, cols)) {
        prntNlStr("hasInvXmtx: Exit 3");
        return false;
    }

    if (isZeroXmtx(trnMtx, cols)) {
        prntNlStr("hasInvXmtx: Exit 4");
        return false;
    }

    if (!idnfXmtx(trnMtx, cols, MTX_OPS.MTX_IS_LIN_INDP)) {
        prntNlStr("hasInvXmtx: Exit 5");
        return false;
    }

    if (!isSqrXmtx(trnMtx, cols)) {
        prntNlStr("hasInvXmtx: Exit 6");
        return false;
    }

    return true;
}

test "XMTX: hasInvXmtx test" {
    var m1: [9]f32 = .{ 3, 2, -3, 4, -3, 6, 1, 0, -1 };
    var m2: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var idtM1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var ret: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var origM1: [9]f32 = .{ 3, 2, -3, 4, -3, 6, 1, 0, -1 };

    const cols: f32 = 3;
    var b: bool = false;

    std.debug.print("\nhasInvXmtx test: initial matrix", .{});
    prntXmtxNl(&m1, 4);

    var sclr: f32 = 0.0;
    b = try rdcXmtxInl(&m1, 3, false, true, &idtM1, 3, false, &sclr);
    try std.testing.expectEqual(true, b);
    cpyXmtx(&m1, &m2);
    cpyXmtx(&origM1, &m1);

    clnXmtx(&m2);
    try std.testing.expectEqual(true, isDiagXmtx(&m2, 3));
    try std.testing.expectEqual(true, isIdtXmtx(&m2, 3));

    std.debug.print("\nClean Short Answer:", .{});
    prntXmtxNl(&m2, 3);
    prntNl();
    try std.testing.expectEqual(true, isRdcFrmXmtx(&m2, 3));

    b = try tmsXmtx(&m2, cols, &m1, cols, &ret, cols);
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&m1, &ret));

    origM1 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    cpyXmtx(&origM1, &m1);
    var trnM1: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    trnXmtx(&m1, cols, &trnM1);

    std.debug.print("\nOrig M1:", .{});
    prntXmtxNl(&origM1, 3);
    prntNl();

    std.debug.print("\nTrn M1:", .{});
    prntXmtxNl(&trnM1, 3);
    prntNl();

    const start = try Instant.now();
    const b1: bool = hasInvXmtx(&origM1, cols, &trnM1);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nhasInvXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("hasInvXmtx", elapsed1);

    try std.testing.expectEqual(true, b1);
    prntNl();
}

///Returns a Boolean value indicating if the provided matrix, with specified column count, is a symmetrical matrix.
///
///  mtx = The matrix to test for symmetry.
///
///  cols = The number of columns in the specified matrix.
///
///  returns = A Boolean value indicating if the matrix is symmetrical.
///
pub fn isSymXmtx(mtx: []f32, cols: usize) bool {
    const rows: usize = mtx.len / cols;
    var r: usize = 0;
    var c: usize = 0;
    var vL: [1]f32 = .{0};
    var vR: [1]f32 = .{0};

    while (r < rows) : (r += 1) {
        c = 0;
        while (c < cols) : (c += 1) {
            if (r != c) {
                vL = .{mtx[(r * cols) + c]};
                vR = .{mtx[(c * cols) + r]};
                if (!equXvec(&vL, &vR)) {
                    return false;
                }
            }
        }
    }
    return true;
}

test "XMTX: isSymXmtx test" {
    var m1: [4]f32 = .{ 1, 0, 0, 1 };
    const cols: usize = 2;

    const start = try Instant.now();
    const b: bool = isSymXmtx(&m1, cols);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisSymXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isSymXmtx", elapsed1);

    try std.testing.expectEqual(true, b);
    prntNl();
}

///Returns a Boolean value indicating if the provided matrix is a zero matrix.
///
///  mtx = The matrix to check for zero rows or columns.
///
///  cols = The number of columns in the matrix.
///
///  returns = A Boolean value indicating if the given matrix is a zero matrix or not.
///
pub fn isZeroXmtx(mtx: []f32, cols: usize) bool {
    const l: usize = mtx.len;
    const rows: usize = l / cols;
    var isZero: bool = true;
    var r: usize = 0;
    var c: usize = 0;

    //check rows
    while (r < rows) : (r += 1) {
        c = 0;
        isZero = true;
        while (c < cols) : (c += 1) {
            if (mtx[((r * cols) + c)] != 0) {
                isZero = false;
                break;
            }
        }

        if (isZero) {
            return true;
        }
    }

    //check cols
    c = 0;
    while (c < cols) : (c += 1) {
        r = 0;
        isZero = true;
        while (r < rows) : (r += 1) {
            if (mtx[((r * cols) + c)] != 0) {
                isZero = false;
                break;
            }
        }

        if (isZero) {
            return true;
        }
    }

    return false;
}

test "XMTX: isZeroXmtx test" {
    var m1: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var m2: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var b: bool = false;

    var start = try Instant.now();
    b = isZeroXmtx(&m1, 3);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisZeroXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isZeroXmtx", elapsed1);

    try std.testing.expectEqual(true, b);

    start = try Instant.now();
    b = isZeroXmtx(&m2, 3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisZeroXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isZeroXmtx", elapsed1);

    try std.testing.expectEqual(false, b);
    prntNl();
}

///Returns a Boolean value indicating if the provided vector is a zero vector.
///
///  vec = The vector to check for zero rows or columns.
///
///  returns = A Boolean value indicating if the given vector is a zero vector or not.
///
pub fn isZeroXvec(vec: []f32) bool {
    const l: usize = vec.len;
    var i: usize = 0;
    while (i < l) {
        if (vec[i] != 0) {
            return false;
        }
        i += 1;
    }
    return true;
}

test "XMTX: isZeroXvec test" {
    var v1: [3]f32 = .{ 0, 0, 0 };
    var v2: [3]f32 = .{ 1, 2, 3 };
    var b: bool = false;

    var start = try Instant.now();
    b = isZeroXvec(&v1);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisZeroXvec #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isZeroXvec", elapsed1);

    try std.testing.expectEqual(true, b);

    start = try Instant.now();
    b = isZeroXvec(&v2);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisZeroXvec #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isZeroXvec", elapsed1);

    try std.testing.expectEqual(false, b);
    prntNl();
}

///Returns a Boolean indicating if the provided vectors are linearly independent or not.
///
///  vecL = The left-hand side vector to compare to vecR.
///
///  vecR = The right-hand side vector to compare to vecL.
///
///  returns = A Boolean value indicating if the two vectors, vecL and vecR, are linearly independent.
///
pub fn isLinIndXvec(vecL: []f32, vecR: []f32) bool {
    if (dotPrdXvec(vecL, vecR) == 0) {
        return true;
    }
    return false;
}

test "XMTX: isLinIndXvec test" {
    var v1: [3]f32 = .{ 1, 0, 0 };
    var v2: [3]f32 = .{ 0, 1, 0 };
    var v3: [3]f32 = .{ 0, 3, 0 };
    var b: bool = false;

    var start = try Instant.now();
    b = isLinIndXvec(&v1, &v2);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisLinIndXvec #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isLinIndXvec", elapsed1);

    try std.testing.expectEqual(true, b);

    start = try Instant.now();
    b = isLinIndXvec(&v2, &v3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisLinIndXvec #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isLinIndXvec", elapsed1);

    try std.testing.expectEqual(false, b);
    prntNl();
}

///Returns a Boolean indicating if the provided vectors are linearly independent, for the given general inner product space, or not.
///
///  vecL = The left-hand side vector to compare to vecR.
///
///  vecR = The right-hand side vector to compare to vecL.
///
///  prdct = A function that takes two vectors as arguments and returns an f32 value.
///
///  returns = A Boolean value indicating if the two vectors, vecL and vecR, are linearly independent, for the given general inner product space.
///
pub fn isLinIndInrPrdctSpcXvec(vecL: []f32, vecR: []f32, prdct: *const fn (l: []f32, r: []f32) f32) bool {
    if (prdct(vecL, vecR) == 0) {
        return true;
    }
    return false;
}

test "XMTX: isLinIndInrPrdctSpcXvec test" {
    var v1: [3]f32 = .{ 1, 0, 0 };
    var v2: [3]f32 = .{ 0, 1, 0 };
    var v3: [3]f32 = .{ 0, 3, 0 };
    var b: bool = false;
    var b1: bool = false;
    const prdct = dotPrdXvec;

    var start = try Instant.now();
    b = isLinIndXvec(&v1, &v2);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisLinIndXvec #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isLinIndXvec", elapsed1);

    try std.testing.expectEqual(true, b);

    start = try Instant.now();
    b1 = isLinIndInrPrdctSpcXvec(&v1, &v2, &prdct);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisLinIndInrPrdctSpcXvec #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isLinIndInrPrdctSpcXvec", elapsed1);

    try std.testing.expectEqual(b, b1);

    start = try Instant.now();
    b = isLinIndXvec(&v2, &v3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisLinIndXvec #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isLinIndXvec", elapsed1);

    try std.testing.expectEqual(false, b);

    start = try Instant.now();
    b1 = isLinIndInrPrdctSpcXvec(&v2, &v3, &prdct);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisLinIndInrPrdctSpcXvec #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isLinIndInrPrdctSpcXvec", elapsed1);

    try std.testing.expectEqual(b, b1);
    prntNl();
}

///Returns true if the vectors of the given matrix, mtx, are orthogonal, linearly independent, using the Euclidean dot product, when tested in series.
///
///  mtx = The matrix to use to test for orthogonality.
///
///  cols = The number of columns in the matrix, mtx.
///
///  alloc = The memory allocator used to create vectors as needed.
///
///  returns = A Boolean indicating if the matrix is made up of linearly independent vectors, orthogonal or an error.
///
pub fn isOrthogonalXmtx(mtx: []f32, cols: usize, alloc: *const std.mem.Allocator) !bool {
    if(mtx.len / cols == 1) {
        return true;
    }
    return try isLinIndXmtx(mtx, cols, alloc);
}

test "XMTX: isOrthogonalXmtx test" {
    const v1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    const v2: [9]f32 = .{ 0, 1, 0, 0, 1, 0, 0, 1, 0 };
    const exp1: bool = true;
    const exp2: bool = false;
    const cols: usize = 3;
    const alloc = std.testing.allocator;
    var b1: bool = false;
    var b2: bool = false;

    var start = try Instant.now();
    b1 = try isOrthogonalXmtx(@constCast(&v1), cols, &alloc);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisOrthogonalXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isOrthogonalXmtx", elapsed1);

    start = try Instant.now();
    b2 = try isOrthogonalXmtx(@constCast(&v2), cols, &alloc);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisOrthogonalXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isOrthogonalXmtx", elapsed1);

    prntNlStrArgs("B1: {} Exp1: {} B2: {} Exp2: {}", .{ b1, exp1, b2, exp2 });
    try std.testing.expectEqual(exp1, b1);
    try std.testing.expectEqual(exp2, b2);
    prntNl();
}

///Returns true if the vectors of the given matrix, mtx, are orthogonal, linearly independent, using the given general inner product space, when tested in series.
///
///  mtx = The matrix to use to test for orthogonality.
///
///  cols = The number of columns in the matrix, mtx.
///
///  prdct = A function that takes two vectors as arguments and returns an f32 value.
///
///  alloc = The memory allocator used to create vectors as needed.
///
///  returns = A Boolean indicating if the matrix is made up of linearly independent vectors, orthogonal or an error.
///
pub fn isOrthogonalInrPrdctSpcXmtx(mtx: []f32, cols: usize, prdct: *const fn (l: []f32, r: []f32) f32, alloc: *const std.mem.Allocator) !bool {
    return try isLinIndInrPrdctSpcXmtx(mtx, cols, prdct, alloc);
}

test "XMTX: isOrthogonalInrPrdctSpcXmtx test" {
    const v1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    const v2: [9]f32 = .{ 0, 1, 0, 0, 1, 0, 0, 1, 0 };
    const exp1: bool = true;
    const exp2: bool = false;
    const cols: usize = 3;
    const alloc = std.testing.allocator;
    const prdct = dotPrdXvec;
    var b1: bool = false;
    var b2: bool = false;
    var b3: bool = false;
    var b4: bool = false;

    var start = try Instant.now();
    b1 = try isOrthogonalXmtx(@constCast(&v1), cols, &alloc);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisOrthogonalXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isOrthogonalXmtx", elapsed1);

    start = try Instant.now();
    b2 = try isOrthogonalXmtx(@constCast(&v2), cols, &alloc);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisOrthogonalXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isOrthogonalXmtx", elapsed1);

    start = try Instant.now();
    b3 = try isOrthogonalInrPrdctSpcXmtx(@constCast(&v1), cols, prdct, &alloc);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisOrthogonalInrPrdctSpcXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isOrthogonalInrPrdctSpcXmtx", elapsed1);

    start = try Instant.now();
    b4 = try isOrthogonalInrPrdctSpcXmtx(@constCast(&v2), cols, prdct, &alloc);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisOrthogonalInrPrdctSpcXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isOrthogonalInrPrdctSpcXmtx", elapsed1);

    prntNlStrArgs("B1: {} Exp1: {} B2: {} Exp2: {} B3: {} B4: {}", .{ b1, exp1, b2, exp2, b3, b4 });
    try std.testing.expectEqual(exp1, b1);
    try std.testing.expectEqual(exp2, b2);
    try std.testing.expectEqual(b3, b1);
    try std.testing.expectEqual(b4, b2);
    prntNl();
}

///Returns true if the vectors of the given matrix, mtx, are orthonormal, linearly independent, using Euclidean dot product with a magnitude of 1, when tested in series.
///
///  mtx = The matrix to use to test for orthonormality.
///
///  cols = The number of columns in the matrix, mtx.
///
///  alloc = The memory allocator used to create vectors as needed.
///
///  returns = A Boolean indicating if the matrix is made up of linearly independent vectors, orthonormal, or an error.
///
pub fn isOrthonormalXmtx(mtx: []f32, cols: usize, alloc: *const std.mem.Allocator) !bool {
    const b: bool = try isLinIndXmtx(mtx, cols, alloc);
    if (b == true) {
        const ret: []f32 = try alloc.*.alloc(f32, cols);
        defer alloc.*.free(ret);

        var i: usize = 0;
        var mag: f32 = 0;
        const l: usize = mtx.len / cols;
        var r: usize = 0;

        while (i < l) {
            clrXmtx(ret);
            cpyLessColRowXmtx(mtx, ret, 0, cols, r, (r + 1), cols, cols);
            mag = magXvec(ret);
            if (!isEquF32(mag, 1.0, false)) {
                return false;
            }
            i += cols;
            r += 1;
        }
        return true;
    }
    return false;
}

test "XMTX: isOrthonormalXmtx test" {
    const v1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    const v2: [9]f32 = .{ 0, 2, 0, 0, 2, 0, 0, 2, 0 };
    const exp1: bool = true;
    const exp2: bool = false;
    const cols: usize = 3;
    const alloc = std.testing.allocator;
    var b1: bool = false;
    var b2: bool = false;

    var start = try Instant.now();
    b1 = try isOrthonormalXmtx(@constCast(&v1), cols, &alloc);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisOrthonormalXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isOrthonormalXmtx", elapsed1);

    start = try Instant.now();
    b2 = try isOrthonormalXmtx(@constCast(&v2), cols, &alloc);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisOrthonormalXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isOrthonormalXmtx", elapsed1);

    prntNlStrArgs("B1: {} Exp1: {} B2: {} Exp2: {}", .{ b1, exp1, b2, exp2 });
    try std.testing.expectEqual(exp1, b1);
    try std.testing.expectEqual(exp2, b2);
    prntNl();
}

///Returns true if the vectors of the given matrix, mtx, are orthonormal, linearly independent, using the given general inner product space with a magnitude of 1, when tested in series.
///
///  mtx = The matrix to use to test for orthonormality.
///
///  cols = The number of columns in the matrix, mtx.
///
///  prdct = A function that takes two vectors as arguments and returns an f32 value.
///
///  alloc = The memory allocator used to create vectors as needed.
///
///  returns = A Boolean indicating if the matrix is made up of linearly independent vectors, orthonormal, or an error.
///
pub fn isOrthonormalInrPrdctSpcXmtx(mtx: []f32, cols: usize, prdct: *const fn (l: []f32, r: []f32) f32, alloc: *const std.mem.Allocator) !bool {
    const b: bool = try isLinIndInrPrdctSpcXmtx(mtx, cols, prdct, alloc);
    if (b == true) {
        const ret: []f32 = try alloc.*.alloc(f32, cols);
        defer alloc.*.free(ret);

        var i: usize = 0;
        var mag: f32 = 0;
        const l: usize = mtx.len / cols;
        var r: usize = 0;

        while (i < l) {
            clrXmtx(ret);
            cpyLessColRowXmtx(mtx, ret, 0, cols, r, (r + 1), cols, cols);
            mag = magXvec(ret);
            if (!isEquF32(mag, 1.0, false)) {
                return false;
            }
            i += cols;
            r += 1;
        }
        return true;
    }
    return false;
}

test "XMTX: isOrthonormalInrPrdctSpcXmtx test" {
    const v1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    const v2: [9]f32 = .{ 0, 2, 0, 0, 2, 0, 0, 2, 0 };
    const exp1: bool = true;
    const exp2: bool = false;
    const cols: usize = 3;
    const alloc = std.testing.allocator;
    const prdct = dotPrdXvec;
    var b1: bool = false;
    var b2: bool = false;
    var b3: bool = false;
    var b4: bool = false;

    var start = try Instant.now();
    b1 = try isOrthonormalXmtx(@constCast(&v1), cols, &alloc);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisOrthonormalXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isOrthonormalXmtx", elapsed1);

    start = try Instant.now();
    b2 = try isOrthonormalXmtx(@constCast(&v2), cols, &alloc);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisOrthonormalXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isOrthonormalXmtx", elapsed1);

    start = try Instant.now();
    b3 = try isOrthonormalInrPrdctSpcXmtx(@constCast(&v1), cols, prdct, &alloc);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisOrthonormalInrPrdctSpcXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isOrthonormalInrPrdctSpcXmtx", elapsed1);

    start = try Instant.now();
    b4 = try isOrthonormalInrPrdctSpcXmtx(@constCast(&v2), cols, prdct, &alloc);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisOrthonormalInrPrdctSpcXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isOrthonormalInrPrdctSpcXmtx", elapsed1);

    prntNlStrArgs("B1: {} Exp1: {} B2: {} Exp2: {} B3: {} B4: {}", .{ b1, exp1, b2, exp2, b3, b4 });
    try std.testing.expectEqual(exp1, b1);
    try std.testing.expectEqual(exp2, b2);
    try std.testing.expectEqual(b3, b1);
    try std.testing.expectEqual(b4, b2);
    prntNl();
}

///Returns true if the vectors of the given matrix, mtx, are linearly independent, when tested in series.
///
///  mtx = The matrix to use to test for linear independence.
///
///  cols = The number of columns in the matrix, mtx, and the vectors vecL and vecR.
///
///  alloc = The memory allocator used to create vectors to hold the left and right comparison of matrix rows.
///
///  returns = A Boolean indicating if the matrix is made up of linearly independent vectors.
///
pub fn isLinIndXmtx(mtx: []f32, cols: usize, alloc: *const std.mem.Allocator) !bool {
    var vecL: []f32 = undefined;
    var vecR: []f32 = undefined;

    vecL = crtXvec(cols, alloc) catch |err| {
        std.debug.print("\nisLinIndXmtx: Error creating new vector vecL: {}", .{err});
        return err;
    };

    vecR = crtXvec(cols, alloc) catch |err| {
        std.debug.print("\nisLinIndXmtx: Error creating new vector vecR: {}", .{err});
        alloc.*.free(vecL);
        return err;
    };

    const b: bool = isLinIndXmtxRef(mtx, vecL, vecR, cols);
    alloc.free(vecL);
    alloc.free(vecR);
    return b;
}

test "XMTX: isLinIndXmtx test" {
    const alloc: std.mem.Allocator = std.testing.allocator;
    var mtx: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var b: bool = false;

    const start = try Instant.now();
    b = try isLinIndXmtx(&mtx, 3, &alloc);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisLinIndXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isLinIndXmtx", elapsed1);

    try std.testing.expectEqual(true, b);
    prntNl();
}

///Returns true if the vectors of the given matrix, mtx, are linearly independent, for the given general inner product space, when tested in series.
///
///  mtx = The matrix to use to test for linear independence.
///
///  cols = The number of columns in the matrix, mtx, and the vectors vecL and vecR.
///
///  prdct = A function that takes two vectors as arguments and returns an f32 value.
///
///  alloc = The memory allocator used to create vectors to hold the left and right comparison of matrix rows.
///
///  returns = A Boolean indicating if the matrix is made up of linearly independent vectors.
///
pub fn isLinIndInrPrdctSpcXmtx(mtx: []f32, cols: usize, prdct: *const fn (l: []f32, r: []f32) f32, alloc: *const std.mem.Allocator) !bool {
    var vecL: []f32 = undefined;
    var vecR: []f32 = undefined;

    vecL = crtXvec(cols, alloc) catch |err| {
        std.debug.print("\nisLinIndInrPrdctSpcXmtx: Error creating new vector vecL: {}", .{err});
        return err;
    };

    vecR = crtXvec(cols, alloc) catch |err| {
        std.debug.print("\nisLinIndInrPrdctSpcXmtx: Error creating new vector vecR: {}", .{err});
        alloc.*.free(vecL);
        return err;
    };

    const b: bool = isLinIndInrPrdctSpcXmtxRef(mtx, vecL, vecR, cols, prdct);
    alloc.free(vecL);
    alloc.free(vecR);
    return b;
}

test "XMTX: isLinIndInrPrdctSpcXmtx test" {
    const alloc: std.mem.Allocator = std.testing.allocator;
    var mtx: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var b: bool = false;
    var b1: bool = false;

    var start = try Instant.now();
    b = try isLinIndXmtx(&mtx, 3, &alloc);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisLinIndXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isLinIndXmtx", elapsed1);

    try std.testing.expectEqual(true, b);
    prntNl();

    start = try Instant.now();
    b1 = try isLinIndInrPrdctSpcXmtx(&mtx, 3, dotPrdXvec, &alloc);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisLinIndInrPrdctSpcXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isLinIndInrPrdctSpcXmtx", elapsed1);

    try std.testing.expectEqual(b1, b);
    prntNl();
}

///Returns true if the vectors of the given matrix, mtx, are linearly independent when tested in series.
///
///  mtx = The matrix to use to test for linear independence.
///
///  vecL = The vector that holds the left comparison row.
///
///  vecR = The vector that holds the right comparison row.
///
///  cols = The number of columns in the matrix, mtx, and the vectors vecL and vecR.
///
///  returns = A Boolean indicating if the matrix is made up of linearly independent vectors.
///
pub fn isLinIndXmtxRef(mtx: []f32, vecL: []f32, vecR: []f32, cols: usize) bool {
    const llen: usize = mtx.len;
    const rows: usize = llen / cols;
    var i: usize = 0;
    var j: i64 = -1;
    var vecLset: bool = false;
    var vecRset: bool = false;
    var lind: bool = false;    

    while (j < (rows - 1)) {
        i = @as(usize, @intCast((j + 1)));
        vecLset = false;
        vecRset = false;

        while (i < (rows - 0)) {
            if (vecLset == false) {
                cpyLessColRowXmtx(mtx, vecL, 0, cols, i, i + 1, cols, cols);
                vecLset = true;
            } else if (vecRset == false) {
                cpyLessColRowXmtx(mtx, vecR, 0, cols, i, i + 1, cols, cols);
                vecRset = true;
            }

            if (vecLset and vecRset) {
                if (VERBOSE) {
                    prntNlStr("\nisLinIndXmtx: Compare L and R:");
                    prntXvecNl(vecL);
                    prntNl();
                    prntXvecNl(vecR);
                    prntNl();
                }

                lind = isLinIndXvec(vecL, vecR);
                if (!lind) {
                    return false;
                } else if (lind) {
                    cpyXvec(vecR, vecL);
                    clrXvec(vecR);
                    vecLset = true;
                    vecRset = false;
                }
            } else {
                lind = false;
            }

            i += 1;
        }
        j += 1;
    }

    return true;
}

test "XMTX: isLinIndXmtxRef test" {
    var mtx: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var vecL: [3]f32 = .{ 0, 0, 0 };
    var vecR: [3]f32 = .{ 0, 0, 0 };
    const cols: usize = 3;
    var b: bool = false;

    const start = try Instant.now();
    b = isLinIndXmtxRef(&mtx, &vecL, &vecR, cols);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisLinIndXmtxRef: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isLinIndXmtxRef", elapsed1);

    try std.testing.expectEqual(true, b);
    prntNl();
}

///Returns true if the vectors of the given matrix, mtx, are linearly independent, for the given general inner product space, when tested in series.
///
///  mtx = The matrix to use to test for linear independence.
///
///  vecL = The vector that holds the left comparison row.
///
///  vecR = The vector that holds the right comparison row.
///
///  cols = The number of columns in the matrix, mtx, and the vectors vecL and vecR.
///
///  prdct = A function that takes two vectors as arguments and returns an f32 value.
///
///  returns = A Boolean indicating if the matrix is made up of linearly independent vectors.
///
pub fn isLinIndInrPrdctSpcXmtxRef(mtx: []f32, vecL: []f32, vecR: []f32, cols: usize, prdct: *const fn (l: []f32, r: []f32) f32) bool {
    const llen: usize = mtx.len;
    const rows: usize = llen / cols;
    var i: usize = 0;
    var j: i64 = -1;
    var vecLset: bool = false;
    var vecRset: bool = false;
    var lind: bool = false;

    while (j < (rows - 1)) {
        i = @as(usize, @intCast((j + 1)));
        vecLset = false;
        vecRset = false;

        while (i < (rows - 1)) {
            if (vecLset == false) {
                cpyLessColRowXmtx(mtx, vecL, 0, cols, i, i + 1, cols, cols);
                vecLset = true;
            } else if (vecRset == false) {
                cpyLessColRowXmtx(mtx, vecR, 0, cols, i, i + 1, cols, cols);
                vecRset = true;
            }

            if (vecLset and vecRset) {
                if (VERBOSE) {
                    prntNlStr("\nisLinIndXmtx: Compare L and R:");
                    prntXvecNl(vecL);
                    prntNl();
                    prntXvecNl(vecR);
                    prntNl();
                }

                lind = isLinIndInrPrdctSpcXvec(vecL, vecR, prdct);
                if (!lind) {
                    return false;
                } else if (lind) {
                    cpyXvec(vecR, vecL);
                    clrXvec(vecR);
                    vecLset = true;
                    vecRset = false;
                }
            } else {
                lind = false;
            }

            i += 1;
        }
        j += 1;
    }

    return true;
}

test "XMTX: isLinIndInrPrdctSpcXmtxRef test" {
    var mtx: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var vecL: [3]f32 = .{ 0, 0, 0 };
    var vecR: [3]f32 = .{ 0, 0, 0 };
    const cols: usize = 3;
    var b: bool = false;
    var b1: bool = false;

    var start = try Instant.now();
    b = isLinIndXmtxRef(&mtx, &vecL, &vecR, cols);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisLinIndXmtxRef: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isLinIndXmtxRef", elapsed1);

    try std.testing.expectEqual(true, b);
    prntNl();

    start = try Instant.now();
    b1 = isLinIndInrPrdctSpcXmtxRef(&mtx, &vecL, &vecR, cols, dotPrdXvec);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisLinIndPrdctSpcXmtxRef: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isLinIndPrdctSpcXmtxRef", elapsed1);

    try std.testing.expectEqual(b1, b);
    prntNl();
}

///Returns the dot product of the two provided vectors, vecL and vecR.
///
///  vecL = The left-hand side vector in the dot product calculation.
///
///  vecR = The right-hand side vector in the dot product calculation.
///
///  returns = The dot product of the two vectors.
///
pub fn dotPrdXvec(vecL: []f32, vecR: []f32) f32 {
    const l: usize = vecL.len;
    var i: usize = 0;
    var val: f32 = 0;
    while (i < l) {
        val += (vecL[i] * vecR[i]);
        i += 1;
    }
    return val;
}

test "XMTX: dotPrdXvec test" {
    var v1: [3]f32 = .{ 1, 0, 0 };
    var v2: [3]f32 = .{ 0, 1, 0 };
    var v3: [3]f32 = .{ 0, 3, 0 };
    var v4: f32 = -1.0;

    var start = try Instant.now();
    v4 = dotPrdXvec(&v1, &v2);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ndotPrdXvec #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("dotPrdXvec", elapsed1);

    try std.testing.expectEqual(true, (v4 == 0));

    start = try Instant.now();
    v4 = dotPrdXvec(&v2, &v3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ndotPrdXvec #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("dotPrdXvec", elapsed1);

    try std.testing.expectEqual(true, (v4 > 0));
    prntNl();
}

///Returns the dot product of the two 3 column vectors provided, vecL and vecR.
///
///  vecL = The left-hand side reference to a vector 3 in the dot product calculation.
///
///  vecR = The right-hand side reference to a vector 3 in the dot product calculation.
///
///  returns = The dot product of the two vectors.
///
pub fn dotPrdXvec3(vecL: *const [3]f32, vecR: *const [3]f32) f32 {
    const l: usize = vecL.len;
    var i: usize = 0;
    var val: f32 = 0;
    while (i < l) {
        val += (vecL[i] * vecR[i]);
        i += 1;
    }
    return val;
}

test "XMTX: dotPrdXvec3 test" {
    var v1: [3]f32 = .{ 1, 0, 0 };
    var v2: [3]f32 = .{ 0, 1, 0 };
    var v3: [3]f32 = .{ 0, 3, 0 };
    var v4: f32 = -1.0;

    var start = try Instant.now();
    v4 = dotPrdXvec3(&v1, &v2);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ndotPrdXvec #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("dotPrdXvec", elapsed1);

    try std.testing.expectEqual(true, (v4 == 0));

    start = try Instant.now();
    v4 = dotPrdXvec3(&v2, &v3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ndotPrdXvec #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("dotPrdXvec", elapsed1);

    try std.testing.expectEqual(true, (v4 > 0));
    prntNl();
}

///Returns the cross-product of the two provided 3x3 vectors, vecL and vecR.
///
///  vecL = The left-hand vector in the cross product calculation.
///
///  vecR = The right-hand vector in the cross product calculation.
///
///  returns = The cross product of the two 3x3 vectors provided.
///
pub fn crsPrdXvec3(vecL: *const [3]f32, vecR: *const [3]f32) [3]f32 {
    return [3]f32{ (vecL[1] * vecR[2]) - (vecL[2] * vecR[1]), (vecL[2] * vecR[0]) - (vecL[0] * vecR[2]), (vecL[0] * vecR[1]) - (vecL[1] * vecR[0]) };
}

test "XMTX: crsPrdXvec3 test" {
    var v1: [3]f32 = .{ 1, 0, 0 };
    var v2: [3]f32 = .{ 0, 1, 0 };
    var v3: [3]f32 = .{ 0, 0, 1 };
    var v4: [3]f32 = .{ 0, 0, 0 };
    var v5: [3]f32 = .{ 0, 0, 0 };

    var start = try Instant.now();
    v4 = crsPrdXvec3(&v1, &v3);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncrsPrdXvec3 #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("crsPrdXvec3", elapsed1);

    start = try Instant.now();
    v5 = crsPrdXvec3(&v1, &v2);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ncrsPrdXvec3 #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("crsPrdXvec3", elapsed1);

    var exp1: [3]f32 = .{ 0, -1, 0 };
    var exp2: [3]f32 = .{ 0, 0, 1 };
    prntXvecNl(&v4);
    prntXvecNl(&v5);

    try std.testing.expectEqual(true, equXvec(&exp1, &v4));
    try std.testing.expectEqual(true, equXvec(&exp2, &v5));
    prntNl();
}

pub fn crsPrdXvec(vecL: []f32, vecR: []f32) ![3]f32 {
    if (vecL.len != 3) {
        prntNlStr("!! Warning, the crsPrdXvec function expect vecL to be a vector 3 !!");
        return Error.InvalidLengths;
    } else if (vecR.len != 3) {
        prntNlStr("!! Warning, the crsPrdXvec function expect vecR to be a vector 3 !!");
        return Error.InvalidLengths;
    }
    //return crsPrdXvec3(vecL[0..3], vecR[0..3]);
    return [3]f32{ (vecL[1] * vecR[2]) - (vecL[2] * vecR[1]), (vecL[2] * vecR[0]) - (vecL[0] * vecR[2]), (vecL[0] * vecR[1]) - (vecL[1] * vecR[0]) };
}

test "XMTX: crsPrdXvec test" {
    var v1: [3]f32 = .{ 1, 0, 0 };
    var v2: [3]f32 = .{ 0, 1, 0 };
    var v3: [3]f32 = .{ 0, 0, 1 };
    var v4: [3]f32 = .{ 0, 0, 0 };
    var v5: [3]f32 = .{ 0, 0, 0 };

    var start = try Instant.now();
    v4 = try crsPrdXvec(&v1, &v3);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncrsPrdXvec #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("crsPrdXvec", elapsed1);

    start = try Instant.now();
    v5 = try crsPrdXvec(&v1, &v2);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ncrsPrdXvec #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("crsPrdXvec", elapsed1);

    var exp1: [3]f32 = .{ 0, -1, 0 };
    var exp2: [3]f32 = .{ 0, 0, 1 };
    prntXvecNl(&v4);
    prntXvecNl(&v5);

    try std.testing.expectEqual(true, equXvec(&exp1, &v4));
    try std.testing.expectEqual(true, equXvec(&exp2, &v5));
    prntNl();
}

///Multiplies each entry in the vec vector by the provided scalar value.
///
///  vec = The vector to apply the multiplication to.
///
///  val = The scalar value to multiply the vector by.
///
pub fn mulXvec(vec: []f32, val: f32) void {
    const l: usize = vec.len;
    var i: usize = 0;
    while (i < l) {
        vec[i] *= val;
        i += 1;
    }
}

test "XMTX: mulXvec test" {
    var v1: [3]f32 = .{ 1, 1, 1 };
    var v2: [3]f32 = .{ 21, 21, 21 };
    const a: f32 = 21;

    const start = try Instant.now();
    mulXvec(&v1, a);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nmulXvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("mulXvec", elapsed1);

    prntXvecNl(&v1);
    try std.testing.expectEqual(true, equXvec(&v1, &v2));
    prntNl();
}

///Divides each entry in the vec vector by the provided value.
///
///  vec = The vector to apply the division to.
///
///  val = The scalar value to divide the vector by.
///
pub fn divXvec(vec: []f32, val: f32) void {
    const l: usize = vec.len;
    var i: usize = 0;
    while (i < l) {
        vec[i] /= val;
        i += 1;
    }
}

test "XMTX: divXvec test" {
    var v1: [3]f32 = .{ 1, 1, 1 };
    var v2: [3]f32 = .{ 4.76190485e-02, 4.76190485e-02, 4.76190485e-02 };
    const a: f32 = 21;

    const start = try Instant.now();
    divXvec(&v1, a);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ndivXvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("divXvec", elapsed1);

    prntXvecNl(&v1);
    try std.testing.expectEqual(true, equXvec(&v1, &v2));
    prntNl();
}

///Adds the value, val, to the data in the vector, vec.
///
///  vec = The vector to apply the addition to.
///
///  val = The scalar value to add to the vector data.
///
pub fn addXvec(vec: []f32, val: f32) void {
    const l: usize = vec.len;
    var i: usize = 0;
    while (i < l) {
        vec[i] += val;
        i += 1;
    }
}

test "XMTX: addXvec test" {
    var v1: [3]f32 = .{ 1, 1, 1 };
    var v2: [3]f32 = .{ 3, 3, 3 };
    const a: f32 = 2;

    const start = try Instant.now();
    addXvec(&v1, a);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\naddXvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("addXvec", elapsed1);

    prntXvecNl(&v1);
    try std.testing.expectEqual(true, equXvec(&v1, &v2));
    prntNl();
}

///Subtracts each entry in the vec vector by the provided value.
///
///  vec = The vector to apply the subtract to.
///
///  val = The scalar value to substract the vector by.
///
pub fn subXvec(vec: []f32, val: f32) void {
    const l: usize = vec.len;
    var i: usize = 0;
    while (i < l) {
        vec[i] -= val;
        i += 1;
    }
}

test "XMTX: subXvec test" {
    var v1: [3]f32 = .{ 3, 3, 3 };
    var v2: [3]f32 = .{ 1, 1, 1 };
    const a: f32 = 2;

    const start = try Instant.now();
    subXvec(&v1, a);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nsubXvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("subXvec", elapsed1);

    prntXvecNl(&v1);
    try std.testing.expectEqual(true, equXvec(&v1, &v2));
    prntNl();
}

///Performs vector multiplication on the provided vectors storing the result in the ret vector.
///
///  vecL = The left-hand vector in the multiplication process.
///
///  vecR = The right-hand vector in the multiplication process.
///
///  ret = The return vector where calculated results are stored.
///
///  returns = A Boolean value indicating if the vector multiplication operation was successful or an error value.
///
pub fn tmsXvec(vecL: []f32, vecR: []f32, ret: []f32) !bool {
    if (vecL.len != vecR.len) {
        prntNlStr("!! Warning tmsXvec vector arguments are not the same legnth !!");
        return Error.InvalidLengths;
    }

    if (ret.len != 1) {
        prntNlStr("!! Warning tmsXvec returns a vector of length 1 !!");
        return Error.InvalidLengths;
    }

    const l: usize = vecL.len;
    var i: usize = 0;
    while (i < l) {
        ret[0] += (vecL[i] * vecR[i]);
        i += 1;
    }
    return true;
}

test "XMTX: tmsXvec test" {
    var v1: [3]f32 = .{ 3, 3, 3 };
    var v2: [3]f32 = .{ 2, 2, 2 };
    var v3: [1]f32 = .{0};
    var exp: [1]f32 = .{18};

    const start = try Instant.now();
    _ = try tmsXvec(&v1, &v2, &v3);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ntmsXvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("tmsXvec", elapsed1);

    prntXvecNl(&v1);
    prntXvecNl(&v2);
    prntXvecNl(&v3);
    try std.testing.expectEqual(true, equXvec(&exp, &v3));
    prntNl();
}

///Performs vector multiplication on the provided vectors storing the result in the ret vector.
///
///  vecL = The left-hand vector in the multiplication process.
///
///  vecR = The right-hand vector in the multiplication process.
///
///  ret = The return vector where calculated results are stored.
///
///  returns = A Boolean value indicating if the vector multiplication operation was successful or an error value.
///
pub fn tmsXvecSclr(vecL: []f32, vecR: []f32, ret: *f32) !bool {
    if (vecL.len != vecR.len) {
        prntNlStr("!! Warning tmsXvec vector arguments are not the same legnth !!");
        return Error.InvalidLengths;
    }

    const l: usize = vecL.len;
    var i: usize = 0;
    while (i < l) {
        ret.* += (vecL[i] * vecR[i]);
        i += 1;
    }
    return true;
}

test "XMTX: tmsXvecSclr test" {
    var v1: [3]f32 = .{ 3, 3, 3 };
    var v2: [3]f32 = .{ 2, 2, 2 };
    var v3: f32 = 0.0;
    const exp: f32 = 18.0;

    const start = try Instant.now();
    _ = try tmsXvecSclr(&v1, &v2, &v3);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ntmsXvecSclr: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("tmsXvecSclr", elapsed1);

    prntXvecNl(&v1);
    prntXvecNl(&v2);
    prntNlStrArgs("v3: {}", .{v3});
    try std.testing.expectEqual(exp, v3);
    prntNl();
}

///Performs vector multiplication on the provided 1 vectors storing the result in the return 1 vector, ret.
///
///  vecL = The left-hand 1 vector in the multiplication process.
///
///  vecR = The right-hand 1 vector in the multiplication process.
///
///  returns = The return 1 vector where calculated results are stored.
///
pub fn tmsXvec1(vecL: *const [1]f32, vecR: *const [1]f32) f32 {
    return (vecL[0] * vecR[0]);
}

test "XMTX: tmsXvec1 test" {
    var v1: [1]f32 = .{3};
    var v2: [1]f32 = .{2};
    var v3: f32 = 0.0;
    const exp: f32 = 6.0;

    const start = try Instant.now();
    v3 = tmsXvec1(&v1, &v2);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ntmsXvec1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("tmsXvec1", elapsed1);

    prntXvecNl(&v1);
    prntXvecNl(&v2);
    prntNlStrArgs("v3: {}", .{v3});
    try std.testing.expectEqual(exp, v3);
    prntNl();
}

///Performs vector multiplication on the provided 1x1 vectors storing the result in the return 2 vector, ret.
///
///  vecL = The left-hand 2 vector in the multiplication process.
///
///  vecR = The right-hand 2 vector in the multiplication process.
///
///  returns = The return 2 vector where calculated results are stored.
///
pub fn tmsXvec2(vecL: *const [2]f32, vecR: *const [2]f32) f32 {
    return ((vecL[0] * vecR[0]) + (vecL[1] * vecR[1]));
}

test "XMTX: tmsXvec2 test" {
    var v1: [2]f32 = .{ 3, 3 };
    var v2: [2]f32 = .{ 2, 2 };
    var v3: f32 = 0.0;
    const exp: f32 = 12.0;

    const start = try Instant.now();
    v3 = tmsXvec2(&v1, &v2);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ntmsXvec2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("tmsXvec2", elapsed1);

    prntXvecNl(&v1);
    prntXvecNl(&v2);
    prntNlStrArgs("v3: {}", .{v3});
    try std.testing.expectEqual(exp, v3);
    prntNl();
}

///Performs vector multiplication on the provided 1x1 vectors storing the result in the return 3 vector, ret.
///
///  vecL = The left-hand 3 vector in the multiplication process.
///
///  vecR = The right-hand 3 vector in the multiplication process.
///
///  returns = The return 3 vector where calculated results are stored.
///
pub fn tmsXvec3(vecL: *const [3]f32, vecR: *const [3]f32) f32 {
    return ((vecL[0] * vecR[0]) + (vecL[1] * vecR[1]) + (vecL[2] * vecR[2]));
}

test "XMTX: tmsXvec3 test" {
    var v1: [3]f32 = .{ 3, 3, 3 };
    var v2: [3]f32 = .{ 2, 2, 2 };
    var v3: f32 = 0.0;
    const exp: f32 = 18.0;

    const start = try Instant.now();
    v3 = tmsXvec3(&v1, &v2);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ntmsXvec3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("tmsXvec3", elapsed1);

    prntXvecNl(&v1);
    prntXvecNl(&v2);
    prntNlStrArgs("v3: {}", .{v3});
    try std.testing.expectEqual(exp, v3);
    prntNl();
}

///Performs vector multiplication on the provided 1x1 vectors storing the result in the return 4 vector, ret.
///
///  vecL = The left-hand 4 vector in the multiplication process.
///
///  vecR = The right-hand 4 vector in the multiplication process.
///
///  returns = The return 4 vector where calculated results are stored.
///
pub fn tmsXvec4(vecL: *const [4]f32, vecR: *const [4]f32) f32 {
    return ((vecL[0] * vecR[0]) + (vecL[1] * vecR[1]) + (vecL[2] * vecR[2]) + (vecL[3] * vecR[3]));
}

test "XMTX: tmsXvec4 test" {
    var v1: [4]f32 = .{ 3, 3, 3, 3 };
    var v2: [4]f32 = .{ 2, 2, 2, 2 };
    var v3: f32 = 0.0;
    const exp: f32 = 24.0;

    const start = try Instant.now();
    v3 = tmsXvec4(&v1, &v2);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ntmsXvec4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("tmsXvec4", elapsed1);

    prntXvecNl(&v1);
    prntXvecNl(&v2);
    prntNlStrArgs("v3: {}", .{v3});
    try std.testing.expectEqual(exp, v3);
    prntNl();
}

///Performs matrix multiplication on the provided matrices with results being stored in the ret matrix.
///
///  mtxL = The matrix on the left-hand side of the matrix multiplication process.
///
///  colsL = The number of columns of the left-hand matrix.
///
///  mtxR = The matrix on the right-hand side of thematrix  mupltiplication process.
///
///  colsR = The number of columns of the right-hand matrix.
///
///  ret = The return matrix that holds all of the calculated values.
///
///  colsRet = The number of columns of the return matrix.
///
///  returns = A Boolean value that indicates if the matrix multiplication succeeded.
///
pub fn tmsXmtx(mtxL: []f32, colsL: usize, mtxR: []f32, colsR: usize, ret: []f32, colsRet: usize) !bool {
    //m X p * p X n = m X n
    //2x2 * 2x3 = 2*3
    //L: row X col * R: row X col = ret: row * col
    const isMtxLSqr: bool = isSqrXmtx(mtxL, colsL);
    const isMtxRSqr: bool = isSqrXmtx(mtxR, colsR);
    const isRetSqr: bool = isSqrXmtx(ret, colsRet);
    const rowsMtxL: usize = mtxL.len / colsL;
    const rowsMtxR: usize = mtxR.len / colsR;
    const rowsMtxRet: usize = ret.len / colsRet;

    if (colsL != rowsMtxR) {
        prntNlStrArgs("rowsL: {} colsL: {} rowsR: {} colsR: {} rowRet: {} colsRet: {}", .{ rowsMtxL, colsL, rowsMtxR, colsR, rowsMtxRet, colsRet });
        prntNlStr("!! Warning tmsXmtx expects the following argument matrix size, m X p * p X n = m X n !!");
        return Error.InvalidLengths;
    }

    if (rowsMtxRet != rowsMtxL) {
        prntNlStrArgs("rowsL: {} colsL: {} rowsR: {} colsR: {} rowRet: {} colsRet: {}", .{ rowsMtxL, colsL, rowsMtxR, colsR, rowsMtxRet, colsRet });
        prntNlStr("!! Warning tmsXmtx expects the following argument matrix size, m X p * p X n = m X n !!");
        return Error.InvalidLengths;
    }

    if (isMtxLSqr and isMtxRSqr) {
        if (!isRetSqr) {
            prntNlStrArgs("rowsL: {} colsL: {} rowsR: {} colsR: {} rowRet: {} colsRet: {}\n", .{ rowsMtxL, colsL, rowsMtxR, colsR, rowsMtxRet, colsRet });
            prntNlStr("!! Warning tmsXmtx expects the return matrix to be square !!");
            return Error.InvalidLengths;
        } else if (colsL != colsR or colsL != colsRet) {
            prntNlStrArgs("rowsL: {} colsL: {} rowsR: {} colsR: {} rowRet: {} colsRet: {}", .{ rowsMtxL, colsL, rowsMtxR, colsR, rowsMtxRet, colsRet });
            prntNlStr("\n!! Warning tmsXmtx expects the matrix columns to match when working with square matrices !!");
            return Error.InvalidLengths;
        } else if (rowsMtxL != rowsMtxR or rowsMtxL != rowsMtxRet) {
            prntNlStrArgs("rowsL: {} colsL: {} rowsR: {} colsR: {} rowRet: {} colsRet: {}", .{ rowsMtxL, colsL, rowsMtxR, colsR, rowsMtxRet, colsRet });
            prntNlStr("\n!! Warning tmsXmtx expects the matrix rows to match when working with square matrices !!");
            return Error.InvalidLengths;
        }
    }

    var r: usize = 0;
    var c: usize = 0;
    var idx: usize = 0;
    while (r < rowsMtxRet) : (r += 1) {
        c = 0;
        while (c < colsRet) : (c += 1) {
            idx = (r * colsRet) + c;
            var k: usize = 0;
            var v: f32 = 0;
            while (k < colsL) : (k += 1) {
                //Crc = sum(k = 1 to n) { Ark * Bkc}
                const aRK = mtxL[((r * colsL) + k)];
                const bKC = mtxR[((k * colsR) + c)];

                //std.debug.print("aRK: {} bKC: {}\n", .{aRK, bKC});
                v += aRK * bKC;
            }
            ret[idx] = v;
        }
    }

    return true;
}

test "XMTX: tmsXmtx test" {
    var m1: [9]f32 = .{ 3, 2, -3, 4, -3, 6, 1, 0, -1 };
    var m2: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var idtM1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var ret: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var origM1: [9]f32 = .{ 3, 2, -3, 4, -3, 6, 1, 0, -1 };
    const cols: f32 = 3;
    var b: bool = false;

    std.debug.print("\ntmsXmtx test: initial matrix", .{});
    prntXmtxNl(&m1, 4);

    var sclr: f32 = 0.0;
    b = try rdcXmtxInl(&m1, 3, false, true, &idtM1, 3, false, &sclr);
    try std.testing.expectEqual(true, b);
    cpyXmtx(&m1, &m2);
    cpyXmtx(&origM1, &m1);

    clnXmtx(&m2);
    try std.testing.expectEqual(true, isDiagXmtx(&m2, 3));
    try std.testing.expectEqual(true, isIdtXmtx(&m2, 3));

    std.debug.print("\nClean Short Answer:", .{});
    prntXmtxNl(&m2, 3);
    prntNl();
    try std.testing.expectEqual(true, isRdcFrmXmtx(&m2, 3));

    const start = try Instant.now();
    b = try tmsXmtx(&m2, cols, &m1, cols, &ret, cols);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ntmsXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("tmsXmtx", elapsed1);

    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&m1, &ret));
    prntNl();
}

///Adds the two provided vectors together storing the adjustment in the vecL argument.
///
///  vecL = The left-hand side of the matrix addition process.
///
///  vecR = The right-hand side of the matrix addition process.
///
pub fn sum1Xvec(vecL: []f32, vecR: []f32) void {
    const l: usize = vecL.len;
    var i: usize = 0;
    while (i < l) {
        vecL[i] += vecR[i];
        i += 1;
    }
}

test "XMTX: sum1Xvec test" {
    var v1: [3]f32 = .{ 1, 1, 1 };
    var v2: [3]f32 = .{ 2, 2, 2 };
    var v3: [3]f32 = .{ 3, 3, 3 };

    const start = try Instant.now();
    sum1Xvec(&v1, &v2);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nsum1Xvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("sum1Xvec", elapsed1);

    prntXvecNl(&v1);
    try std.testing.expectEqual(true, equXvec(&v1, &v3));
    prntNl();
}

///Adds the three provided vectors together storing the adjustment in the vecL argument.
///
///  vecL = The left-hand side of the matrix addition process.
///
///  vec1 = The first matrix to add to vecL.
///
///  vec2 = The second matrix to add to vecL.
///
pub fn sum2Xvec(vecL: []f32, vec1: []f32, vec2: []f32) void {
    const l: usize = vecL.len;
    var i: usize = 0;
    while (i < l) {
        vecL[i] += (vec1[i] + vec2[i]);
        i += 1;
    }
}

test "XMTX: sum2Xvec test" {
    var v1: [3]f32 = .{ 1, 1, 1 };
    var v2: [3]f32 = .{ 2, 2, 2 };
    var v3: [3]f32 = .{ 5, 5, 5 };

    const start = try Instant.now();
    sum2Xvec(&v1, &v2, &v2);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nsum2Xvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("sum2Xvec", elapsed1);

    prntXvecNl(&v1);
    try std.testing.expectEqual(true, equXvec(&v1, &v3));
    prntNl();
}

///Subtracts the two provided vectors from each other storing the adjustment in the vecL argument.
///
///  vecL = The left-hand side of the matrix subtraction process.
///
///  vecR = The right-hand side of the matrix subtraction process.
pub fn diff1Xvec(vecL: []f32, vecR: []f32) void {
    const l: usize = vecL.len;
    var i: usize = 0;
    while (i < l) {
        vecL[i] -= vecR[i];
        i += 1;
    }
}

test "XMTX: diff1Xvec test" {
    var v1: [3]f32 = .{ 3, 3, 3 };
    var v2: [3]f32 = .{ 1, 1, 1 };
    var v3: [3]f32 = .{ 2, 2, 2 };

    const start = try Instant.now();
    diff1Xvec(&v1, &v2);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ndiff1Xvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("diff1Xvec", elapsed1);

    prntXvecNl(&v1);
    try std.testing.expectEqual(true, equXvec(&v1, &v3));
    prntNl();
}

///Subtracts the three provided vectors from eachother storing the adjustment in the vecL argument.
///
///  vecL = The left-hand side of the matrix subtraction process.
///
///  vec1 = The first matrix to subtract to vecL.
///
///  vec2 = The second matrix to subtract to vecL.
///
pub fn diff2Xvec(vecL: []f32, vec1: []f32, vec2: []f32) void {
    const l: usize = vecL.len;
    var i: usize = 0;
    while (i < l) {
        vecL[i] += (vec1[i] - vec2[i]);
        i += 1;
    }
}

test "XMTX: diff2Xvec test" {
    var v1: [3]f32 = .{ 3, 3, 3 };
    var v2: [3]f32 = .{ 1, 1, 1 };
    var v3: [3]f32 = .{ 3, 3, 3 };

    const start = try Instant.now();
    diff2Xvec(&v1, &v2, &v2);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ndiff2Xvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("diff2Xvec", elapsed1);

    prntXvecNl(&v1);
    try std.testing.expectEqual(true, equXvec(&v1, &v3));
    prntNl();
}

///Prints the matrix, followed by a new line, with the specified number of columns.
///
///  mtx = The matrix to print.
///
///  cols = The number of columns in the matrix.
///
pub fn prntXmtx(mtx: []f32, cols: usize) void {
    //std.debug.print("", .{});
    const l: usize = mtx.len / cols;
    var i: usize = 0;
    while (i < l) {
        const vec: []f32 = mtx[(i * cols)..((i * cols) + cols)];
        std.debug.print("{}: ", .{i});
        prntXvec(vec);
        i += 1;
    }
    std.debug.print("\n", .{});
}

test "XMTX: prntXmtx test" {
    var v1: [9]f32 = .{ 3, 3, 3, 0, 0, 0, 1, 1, 1 };

    const start = try Instant.now();
    prntXmtx(&v1, 3);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nprntXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("prntXmtx", elapsed1);

    prntNl();
}

///Prints the matrix, after printing a new line and followed by a new line, with the specified number of columns.
///
///  mtx = The matrix to print.
///
///  cols = The number of columns in the matrix.
///
pub fn prntXmtxNl(mtx: []f32, cols: usize) void {
    std.debug.print("\n", .{});
    const l: usize = mtx.len / cols;
    var i: usize = 0;
    while (i < l) {
        const vec: []f32 = mtx[(i * cols)..((i * cols) + cols)];
        std.debug.print("{}: ", .{i});
        prntXvec(vec);
        i += 1;
    }
    //std.debug.print("\n", .{});
}

test "XMTX: prntXmtxNl test" {
    var v1: [9]f32 = .{ 3, 3, 3, 0, 0, 0, 1, 1, 1 };

    const start = try Instant.now();
    prntXmtxNl(&v1, 3);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nprntXmtxNl: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("prntXmtxNl", elapsed1);

    prntNl();
}

///Prints the specified vector.
///
///  vec = The vector to print.
///
pub fn prntXvec(vec: []f32) void {
    const l: usize = vec.len;
    var i: usize = 0;
    var c: u8 = '-';
    while (i < l) {
        if (i == 0) {
            c = 'x';
        } else if (i == 1) {
            c = 'y';
        } else if (i == 2) {
            c = 'z';
        } else {
            c = 'w';
        }
        std.debug.print("{c}: {} ", .{ c, vec[i] });
        i += 1;
    }
    std.debug.print("\n", .{});
}

test "XMTX: prntXvec test" {
    var v1: [3]f32 = .{ 3, 3, 3 };

    const start = try Instant.now();
    prntXvecNl(&v1);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nprntXvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("prntXvec", elapsed1);

    prntNl();
}

///Prints the specified vector starting on a new line.
///
///  vec = The vector to print.
///
pub fn prntXvecNl(vec: []f32) void {
    std.debug.print("\n", .{});
    const l: usize = vec.len;
    var i: usize = 0;
    var c: u8 = '-';
    while (i < l) {
        if (i == 0) {
            c = 'x';
        } else if (i == 1) {
            c = 'y';
        } else if (i == 2) {
            c = 'z';
        } else {
            c = 'w';
        }
        std.debug.print("{c}: {} ", .{ c, vec[i] });
        i += 1;
    }
    //std.debug.print("\n", .{});
}

test "XMTX: prntXvecNl test" {
    var v1: [3]f32 = .{ 3, 3, 3 };

    const start = try Instant.now();
    prntXvecNl(&v1);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nprntXvecNl: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("prntXvecNl", elapsed1);

    prntNl();
}

///Calculates the angle between the two provided vectors.
///
///  vecL = The left-hand vector in the calculation.
///
///  vecR = The right-hand vector in the calculation.
///
///  returns = The angle in radians between the two vectors.
///
pub fn aglBtwnXvec(vecL: []f32, vecR: []f32) f32 {
    const dotP = dotPrdXvec(vecL, vecR);
    const mag1 = magXvec(vecL);
    const mag2 = magXvec(vecR);
    const cosA: f32 = (dotP / (mag1 * mag2));
    const arcCosA = std.math.acos(cosA);
    return arcCosA;
}

test "XMTX: aglBtwnXvec test" {
    var v1: [3]f32 = .{ 1, 0, 0 };
    var v2: [3]f32 = .{ 0, 0, 1 };
    var angle: f32 = -1.0;

    const start = try Instant.now();
    angle = aglBtwnXvec(&v1, &v2);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\naglBtwnXvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("aglBtwnXvec", elapsed1);

    const exp: f32 = 1.57079625e+00;
    prntXvecNl(&v1);
    std.debug.print("\nAngle RAD: {}", .{angle});
    std.debug.print("\nAngle DEG: {}", .{rad2Deg(angle)});
    try std.testing.expectEqual(exp, angle);
    prntNl();
}

///A function for converting radians to degrees.
///
///  rad = The number of radians to convert.
///
///  returns = The number of degrees.
///
pub fn rad2Deg(rad: f32) f32 {
    return (rad * (180.0 / std.math.pi));
}

test "XMTX: rad2Deg test" {
    var v1: [3]f32 = .{ 1, 0, 0 };
    var v2: [3]f32 = .{ 0, 0, 1 };

    const start = try Instant.now();
    const angleRad: f32 = aglBtwnXvec(&v1, &v2);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nrad2Deg: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rad2Deg", elapsed1);

    const exp1: f32 = 1.57079625e+00;
    const exp2: f32 = 90;
    const angleDeg: f32 = rad2Deg(angleRad);

    prntXvecNl(&v1);
    std.debug.print("\nAngle RAD: {}", .{angleRad});
    std.debug.print("\nAngle DEG: {}", .{angleDeg});
    try std.testing.expectEqual(exp1, angleRad);
    try std.testing.expectEqual(exp2, angleDeg);
    prntNl();
}

///A function for converting degrees to radians.
///
///  deg = The number of degrees to convert.
///
///  returns = The numbers of radians.
///
pub fn deg2Rad(deg: f32) f32 {
    return (deg * (std.math.pi / 180.0));
}

test "XMTX: deg2Rad test" {
    var v1: [3]f32 = .{ 1, 0, 0 };
    var v2: [3]f32 = .{ 0, 0, 1 };

    const angleRad: f32 = aglBtwnXvec(&v1, &v2);
    const exp1: f32 = 1.57079625e+00;
    const exp2: f32 = 90;
    const angleDeg: f32 = rad2Deg(angleRad);

    prntXvecNl(&v1);
    std.debug.print("\nAngle RAD: {}", .{angleRad});
    std.debug.print("\nAngle DEG: {}", .{angleDeg});
    try std.testing.expectEqual(exp1, angleRad);
    try std.testing.expectEqual(exp2, angleDeg);

    const start = try Instant.now();
    const angleRad2: f32 = deg2Rad(angleDeg);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ndeg2Rad: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("deg2Rad", elapsed1);

    std.debug.print("\nAngle RAD: {}", .{angleRad});
    std.debug.print("\nAngle RAD 2: {}", .{angleRad2});

    try std.testing.expectEqual(true, isEquF32(angleRad2, angleRad, false));
    prntNl();
}

///Determines if the two provided matrices are equal.
///
///  mtxL = The left-hand side of the matrix comparison.
///
///  mtxR = The right-hand side of the matrix comparison.
///
///  returns = A Boolean value indicating if the matrices are equal.
///
pub fn equXmtx(mtxL: []f32, mtxR: []f32) bool {
    return equXvec(mtxL, mtxR);
}

test "XMTX: equXmtx test" {
    var m1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var m2: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var b: bool = false;

    const start = try Instant.now();
    b = equXmtx(&m1, &m2);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nequXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("equXmtx", elapsed1);

    prntXmtxNl(&m1, 3);
    try std.testing.expectEqual(true, b);
    prntNl();
}

///Determines if the two provided vectors are equal.
///
///  vecL = The left-hand side of the vector comparison.
///
///  vecR = The right-hand side of the vector comparison.
///
///  returns = A Boolean value indicating if the vectors are equal.
///
pub fn equXvec(vecL: []f32, vecR: []f32) bool {
    return equXvecWrkr(vecL, vecR, COMPARE_MODE_EXACT);
}

test "XMTX: equXvec test" {
    var v1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var v2: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var b: bool = false;

    const start = try Instant.now();
    b = equXvec(&v1, &v2);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nequXvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("equXvec", elapsed1);

    prntXvecNl(&v1);
    try std.testing.expectEqual(true, b);
    prntNl();
}

///A private function that handles matrix and vector value comparison supporting two modes, direct compare and precision compare.
///
///  vecL = The left-hand side of the vector comparison.
///
///  vecR = The right-hand side of the vector comparison.
///
///  compareModeExact = A Boolean value indicating if an exact comparison mode should be used.
///
///  returns = A Boolean value indicating if the vectors are equal.
///
pub fn equXvecWrkr(vecL: []f32, vecR: []f32, compareModeExact: bool) bool {
    const l: usize = vecL.len;
    var i: usize = 0;
    while (i < l) {
        if (!compareModeExact) {
            if (absF32(vecL[i] - vecR[i]) > ZERO_F32) {
                return false;
            }
        } else {
            if (vecL[i] != vecR[i]) {
                return false;
            }
        }
        i += 1;
    }
    return true;
}

test "XMTX: equXvecWrkr test" {
    var v1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var v2: [9]f32 = .{ 1.0001, 2.0001, 3.0001, 4, 5, 6, 7, 8, 9 };
    var b: bool = false;

    var start = try Instant.now();
    b = equXvecWrkr(&v1, &v2, true);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nequXvecWrkr: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("equXvecWrkr", elapsed1);

    prntXvecNl(&v1);
    try std.testing.expectEqual(false, b);

    start = try Instant.now();
    b = equXvecWrkr(&v1, &v2, false);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nequXvecWrkr: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("equXvecWrkr", elapsed1);

    prntXvecNl(&v1);
    try std.testing.expectEqual(true, b);
    prntNl();
}

///Converts the provided matrix data into absolute values.
///
///  mtx = The matrix data to convert to absolute values.
///
pub fn absXmtx(mtx: []f32) void {
    return absXvec(mtx);
}

test "XMTX: absXmtx test" {
    var v1: [9]f32 = .{ -1, -2, -3, 4, 5, 6, -7, -8, -9 };
    var v2: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    const start = try Instant.now();
    absXmtx(&v1);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nabsXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("absXmtx", elapsed1);

    try std.testing.expectEqual(true, equXvec(&v1, &v2));
    prntNl();
}

///Converts the provided vector data into absolute values.
///
///  vec = The vector data to convert to absolute values.
///
pub fn absXvec(vec: []f32) void {
    const l: usize = vec.len;
    var i: usize = 0;
    while (i < l) {
        vec[i] = absF32(vec[i]);
        i += 1;
    }
}

test "XMTX: absXvec test" {
    var v1: [9]f32 = .{ -1, -2, -3, 4, 5, 6, -7, -8, -9 };
    var v2: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    const start = try Instant.now();
    absXvec(&v1);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nabsXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("absXmtx", elapsed1);

    try std.testing.expectEqual(true, equXvec(&v1, &v2));
    prntNl();
}

///Normalizes the provided vector based on a calculated magnitude.
///
///  vec = The vector to normalize the length of.
///
pub fn nrmXvec(vec: []f32) void {
    const magVec: f32 = magXvec(vec);
    const invMagVec: f32 = 1 / magVec;
    mulXvec(vec, invMagVec);
}

test "XMTX: nrmXvec test" {
    var v1: [3]f32 = .{ 5, 10, 15 };

    const start = try Instant.now();
    nrmXvec(&v1);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nnrmXvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("nrmXvec", elapsed1);

    prntXvecNl(&v1);
    try std.testing.expectEqual(@as(f32, 2.67261236e-01), v1[0]);
    try std.testing.expectEqual(@as(f32, 5.34522473e-01), v1[1]);
    try std.testing.expectEqual(@as(f32, 8.01783740e-01), v1[2]);
    prntNl();
}

///Projects vector P onto vector Q. Alters the vecQ argument.
///
///  vecP = The vector to project onto.
///
///  vecQ = The vector to be projected on.
///
///  returns = The parameters of the vector P that projects onto Q, stored in vecQ.
///
pub fn projXvec_VecP_Onto_VecQ(vecP: []f32, vecQ: []f32) []f32 {
    const dotProd: f32 = dotPrdXvec(vecP, vecQ);
    const mag: f32 = magXvec(vecQ);
    const magSqr: f32 = (mag * mag);
    const mul: f32 = (dotProd / magSqr);
    mulXvec(vecQ, mul);
    return vecQ;
}

test "XMTX: projXvec_VecP_Onto_VecQ test" {
    //Q = 2 2 1
    //P = 1 -2 0
    var Q: [3]f32 = .{ 2, 2, 1 };
    var P: [3]f32 = .{ 1, -2, 0 };
    var exp: [3]f32 = .{ (-4.0 / 9.0), (-4.0 / 9.0), (-2.0 / 9.0) };

    const start = try Instant.now();
    const ret: []f32 = projXvec_VecP_Onto_VecQ(&P, &Q);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nprojXvec_VecP_Onto_VecQ: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("projXvec_VecP_Onto_VecQ", elapsed1);

    clnXmtx(ret);
    clnXmtx(&exp);

    prntXvecNl(ret);
    prntNl();

    prntXvecNl(&exp);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(ret, &exp));
    prntNl();
}

///Projects vector P onto vector unit vector Q. Alters the vecQ argument.
///
///  vecP = The vector to project onto.
///
///  vecQ = The unit vector to be projected on.
///
///  returns = The parameters of the vector P that projects onto Q, stored in vecQ.
///
///  NOTE: progUontoV = progPontoQ; if v is a unit vector then <v,v> = ||v||^2 = 1, which simplifies the formula to <u,v>v.
pub fn projUnitXvec_VecP_Onto_VecQ(vecP: []f32, vecQ: []f32) []f32 {
    //Larson, Edwards: Chapter 5: remark: pg 272
    const dotProd: f32 = dotPrdXvec(vecP, vecQ);
    mulXvec(vecQ, dotProd);
    return vecQ;
}

test "XMTX: projUnitXvec_VecP_Onto_VecQ test" {
    //Q = 2 2 1
    //P = 1 -2 0
    var Q: [3]f32 = .{ 2, 2, 1 };
    nrmXvec(&Q);
    var P: [3]f32 = .{ 1, -2, 0 };
    nrmXvec(&P);
    var exp: [3]f32 = .{ -1.76676996e-02, -1.76676996e-02, -8.83384980e-03 };

    var start = try Instant.now();
    const ret1: []f32 = projXvec_VecP_Onto_VecQ(&P, &Q);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nprojXvec_VecP_Onto_VecQ: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("projXvec_VecP_Onto_VecQ", elapsed1);

    start = try Instant.now();
    const ret2: []f32 = projUnitXvec_VecP_Onto_VecQ(&P, &Q);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nprojUnitXvec_VecP_Onto_VecQ: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("projUnitXvec_VecP_Onto_VecQ", elapsed1);

    clnXvec(ret1);
    clnXvec(ret2);
    clnXvec(&exp);

    prntXvecNl(ret1);
    prntNl();
    prntXvecNl(ret2);
    prntNl();
    prntXvecNl(&exp);
    prntNl();

    try std.testing.expectEqual(true, equXvecWrkr(ret1, ret2, false));
    try std.testing.expectEqual(true, equXvecWrkr(ret1, &exp, false));
    try std.testing.expectEqual(true, equXvecWrkr(ret2, &exp, false));
    prntNl();
}

///Finds the vector perpendicular to vectors P and Q.
///
///  vecP = The left-hand vector in the calculation to find a perpendicular vector.
///
///  vecQ = The right-hand vector in the calculation to find a perpendicular vector.
///
///  returns = A vector perpendicular to vectors P and Q.
///
pub fn perpXvec_VecP_To_VecQ(vecP: []f32, vecQ: []f32) []f32 {
    diff1Xvec(vecP, projXvec_VecP_Onto_VecQ(vecP, vecQ));
    return vecP;
}

test "XMTX: perpXvec_VecP_To_VecQ test" {
    //Q = 2 2 1
    //P = 1 -2 0
    var Q: [3]f32 = .{ 2, 2, 1 };
    var P: [3]f32 = .{ 1, -2, 0 };
    var exp: [3]f32 = .{ 1.44444441e+00, -1.55555558e+00, 2.22222223e-01 };

    const start = try Instant.now();
    const ret: []f32 = perpXvec_VecP_To_VecQ(&P, &Q);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nperpXvec_VecP_To_VecQ: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("perpXvec_VecP_To_VecQ", elapsed1);

    clnXmtx(ret);
    clnXmtx(&exp);

    prntXvecNl(ret);
    prntNl();

    prntXvecNl(&exp);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(ret, &exp));
    prntNl();
}

///A function that handles multiple operations on a provided matrix. This function returns information about the matrix.
///
///  mtx = The matrix to process.
///
///  cols = The numbers of columns in the matrix.
///
///  op = The MTX_OPS operation to perform MTX_IS_LIN_INDP or MTX_IS_INVERTIBLE.
///
///  returns = A Boolean indicating the result of the information operation.
///
pub fn idnfXmtx(mtx: []f32, cols: usize, op: MTX_OPS) bool {
    const l: usize = mtx.len / cols;
    var i: usize = 0;
    var res: bool = true;
    while (i < l) {
        const vec: []f32 = mtx[(i * cols)..((i * cols) + cols)];
        if (op == MTX_OPS.MTX_IS_LIN_INDP or op == MTX_OPS.MTX_IS_INVERTIBLE) {
            if (i + 1 < l) {
                i += 1;
                const vecR: []f32 = mtx[(i * cols)..((i * cols) + cols)];
                if (!isLinIndXvec(vec, vecR)) {
                    res = false;
                    break;
                }
                i -= 1;
            } else {
                const vecR: []f32 = mtx[0..cols];
                if (!isLinIndXvec(vec, vecR)) {
                    res = false;
                    break;
                }
            }
        } else if (op == MTX_OPS.MTX_IS_ZERO) {
            if (!isZeroXvec(vec)) {
                res = false;
                break;
            }
        }
        i += 1;
    }
    return res;
}

test "XMTX: idnfXmtx tests" {
    var mtx: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    std.debug.print("\n", .{});

    var start = try Instant.now();
    var b = idnfXmtx(&mtx, 3.0, MTX_OPS.MTX_IS_ZERO);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nidnfXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("idnfXmtx", elapsed1);

    try std.testing.expectEqual(true, b);
    mtx = .{ 1, 0, 0, 0, 1, 0 };

    start = try Instant.now();
    b = idnfXmtx(&mtx, 3.0, MTX_OPS.MTX_IS_LIN_INDP);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nidnfXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("idnfXmtx", elapsed1);

    try std.testing.expectEqual(true, b);
    prntNl();
}

///A function that handles multiple operations on a provided matrix. This function returns an adjusted version of the matrix.
///
///  mtx = The matrix to process.
///
///  val = The adjustment value used in the specified process.
///
///  cols = The numbers of columns in the matrix.
///
///  op = The MTX_OPS operation to perform, MTX_MUL, MTX_DIV, MTX_ADD, MTX_SUB, MTX_PRNT, MTX_NRM, MTX_ABS.
///
pub fn procXmtx(mtx: []f32, val: f32, cols: usize, op: MTX_OPS) void {
    const l: usize = mtx.len / cols;
    var i: usize = 0;

    if (op == MTX_OPS.MTX_PRNT) {
        std.debug.print("\n", .{});
    }

    while (i < l) {
        const vec: []f32 = mtx[(i * cols)..((i * cols) + cols)];
        if (op == MTX_OPS.MTX_MUL) {
            mulXvec(vec, val);
        } else if (op == MTX_OPS.MTX_DIV) {
            divXvec(vec, val);
        } else if (op == MTX_OPS.MTX_ADD) {
            addXvec(vec, val);
        } else if (op == MTX_OPS.MTX_SUB) {
            subXvec(vec, val);
        } else if (op == MTX_OPS.MTX_PRNT) {
            std.debug.print("{}: ", .{i});
            prntXvec(vec);
        } else if (op == MTX_OPS.MTX_NRM) {
            nrmXvec(vec);
        } else if (op == MTX_OPS.MTX_ABS) {
            absXvec(vec);
        }
        i += 1;
    }

    if (op == MTX_OPS.MTX_PRNT) {
        std.debug.print("\n", .{});
    }
}

test "XMTX: procXmtx tests" {
    var mtx: [6]f32 = .{ 1, 3, 5, 10, 3, 6 };
    var exp: [6]f32 = .{ 2, 6, 10, 20, 6, 12 };
    prntNl();

    var start = try Instant.now();
    procXmtx(&mtx, 0, 3.0, MTX_OPS.MTX_PRNT);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nprocXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("procXmtx", elapsed1);

    start = try Instant.now();
    procXmtx(&mtx, 2.0, 3.0, MTX_OPS.MTX_MUL);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nprocXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("procXmtx", elapsed1);

    start = try Instant.now();
    procXmtx(&mtx, 0, 3.0, MTX_OPS.MTX_PRNT);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nprocXmtx #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("procXmtx", elapsed1);

    try std.testing.expectEqual(true, equXvec(&exp, &mtx));
    exp = .{ 1, 3, 5, 10, 3, 6 };

    start = try Instant.now();
    procXmtx(&mtx, 2.0, 3.0, MTX_OPS.MTX_DIV);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nprocXmtx #4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("procXmtx", elapsed1);

    start = try Instant.now();
    procXmtx(&mtx, 0, 3.0, MTX_OPS.MTX_PRNT);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nprocXmtx #5: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("procXmtx", elapsed1);

    try std.testing.expectEqual(true, equXvec(&exp, &mtx));
    exp = .{ 2, 4, 6, 11, 4, 7 };

    start = try Instant.now();
    procXmtx(&mtx, 1.0, 3.0, MTX_OPS.MTX_ADD);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nprocXmtx #6: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("procXmtx", elapsed1);

    start = try Instant.now();
    procXmtx(&mtx, 0, 3.0, MTX_OPS.MTX_PRNT);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nprocXmtx #7: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("procXmtx", elapsed1);

    try std.testing.expectEqual(true, equXvec(&exp, &mtx));
    exp = .{ 1, 3, 5, 10, 3, 6 };

    start = try Instant.now();
    procXmtx(&mtx, 1.0, 3.0, MTX_OPS.MTX_SUB);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nprocXmtx #8: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("procXmtx", elapsed1);

    start = try Instant.now();
    procXmtx(&mtx, 0, 3.0, MTX_OPS.MTX_PRNT);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nprocXmtx #9: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("procXmtx", elapsed1);

    try std.testing.expectEqual(true, equXvec(&exp, &mtx));
    mtx = .{ -1, -3, -5, -10, -3, -6 };
    exp = .{ 1, 3, 5, 10, 3, 6 };

    start = try Instant.now();
    procXmtx(&mtx, 2.0, 3.0, MTX_OPS.MTX_ABS);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nprocXmtx #10: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("procXmtx", elapsed1);

    start = try Instant.now();
    procXmtx(&mtx, 0, 3.0, MTX_OPS.MTX_PRNT);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nprocXmtx #11: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("procXmtx", elapsed1);

    try std.testing.expectEqual(true, equXvec(&exp, &mtx));
    exp = .{ 1.69030845e-01, 5.07092535e-01, 8.45154225e-01, 8.30454826e-01, 2.49136447e-01, 4.98272895e-01 };

    start = try Instant.now();
    procXmtx(&mtx, 2.0, 3.0, MTX_OPS.MTX_NRM);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nprocXmtx #12: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("procXmtx", elapsed1);

    start = try Instant.now();
    procXmtx(&mtx, 0, 3.0, MTX_OPS.MTX_PRNT);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nprocXmtx #13: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("procXmtx", elapsed1);

    try std.testing.expectEqual(true, equXvec(&exp, &mtx));
    prntNl();
}

///Converts the provided matrix to reduced row eschelon form using Gauss-Jordan Elimination and optionaly calculate the matrix inverse. Alters the provided matrix inline.
///
///  mtx = The matrix to reduce.
///
///  cols = The number of columns in the matrix.
///
///  hasAug = A Boolean indicating if the matrix is an augmented matrix.
///
///  hasIdtMtx = A Boolean indicating if the identity matrix has been provided.
///
///  idtMtx = The identity matrix associated with the mtx matrix provided.
///
///  dim = The number of matrix columns that must be zero for a zero row to exist.
///
///  triagRdcOnly = A Boolean value indicating if the reduction operation should stop when the matrix is triangular.
///
///  sclr = A pointer to a floating point variable that keeps track of the scalar multiplication performed against the matrix, mtx.
///
///  returns = A Boolean value indicating  if the matrix was reduced successfuly.
///
pub fn rdcXmtxInl(mtx: []f32, cols: usize, hasAug: bool, hasIdtMtx: bool, idtMtx: []f32, dim: usize, triagRdcOnly: bool, sclr: *f32) !bool {
    const rows: usize = mtx.len / cols;
    var r: usize = 0;
    var c: usize = 0;
    const errRow: usize = mtx.len;
    var diff: usize = 0;
    var idtCols: usize = cols;
    const verbose: bool = false;
    var sclMul: f32 = 1.0;

    if (hasAug) {
        diff = 1;
        if (dim < (cols - diff)) {
            diff = (cols - dim);
        }
    }

    idtCols -= diff;
    const idtRows: usize = idtMtx.len / idtCols;
    var startingCol: usize = 0;

    if (hasAug) {
        if (rows != (cols - diff)) {
            prntNlStr("!!Warning rdcXmtxInl matrix should be square or have rows + 1 columns !!");
            return Error.InvalidLengths;
        }
    } else {
        if (rows != cols) {
            prntNlStr("!!Warning rdcXmtxInl matrix should be square or have rows + 1 columns !!");
            return Error.InvalidLengths;
        }
    }

    if (hasIdtMtx) {
        if (idtRows != idtCols) {
            prntNlStr("!!Warning rdcXmtxInl identity matrix should be square !!");
            return Error.InvalidLengths;
        }
    }

    if (verbose) {
        prntNlStr("Start Mtx");
        prntXmtxNl(mtx, cols);
        prntNl();
    }

    while (r < rows) : (r += 1) {
        c = startingCol;

        var fndRow: usize = 0;
        if (mtx[(r * cols) + c] == 0) {
            fndRow = fndLgstRowAbsXmtx(mtx, cols, c, r);
            if (fndRow != r and fndRow != errRow) {
                if (verbose) {
                    prntNlStrArgs("Alternate row {} with {}", .{ fndRow, r });
                }

                altXmtxRowsInl(fndRow, r, cols, mtx);
                if (hasIdtMtx) {
                    altXmtxRowsInl(fndRow, r, idtCols, idtMtx);
                }

                if (verbose) {
                    prntXmtxNl(mtx, cols);
                    prntNl();
                }
            }
        }

        var v: f32 = 0;
        var tmp: f32 = 0;
        if (mtx[(r * cols) + c] != 1) {
            v = mtx[(r * cols) + c];
            if (v != 0.0) {
                tmp = (1.0 / v);
                sclMul *= v;

                if (verbose) {
                    prntNlStrArgs("{} Scalar multiply row: {} by {}", .{ c, r, tmp });
                }

                sclMulXmtxRowsInl(r, tmp, cols, mtx);
                if (hasIdtMtx) {
                    sclMulXmtxRowsInl(r, tmp, idtCols, idtMtx);
                }
            }

            if (verbose) {
                prntXmtxNl(mtx, cols);
                prntNl();
            }
        }

        //for other rows
        var z: usize = 0;

        if (triagRdcOnly) {
            z = r + 1;
        }

        while (z < rows) : (z += 1) {
            if (z != r) {
                const tmp2: f32 = -1.0;
                const tmp3: f32 = mtx[((z * cols) + c)] * tmp2;
                if (verbose) {
                    prntNlStrArgs("{} Add scalar multiple row: {} to row {} val {}, {}, {}", .{ c, r, z, tmp3, mtx[((z * cols) + c)], ((z * cols) + c) });
                }

                addSclMulXmtxRowsInl(r, z, tmp3, cols, mtx);
                if (hasIdtMtx) {
                    addSclMulXmtxRowsInl(r, z, tmp3, idtCols, idtMtx);
                }

                if (isZeroXvec(mtx[(z * cols)..((z * cols) + dim)])) {
                    if (hasAug == false or (hasAug == true and isEquF32(mtx[(z * cols) + (cols - 1)], 0, true) == false)) {
                        if (verbose) {
                            prntNlStrArgs("Found zero vector index start {} and stop {}", .{ (z * cols), ((z * cols) + cols - 1) });
                            prntXvecNl(mtx[(z * cols)..((z * cols) + cols - 1)]);
                            prntNl();

                            prntNlStrArgs("nMatrix state:", .{});
                            prntXmtxNl(mtx, cols);
                            prntNl();
                        }
                        return false;
                    }
                }

                if (verbose) {
                    prntXmtxNl(mtx, cols);
                    prntNl();
                }
            }
        }

        startingCol += 1;
        if (startingCol >= (cols - diff)) {
            break;
        }

        if (verbose) {
            prntNlStrArgs("Row {} Mtx", .{r});
            prntXmtxNl(mtx, cols);
            prntNl();
        }
    }

    if (verbose) {
        prntNlStrArgs("Scalar Multiple: {}", .{sclMul});
        prntXmtxNl(mtx, cols);
        prntNl();
    }

    sclr.* = sclMul;
    return true;
}

test "XMTX: rdcXmtxInl test" {
    const origM1: [9]f32 = .{ 3, 2, -3, 4, -3, 6, 1, 0, -1 };
    _ = origM1;

    var m1: [12]f32 = .{ 3, 2, -3, -13, 4, -3, 6, 7, 1, 0, -1, -5 };
    var m2: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var m3: [9]f32 = .{ 2, 3, 8, 6, 0, -3, -1, 3, 2 };
    var origM3: [9]f32 = .{ 2, 3, 8, 6, 0, -3, -1, 3, 2 };

    var m4: [9]f32 = .{ 1, -1, 0, 1, 0, -1, -6, 2, 3 };
    var origM4: [9]f32 = .{ 1, -1, 0, 1, 0, -1, -6, 2, 3 };
    var idtM1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var idtM3: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    var idtM4: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var ret: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var retVal: bool = false;

    std.debug.print("\nrdcXmtxInl test: initial matrix", .{});
    prntXmtxNl(&m1, 4);
    var sclr: f32 = 0.0;

    var start = try Instant.now();
    _ = try rdcXmtxInl(&m1, 4, true, true, &idtM1, 3, false, &sclr);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nrdcXmtxInl #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rdcXmtxInl", elapsed1);

    cpyLessXmtx(&m1, &m2, 4, 3);
    clnXmtx(&m2);
    try std.testing.expectEqual(true, isDiagXmtx(&m2, 3));
    try std.testing.expectEqual(true, isIdtXmtx(&m2, 3));

    std.debug.print("\nFull Answer:", .{});
    prntXmtxNl(&m1, 4);
    std.debug.print("\nClean Short Answer:", .{});
    prntXmtxNl(&m2, 3);
    prntNl();
    try std.testing.expectEqual(true, isRdcFrmXmtx(&m2, 3));
    sclr = 0.0;

    start = try Instant.now();
    _ = try rdcXmtxInl(&m3, 3, false, true, &idtM3, 3, false, &sclr);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrdcXmtxInl #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rdcXmtxInl", elapsed1);

    clnXmtx(&m3);
    clnXmtx(&idtM3);
    cpyXmtx(&m3, &m2);
    std.debug.print("\nInverse:", .{});
    prntXmtxNl(&idtM3, 3);
    prntNl();

    std.debug.print("\nTest:", .{});
    retVal = try tmsXmtx(&origM3, 3, &idtM3, 3, &ret, 3);
    clnXmtx(&ret);
    prntXmtxNl(&ret, 3);
    prntNl();
    sclr = 0.0;

    start = try Instant.now();
    _ = try rdcXmtxInl(&m4, 3, false, true, &idtM4, 3, false, &sclr);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nrdcXmtxInl #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rdcXmtxInl", elapsed1);

    std.debug.print("\nInverse 2:", .{});
    prntXmtxNl(&idtM4, 3);
    prntNl();

    ret = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    std.debug.print("\nTest 2:", .{});
    retVal = try tmsXmtx(&origM4, 3, &idtM4, 3, &ret, 3);
    clnXmtx(&ret);
    prntXmtxNl(&ret, 3);
    prntNl();
}

///Converts the provided matrix to reduced row eschelon form using Gauss-Jordan Elimination and optionaly calculate the matrix inverse.
///
///  mtx = The matrix to reduce.
///
///  cols = The number of columns in the matrix.
///
///  hasAug = A Boolean indicating if the matrix is an augmented matrix.
///
///  ret = The matrix that holds the reduced matrix.
///
///  hasIdtMtx = A Boolean indicating if the identity matrix has been provided.
///
///  idtMtx = The identity matrix associated with the mtx matrix provided.
///
///  dim = The number of matrix columns that must be zero for a zero row to exist.
///
///  triagRdcOnly = A Boolean value indicating if the reduction operation should stop when the matrix is triangular.
///
///  sclr = A pointer to a floating point variable that keeps track of the scalar multiplication performed against the matrix, mtx.
///
///  returns = A Boolean value indicating  if the matrix was reduced successfuly.
///
pub fn rdcXmtx(mtx: []f32, cols: usize, hasAug: bool, ret: []f32, hasIdtMtx: bool, idtMtx: []f32, dim: usize, triagRdcOnly: bool, sclr: *f32) !bool {
    cpyXmtx(mtx, ret);
    return try rdcXmtxInl(ret, cols, hasAug, hasIdtMtx, idtMtx, dim, triagRdcOnly, sclr);
}

test "XMTX: rdcXmtx test" {
    var m1: [9]f32 = .{ 3, 2, -3, 4, -3, 6, 1, 0, -1 };
    var m2: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var idtM1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var origM1: [9]f32 = .{ 3, 2, -3, 4, -3, 6, 1, 0, -1 };

    const cols: f32 = 3;
    _ = cols;
    var b: bool = false;
    std.debug.print("\nrdcXmtx test: initial matrix", .{});
    prntXmtxNl(&m1, 4);
    var sclr: f32 = 0.0;

    const start = try Instant.now();
    b = try rdcXmtxInl(&m1, 3, false, true, &idtM1, 3, false, &sclr);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nrdcXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rdcXmtx", elapsed1);

    try std.testing.expectEqual(true, b);
    cpyXmtx(&m1, &m2);
    cpyXmtx(&origM1, &m1);

    clnXmtx(&m2);
    try std.testing.expectEqual(true, isDiagXmtx(&m2, 3));
    try std.testing.expectEqual(true, isIdtXmtx(&m2, 3));

    std.debug.print("\nClean Short Answer:", .{});
    prntXmtxNl(&m2, 3);
    prntNl();
    try std.testing.expectEqual(true, isRdcFrmXmtx(&m2, 3));
    prntNl();
}

///Prints a new line followed by the given string and an ending new line.
pub fn prntNlStr(comptime s: []const u8) void {
    std.debug.print("\n {s}", .{s});
}

test "XMTX: prntNlStr test" {
    const start = try Instant.now();
    prntNlStr("prntNlStr");
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nprntNlStr: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("prntNlStr", elapsed1);
    prntNl();
}

///Prints a new line followed by the given string formatted with the given args.
pub fn prntNlStrArgs(comptime s: []const u8, args: anytype) void {
    const test_allocator = std.testing.allocator;
    const str = std.fmt.allocPrint(test_allocator, s, args) catch |err| {
        std.debug.print("prntNlStrArgs: Exception caught: {}", .{err});
        return;
    };
    std.debug.print("\n {s}", .{str});
    test_allocator.free(str);
}

test "XMTX: prntNlStrArgs test" {
    const start = try Instant.now();
    prntNlStrArgs("prntNlStrArgs: {}", .{0});
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nprntNlStrArgs: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("prntNlStrArgs", elapsed1);
    prntNl();
}

///Prints a newline.
pub fn prntNl() void {
    std.debug.print("\n", .{});
}

test "XMTX: prntNl test" {
    const start = try Instant.now();
    prntNl();
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nprntNl: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("prntNl", elapsed1);
    prntNl();
}

///Finds the largest row, by absolute value, in the given matrix, at the specified column, starting on the given row.
///
///  mtx = The matrix to search through.
///
///  cols = The number of columns in the matrix.
///
///  targetCol = The target column to look for values.
///
///  startingRow = The row to start the search on.
///
///  returns = The index of the row found or an illegal row index, the number of rows.
///
pub fn fndLgstRowAbsXmtx(mtx: []f32, cols: usize, targetCol: usize, startingRow: usize) usize {
    const l: usize = mtx.len;
    var i: usize = (cols * startingRow);
    var min: f32 = 0;
    var targetRow: usize = l;
    var found: bool = false;
    while (i < l) : (i += 1) {
        const row = i / cols;
        const col = i % cols;
        const val: f32 = mtx[i];
        const absVal: f32 = absF32(val);
        if (row >= startingRow and col == targetCol and absVal > min and absVal != 0) {
            min = absVal;
            targetRow = row;
            found = true;
        }
    }
    return targetRow;
}

test "XMTX: fndLgstRowAbsXmtx test" {
    var mtx: [12]f32 = .{ 1, -3, 21, -10, 3, 6, 7, -13, 5, 0, 0, 0 };
    const exp: [9]usize = .{ 1, 1, 2, 2, 0, 1, 12, 12, 12 };

    var start = try Instant.now();
    const v0: usize = fndLgstRowAbsXmtx(&mtx, 3, 0, 0);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nfndLgstRowAbsXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("fndLgstRowAbsXmtx", elapsed1);

    if (v0 >= mtx.len) {
        std.debug.print("\nNo matching row found, v0", .{});
    }

    start = try Instant.now();
    const v1: usize = fndLgstRowAbsXmtx(&mtx, 3, 0, 1);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nfndLgstRowAbsXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("fndLgstRowAbsXmtx", elapsed1);

    if (v1 >= mtx.len) {
        std.debug.print("\nNo matching row found, v1", .{});
    }

    start = try Instant.now();
    const v2: usize = fndLgstRowAbsXmtx(&mtx, 3, 1, 0);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nfndLgstRowAbsXmtx #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("fndLgstRowAbsXmtx", elapsed1);

    if (v2 >= mtx.len) {
        std.debug.print("\nNo matching row found, v2", .{});
    }

    start = try Instant.now();
    const v3: usize = fndLgstRowAbsXmtx(&mtx, 3, 1, 1);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nfndLgstRowAbsXmtx #4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("fndLgstRowAbsXmtx", elapsed1);

    if (v3 >= mtx.len) {
        std.debug.print("\nNo matching row found, v3", .{});
    }

    start = try Instant.now();
    const v4: usize = fndLgstRowAbsXmtx(&mtx, 3, 2, 0);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nfndLgstRowAbsXmtx #5: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("fndLgstRowAbsXmtx", elapsed1);

    if (v4 >= mtx.len) {
        std.debug.print("\nNo matching row found, v4", .{});
    }

    start = try Instant.now();
    const v5: usize = fndLgstRowAbsXmtx(&mtx, 3, 2, 1);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nfndLgstRowAbsXmtx #6: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("fndLgstRowAbsXmtx", elapsed1);

    if (v5 >= mtx.len) {
        std.debug.print("\nNo matching row found, v5", .{});
    }

    start = try Instant.now();
    const v6: usize = fndLgstRowAbsXmtx(&mtx, 3, 3, 0);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nfndLgstRowAbsXmtx #7: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("fndLgstRowAbsXmtx", elapsed1);

    if (v6 >= mtx.len) {
        std.debug.print("\nNo matching row found, v6", .{});
    }

    start = try Instant.now();
    const v7: usize = fndLgstRowAbsXmtx(&mtx, 3, 3, 1);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nfndLgstRowAbsXmtx #8: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("fndLgstRowAbsXmtx", elapsed1);

    if (v7 >= mtx.len) {
        std.debug.print("\nNo matching row found, v7", .{});
    }

    start = try Instant.now();
    const v8: usize = fndLgstRowAbsXmtx(&mtx, 3, 3, 2);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nfndLgstRowAbsXmtx #9: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("fndLgstRowAbsXmtx", elapsed1);

    if (v8 >= mtx.len) {
        std.debug.print("\nNo matching row found, v8", .{});
    }

    std.debug.print("\nv0: {} v1: {} v2: {} v3: {} v4: {} v5: {} v6: {} v7: {} v8: {}", .{ v0, v1, v2, v3, v4, v5, v6, v7, v8 });
    try std.testing.expectEqual(exp[0], v0);
    try std.testing.expectEqual(exp[1], v1);
    try std.testing.expectEqual(exp[2], v2);
    try std.testing.expectEqual(exp[3], v3);
    try std.testing.expectEqual(exp[4], v4);
    try std.testing.expectEqual(exp[5], v5);
    try std.testing.expectEqual(exp[6], v6);
    try std.testing.expectEqual(exp[7], v7);
    try std.testing.expectEqual(exp[8], v8);

    prntNl();
}

///Returns a Boolean indicating if the provided matrix is reduced.
///
///  mtx = The matrix to process.
///
///  cols = The number of columns in the matrix.
///
///  returns = A Boolean value indicating if the given matrix, mtx, is in reduced form.
///
pub fn isRdcFrmXmtx(mtx: []f32, cols: usize) bool {
    const l: usize = mtx.len;
    var i: usize = 0;
    var foundZero: bool = false;

    while (i < l) : (i += 1) {
        const row = i / cols;
        const col = i % cols;

        if (row > col) {
            if (mtx[i] != 0) {
                return false;
            }
        } else if (row == col) {
            if (mtx[i] != 1) {
                if (mtx[i] == 0 and !foundZero) {
                    foundZero = true;
                } else {
                    return false;
                }
            } else {
                if (foundZero) {
                    return false;
                }
            }
        }
    }
    return true;
}

test "XMTX: isRdcFrmXmtx test" {
    var mtx: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var b: bool = false;
    prntXmtxNl(&mtx, 3);

    var start = try Instant.now();
    b = isRdcFrmXmtx(&mtx, 3);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisRdcFrmXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isRdcFrmXmtx", elapsed1);

    try std.testing.expectEqual(true, b);
    mtx = .{ 1, 0, 0, 0, 1, 0, 0, 0, 0 };
    prntXmtxNl(&mtx, 3);

    start = try Instant.now();
    b = isRdcFrmXmtx(&mtx, 3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisRdcFrmXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isRdcFrmXmtx", elapsed1);

    try std.testing.expectEqual(true, b);
    mtx = .{ 1, 0, 0, 0, 0, 0, 0, 0, 1 };
    prntXmtxNl(&mtx, 3);

    start = try Instant.now();
    b = isRdcFrmXmtx(&mtx, 3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisRdcFrmXmtx #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isRdcFrmXmtx", elapsed1);

    try std.testing.expectEqual(false, b);
    prntNl();
}

///Adds a scalar multiple of one matrix row to another.
///Processes the operation inline and stores the results in the provided matrix.
///
///  srcRow = The source row, multiplied by mul, to be added to the destination row.
///
///  dstRow = The destination row to add the scalar multiplied source row to.
///
///  mul = The multiplier to apply to the source row.
///
///  cols = The number of columns in the matrix.
///
///  mtx = The matrix used in the row operation.
///
pub fn addSclMulXmtxRowsInl(srcRow: usize, dstRow: usize, mul: f32, cols: usize, mtx: []f32) void {
    const l: usize = mtx.len;
    var i: usize = (dstRow * cols);
    while (i < l) {
        if (i >= (dstRow * cols) and i < ((dstRow * cols) + cols)) {
            const row: usize = i / cols;
            _ = row;
            const col: usize = i % cols;
            //std.debug.print("\n{} {} {} {} Src row val {} multiplier {} dst row val {} final dst row val {}", .{i, col, cols, row, mtx[((srcRow * cols) + col)], mul, mtx[i], (mtx[i] + (mtx[((srcRow * cols) + col)] * mul)) });
            mtx[i] += (mtx[((srcRow * cols) + col)] * mul);
            if (i == ((dstRow * cols) + cols - 1)) {
                break;
            }
        }
        i += 1;
    }
}

test "XMTX: addSclMulXmtxRowsInl test" {
    var mtx: [9]f32 = .{ 1, 1, 1, 2, 2, 2, 3, 3, 3 };
    var exp: [9]f32 = .{ 1, 1, 1, 2, 2, 2, 7, 7, 7 };
    prntXmtxNl(&mtx, 3);

    const start = try Instant.now();
    addSclMulXmtxRowsInl(1, 2, 2, 3, &mtx);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\naddSclMulXmtxRowsInl: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("addSclMulXmtxRowsInl", elapsed1);

    prntXmtxNl(&mtx, 3);
    try std.testing.expectEqual(true, equXvec(&exp, &mtx));
    prntNl();
}

///Adds a scalar multiple of one matrix row to another.
///
///  srcRow = The source row, multiplied by mul, to be added to the destination row.
///
///  dstRow = The destination row to add the scalar multiplied source row to.
///
///  mul = The multiplier to apply to the source row.
///
///  cols = The number of columns in the matrix.
///
///  mtx = The matrix used in the row operation.
///
///  ret = The return matrix used in the row operation.
///
pub fn addSclMulXmtxRows(srcRow: usize, dstRow: usize, mul: f32, cols: usize, mtx: []f32, ret: []f32) void {
    cpyXmtx(mtx, ret);
    addSclMulXmtxRowsInl(srcRow, dstRow, mul, cols, ret);
}

test "XMTX: addSclMulXmtxRows test" {
    var mtx: [9]f32 = .{ 1, 1, 1, 2, 2, 2, 3, 3, 3 };
    var res: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var exp: [9]f32 = .{ 1, 1, 1, 2, 2, 2, 7, 7, 7 };
    prntXmtxNl(&mtx, 3);

    const start = try Instant.now();
    addSclMulXmtxRows(1, 2, 2, 3, &mtx, &res);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\naddSclMulXmtxRows: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("addSclMulXmtxRows", elapsed1);

    prntXmtxNl(&res, 3);
    try std.testing.expectEqual(true, equXvec(&exp, &res));
    prntNl();
}

///Multiplies a matrix row by a scalar. Performs the operation inline in the provided matrix.
///
///  srcRow = The source row to apply the multiplication to.
///
///  mul = The value to apply in the multiplication.
///
///  cols = The number of columns in the matrix.
///
///  mtx = The matrix to use in this operation.
///
pub fn sclMulXmtxRowsInl(srcRow: usize, mul: f32, cols: usize, mtx: []f32) void {
    const l: usize = mtx.len;
    var i: usize = (srcRow * cols);
    while (i < l) {
        if (i >= (srcRow * cols) and i < ((srcRow * cols) + cols)) {
            //std.debug.print("\n{} times {} = {}", .{mtx[i], mul, (mtx[i] * mul)});
            mtx[i] *= mul;
            if (i == ((srcRow * cols) + cols - 1)) {
                break;
            }
        }
        i += 1;
    }
}

test "XMTX: sclMulXmtxRowsInl test" {
    var mtx: [9]f32 = .{ 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    var exp: [9]f32 = .{ 3, 3, 3, 1, 1, 1, 1, 1, 1 };
    prntXmtxNl(&mtx, 3);

    const start = try Instant.now();
    sclMulXmtxRowsInl(0, 3, 3, &mtx);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nsclMulXmtxRowsInl: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("sclMulXmtxRowsInl", elapsed1);

    prntXmtxNl(&mtx, 3);
    try std.testing.expectEqual(true, equXvec(&exp, &mtx));
    prntNl();
}

///Multiplies a matrix row by a scalar.
///
///  srcRow = The source row to apply the multiplication to.
///
///  mul = The value to apply in the multiplication.
///
///  cols = The number of columns in the matrix.
///
///  mtx = The matrix to use in this operation.
///
///  ret = The return matrix that holds the newly calculated variables.
///
pub fn sclMulXmtxRows(srcRow: usize, mul: f32, cols: usize, mtx: []f32, ret: []f32) void {
    cpyXmtx(mtx, ret);
    sclMulXmtxRowsInl(srcRow, mul, cols, ret);
}

test "XMTX: sclMulXmtxRows test" {
    var mtx: [9]f32 = .{ 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    var res: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var exp: [9]f32 = .{ 3, 3, 3, 1, 1, 1, 1, 1, 1 };
    prntXmtxNl(&mtx, 3);

    const start = try Instant.now();
    sclMulXmtxRows(0, 3, 3, &mtx, &res);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nsclMulXmtxRows: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("sclMulXmtxRows", elapsed1);

    prntXmtxNl(&res, 3);
    try std.testing.expectEqual(true, equXvec(&exp, &res));
    prntNl();
}

///Adds a scalar to a matrix row. Performs the operation inline in the provided matrix.
///
///  srcRow = The source row to apply the addition to.
///
///  amt = The value to apply in the addition.
///
///  cols = The number of columns in the matrix.
///
///  mtx = The matrix to use in this operation.
///
pub fn addSclXmtxRowsInl(srcRow: usize, amt: f32, cols: usize, mtx: []f32) void {
    const l: usize = mtx.len;
    var i: usize = (srcRow * cols);
    while (i < l) {
        if (i >= (srcRow * cols) and i < ((srcRow * cols) + cols)) {
            //std.debug.print("\n{} times {} = {}", .{mtx[i], mul, (mtx[i] * mul)});
            mtx[i] += amt;
            if (i == ((srcRow * cols) + cols - 1)) {
                break;
            }
        }
        i += 1;
    }
}

test "XMTX: addSclXmtxRowsInl test" {
    var mtx: [9]f32 = .{ 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    var exp: [9]f32 = .{ 4, 4, 4, 1, 1, 1, 1, 1, 1 };
    prntXmtxNl(&mtx, 3);

    const start = try Instant.now();
    addSclXmtxRowsInl(0, 3, 3, &mtx);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\naddSclXmtxRowsInl: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("addSclXmtxRowsInl", elapsed1);

    prntXmtxNl(&mtx, 3);
    try std.testing.expectEqual(true, equXvec(&exp, &mtx));
    prntNl();
}

///Adds a scalar to a matrix row.
///
///  srcRow = The source row to apply the addition to.
///
///  amt = The value to apply in the addition.
///
///  cols = The number of columns in the matrix.
///
///  mtx = The matrix to use in this operation.
///
///  ret = The return matrix that holds the newly calculated variables.
///
pub fn addSclXmtxRows(srcRow: usize, amt: f32, cols: usize, mtx: []f32, ret: []f32) void {
    cpyXmtx(mtx, ret);
    addSclXmtxRowsInl(srcRow, amt, cols, ret);
}

test "XMTX: addSclXmtxRows test" {
    var mtx: [9]f32 = .{ 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    var res: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var exp: [9]f32 = .{ 4, 4, 4, 1, 1, 1, 1, 1, 1 };
    prntXmtxNl(&mtx, 3);

    const start = try Instant.now();
    addSclXmtxRows(0, 3, 3, &mtx, &res);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\naddSclXmtxRows: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("addSclXmtxRows", elapsed1);

    prntXmtxNl(&res, 3);
    try std.testing.expectEqual(true, equXvec(&exp, &res));
    prntNl();
}

///Alternate two rows of a the given matrix, mtx. Performs the operation inline.
///
///  srcRow = The source row to alternate with the destination row.
///
///  dstRow = The destination row to alternate with the source row.
///
///  cols = The number of columns in the matrix.
///
///  mtx = The matrix to perform the operation on.
///
pub fn altXmtxRowsInl(srcRow: usize, dstRow: usize, cols: usize, mtx: []f32) void {
    const l: usize = mtx.len;
    var i: usize = 0;
    while (i < l) {
        if (i >= (srcRow * cols) and i < ((srcRow * cols) + cols)) {
            const a = mtx[i];
            const b = i % cols;
            const c = (dstRow * cols) + b;
            const d = mtx[c];
            mtx[c] = a;
            mtx[i] = d;
            if (i == ((srcRow * cols) + cols - 1)) {
                break;
            }
        }
        i += 1;
    }
}

test "XMTX: altXmtxRowsInl test" {
    var mtx: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var exp: [9]f32 = .{ 7, 8, 9, 4, 5, 6, 1, 2, 3 };
    prntXmtxNl(&mtx, 3);

    const start = try Instant.now();
    altXmtxRowsInl(0, 2, 3, &mtx);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\naltXmtxRowsInl: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("altXmtxRowsInl", elapsed1);

    prntXmtxNl(&mtx, 3);
    try std.testing.expectEqual(true, equXvec(&exp, &mtx));
    prntNl();
}

///Alternate two rows of a matrix.
///
///  srcRow = The source row to alternate with the destination row.
///
///  dstRow = The destination row to alternate with the source row.
///
///  cols = The number of columns in the matrix.
///
///  mtx = The matrix to perform the operation on.
///
///  ret = The return matrix that holds the newly calculated variables.
///
pub fn altXmtxRows(srcRow: usize, dstRow: usize, cols: usize, mtx: []f32, ret: []f32) void {
    cpyXmtx(mtx, ret);
    altXmtxRowsInl(srcRow, dstRow, cols, ret);
}

test "XMTX: altXmtxRows test" {
    var mtx: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var exp: [9]f32 = .{ 7, 8, 9, 4, 5, 6, 1, 2, 3 };
    var alt: [9]f32 = .{ 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    prntXmtxNl(&mtx, 3);

    const start = try Instant.now();
    altXmtxRows(0, 2, 3, &mtx, &alt);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\naltXmtxRows: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("altXmtxRows", elapsed1);

    prntXmtxNl(&alt, 3);
    try std.testing.expectEqual(true, equXvec(&exp, &alt));
    prntNl();
}

///Transposes the provided matrix into the return matrix. Operation is performed inline. Doesn't work with rectangular matrices, square matrices only.
///
///  mtx = The source matrix to rn the transpose operation on.
///
///  cols = The number of columns in the matrix.
///
///  alloc = A allocator used to create the return matrix.
///
pub fn trnXmtxInl(mtx: []f32, cols: usize, alloc: *const std.mem.Allocator) !void {
    const ret: []f32 = try crtXmtxEz(mtx.len, alloc);
    trnXmtxRect(mtx, cols, ret, cols);
    cpyXmtx(ret, mtx);
    alloc.free(ret);
}

test "XMTX: trnXmtxInl test" {
    var mtx: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var exp: [9]f32 = .{ 1, 4, 7, 2, 5, 8, 3, 6, 9 };
    prntXmtxNl(&mtx, 3);
    const alloc: std.mem.Allocator = std.testing.allocator;

    const start = try Instant.now();
    try trnXmtxInl(&mtx, 3, &alloc);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ntrnXmtxInl: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("trnXmtxInl", elapsed1);

    prntXmtxNl(&mtx, 3);
    try std.testing.expectEqual(true, equXvec(&exp, &mtx));
    prntNl();
}

///Transposes the provided matrix into the return matrix.
///
///  mtx = The source matrix to rn the transpose operation on.
///
///  cols = The number of columns in the matrix.
///
///  ret = The return matrix with the new transposed values.
///
pub fn trnXmtx(mtx: []f32, cols: usize, ret: []f32) void {
    trnXmtxRect(mtx, cols, ret, cols);
}

test "XMTX: trnXmtx test" {
    var mtx: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var exp: [9]f32 = .{ 1, 4, 7, 2, 5, 8, 3, 6, 9 };
    var trn: [9]f32 = .{ 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    prntXmtxNl(&mtx, 3);

    const start = try Instant.now();
    trnXmtx(&mtx, 3, &trn);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ntrnXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("trnXmtx", elapsed1);

    prntXmtxNl(&trn, 3);
    try std.testing.expectEqual(true, equXvec(&exp, &trn));
    prntNl();
}

///Transposes the provided matrix into the return matrix. Supports rectangular transposition using the return column argument.
///
///  mtx = The source matrix to rn the transpose operation on.
///
///  cols = The number of columns in the matrix.
///
///  ret = The return matrix with the new transposed values.
///
///  retCols = The number of columns in the return matrix.
///
pub fn trnXmtxRect(mtx: []f32, cols: usize, ret: []f32, retCols: usize) void {
    const l = mtx.len;
    const rows: usize = (l / cols);
    const lcols: usize = cols;
    var r: usize = 0;
    var c: usize = 0;

    while (r < rows) : (r += 1) {
        c = 0;
        while (c < lcols) : (c += 1) {
            ret[(c * retCols) + r] = mtx[(r * cols) + c];
        }
    }
}

test "XMTX: trnXmtxRect test" {
    var a: [9]f32 = .{ 2, 1, -2, -1, 0, 3, 0, -2, 1 };
    var b: [6]f32 = .{ 3, 1, 2, -1, 3, 0 };
    var aT: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var expAt: [9]f32 = .{ 2, -1, 0, 1, 0, -2, -2, 3, 1 };
    var bT: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    var expBt: [6]f32 = .{ 3, 2, 3, 1, -1, 0 };

    prntXmtxNl(&a, 3);
    prntNl();

    trnXmtx(&a, 3, &aT);
    prntXmtxNl(&aT, 3);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&aT, &expAt));

    prntXmtxNl(&b, 2);
    prntNl();

    const start = try Instant.now();
    trnXmtxRect(&b, 2, &bT, 3);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ntrnXmtxRect: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("trnXmtxRect", elapsed1);

    prntXmtxNl(&bT, 3);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&bT, &expBt));
    prntNl();
}

///Transposes the provided matrix into the mtx argument. Supports rectangular transposition using the return column argument.
///
///  mtx = The source matrix to run the transpose operation on and store the results to.
///
///  cols = The number of columns in the matrix.
///
///  retCols = The number of columns in the return matrix.
///
pub fn trnXmtxRectInl(mtx: []f32, cols: usize, retCols: usize, alloc: *const std.mem.Allocator) !void {
    const l = mtx.len;
    const rows: usize = (l / cols);
    const lcols: usize = cols;
    var r: usize = 0;
    var c: usize = 0;
    var ret: []f32 = try alloc.*.alloc(f32, mtx.len);
    defer alloc.*.free(ret);

    while (r < rows) : (r += 1) {
        c = 0;
        while (c < lcols) : (c += 1) {
            ret[(c * retCols) + r] = mtx[(r * cols) + c];
        }
    }
    cpyXmtx(ret, mtx);
}

test "XMTX: trnXmtxRectInl test" {
    var a: [9]f32 = .{ 2, 1, -2, -1, 0, 3, 0, -2, 1 };
    var b: [6]f32 = .{ 3, 1, 2, -1, 3, 0 };
    var aT: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var expAt: [9]f32 = .{ 2, -1, 0, 1, 0, -2, -2, 3, 1 };
    var bT: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    var bT2: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    var expBt: [6]f32 = .{ 3, 2, 3, 1, -1, 0 };

    prntXmtxNl(&a, 3);
    prntNl();

    trnXmtx(&a, 3, &aT);
    prntXmtxNl(&aT, 3);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&aT, &expAt));

    prntXmtxNl(&b, 2);
    prntNl();

    var start = try Instant.now();
    trnXmtxRect(&b, 2, &bT, 3);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ntrnXmtxRect: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("trnXmtxRect", elapsed1);

    prntXmtxNl(&bT, 3);
    try std.testing.expectEqual(true, equXmtx(&bT, &expBt));
    prntNl();

    cpyXmtx(&b, &bT2);

    start = try Instant.now();
    try trnXmtxRectInl(&bT2, 2, 3, &std.testing.allocator);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ntrnXmtxRectInl: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("trnXmtxRectInl", elapsed1);

    prntXmtxNl(&bT2, 3);
    try std.testing.expectEqual(true, equXmtx(&bT, &bT2));
    prntNl();
}

///Indicates if the provided matrix is square.
///
///  mtx = The matrix to scan.
///
///  cols = The number of columns in the matrix.
///
///  returns = A Boolean value indicating if the matrix is a square matrix.
///
pub fn isSqrXmtx(mtx: []f32, cols: usize) bool {
    const l: usize = mtx.len;
    const rows = l / cols;
    return (rows == cols);
}

test "XMTX: isSqrXmtx test" {
    var mtx1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var mtx2: [12]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1 };
    prntXmtxNl(&mtx1, 3);
    prntXmtxNl(&mtx2, 3);

    var start = try Instant.now();
    const b1: bool = isSqrXmtx(&mtx1, 3);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisSqrXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isSqrXmtx", elapsed1);

    start = try Instant.now();
    const b2: bool = isSqrXmtx(&mtx2, 3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisSqrXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isSqrXmtx", elapsed1);

    try std.testing.expectEqual(true, b1);
    try std.testing.expectEqual(false, b2);
    prntNl();
}

///Indicates if the provided matrix is diagonal.
///
///  mtx = The matrix to scan.
///
///  cols = The number of columns in the matrix.
///
///  returns = A Boolean value indicating if the matrix is a diagonal matrix.
///
pub fn isDiagXmtx(mtx: []f32, cols: usize) bool {
    const l: usize = mtx.len;
    var i: usize = 0;
    while (i < l) : (i += 1) {
        const r = i % cols;
        const c = i / cols;
        if (r != c) {
            if (mtx[i] != 0) {
                return false;
            }
        }
    }
    return true;
}

test "XMTX: isDiagXmtx test" {
    var mtx1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var mtx2: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    prntXmtxNl(&mtx1, 3);
    prntXmtxNl(&mtx2, 3);

    var start = try Instant.now();
    const b1: bool = isDiagXmtx(&mtx1, 3);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisDiagXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isDiagXmtx", elapsed1);

    start = try Instant.now();
    const b2: bool = isDiagXmtx(&mtx2, 3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisDiagXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isDiagXmtx", elapsed1);

    try std.testing.expectEqual(false, b1);
    try std.testing.expectEqual(true, b2);
    prntNl();
}

///Indicates if the provided matrix is an identity matrix.
///
///  mtx = The matrix to scan.
///
///  cols = The number of columns in the matrix.
///
///  returns = A Boolean value indicating if the matrix is an identity matrix.
///
pub fn isIdtXmtx(mtx: []f32, cols: usize) bool {
    const l: usize = mtx.len;
    var i: usize = 0;
    while (i < l) : (i += 1) {
        const r = i % cols;
        const c = i / cols;
        if (r == c) {
            if (mtx[i] != 1) {
                return false;
            }
        } else {
            if (mtx[i] != 0) {
                return false;
            }
        }
    }
    return true;
}

test "XMTX: isIdtXmtx test" {
    var mtx1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var mtx2: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    prntXmtxNl(&mtx1, 3);
    prntXmtxNl(&mtx2, 3);

    var start = try Instant.now();
    const b1: bool = isIdtXmtx(&mtx1, 3);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisIdtXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isIdtXmtx", elapsed1);

    start = try Instant.now();
    const b2: bool = isIdtXmtx(&mtx2, 3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisIdtXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isIdtXmtx", elapsed1);

    try std.testing.expectEqual(false, b1);
    try std.testing.expectEqual(true, b2);
    prntNl();
}

///Creates a smaller matrix for cofactor calculations based on the provided, larger, initial matrix mtx.
///Returns a Boolean value indicating if the operation to create the cofactor matrix and store it in the ret matrix argument.
///
///  mtx = The matrix to find a cofactor matrix for.
///
///  cols = The number of cols in the matrix.
///
///  cofR = The row to find the cofactor matrix for.
///
///  cofC = The col to find the cofactor matrix for.
///
///  ret = The matrix to store the resulting cofactor matrix.
///
///  colsRet = The the number of columns in the return matrix.
///
///  returns = A Boolean value indicating if the operation was a success.
///
pub fn cofXmtx(mtx: []f32, cols: usize, cofR: usize, cofC: usize, ret: []f32, colsRet: usize) bool {
    const rows: usize = mtx.len / cols;
    const rowsRet: usize = ret.len / colsRet;

    if (rowsRet != (rows - 1)) {
        std.debug.print("\n!!Warning cofXmtx matrix should have 1 less than the matrix rows !!", .{});
        return false;
    }

    if (colsRet != (cols - 1)) {
        std.debug.print("\n!!Warning cofXmtx matrix should have 1 less than the matrix columns !!", .{});
        return false;
    }

    var r: usize = 0;
    var c: usize = 0;
    var tr: usize = 0;
    var tc: usize = 0;
    while (r < rows) : (r += 1) {
        c = 0;
        tc = 0;
        if (r != cofR) {
            while (c < cols) : (c += 1) {
                if (c != cofC) {
                    ret[((tr * colsRet) + tc)] = mtx[((r * cols) + c)];
                    tc += 1;
                }
            }
            tr += 1;
        }
    }
    return true;
}

test "XMTX: cofXmtx test" {
    //A = 1 2 3
    //    4 5 6
    //    7 2 9

    //detA = (1)(33)          - (2)(-6)          + (3)(-27)
    //detA = (1)(det 5 6 2 9) - (2)(det 4 6 7 9) + (3)(det 4 5 7 2)
    //detA = -36;

    //cofA = 33   6    -27
    //       -12  -12  12
    //       -3   6    -3

    //A^-1 = -11/12  1/3  1/12
    //       -1/6    1/3  -1/6
    //       3/4     -1/3 1/12

    var A: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 2, 9 };
    const cofA: [9]f32 = .{ 33, 6, -27, -12, -12, 12, -3, 6, -3 };
    var res: [4]f32 = .{ 0, 0, 0, 0 };
    const resCols: f32 = 2;
    const cols: f32 = 3;
    var b: bool = false;

    var start = try Instant.now();
    b = cofXmtx(&A, cols, 0, 0, &res, resCols);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncofXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cofXmtx", elapsed1);

    prntXmtxNl(&res, resCols);
    prntNl();
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(cofA[0], (detXmtx2(&res) * cofXmtxSign(0, 0, true)));

    start = try Instant.now();
    b = cofXmtx(&A, cols, 0, 1, &res, resCols);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ncofXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cofXmtx", elapsed1);

    prntXmtxNl(&res, resCols);
    prntNl();
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(cofA[1], (detXmtx2(&res) * cofXmtxSign(0, 1, true)));
    prntNl();
}

///Returns a 1 or -1 depending on the cofactor's row and column. An offset if automatically added for zero based row and column indices.
///
///  cofR = The row of the desired cofactor.
///
///  cofC = The column of the desired cofactor.
///
///  zeroBased = A Boolean value indicating if the row and column values provided are zero based indices or 1 based.
///
///  returns = A 1 or -1 which is the calculated sign of the cofactor.
///
pub fn cofXmtxSign(cofR: usize, cofC: usize, zeroBased: bool) f32 {
    var na: f32 = 1;
    if (zeroBased) {
        na = @floatFromInt((cofR + 1 + cofC + 1)); //@as(f32, @floatFromInt((cofR + 1 + cofC + 1)));
    } else {
        na = @floatFromInt((cofR + cofC)); //@as(f32, @floatFromInt((cofR + cofC)));
    }
    return std.math.pow(f32, -1.0, na);
}

test "XMTX: cofXmtxSign test" {
    var c: usize = 0;
    var r: usize = 0;
    var ret: f32 = -1.0;

    var start = try Instant.now();
    ret = cofXmtxSign(r, c, true);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncofXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cofXmtx", elapsed1);

    try std.testing.expectEqual(@as(f32, 1.0), ret);
    c = 1;
    r = 0;

    start = try Instant.now();
    ret = cofXmtxSign(r, c, true);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ncofXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cofXmtx", elapsed1);

    try std.testing.expectEqual(@as(f32, -1.0), ret);
    prntNl();
}

///Returns a new matrix based on the 4x4 matrix, mtx, with the specified row and column data removed.
///
///  mtx = The 4x4 matrix to use as the basis of the new matrix.
///
///  row = The row to skip in the source matrix.
///
///  col = The col to skip in the source matrix.
///
///  returns = A new 3x3 matrix based on, mtx, with the specified row and column removed.
///
pub fn rmvRowColXmtx4(mtx: *[16]f32, row: f32, col: f32) [9]f32 {
    const cols: usize = 4;
    const destCols: usize = 3;
    const r: usize = 4;
    const c: usize = 4;
    var i: usize = 0;
    var j: usize = 0;
    var actI: usize = 0;
    var actJ: usize = 0;
    var dest: [9]f32 = std.mem.zeroes([9]f32); //.{};
    const actRow: usize = @intFromFloat(row); //@as(usize, @intFromFloat(row));
    const actCol: usize = @intFromFloat(col); //@as(usize, @intFromFloat(col));

    while (i < r) {
        j = 0;
        actJ = 0;
        if (i != actRow) {
            while (j < c) {
                if (j != actCol) {
                    dest[(actI * destCols) + actJ] = mtx[(i * cols) + j];
                    actJ += 1;
                }
                j += 1;
            }
            actI += 1;
        }
        i += 1;
    }
    return dest;
}

test "XMTX: rmvRowColXmtx4 test" {
    var A: [16]f32 = .{ 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 };
    //A = 1  1  1  1
    //    2  2  2  2
    //    3  3  3  3
    //    4  4  4  4
    var expA: [9]f32 = .{ 2, 2, 2, 3, 3, 3, 4, 4, 4 };
    const r: f32 = 0.0;
    const c: f32 = 0.0;
    var res: [9]f32 = std.mem.zeroes([9]f32); //.{};

    const start = try Instant.now();
    res = rmvRowColXmtx4(&A, r, c);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nrmvRowColXmtx4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rmvRowColXmtx4", elapsed1);

    std.debug.print("\nrmvRowColXmtx4:", .{});
    prntXmtxNl(&res, 2);
    try std.testing.expectEqual(true, equXmtx(&expA, &res));
    prntNl();
}

///Returns a new matrix based on the 3x3 matrix, mtx, with the specified row and column data removed.
///
///  mtx = The 3x3 matrix to use as the basis of the new matrix.
///
///  row = The row to skip in the source matrix.
///
///  col = The col to skip in the source matrix.
///
///  returns = A new 2x2 matrix based on, mtx, with the specified row and column removed.
///
pub fn rmvRowColXmtx3(mtx: *[9]f32, row: f32, col: f32) [4]f32 {
    const cols: usize = 3;
    const destCols: usize = 2;
    const r: usize = 3;
    const c: usize = 3;
    var i: usize = 0;
    var j: usize = 0;
    var actI: usize = 0;
    var actJ: usize = 0;
    var dest: [4]f32 = std.mem.zeroes([4]f32);
    const actRow: usize = @intFromFloat(row);
    const actCol: usize = @intFromFloat(col);

    while (i < r) {
        j = 0;
        actJ = 0;
        if (i != actRow) {
            while (j < c) {
                if (j != actCol) {
                    dest[(actI * destCols) + actJ] = mtx[(i * cols) + j];
                    actJ += 1;
                }
                j += 1;
            }
            actI += 1;
        }
        i += 1;
    }
    return dest;
}

test "XMTX: rmvRowColXmtx3 test" {
    var A: [9]f32 = .{ 1, 1, 1, 2, 2, 2, 3, 3, 3 };
    //A = 1  1  1
    //    2  2  2
    //    3  3  3
    var expA: [4]f32 = .{ 2, 2, 3, 3 };
    const r: f32 = 0.0;
    const c: f32 = 0.0;
    var res: [4]f32 = std.mem.zeroes([4]f32); //.{};

    const start = try Instant.now();
    res = rmvRowColXmtx3(&A, r, c);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nrmvRowColXmtx3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rmvRowColXmtx3", elapsed1);

    std.debug.print("\nrmvRowColXmtx3:", .{});
    prntXmtxNl(&res, 2);
    try std.testing.expectEqual(true, equXmtx(&expA, &res));
    prntNl();
}

///Returns a new matrix based on the 2x2 matrix, mtx, with the specified row and column data removed.
///
///  mtx = The 2x2 matrix to use as the basis of the new matrix.
///
///  row = The row to skip in the source matrix.
///
///  col = The col to skip in the source matrix.
///
///  returns = A new 1x1 matrix based on, mtx, with the specified row and column removed.
///
pub fn rmvRowColXmtx2(mtx: *[4]f32, row: f32, col: f32) [1]f32 {
    const cols: usize = 2;
    const destCols: usize = 1;
    const r: usize = 2;
    const c: usize = 2;
    var i: usize = 0;
    var j: usize = 0;
    var actI: usize = 0;
    var actJ: usize = 0;
    var dest: [1]f32 = std.mem.zeroes([1]f32);
    const actRow: usize = @intFromFloat(row);
    const actCol: usize = @intFromFloat(col);

    while (i < r) {
        j = 0;
        actJ = 0;
        if (i != actRow) {
            while (j < c) {
                if (j != actCol) {
                    dest[(actI * destCols) + actJ] = mtx[(i * cols) + j];
                    actJ += 1;
                }
                j += 1;
            }
            actI += 1;
        }
        i += 1;
    }
    return dest;
}

test "XMTX: rmvRowColXmtx2 test" {
    var A: [4]f32 = .{ 1, 1, 2, 2 };
    //A = 1  1
    //    2  2
    var expA: [1]f32 = .{2};
    const r: f32 = 0.0;
    const c: f32 = 0.0;
    var res: [1]f32 = std.mem.zeroes([1]f32);

    const start = try Instant.now();
    res = rmvRowColXmtx2(&A, r, c);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nrmvRowColXmtx2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rmvRowColXmtx2", elapsed1);

    std.debug.print("\nrmvRowColXmtx2:", .{});
    prntXmtxNl(&res, 2);
    try std.testing.expectEqual(true, equXmtx(&expA, &res));
    prntNl();
}

///Returns the adjoint matrix for a given 4x4 matrix, mtx.
///
///  mtx = The 4x4 matrix to determine the adjoint matrix for.
///
///  returns = A new 4x4 matrix that contains the adjoint matrix of mtx.
///
pub fn adjXmtx4(mtx: *[16]f32) [16]f32 {
    var ret: [16]f32 = std.mem.zeroes([16]f32); //.{};
    var cofA: [16]f32 = cofXmtx4(mtx);
    trnXmtx(&cofA, 4, &ret);
    return ret;
}

test "XMTX: adjXmtx4 test" {
    var A: [16]f32 = .{ 1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1 };
    var retA: [16]f32 = std.mem.zeroes([16]f32); //.{};

    const start = try Instant.now();
    retA = adjXmtx4(&A);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nadjXmtx4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("adjXmtx4", elapsed1);

    var expA: [16]f32 = .{ -4, -4, -4, 4, -4, -4, 4, -4, -4, 4, -4, -4, 4, -4, -4, -4 };
    try std.testing.expectEqual(true, equXmtx(&expA, &retA));
    prntNl();
}

///Returns the adjoint of the given 3x3 matrix, mtx.
///
///  mtx = The 3x3 matrix to find the adjoint matrix for.
///
///  returns = A new 3x3 matrix with the transpose of the cofactors of the original matrix, mtx.
///
pub fn adjXmtx3(mtx: *[9]f32) [9]f32 {
    var ret: [9]f32 = std.mem.zeroes([9]f32);
    var cofA: [9]f32 = cofXmtx3(mtx);
    trnXmtx(&cofA, 3, &ret);
    return ret;
}

test "XMTX: adjXmtx3 test" {
    var A: [9]f32 = .{ -1, 3, 2, 0, -2, 1, 1, 0, -2 };
    var retA: [9]f32 = std.mem.zeroes([9]f32); //.{};

    const start = try Instant.now();
    retA = adjXmtx3(&A);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nadjXmtx3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("adjXmtx3", elapsed1);

    var expA: [9]f32 = .{ 4, 6, 7, 1, 0, 1, 2, 3, 2 };
    try std.testing.expectEqual(true, equXmtx(&expA, &retA));
    prntNl();
}

///Returns the adjoint of the given 2x2 matrix, mtx.
///
///  mtx = The 2x2 matrix to find the adjoint matrix for.
///
///  returns = A new 2x2 matrix with the transpose of the cofactors of the original matrix, mtx.
///
pub fn adjXmtx2(mtx: *[4]f32) [4]f32 {
    //|  d -b | or |  mtx[3] -mtx[1] |
    //| -c  a |    | -mtx[2]  mtx[0] |
    return [4]f32{ mtx[3], (-1.0 * mtx[1]), (-1.0 * mtx[2]), mtx[0] };
}

test "XMTX: adjXmtx2 test" {
    var A: [4]f32 = .{ -1, 3, 2, 0 };
    //| -1  3 |
    //|  2  0 |
    var retA: [4]f32 = std.mem.zeroes([4]f32);

    const start = try Instant.now();
    retA = adjXmtx2(&A);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nadjXmtx2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("adjXmtx2", elapsed1);

    var expA: [4]f32 = .{ 0, -3, -2, -1 };
    try std.testing.expectEqual(true, equXmtx(&expA, &retA));
    prntNl();
}

///Returns the adjoint of the given 1x1 matrix, mtx.
///
///  mtx = The 1x1 matrix to find the adjoint matrix for.
///
///  returns = A new 1x1 matrix with the transpose of the cofactors of the original matrix, mtx.
///
pub fn adjXmtx1(mtx: *[1]f32) [1]f32 {
    return [1]f32{mtx[0]};
}

test "XMTX: adjXmtx1 test" {
    var A: [1]f32 = .{-1};
    var retA: [1]f32 = std.mem.zeroes([1]f32);

    const start = try Instant.now();
    retA = adjXmtx1(&A);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nadjXmtx1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("adjXmtx1", elapsed1);

    var expA: [1]f32 = .{-1};
    try std.testing.expectEqual(true, equXmtx(&expA, &retA));
    prntNl();
}

///Returns a matrix of full co-factors, sign * minor for the given 4x4 matrix.
///
///  mtx = The 4x4 matrix to find cofactors for.
///
///  returns = A new 4x4 matrix of cofactors for the given matrix, mtx.
///
pub fn cofXmtx4(mtx: *[16]f32) [16]f32 {
    var mnr = mnrXmtx4(mtx);

    //std.debug.print("\n111 cofXmtx4 test", .{});
    //prntXmtx(&mnr, 4);

    const cofSgn = cofXmtxSign4();

    //std.debug.print("\n222 cofXmtx4 test", .{});
    //prntXmtx(&cofSgn, 4);

    var i: usize = 0;
    const l: usize = mtx.len;

    while (i < l) : (i += 1) {
        mnr[i] = (mnr[i] * cofSgn[i]);
    }

    return mnr;
}

test "XMTX: cofXmtx4 test" {
    //A = 1  2  3  4
    //    0  1  2  3
    //    -1 3  5  6
    //    1  1  1  1
    var A: [16]f32 = .{ 1, 3, -2, 1, 5, 1, 0, -1, 0, 1, 0, -2, 2, -1, 0, 3 };
    var exp: [16]f32 = .{ 0, 0, 3, 0, -2, 8, 13, 4, 4, -34, -56, -14, 2, -20, -34, -10 };
    var cof: [16]f32 = undefined;

    const start = try Instant.now();
    cof = cofXmtx4(&A);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncofXmtx4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cofXmtx4", elapsed1);

    std.debug.print("\nAAA cofXmtx4 test", .{});
    clnXmtx(&cof);
    prntXmtxNl(&cof, 4);

    std.debug.print("\nBBB cofXmtx4 test", .{});
    prntXmtxNl(&exp, 4);
    const b: bool = equXmtx(&cof, &exp);

    std.debug.print("\nCCC cofXmtx4: {}", .{b});
    std.debug.print("\nDDD cofXmtx4 test", .{});
    prntXmtxNl(&A, 4);

    try std.testing.expectEqual(true, b);
    prntNl();
}

///Returns a matrix of full cofactors, sign * minor, for the given 3x3 matrix.
///
///  mtx = The 3x3 matrix to find cofactors for.
///
///  returns = A new 3x3 matrix of cofactors for the given matrix, mtx.
///
pub fn cofXmtx3(mtx: *[9]f32) [9]f32 {
    var mnr = mnrXmtx3(mtx);
    const cofSgn = cofXmtxSign3();
    var i: usize = 0;
    const l: usize = mtx.len;

    while (i < l) : (i += 1) {
        mnr[i] = (mnr[i] * cofSgn[i]);
    }

    return mnr;
}

test "XMTX: cofXmtx3 test" {
    //A = 0  2  1
    //    3  -1 2
    //    4  0  1
    var A: [9]f32 = .{ 0, 2, 1, 3, -1, 2, 4, 0, 1 };
    var exp: [9]f32 = .{ -1, 5, 4, -2, -4, 8, 5, 3, -6 };
    var cof: [9]f32 = undefined;

    const start = try Instant.now();
    cof = cofXmtx3(&A);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncofXmtx3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cofXmtx3", elapsed1);

    std.debug.print("\nAAA cofXmtx3 test", .{});
    prntXmtxNl(&cof, 3);
    std.debug.print("\nBBB cofXmtx3 test", .{});
    prntXmtxNl(&exp, 3);

    const b: bool = equXmtx(&cof, &exp);
    std.debug.print("\nCCC cofXmtx3: {}", .{b});
    try std.testing.expectEqual(true, b);
    prntNl();
}

///Returns a matrix of full cofactors, sign * minor, for the given 2x2 matrix.
///
///  mtx = The 2x2 matrix to find cofactors for.
///
///  returns = A new 2x2 matrix of cofactors for the given matrix, mtx.
///
pub fn cofXmtx2(mtx: *[4]f32) [4]f32 {
    var mnr = mnrXmtx2(mtx);
    const cofSgn = cofXmtxSign2();
    var i: usize = 0;
    const l: usize = mtx.len;

    while (i < l) : (i += 1) {
        mnr[i] = (mnr[i] * cofSgn[i]);
    }

    return mnr;
}

test "XMTX: cofXmtx2 test" {
    //A = 0  2
    //    3  -1
    var A: [4]f32 = .{ 0, 2, 3, -1 };
    var exp: [4]f32 = .{ -1, -3, -2, 0 };
    var cof: [4]f32 = undefined;

    const start = try Instant.now();
    cof = cofXmtx2(&A);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncofXmtx2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cofXmtx2", elapsed1);

    std.debug.print("\nAAA cofXmtx2 test", .{});
    prntXmtxNl(&cof, 2);
    std.debug.print("\nBBB cofXmtx2 test", .{});
    prntXmtxNl(&exp, 2);

    const b: bool = equXmtx(&cof, &exp);
    std.debug.print("\nCCC cofXmtx2: {}", .{b});
    try std.testing.expectEqual(true, b);
    prntNl();
}

///Returns a matrix of full cofactors, sign * minor, for the given 1x1 matrix.
///
///  mtx = The 1x1 matrix to find cofactor for.
///
///  returns = A new 1x1 matrix with cofactor for the given matrix, mtx.
///
pub fn cofXmtx1(mtx: *[1]f32) [1]f32 {
    _ = mtx;
    return [1]f32{1.0};
}

test "XMTX: cofXmtx1 test" {
    //A = 0
    var A: [1]f32 = .{0};
    var exp: [1]f32 = .{1};
    var cof: [1]f32 = undefined;

    const start = try Instant.now();
    cof = cofXmtx1(&A);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncofXmtx1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cofXmtx1", elapsed1);

    std.debug.print("\nAAA cofXmtx1 test", .{});
    prntXmtxNl(&cof, 1);
    std.debug.print("\nBBB cofXmtx1 test", .{});
    prntXmtxNl(&exp, 1);

    const b: bool = equXmtx(&cof, &exp);
    std.debug.print("\nCCC cofXmtx1: {}", .{b});
    try std.testing.expectEqual(true, b);
    prntNl();
}

///Calculates a matrix of minors for the given 4x4 matrix.
///
///  mtx = The 4x4 matrix to find minors for.
///
///  ret = Stores a new 4x4 matrix of the calculated minors for the provided matrix, mtx.
///
pub fn mnrXmtx4Ret(mtx: *[16]f32, ret: *[16]f32) void {
    var T1: [9]f32 = .{ mtx[5], mtx[6], mtx[7], mtx[9], mtx[10], mtx[11], mtx[13], mtx[14], mtx[15] };
    var T2: [9]f32 = .{ mtx[4], mtx[6], mtx[7], mtx[8], mtx[10], mtx[11], mtx[12], mtx[14], mtx[15] };
    var T3: [9]f32 = .{ mtx[4], mtx[5], mtx[7], mtx[8], mtx[9], mtx[11], mtx[12], mtx[13], mtx[15] };
    var T4: [9]f32 = .{ mtx[4], mtx[5], mtx[6], mtx[8], mtx[9], mtx[10], mtx[12], mtx[13], mtx[14] };

    var T5: [9]f32 = .{ mtx[1], mtx[2], mtx[3], mtx[9], mtx[10], mtx[11], mtx[13], mtx[14], mtx[15] };
    var T6: [9]f32 = .{ mtx[0], mtx[2], mtx[3], mtx[8], mtx[10], mtx[11], mtx[12], mtx[14], mtx[15] };
    var T7: [9]f32 = .{ mtx[0], mtx[1], mtx[3], mtx[8], mtx[9], mtx[11], mtx[12], mtx[13], mtx[15] };
    var T8: [9]f32 = .{ mtx[0], mtx[1], mtx[2], mtx[8], mtx[9], mtx[10], mtx[12], mtx[13], mtx[14] };

    var T9: [9]f32 = .{ mtx[1], mtx[2], mtx[3], mtx[5], mtx[6], mtx[7], mtx[13], mtx[14], mtx[15] };
    var T10: [9]f32 = .{ mtx[0], mtx[2], mtx[3], mtx[4], mtx[6], mtx[7], mtx[12], mtx[14], mtx[15] };
    var T11: [9]f32 = .{ mtx[0], mtx[1], mtx[3], mtx[4], mtx[5], mtx[7], mtx[12], mtx[13], mtx[15] };
    var T12: [9]f32 = .{ mtx[0], mtx[1], mtx[2], mtx[4], mtx[5], mtx[6], mtx[12], mtx[13], mtx[14] };

    var T13: [9]f32 = .{ mtx[1], mtx[2], mtx[3], mtx[5], mtx[6], mtx[7], mtx[9], mtx[10], mtx[11] };
    var T14: [9]f32 = .{ mtx[0], mtx[2], mtx[3], mtx[4], mtx[6], mtx[7], mtx[8], mtx[10], mtx[11] };
    var T15: [9]f32 = .{ mtx[0], mtx[1], mtx[3], mtx[4], mtx[5], mtx[7], mtx[8], mtx[9], mtx[11] };
    var T16: [9]f32 = .{ mtx[0], mtx[1], mtx[2], mtx[4], mtx[5], mtx[6], mtx[8], mtx[9], mtx[10] };

    ret[0] = detXmtx3(&T1);
    ret[1] = detXmtx3(&T2);
    ret[2] = detXmtx3(&T3);
    ret[3] = detXmtx3(&T4);

    ret[4] = detXmtx3(&T5);
    ret[5] = detXmtx3(&T6);
    ret[6] = detXmtx3(&T7);
    ret[7] = detXmtx3(&T8);

    ret[8] = detXmtx3(&T9);
    ret[9] = detXmtx3(&T10);
    ret[10] = detXmtx3(&T11);
    ret[11] = detXmtx3(&T12);

    ret[12] = detXmtx3(&T13);
    ret[13] = detXmtx3(&T14);
    ret[14] = detXmtx3(&T15);
    ret[15] = detXmtx3(&T16);
}

test "XMTX: mnrXmtx4Ret test" {
    var A: [16]f32 = .{ 1, 0, 3, 5, 2, 1, 4, 6, 1, 2, 3, 4, 4, 3, 2, 1 };
    var mnrA: [16]f32 = std.mem.zeroes([16]f32); //.{};
    var exp: [16]f32 = .{ 5, -10, -35, -20, 5, -10, -35, -20, 1, -2, -7, -4, -1, 2, 7, 4 };

    const start = try Instant.now();
    mnrXmtx4Ret(&A, &mnrA);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nmnrXmtx4Ret: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("mnrXmtx4Ret", elapsed1);

    try std.testing.expectEqual(true, equXmtx(&mnrA, &exp));
    prntNl();
}

///Returns a matrix of minors for the given 4x4 matrix.
///
///  mtx = The 4x4 matrix to find minors for.
///
///  returns = A new 4x4 matrix of the calculated minors for the provided matrix, mtx.
///
pub fn mnrXmtx4(mtx: *[16]f32) [16]f32 {
    //A = a b c d     0  1  2  3
    //    e f g h     4  5  6  7
    //    i j k l     8  9  10 11
    //    m n o p     12 13 14 15

    //5  6  7
    //9  10 11
    //13 14 15
    var T1: [9]f32 = .{ mtx[5], mtx[6], mtx[7], mtx[9], mtx[10], mtx[11], mtx[13], mtx[14], mtx[15] };

    //4  6  7
    //8  10 11
    //12 14 15
    var T2: [9]f32 = .{ mtx[4], mtx[6], mtx[7], mtx[8], mtx[10], mtx[11], mtx[12], mtx[14], mtx[15] };

    //4  5  7
    //8  9  11
    //12 13 15
    var T3: [9]f32 = .{ mtx[4], mtx[5], mtx[7], mtx[8], mtx[9], mtx[11], mtx[12], mtx[13], mtx[15] };

    //4  5  6
    //8  9  10
    //12 13 14
    var T4: [9]f32 = .{ mtx[4], mtx[5], mtx[6], mtx[8], mtx[9], mtx[10], mtx[12], mtx[13], mtx[14] };

    //1  2  3
    //9  10 11
    //13 14 15
    var T5: [9]f32 = .{ mtx[1], mtx[2], mtx[3], mtx[9], mtx[10], mtx[11], mtx[13], mtx[14], mtx[15] };

    //0  2  3
    //8  10 11
    //12 14 15
    var T6: [9]f32 = .{ mtx[0], mtx[2], mtx[3], mtx[8], mtx[10], mtx[11], mtx[12], mtx[14], mtx[15] };

    //0  1  3
    //8  9  11
    //12 13 15
    var T7: [9]f32 = .{ mtx[0], mtx[1], mtx[3], mtx[8], mtx[9], mtx[11], mtx[12], mtx[13], mtx[15] };

    //0  1  2
    //8  9  10
    //12 13 14
    var T8: [9]f32 = .{ mtx[0], mtx[1], mtx[2], mtx[8], mtx[9], mtx[10], mtx[12], mtx[13], mtx[14] };

    //1  2  3
    //5  6  7
    //13 14 15
    var T9: [9]f32 = .{ mtx[1], mtx[2], mtx[3], mtx[5], mtx[6], mtx[7], mtx[13], mtx[14], mtx[15] };

    //0  2  3
    //4  6  7
    //12 14 15
    var T10: [9]f32 = .{ mtx[0], mtx[2], mtx[3], mtx[4], mtx[6], mtx[7], mtx[12], mtx[14], mtx[15] };

    //0  1  3
    //4  5  7
    //12 13 15
    var T11: [9]f32 = .{ mtx[0], mtx[1], mtx[3], mtx[4], mtx[5], mtx[7], mtx[12], mtx[13], mtx[15] };

    //0  1  2
    //4  5  6
    //12 13 14
    var T12: [9]f32 = .{ mtx[0], mtx[1], mtx[2], mtx[4], mtx[5], mtx[6], mtx[12], mtx[13], mtx[14] };

    //1  2  3
    //5  6  7
    //9  10 11
    var T13: [9]f32 = .{ mtx[1], mtx[2], mtx[3], mtx[5], mtx[6], mtx[7], mtx[9], mtx[10], mtx[11] };

    //0  2  3
    //4  6  7
    //8  10 11
    var T14: [9]f32 = .{ mtx[0], mtx[2], mtx[3], mtx[4], mtx[6], mtx[7], mtx[8], mtx[10], mtx[11] };

    //0  1  3
    //4  5  7
    //8  9  11
    var T15: [9]f32 = .{ mtx[0], mtx[1], mtx[3], mtx[4], mtx[5], mtx[7], mtx[8], mtx[9], mtx[11] };

    //0  1  2
    //4  5  6
    //8  9  10
    var T16: [9]f32 = .{ mtx[0], mtx[1], mtx[2], mtx[4], mtx[5], mtx[6], mtx[8], mtx[9], mtx[10] };

    return .{
        detXmtx3(&T1),
        detXmtx3(&T2),
        detXmtx3(&T3),
        detXmtx3(&T4),

        detXmtx3(&T5),
        detXmtx3(&T6),
        detXmtx3(&T7),
        detXmtx3(&T8),

        detXmtx3(&T9),
        detXmtx3(&T10),
        detXmtx3(&T11),
        detXmtx3(&T12),

        detXmtx3(&T13),
        detXmtx3(&T14),
        detXmtx3(&T15),
        detXmtx3(&T16),
    };
}

test "XMTX: mnrXmtx4 test" {
    var A: [16]f32 = .{ 1, 0, 3, 5, 2, 1, 4, 6, 1, 2, 3, 4, 4, 3, 2, 1 };
    var mnrA: [16]f32 = std.mem.zeroes([16]f32); //.{};
    var exp: [16]f32 = .{ 5, -10, -35, -20, 5, -10, -35, -20, 1, -2, -7, -4, -1, 2, 7, 4 };

    const start = try Instant.now();
    mnrA = mnrXmtx4(&A);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nmnrXmtx4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("mnrXmtx4", elapsed1);

    try std.testing.expectEqual(true, equXmtx(&mnrA, &exp));
    prntNl();
}

///Calculates a matrix of minors for the given 3x3 matrix.
///
///  mtx = The 3x3 matrix to find minors for.
///
///  ret = Stores a new 3x3 matrix of the calculated minors for the provided matrix, mtx.
///
pub fn mnrXmtx3Ret(mtx: *[9]f32, ret: *[9]f32) void {
    var T1: [4]f32 = .{ mtx[4], mtx[5], mtx[7], mtx[8] };
    var T2: [4]f32 = .{ mtx[3], mtx[5], mtx[6], mtx[8] };
    var T3: [4]f32 = .{ mtx[3], mtx[4], mtx[6], mtx[7] };

    var T4: [4]f32 = .{ mtx[1], mtx[2], mtx[7], mtx[8] };
    var T5: [4]f32 = .{ mtx[0], mtx[2], mtx[6], mtx[8] };
    var T6: [4]f32 = .{ mtx[0], mtx[1], mtx[6], mtx[7] };

    var T7: [4]f32 = .{ mtx[1], mtx[2], mtx[4], mtx[5] };
    var T8: [4]f32 = .{ mtx[0], mtx[2], mtx[3], mtx[5] };
    var T9: [4]f32 = .{ mtx[0], mtx[1], mtx[3], mtx[4] };

    ret[0] = detXmtx2(&T1);
    ret[1] = detXmtx2(&T2);
    ret[2] = detXmtx2(&T3);

    ret[3] = detXmtx2(&T4);
    ret[4] = detXmtx2(&T5);
    ret[5] = detXmtx2(&T6);

    ret[6] = detXmtx2(&T7);
    ret[7] = detXmtx2(&T8);
    ret[8] = detXmtx2(&T9);
}

test "XMTX: mnrXmtx3Ret test" {
    var A: [9]f32 = .{ 1, 0, 3, 2, 1, 4, 1, 2, 3 };
    var mnrA: [9]f32 = std.mem.zeroes([9]f32); //.{};
    var exp: [9]f32 = .{ -5, 2, 3, -6, 0, 2, -3, -2, 1 };

    const start = try Instant.now();
    mnrXmtx3Ret(&A, &mnrA);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nmnrXmtx3Ret: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("mnrXmtx3Ret", elapsed1);

    try std.testing.expectEqual(true, equXmtx(&mnrA, &exp));
    prntNl();
}

///Returns a matrix of minors for the given 3x3 matrix.
///
///  mtx = The 3x3 matrix to find minors for.
///
///  returns = A new 3x3 matrix of the calculated minors for the provided matrix, mtx.
///
pub fn mnrXmtx3(mtx: *[9]f32) [9]f32 {
    //Given 3x3 matrix A
    //A = a b c     0 1 2
    //    d e f     3 4 5
    //    g h i     6 7 8
    var T1: [4]f32 = .{ mtx[4], mtx[5], mtx[7], mtx[8] };
    var T2: [4]f32 = .{ mtx[3], mtx[5], mtx[6], mtx[8] };
    var T3: [4]f32 = .{ mtx[3], mtx[4], mtx[6], mtx[7] };

    var T4: [4]f32 = .{ mtx[1], mtx[2], mtx[7], mtx[8] };
    var T5: [4]f32 = .{ mtx[0], mtx[2], mtx[6], mtx[8] };
    var T6: [4]f32 = .{ mtx[0], mtx[1], mtx[6], mtx[7] };

    var T7: [4]f32 = .{ mtx[1], mtx[2], mtx[4], mtx[5] };
    var T8: [4]f32 = .{ mtx[0], mtx[2], mtx[3], mtx[5] };
    var T9: [4]f32 = .{ mtx[0], mtx[1], mtx[3], mtx[4] };

    return .{ detXmtx2(&T1), detXmtx2(&T2), detXmtx2(&T3), detXmtx2(&T4), detXmtx2(&T5), detXmtx2(&T6), detXmtx2(&T7), detXmtx2(&T8), detXmtx2(&T9) };
}

test "XMTX: mnrXmtx3 test" {
    var A: [9]f32 = .{ 1, 0, 3, 2, 1, 4, 1, 2, 3 };
    var mnrA: [9]f32 = std.mem.zeroes([9]f32); //.{};
    var exp: [9]f32 = .{ -5, 2, 3, -6, 0, 2, -3, -2, 1 };

    const start = try Instant.now();
    mnrA = mnrXmtx3(&A);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nmnrXmtx3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("mnrXmtx3", elapsed1);

    try std.testing.expectEqual(true, equXmtx(&mnrA, &exp));
    prntNl();
}

///Calculates a matrix of minors for the given 2x2 matrix.
///
///  mtx = The 2x2 matrix to find minors for.
///
///  ret = Stores a new 2x2 matrix of the calculated minors for the provided matrix, mtx.
///
pub fn mnrXmtx2Ret(mtx: *[4]f32, ret: *[4]f32) void {
    //Given 2x2 matrix A
    //A = a b     0 1
    //    c d     2 3
    var T1: [1]f32 = .{mtx[3]};
    var T2: [1]f32 = .{mtx[2]};
    var T3: [1]f32 = .{mtx[1]};
    var T4: [1]f32 = .{mtx[0]};

    ret[0] = detXmtx1(&T1);
    ret[1] = detXmtx1(&T2);
    ret[2] = detXmtx1(&T3);
    ret[3] = detXmtx1(&T4);
}

test "XMTX: mnrXmtx2Ret test" {
    var A: [4]f32 = .{ 1, -8, 5, 4 };
    var mnrA: [4]f32 = std.mem.zeroes([4]f32);
    var exp: [4]f32 = .{ 4, 5, -8, 1 };

    const start = try Instant.now();
    mnrXmtx2Ret(&A, &mnrA);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nmnrXmtx2Ret: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("mnrXmtx2Ret", elapsed1);

    try std.testing.expectEqual(true, equXmtx(&mnrA, &exp));
    prntNl();
}

///Returns a matrix of minors for the given 2x2 matrix.
///
///  mtx = The 2x2 matrix to find minors for.
///
///  returns = A new 2x2 matrix of the calculated minors for the provided matrix, mtx.
///
pub fn mnrXmtx2(mtx: *[4]f32) [4]f32 {
    //Given 2x2 matrix A
    //A = a b     0 1
    //    c d     2 3
    var T1: [1]f32 = .{mtx[3]};
    var T2: [1]f32 = .{mtx[2]};
    var T3: [1]f32 = .{mtx[1]};
    var T4: [1]f32 = .{mtx[0]};

    return .{ detXmtx1(&T1), detXmtx1(&T2), detXmtx1(&T3), detXmtx1(&T4) };
}

test "XMTX: mnrXmtx2 test" {
    var A: [4]f32 = .{ 1, -8, 5, 4 };
    var mnrA: [4]f32 = std.mem.zeroes([4]f32);
    var exp: [4]f32 = .{ 4, 5, -8, 1 };

    const start = try Instant.now();
    mnrA = mnrXmtx2(&A);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nmnrXmtx2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("mnrXmtx2", elapsed1);

    try std.testing.expectEqual(true, equXmtx(&mnrA, &exp));
    prntNl();
}

///Calculates a matrix of minors for the given 1x1 matrix.
///
///  mtx = The 1x1 matrix to find minors for.
///
///  ret = Stores a new 1x1 matrix of the calculated minors for the provided matrix, mtx.
///
pub fn mnrXmtx1Ret(mtx: *[1]f32, ret: *[1]f32) void {
    _ = mtx[0];
    ret[0] = 1.0; //det[A] = a * det[] = a * 1
}

test "XMTX: mnrXmtx1Ret test" {
    var A: [1]f32 = .{8};
    var mnrA: [1]f32 = std.mem.zeroes([1]f32);
    var exp: [1]f32 = .{1};

    const start = try Instant.now();
    mnrXmtx1Ret(&A, &mnrA);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nmnrXmtx1Ret: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("mnrXmtx1Ret", elapsed1);

    try std.testing.expectEqual(true, equXmtx(&mnrA, &exp));
    prntNl();
}

///Returns a matrix of minors for the given 1x1 matrix.
///
///  mtx = The 1x1 matrix to find minors for.
///
///  returns = A new 1x1 matrix of the calculated minors for the provided matrix, mtx.
///
pub fn mnrXmtx1(mtx: *[1]f32) [1]f32 {
    _ = mtx;
    return .{1.0};
}

test "XMTX: mnrXmtx1 test" {
    var A: [1]f32 = .{8};
    var mnrA: [1]f32 = std.mem.zeroes([1]f32);
    var exp: [1]f32 = .{1};

    const start = try Instant.now();
    mnrA = mnrXmtx1(&A);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nmnrXmtx1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("mnrXmtx1", elapsed1);

    try std.testing.expectEqual(true, equXmtx(&mnrA, &exp));
    prntNl();
}

///Returns a matrix of signed values, 1 or -1, for the associated matrix of minors to determine the cofactor matrix for a given 4x4 matrix.
pub fn cofXmtxSign4() [16]f32 {
    return .{ 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1 };
}

test "XMTX: cofXmtxSign4 test" {
    var cofSgn: [16]f32 = undefined;

    const start = try Instant.now();
    cofSgn = cofXmtxSign4();
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncofXmtxSign4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cofXmtxSign4", elapsed1);

    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[0]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[1]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[2]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[3]);

    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[4]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[5]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[6]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[7]);

    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[8]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[9]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[10]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[11]);

    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[12]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[13]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[14]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[15]);
    prntNl();
}

///Calculates and stores, in ret, a matrix of signed values, 1 or -1, for the associated matrix of minors to determine the cofactor matrix for a given 4x4 matrix.
pub fn cofXmtxSign4Ret(ret: *[16]f32) void {
    ret[0] = 1;
    ret[1] = -1;
    ret[2] = 1;
    ret[3] = -1;

    ret[4] = -1;
    ret[5] = 1;
    ret[6] = -1;
    ret[7] = 1;

    ret[8] = 1;
    ret[9] = -1;
    ret[10] = 1;
    ret[11] = -1;

    ret[12] = -1;
    ret[13] = 1;
    ret[14] = -1;
    ret[15] = 1;

    //return .{ 1,-1, 1,-1,
    //         -1, 1,-1, 1,
    //          1,-1, 1,-1,
    //         -1, 1,-1, 1
    //};
}

test "XMTX: cofXmtxSign4Ret test" {
    var cofSgn: [16]f32 = undefined;

    const start = try Instant.now();
    cofXmtxSign4Ret(&cofSgn);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncofXmtxSign4Ret: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cofXmtxSign4Ret", elapsed1);

    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[0]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[1]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[2]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[3]);

    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[4]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[5]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[6]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[7]);

    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[8]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[9]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[10]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[11]);

    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[12]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[13]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[14]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[15]);
    prntNl();
}

///Returns a matrix of signed values, 1 or -1, for the associated matrix of minors to determine the cofactor matrix for a given 3x3 matrix.
pub fn cofXmtxSign3() [9]f32 {
    return .{ 1, -1, 1, -1, 1, -1, 1, -1, 1 };
}

test "XMTX: cofXmtxSign3 test" {
    var cofSgn: [9]f32 = undefined;

    const start = try Instant.now();
    cofSgn = cofXmtxSign3();
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncofXmtxSign3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cofXmtxSign3", elapsed1);

    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[0]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[1]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[2]);

    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[3]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[4]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[5]);

    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[6]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[7]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[8]);
    prntNl();
}

///Calculates and stores, in ret, a matrix of signed values, 1 or -1, for the associated matrix of minors to determine the cofactor matrix for a given 3x3 matrix.
pub fn cofXmtxSign3Ret(ret: *[9]f32) void {
    ret[0] = 1;
    ret[1] = -1;
    ret[2] = 1;

    ret[3] = -1;
    ret[4] = 1;
    ret[5] = -1;

    ret[6] = 1;
    ret[7] = -1;
    ret[8] = 1;

    //return .{ 1,-1, 1,
    //           -1, 1,-1,
    //            1,-1, 1 };
}

test "XMTX: cofXmtxSign3Ret test" {
    var cofSgn: [9]f32 = undefined;

    const start = try Instant.now();
    cofXmtxSign3Ret(&cofSgn);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncofXmtxSign3Ret: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cofXmtxSign3Ret", elapsed1);

    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[0]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[1]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[2]);

    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[3]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[4]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[5]);

    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[6]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[7]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[8]);
    prntNl();
}

///Returns a matrix of signed values, 1 or -1, for the associated matrix of minors to determine the cofactor matrix for a given 2x2 matrix.
pub fn cofXmtxSign2() [4]f32 {
    return .{ 1, -1, -1, 1 };
}

test "XMTX: cofXmtxSign2 test" {
    var cofSgn: [4]f32 = undefined;

    const start = try Instant.now();
    cofSgn = cofXmtxSign2();
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncofXmtxSign2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cofXmtxSign2", elapsed1);

    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[0]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[1]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[2]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[3]);
    prntNl();
}

///Calculates and stores, in ret, a matrix of signed values, 1 or -1, for the associated matrix of minors to determine the cofactor matrix for a given 2x2 matrix.
pub fn cofXmtxSign2Ret(ret: *[4]f32) void {
    ret[0] = 1;
    ret[1] = -1;
    ret[2] = -1;
    ret[3] = 1;
    //return .{ 1,-1,
    //         -1, 1 };
}

test "XMTX: cofXmtxSign2Ret test" {
    var cofSgn: [4]f32 = undefined;

    const start = try Instant.now();
    cofXmtxSign2Ret(&cofSgn);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncofXmtxSign2Ret: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cofXmtxSign2Ret", elapsed1);

    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[0]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[1]);
    try std.testing.expectEqual(@as(f32, -1.0), cofSgn[2]);
    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[3]);
    prntNl();
}

///Returns a matrix of signed values, 1 or -1, for the associated matrix of minors to determine the cofactor matrix for a given 1x1 matrix.
pub fn cofXmtxSign1() [1]f32 {
    return .{1};
}

test "XMTX: cofXmtxSign1 test" {
    var cofSgn: [1]f32 = undefined;

    const start = try Instant.now();
    cofSgn = cofXmtxSign1();
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncofXmtxSign1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cofXmtxSign1", elapsed1);

    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[0]);
    prntNl();
}

///Calculates and stores, in ret, a matrix of signed values, 1 or -1, for the associated matrix of minors to determine the cofactor matrix for a given 1x1 matrix.
pub fn cofXmtxSign1Ret(ret: *[1]f32) void {
    ret[0] = 1;
}

test "XMTX: cofXmtxSign1Ret test" {
    var cofSgn: [1]f32 = undefined;

    const start = try Instant.now();
    cofXmtxSign1Ret(&cofSgn);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncofXmtxSign1Ret: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cofXmtxSign1Ret", elapsed1);

    try std.testing.expectEqual(@as(f32, 1.0), cofSgn[0]);
    prntNl();
}

///Returns the determinant of a 1x1 matrix.
///
///  mtx = The 1x1 matrix to find the determinant for.
///
///  returns = The value of the determinant.
///
pub fn detXmtx1(mtx: *[1]f32) f32 {
    return mtx[0];
}

test "XMTX: detXmtx1 test" {
    var m1: [1]f32 = .{2};
    var detM1: f32 = undefined;

    const start = try Instant.now();
    detM1 = detXmtx1(&m1);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ndetXmtx1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("detXmtx1", elapsed1);

    try std.testing.expectEqual(m1[0], detM1);
    prntNl();
}

///Returns the determinant of a 2x2 matrix.
///
///  mtx = The 2x2 matrix to find the determinant for.
///
///  returns = The determinant of the diagonal matrix, mtx, provided.
///
pub fn detXmtx2(mtx: *[4]f32) f32 {
    return ((mtx[0] * mtx[3]) - (mtx[1] * mtx[2]));
}

test "XMTX: detXmtx2 test" {
    var m1: [4]f32 = .{ 2, 1, 7, 4 };
    var detM1: f32 = undefined;

    const start = try Instant.now();
    detM1 = detXmtx2(&m1);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ndetXmtx2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("detXmtx2", elapsed1);

    try std.testing.expectEqual(((m1[0] * m1[3]) - (m1[1] * m1[2])), detM1);
    prntNl();
}

///Returns the determinant of a 3x3 matrix.
///
///  mtx = The 3x3 matrix to find the determinant for.
///
///  returns = The determinant of the diagonal matrix, mtx, provided.
///
pub fn detXmtx3(mtx: *[9]f32) f32 {
    //11 12 13 11 12
    //21 22 23 21 22
    //31 32 33 31 32
    //0  1  2
    //3  4  5
    //6  7  8
    return ((mtx[0] * mtx[4] * mtx[8]) //11 22 33
    + (mtx[1] * mtx[5] * mtx[6]) //12 23 31
    + (mtx[2] * mtx[3] * mtx[7]) //13 21 32
    - (mtx[6] * mtx[4] * mtx[2]) //31 22 13
    - (mtx[7] * mtx[5] * mtx[0]) //32 23 11
    - (mtx[8] * mtx[3] * mtx[1])); //33 21 12
}

test "XMTX: detXmtx3 test" {
    prntNl();
    var A: [9]f32 = .{ 6, 1, 1, 4, -2, 5, 2, 8, 7 };
    const exp: f32 = -306;
    var detA: f32 = undefined;

    const start = try Instant.now();
    detA = detXmtx3(&A);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ndetXmtx3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("detXmtx3", elapsed1);

    try std.testing.expectEqual(exp, detA);
    prntNl();
}

///Returns the determinant of a 4x4 matrix.
///
///  mtx = The 4x4 matrix to find the determinant for.
///
///  returns = The determinant of the diagonal matrix, mtx, provided.
///
pub fn detXmtx4(mtx: *[16]f32) f32 {
    var d1: [9]f32 = rmvRowColXmtx4(mtx, 0, 0);
    const det1: f32 = detXmtx3(&d1);
    var d2: [9]f32 = rmvRowColXmtx4(mtx, 0, 1);
    const det2: f32 = detXmtx3(&d2);
    var d3: [9]f32 = rmvRowColXmtx4(mtx, 0, 2);
    const det3: f32 = detXmtx3(&d3);
    var d4: [9]f32 = rmvRowColXmtx4(mtx, 0, 3);
    const det4: f32 = detXmtx3(&d4);
    return ((mtx[0] * det1) - (mtx[1] * det2) + (mtx[2] * det3) - (mtx[3] * det4));
}

test "XMTX: detXmtx4 test" {
    var A: [16]f32 = .{ 4, 3, 2, 2, 0, 1, -3, 3, 0, -1, 3, 3, 0, 3, 1, 1 };
    const expA: f32 = -240.0;
    var detA: f32 = undefined;

    const start = try Instant.now();
    detA = detXmtx4(&A);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ndetXmtx4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("detXmtx4", elapsed1);

    try std.testing.expectEqual(true, isEquF32(expA, detA, true));
    prntNl();
}

///Returns the determinant of a diagonal or triangular matrix.
///
///  mtx = A matrix to calculate the determinant for.
///
///  cols = The number of columns and rows in the diagonal matrix.
///
///  returns = The determinant of the diagonal matrix, mtx, provided.
///
pub fn detTriangXmtx(mtx: []f32, cols: usize) f32 {
    return detDiagXmtx(mtx, cols);
}

test "XMTX: detTriangXmtx test" {
    var A: [9]f32 = .{ 2, 0, 0, 0, 2, 0, 0, 0, 2 };
    var B: [16]f32 = .{ 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 3 };
    var detA: f32 = 0;
    var detB: f32 = 0;

    var start = try Instant.now();
    detA = detTriangXmtx(&A, 3);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ndetTriangXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("detTriangXmtx", elapsed1);

    start = try Instant.now();
    detB = detTriangXmtx(&B, 4);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ndetTriangXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("detTriangXmtx", elapsed1);

    std.debug.print("\ndetA: {} detB: {}", .{ detA, detB });
    try std.testing.expectEqual(true, isEquF32(detA, 8.0, true));
    try std.testing.expectEqual(true, isEquF32(detB, 81.0, true));
    prntNl();
}

///Returns the determinant of a diagonal or triangular matrix.
///
///  mtx = A matrix to calculate the determinant for.
///
///  cols = The number of columns and rows in the diagonal matrix.
///
///  returns = The determinant of the diagonal matrix, mtx, provided.
///
pub fn detDiagXmtx(mtx: []f32, cols: usize) f32 {
    var row: usize = 0;
    var col: usize = 0;
    var ret: f32 = 1;

    while (row < cols) : (row += 1) {
        col = 0;
        while (col < cols) : (col += 1) {
            if (row == col) {
                ret *= mtx[(row * cols) + col];
            }
        }
    }
    return ret;
}

test "XMTX: detDiagXmtx test" {
    var A: [9]f32 = .{ 2, 0, 0, 0, 2, 0, 0, 0, 2 };
    var B: [16]f32 = .{ 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 3 };
    var detA: f32 = 0;
    var detB: f32 = 0;

    var start = try Instant.now();
    detA = detDiagXmtx(&A, 3);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ndetTriangXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("detTriangXmtx", elapsed1);

    start = try Instant.now();
    detB = detDiagXmtx(&B, 4);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ndetTriangXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("detTriangXmtx", elapsed1);

    std.debug.print("\ndetA: {} detB: {}", .{ detA, detB });
    try std.testing.expectEqual(true, isEquF32(detA, 8.0, true));
    try std.testing.expectEqual(true, isEquF32(detB, 81.0, true));
    prntNl();
}

///Returns the determinant of the specified matrix for the given matrix row.
///
///  mtx = The matrix to find the determinant for.
///
///  cols = The number of columns in the matrix.
///
///  alloc = The allocator used by the function to create a return matrix for the determinants as the recursion makes them smaller each iteration.
///
///  cofR = The row to use as the basis for calculating the determinant.
///
///  returns = The determinant value for the given matrix,mtx.
///
pub fn detXmtx(mtx: []f32, cols: usize, alloc: *const std.mem.Allocator, cofR: usize) !f32 {
    var c: usize = 0;
    var a: f32 = 0;
    const rows: usize = mtx.len / cols;
    const resRows = rows - 1;
    const resCols = cols - 1;
    var b: bool = false;

    if (resCols == 0) {
        return mtx[0];
    } else if (resCols == 1) {
        //a b
        //c d
        //ad - bc
        //std.debug.print("\n({} * {}) - ({} * {}) times = {}:", .{mtx[0], mtx[3], mtx[1], mtx[2], (((mtx[0] * mtx[3]) - (mtx[1] * mtx[2]))) });
        return ((mtx[0] * mtx[3]) - (mtx[1] * mtx[2]));
    } else {
        var tot: f32 = 0;
        const res: []f32 = try alloc.*.alloc(f32, resRows * resCols);
        while (c < cols) : (c += 1) {
            b = cofXmtx(mtx, cols, cofR, c, res, resCols);
            a = mtx[c];

            //const na: f32 = @intToFloat(f32, (c + 2));
            //a *= std.math.pow(f32, -1.0, na);
            a *= cofXmtxSign(cofR, c, true);
            var nv: f32 = 0;

            if (!b) {
                //warning
                return error.OutOfMemory;
            } else {
                nv = try detXmtx(res, resCols, alloc, 0);
                tot += (a * nv);
            }

            //std.debug.print("\na {} times nv {} {} res[0] = {}:", .{a, nv, (a * nv), tot});
            //prntXmtxNl(res, resCols);
            //prntNl();
        }

        alloc.free(res);
        return tot;
    }
}

test "XMTX: detXmtx test" {
    const detRow: usize = 0;
    var m2: [4]f32 = .{ 5, 6, 8, 9 };
    var cols: usize = 2;

    //a b
    //c d
    //ad - bc
    //5(9) - 6(8) = -3
    const alloc: std.mem.Allocator = std.testing.allocator;
    var v: f32 = -1.0;

    var start = try Instant.now();
    v = try detXmtx(&m2, cols, &alloc, detRow);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ndetXmtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("detXmtx", elapsed1);

    var exp: f32 = -3.0;
    std.debug.print("\nFound detXmtx 2x2 value: {}", .{v});
    try std.testing.expectEqual(exp, v);
    var m3: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    cols = 3;

    start = try Instant.now();
    v = try detXmtx(&m3, cols, &alloc, detRow);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ndetXmtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("detXmtx", elapsed1);

    exp = 0;
    std.debug.print("\nFound detXmtx 3x3 value: {}", .{v});
    try std.testing.expectEqual(exp, v);
    var m4: [16]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };
    cols = 4;

    start = try Instant.now();
    v = try detXmtx(&m4, cols, &alloc, detRow);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ndetXmtx #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("detXmtx", elapsed1);

    exp = 0;
    std.debug.print("\nFound detXmtx 4x4 value: {}", .{v});
    try std.testing.expectEqual(exp, v);
    prntNl();
}

///Calculates the inverse of the provided 1x1 matrix using the determinant and stores the result in the return matrix.
///
///  mtx = The 1x1 matrix to calculate the inverse for.
///
///  det = The determinant of the provided matrix.
///
///  ret = The return 2x2 matrix to hold the inverted matrix values.
///
///  returns = A Boolean value indicating the operation was a success.
///
pub fn getInvFromDet1(mtx: *[1]f32, det: f32, ret: *[1]f32) void {
    _ = mtx;
    ret[0] = (1.0 / det);
}

test "XMTX: getInvFromDet1 test" {
    var m1: [1]f32 = .{7.0};
    var invM1: [1]f32 = .{(1.0 / 7.0)};
    var res: [1]f32 = std.mem.zeroes([1]f32);
    const detM1: f32 = detXmtx1(&m1);

    const start = try Instant.now();
    getInvFromDet1(&m1, detM1, &res);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ngetInvFromDet1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getInvFromDet1", elapsed1);

    clnXmtx(&res);
    clnXmtx(&invM1);

    prntNlStr("Res:");
    prntXmtxNl(&res, 1);
    prntNl();

    prntNlStr("invM1:");
    prntXmtxNl(&invM1, 1);
    prntNl();

    try std.testing.expectEqual(true, equXmtx(&res, &invM1));
    prntNl();
}

///Returns the inverse of a 1x1 matrix which is the reciprocal of the single matrix entry.
///
///  mtx = The 1x1 matrix to calculate the inverse for.
///
///  returns = The inverse of the provided 1x1 matrix, mtx.
///
pub fn getInvXmtx1(mtx: *[1]f32) [1]f32 {
    return .{(1.0 / mtx[0])};
}

test "XMTX: getInvXmtx1 test" {
    var m1: [1]f32 = .{7.0};
    var invM1: [1]f32 = .{(1.0 / 7.0)};
    var res: [1]f32 = std.mem.zeroes([1]f32);

    const start = try Instant.now();
    res = getInvXmtx1(&m1);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ngetInvXmtx1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getInvXmtx1", elapsed1);

    clnXmtx(&res);
    clnXmtx(&invM1);

    prntNlStr("Res:");
    prntXmtxNl(&res, 1);
    prntNl();

    prntNlStr("invM1:");
    prntXmtxNl(&invM1, 1);
    prntNl();

    try std.testing.expectEqual(true, equXmtx(&res, &invM1));
    prntNl();
}

///Calculates the inverse of the provided 2x2 matrix using the determinant and stores the result in the return matrix.
///
///  mtx = The 2x2 matrix to calculate the inverse for.
///
///  det = The determinant of the provided matrix.
///
///  ret = The return 2x2 matrix to hold the inverted matrix values.
///
///  returns = A Boolean value indicating the operation was a success.
///
pub fn getInvFromDet2(mtx: *[4]f32, det: f32, ret: *[4]f32) void {
    const cols: usize = 2;
    ret[(0 * cols) + 0] = mtx[(1 * cols) + 1] / det; //set 0x0 from 1x1
    ret[(1 * cols) + 1] = mtx[(0 * cols) + 0] / det; //set 1x1 from 0x0
    ret[(0 * cols) + 1] = -1.0 * mtx[(0 * cols) + 1] / det; //set 0x1
    ret[(1 * cols) + 0] = -1.0 * mtx[(1 * cols) + 0] / det; //set 1x0
}

test "XMTX: getInvFromDet2 test" {
    var m1: [4]f32 = .{ 2, 1, 7, 4 };
    var invM1: [4]f32 = .{ 4, -1, -7, 2 };
    var res: [4]f32 = std.mem.zeroes([4]f32); //.{};
    const detM1: f32 = detXmtx2(&m1);

    const start = try Instant.now();
    getInvFromDet2(&m1, detM1, &res);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ngetInvFromDet2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getInvFromDet2", elapsed1);

    try std.testing.expectEqual(true, equXmtx(&res, &invM1));
    prntNl();
}

///Returns the inverse of a 2x2 matrix with direct determinant calculation.
///
///  mtx = The 2x2 matrix to calculate the inverse for.
///
///  returns = The inverse of the provided 2x2 matrix, mtx.
///
pub fn getInvXmtx2(mtx: *[4]f32) ![4]f32 {
    if (((mtx[0] * mtx[3]) - (mtx[1] * mtx[2])) == 0) {
        prntNlStr("!! Warning the getInvXmtx2 function expects the mtx argument to be invertible, non singular !!");
        return Error.NonInvertibleMatrix;
    }

    const nmtx: [4]f32 = .{ mtx[3], (-1.0 * mtx[1]), (-1.0 * mtx[2]), mtx[0] };
    mulXvec(@constCast(&nmtx), (1.0 / ((mtx[0] * mtx[3]) - (mtx[1] * mtx[2]))));
    return nmtx;
}

test "XMTX: getInvXmtx2 test" {
    var m1: [4]f32 = .{ 2, 1, 7, 4 };
    var invM1: [4]f32 = .{ 4, -1, -7, 2 };
    var res: [4]f32 = std.mem.zeroes([4]f32);
    var res2: [4]f32 = std.mem.zeroes([4]f32);
    const detM1: f32 = detXmtx2(&m1);

    var start = try Instant.now();
    getInvFromDet2(&m1, detM1, &res);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ngetInvFromDet2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getInvFromDet2", elapsed1);

    try std.testing.expectEqual(true, equXmtx(&res, &invM1));
    prntNl();

    start = try Instant.now();
    res2 = try getInvXmtx2(&m1);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ngetInvXmtx2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getInvXmtx2", elapsed1);

    try std.testing.expectEqual(true, equXmtx(&res2, &res));
    prntNl();
}

///Calculates the inverse of the provided 3x3 matrix using the determinant and stores the result in the return matrix.
///
///  mtx = The 3x3 matrix to calculate the inverse for.
///
///  det = The determinant of the provided matrix.
///
///  ret = The return 3x3 matrix to hold the inverted matrix values.
///
///  returns = A Boolean value indicating if the operation was a success.
///
pub fn getInvFromDet3(mtx: *[9]f32, det: f32, ret: *[9]f32) void {
    const cols: usize = 3;
    //COL 0
    //0x0                  1x1                   2x2                     1x2                   2x1
    ret[(0 * cols) + 0] = (mtx[(1 * cols) + 1] * mtx[(2 * cols) + 2]) - (mtx[(1 * cols) + 2] * mtx[(2 * cols) + 1]); //set 0x0

    //1x0                  1x2                   2x0                     1x0                   2x2
    ret[(1 * cols) + 0] = (mtx[(1 * cols) + 2] * mtx[(2 * cols) + 0]) - (mtx[(1 * cols) + 0] * mtx[(2 * cols) + 2]); //set 1x0

    //2x0                  1x0                   2x1                     1x1                   2x0
    ret[(2 * cols) + 0] = (mtx[(1 * cols) + 0] * mtx[(2 * cols) + 1]) - (mtx[(1 * cols) + 1] * mtx[(2 * cols) + 0]); //set 2x0

    //COL 1
    //0x1                  0x2                   2x1                     0x1                   2x2
    ret[(0 * cols) + 1] = (mtx[(0 * cols) + 2] * mtx[(2 * cols) + 1]) - (mtx[(0 * cols) + 1] * mtx[(2 * cols) + 2]); //set 0x1

    //1x1                  0x0                   2x2                     0x2                   2x0
    ret[(1 * cols) + 1] = (mtx[(0 * cols) + 0] * mtx[(2 * cols) + 2]) - (mtx[(0 * cols) + 2] * mtx[(2 * cols) + 0]); //set 1x1

    //2x1                  0x1                   2x0                     0x0                   2x1
    ret[(2 * cols) + 1] = (mtx[(0 * cols) + 1] * mtx[(2 * cols) + 0]) - (mtx[(0 * cols) + 0] * mtx[(2 * cols) + 1]); //set 2x1

    //COL 2
    //0x2                  0x1                   1x2                     0x2                   1x1
    ret[(0 * cols) + 2] = (mtx[(0 * cols) + 1] * mtx[(1 * cols) + 2]) - (mtx[(0 * cols) + 2] * mtx[(1 * cols) + 1]); //set 0x2

    //1x2                  0x2                   1x0                     0x0                   1x2
    ret[(1 * cols) + 2] = (mtx[(0 * cols) + 2] * mtx[(1 * cols) + 0]) - (mtx[(0 * cols) + 0] * mtx[(1 * cols) + 2]); //set 1x2

    //2x2                  0x0                   1x1                     0x1                   1x0
    ret[(2 * cols) + 2] = (mtx[(0 * cols) + 0] * mtx[(1 * cols) + 1]) - (mtx[(0 * cols) + 1] * mtx[(1 * cols) + 0]); //set 2x2

    mulXvec(ret, (1.0 / det));
}

test "XMTX: getInvFromDet3 test" {
    //A = 1 2 3
    //    4 5 6
    //    7 2 9

    //detA = (1)(33)          - (2)(-6)          + (3)(-27)
    //detA = (1)(det 5 6 2 9) - (2)(det 4 6 7 9) + (3)(det 4 5 7 2)
    //detA = -36;

    //A^-1 = -11/12  1/3  1/12
    //       -1/6    1/3  -1/6
    //       3/4     -1/3 1/12

    var A: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 2, 9 };
    var res: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const detA: f32 = -36;
    var invA: [9]f32 = .{ (-11.0 / 12.0), (1.0 / 3.0), (1.0 / 12.0), (-1.0 / 6.0), (1.0 / 3.0), (-1.0 / 6.0), (3.0 / 4.0), (-1.0 / 3.0), (1.0 / 12.0) };

    const start = try Instant.now();
    getInvFromDet3(&A, detA, &res);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ngetInvFromDet3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getInvFromDet3", elapsed1);

    try std.testing.expectEqual(true, equXmtx(&res, &invA));
    prntNl();
}

///Calculates the inverse of the provided 4x4 matrix using the determinant and stores the result in the return matrix.
///
///  mtx = The 4x4 matrix to calculate the inverse for.
///
///  det = The determinant of the provided matrix.
///
///  ret = The return 4x4 matrix to hold the inverted matrix values.
///
///  returns = A Boolean value indicating if the operation was a success.
///
pub fn getInvFromDet4(mtx: *[16]f32, det: f32, ret: *[16]f32) void {
    const cols: usize = 4;
    //ROW 0
    //COL 0
    //0x0
    ret[(0 * cols) + 0] = cofXmtxSign(0, 0, true) * (
    //1x1                    2x2                   3x3
        (mtx[(1 * cols) + 1] * mtx[(2 * cols) + 2] * mtx[(3 * cols) + 3])
    //1x2                    2x3                   3x1
    + (mtx[(1 * cols) + 2] * mtx[(2 * cols) + 3] * mtx[(3 * cols) + 1])
    //1x3                    2x1                   3x2
    + (mtx[(1 * cols) + 3] * mtx[(2 * cols) + 1] * mtx[(3 * cols) + 2])
    //1x3                    2x2                   3x1
    - (mtx[(1 * cols) + 3] * mtx[(2 * cols) + 2] * mtx[(3 * cols) + 1])
    //1x2                    2x1                   3x3
    - (mtx[(1 * cols) + 2] * mtx[(2 * cols) + 1] * mtx[(3 * cols) + 3])
    //1x1                    2x3                   3x2
    - (mtx[(1 * cols) + 1] * mtx[(2 * cols) + 3] * mtx[(3 * cols) + 2]));

    //COL 1
    //0x1
    ret[(0 * cols) + 1] = cofXmtxSign(0, 1, true) * (
    //1x0                    2x2                   3x3
        (mtx[(1 * cols) + 0] * mtx[(2 * cols) + 2] * mtx[(3 * cols) + 3])
    //1x2                    2x3                   3x0
    + (mtx[(1 * cols) + 2] * mtx[(2 * cols) + 3] * mtx[(3 * cols) + 0])
    //1x3                    2x0                   3x2
    + (mtx[(1 * cols) + 3] * mtx[(2 * cols) + 0] * mtx[(3 * cols) + 2])
    //1x3                    2x2                   3x0
    - (mtx[(1 * cols) + 3] * mtx[(2 * cols) + 2] * mtx[(3 * cols) + 0])
    //1x2                    2x0                   3x3
    - (mtx[(1 * cols) + 2] * mtx[(2 * cols) + 0] * mtx[(3 * cols) + 3])
    //1x0                    2x3                   3x2
    - (mtx[(1 * cols) + 0] * mtx[(2 * cols) + 3] * mtx[(3 * cols) + 2]));

    //COL 2
    //0x2
    ret[(0 * cols) + 2] = cofXmtxSign(0, 2, true) * (
    //1x0                    2x1                   3x3
        (mtx[(1 * cols) + 0] * mtx[(2 * cols) + 1] * mtx[(3 * cols) + 3])
    //1x1                    2x3                   3x0
    + (mtx[(1 * cols) + 1] * mtx[(2 * cols) + 3] * mtx[(3 * cols) + 0])
    //1x3                    2x0                   3x1
    + (mtx[(1 * cols) + 3] * mtx[(2 * cols) + 0] * mtx[(3 * cols) + 1])
    //1x3                    2x1                   3x0
    - (mtx[(1 * cols) + 3] * mtx[(2 * cols) + 1] * mtx[(3 * cols) + 0])
    //1x1                    2x0                   3x3
    - (mtx[(1 * cols) + 1] * mtx[(2 * cols) + 0] * mtx[(3 * cols) + 3])
    //1x0                    2x3                   3x1
    - (mtx[(1 * cols) + 0] * mtx[(2 * cols) + 3] * mtx[(3 * cols) + 1]));

    //COL 3
    //0x3
    ret[(0 * cols) + 3] = cofXmtxSign(0, 3, true) * (
    //1x0                    2x1                   3x2
        (mtx[(1 * cols) + 0] * mtx[(2 * cols) + 1] * mtx[(3 * cols) + 2])
    //1x1                    2x2                   3x0
    + (mtx[(1 * cols) + 1] * mtx[(2 * cols) + 2] * mtx[(3 * cols) + 0])
    //1x2                    2x0                   3x1
    + (mtx[(1 * cols) + 2] * mtx[(2 * cols) + 0] * mtx[(3 * cols) + 1])
    //1x2                    2x1                   3x0
    - (mtx[(1 * cols) + 2] * mtx[(2 * cols) + 1] * mtx[(3 * cols) + 0])
    //1x1                    2x0                   3x2
    - (mtx[(1 * cols) + 1] * mtx[(2 * cols) + 0] * mtx[(3 * cols) + 2])
    //1x0                    2x2                   3x1
    - (mtx[(1 * cols) + 0] * mtx[(2 * cols) + 2] * mtx[(3 * cols) + 1]));

    //ROW 1
    //COL 0
    //1x0
    ret[(1 * cols) + 0] = cofXmtxSign(1, 0, true) * (
    //0x1                    2x2                   3x3
        (mtx[(0 * cols) + 1] * mtx[(2 * cols) + 2] * mtx[(3 * cols) + 3])
    //0x2                    2x3                   3x1
    + (mtx[(0 * cols) + 2] * mtx[(2 * cols) + 3] * mtx[(3 * cols) + 1])
    //0x3                    2x1                   3x2
    + (mtx[(0 * cols) + 3] * mtx[(2 * cols) + 1] * mtx[(3 * cols) + 2])
    //0x3                    2x2                   3x1
    - (mtx[(0 * cols) + 3] * mtx[(2 * cols) + 2] * mtx[(3 * cols) + 1])
    //0x2                    2x1                   3x3
    - (mtx[(0 * cols) + 2] * mtx[(2 * cols) + 1] * mtx[(3 * cols) + 3])
    //0x1                    2x3                   3x2
    - (mtx[(0 * cols) + 1] * mtx[(2 * cols) + 3] * mtx[(3 * cols) + 2]));

    //COL 1
    //1x1
    ret[(1 * cols) + 1] = cofXmtxSign(1, 1, true) * (
    //0x0                    2x2                   3x3
        (mtx[(0 * cols) + 0] * mtx[(2 * cols) + 2] * mtx[(3 * cols) + 3])
    //0x2                    2x3                   3x0
    + (mtx[(0 * cols) + 2] * mtx[(2 * cols) + 3] * mtx[(3 * cols) + 0])
    //0x3                    2x0                   3x2
    + (mtx[(0 * cols) + 3] * mtx[(2 * cols) + 0] * mtx[(3 * cols) + 2])
    //0x3                    2x2                   3x0
    - (mtx[(0 * cols) + 3] * mtx[(2 * cols) + 2] * mtx[(3 * cols) + 0])
    //0x2                    2x0                   3x3
    - (mtx[(0 * cols) + 2] * mtx[(2 * cols) + 0] * mtx[(3 * cols) + 3])
    //0x0                    2x3                   3x2
    - (mtx[(0 * cols) + 0] * mtx[(2 * cols) + 3] * mtx[(3 * cols) + 2]));

    //COL 2
    //1x2
    ret[(1 * cols) + 2] = cofXmtxSign(1, 2, true) * (
    //0x0                    2x1                   3x3
        (mtx[(0 * cols) + 0] * mtx[(2 * cols) + 1] * mtx[(3 * cols) + 3])
    //0x1                    2x3                   3x0
    + (mtx[(0 * cols) + 1] * mtx[(2 * cols) + 3] * mtx[(3 * cols) + 0])
    //0x3                    2x0                   3x1
    + (mtx[(0 * cols) + 3] * mtx[(2 * cols) + 0] * mtx[(3 * cols) + 1])
    //0x3                    2x1                   3x0
    - (mtx[(0 * cols) + 3] * mtx[(2 * cols) + 1] * mtx[(3 * cols) + 0])
    //0x1                    2x0                   3x3
    - (mtx[(0 * cols) + 1] * mtx[(2 * cols) + 0] * mtx[(3 * cols) + 3])
    //0x0                    2x3                   3x1
    - (mtx[(0 * cols) + 0] * mtx[(2 * cols) + 3] * mtx[(3 * cols) + 1]));

    //COL 3
    //1x3
    ret[(1 * cols) + 3] = cofXmtxSign(1, 3, true) * (
    //0x0                    2x1                   3x2
        (mtx[(0 * cols) + 0] * mtx[(2 * cols) + 1] * mtx[(3 * cols) + 2])
    //0x1                    2x2                   3x0
    + (mtx[(0 * cols) + 1] * mtx[(2 * cols) + 2] * mtx[(3 * cols) + 0])
    //0x2                    2x0                   3x1
    + (mtx[(0 * cols) + 2] * mtx[(2 * cols) + 0] * mtx[(3 * cols) + 1])
    //0x2                    2x1                   3x0
    - (mtx[(0 * cols) + 2] * mtx[(2 * cols) + 1] * mtx[(3 * cols) + 0])
    //0x1                    2x0                   3x2
    - (mtx[(0 * cols) + 1] * mtx[(2 * cols) + 0] * mtx[(3 * cols) + 2])
    //0x0                    2x2                   3x1
    - (mtx[(0 * cols) + 0] * mtx[(2 * cols) + 2] * mtx[(3 * cols) + 1]));

    //ROW 2
    //COL 0
    //2x0
    ret[(2 * cols) + 0] = cofXmtxSign(2, 0, true) * (
    //0x1                    1x2                   3x3
        (mtx[(0 * cols) + 1] * mtx[(1 * cols) + 2] * mtx[(3 * cols) + 3])
    //0x2                    1x3                   3x1
    + (mtx[(0 * cols) + 2] * mtx[(1 * cols) + 3] * mtx[(3 * cols) + 1])
    //0x3                    1x1                   3x2
    + (mtx[(0 * cols) + 3] * mtx[(1 * cols) + 1] * mtx[(3 * cols) + 2])
    //0x3                    1x2                   3x1
    - (mtx[(0 * cols) + 3] * mtx[(1 * cols) + 2] * mtx[(3 * cols) + 1])
    //0x2                    1x1                   3x3
    - (mtx[(0 * cols) + 2] * mtx[(1 * cols) + 1] * mtx[(3 * cols) + 3])
    //0x1                    1x3                   3x2
    - (mtx[(0 * cols) + 1] * mtx[(1 * cols) + 3] * mtx[(3 * cols) + 2]));

    //COL 1
    //2x1
    ret[(2 * cols) + 1] = cofXmtxSign(2, 1, true) * (
    //0x0                    1x2                   3x3
        (mtx[(0 * cols) + 0] * mtx[(1 * cols) + 2] * mtx[(3 * cols) + 3])
    //0x2                    1x3                   3x0
    + (mtx[(0 * cols) + 2] * mtx[(1 * cols) + 3] * mtx[(3 * cols) + 0])
    //0x3                    1x0                   3x2
    + (mtx[(0 * cols) + 3] * mtx[(1 * cols) + 0] * mtx[(3 * cols) + 2])
    //0x3                    1x2                   3x0
    - (mtx[(0 * cols) + 3] * mtx[(1 * cols) + 2] * mtx[(3 * cols) + 0])
    //0x2                    1x0                   3x3
    - (mtx[(0 * cols) + 2] * mtx[(1 * cols) + 0] * mtx[(3 * cols) + 3])
    //0x0                    1x3                   3x2
    - (mtx[(0 * cols) + 0] * mtx[(1 * cols) + 3] * mtx[(3 * cols) + 2]));

    //COL 2
    //2x2
    ret[(2 * cols) + 2] = cofXmtxSign(2, 2, true) * (
    //0x0                    1x1                   3x3
        (mtx[(0 * cols) + 0] * mtx[(1 * cols) + 1] * mtx[(3 * cols) + 3])
    //0x1                    1x3                   3x0
    + (mtx[(0 * cols) + 1] * mtx[(1 * cols) + 3] * mtx[(3 * cols) + 0])
    //0x3                    1x0                   3x1
    + (mtx[(0 * cols) + 3] * mtx[(1 * cols) + 0] * mtx[(3 * cols) + 1])
    //0x3                    1x1                   3x0
    - (mtx[(0 * cols) + 3] * mtx[(1 * cols) + 1] * mtx[(3 * cols) + 0])
    //0x1                    1x0                   3x3
    - (mtx[(0 * cols) + 1] * mtx[(1 * cols) + 0] * mtx[(3 * cols) + 3])
    //0x0                    1x3                   3x1
    - (mtx[(0 * cols) + 0] * mtx[(1 * cols) + 3] * mtx[(3 * cols) + 1]));

    //COL 3
    //2x3
    ret[(2 * cols) + 3] = cofXmtxSign(2, 3, true) * (
    //0x0                    1x1                   3x2
        (mtx[(0 * cols) + 0] * mtx[(1 * cols) + 1] * mtx[(3 * cols) + 2])
    //0x1                    1x2                   3x0
    + (mtx[(0 * cols) + 1] * mtx[(1 * cols) + 2] * mtx[(3 * cols) + 0])
    //0x2                    1x0                   3x1
    + (mtx[(0 * cols) + 2] * mtx[(1 * cols) + 0] * mtx[(3 * cols) + 1])
    //0x2                    1x1                   3x0
    - (mtx[(0 * cols) + 2] * mtx[(1 * cols) + 1] * mtx[(3 * cols) + 0])
    //0x1                    1x0                   3x2
    - (mtx[(0 * cols) + 1] * mtx[(1 * cols) + 0] * mtx[(3 * cols) + 2])
    //0x0                    1x2                   3x1
    - (mtx[(0 * cols) + 0] * mtx[(1 * cols) + 2] * mtx[(3 * cols) + 1]));

    //ROW 3
    //COL 0
    //3x0
    ret[(3 * cols) + 0] = cofXmtxSign(3, 0, true) * (
    // 0x1                    1x2                   2x3
        (mtx[(0 * cols) + 1] * mtx[(1 * cols) + 2] * mtx[(2 * cols) + 3])
    //0x2                    1x3                   2x1
    + (mtx[(0 * cols) + 2] * mtx[(1 * cols) + 3] * mtx[(2 * cols) + 1])
    //0x3                    1x1                   2x2
    + (mtx[(0 * cols) + 3] * mtx[(1 * cols) + 1] * mtx[(2 * cols) + 2])
    //0x3                    1x2                   2x1
    - (mtx[(0 * cols) + 3] * mtx[(1 * cols) + 2] * mtx[(2 * cols) + 1])
    //0x2                    1x1                   2x3
    - (mtx[(0 * cols) + 2] * mtx[(1 * cols) + 1] * mtx[(2 * cols) + 3])
    //0x1                    1x3                   2x2
    - (mtx[(0 * cols) + 1] * mtx[(1 * cols) + 3] * mtx[(2 * cols) + 2]));

    //COL 1
    //3x1
    ret[(3 * cols) + 1] = cofXmtxSign(3, 1, true) * (
    //0x0                    1x2                   2x3
        (mtx[(0 * cols) + 0] * mtx[(1 * cols) + 2] * mtx[(2 * cols) + 3])
    //0x2                    1x3                   2x0
    + (mtx[(0 * cols) + 2] * mtx[(1 * cols) + 3] * mtx[(2 * cols) + 0])
    //0x3                    1x0                   2x2
    + (mtx[(0 * cols) + 3] * mtx[(1 * cols) + 0] * mtx[(2 * cols) + 2])
    //0x3                    1x2                   2x0
    - (mtx[(0 * cols) + 3] * mtx[(1 * cols) + 2] * mtx[(2 * cols) + 0])
    //0x2                    1x0                   2x3
    - (mtx[(0 * cols) + 2] * mtx[(1 * cols) + 0] * mtx[(2 * cols) + 3])
    //0x0                    1x3                   2x2
    - (mtx[(0 * cols) + 0] * mtx[(1 * cols) + 3] * mtx[(2 * cols) + 2]));

    //COL 2
    //3x2
    ret[(3 * cols) + 2] = cofXmtxSign(3, 2, true) * (
    //0x0                    1x1                   2x3
        (mtx[(0 * cols) + 0] * mtx[(1 * cols) + 1] * mtx[(2 * cols) + 3])
    //0x1                    1x3                   2x0
    + (mtx[(0 * cols) + 1] * mtx[(1 * cols) + 3] * mtx[(2 * cols) + 0])
    //0x3                    1x0                   2x1
    + (mtx[(0 * cols) + 3] * mtx[(1 * cols) + 0] * mtx[(2 * cols) + 1])
    //0x3                    1x1                   2x0
    - (mtx[(0 * cols) + 3] * mtx[(1 * cols) + 1] * mtx[(2 * cols) + 0])
    //0x1                    1x0                   2x3
    - (mtx[(0 * cols) + 1] * mtx[(1 * cols) + 0] * mtx[(2 * cols) + 3])
    //0x0                    1x3                   2x1
    - (mtx[(0 * cols) + 0] * mtx[(1 * cols) + 3] * mtx[(2 * cols) + 1]));

    //COL 3
    //3x3
    ret[(3 * cols) + 3] = cofXmtxSign(3, 3, true) * (
    //0x0                    1x1                   2x2
        (mtx[(0 * cols) + 0] * mtx[(1 * cols) + 1] * mtx[(2 * cols) + 2])
    //0x1                    1x2                   2x0
    + (mtx[(0 * cols) + 1] * mtx[(1 * cols) + 2] * mtx[(2 * cols) + 0])
    //0x2                    1x0                   2x1
    + (mtx[(0 * cols) + 2] * mtx[(1 * cols) + 0] * mtx[(2 * cols) + 1])
    //0x2                    1x1                   2x0
    - (mtx[(0 * cols) + 2] * mtx[(1 * cols) + 1] * mtx[(2 * cols) + 0])
    //0x1                    1x0                   2x2
    - (mtx[(0 * cols) + 1] * mtx[(1 * cols) + 0] * mtx[(2 * cols) + 2])
    //0x0                    1x2                   2x1
    - (mtx[(0 * cols) + 0] * mtx[(1 * cols) + 2] * mtx[(2 * cols) + 1]));

    if (VERBOSE) {
        prntNl();
        prntXmtxNl(ret, 4);
        prntNl();
    }

    mulXvec(ret, (1.0 / det));
}

test "XMTX: getInvFromDet4 test" {
    var ret: [16]f32 = std.mem.zeroes([16]f32); //.{};
    var mtx: [16]f32 = .{ 1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1 };
    var expMtx: [16]f32 = .{ (1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0), (-1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0), (-1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0), (-1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0), (-1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0), (1.0 / 4.0) };
    const det: f32 = detXmtx4(&mtx);
    const expDet: f32 = -16;

    const start = try Instant.now();
    getInvFromDet4(&mtx, det, &ret);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ngetInvFromDet4: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getInvFromDet4", elapsed1);

    const cols: usize = 4;
    std.debug.print("\ngetInvFromDet4 det: {} expDet: {}", .{ det, expDet });
    prntXmtxNl(&ret, cols);
    prntNl();

    std.debug.print("\ngetInvFromDet4 det:", .{});
    prntXmtxNl(&expMtx, cols);
    prntNl();

    try std.testing.expectEqual(true, isEquF32(expDet, det, true));
    try std.testing.expectEqual(true, equXmtx(&expMtx, &ret));
    prntNl();

    //TODO: add more tests
}

///Creates a Cramer's Rule support matrix based on the provided augmented matrix, mtx, and other arguments.
///
///  mtx = The augmented matrix to create a support matrix for.
///
///  cols = The number of columns in the mtx matrix.
///
///  srcCol = The last column in the mtx matrix that contains constant values.
///
///  dstCol = The target column in the return matrix, ret, to store data.
///
///  ret = The return matrix that has it's dstCol, destination column, replaced with the mtx matrix's constant values.
///
///  returns = A Boolean value indicating if the operation was a success or an error.
///
pub fn getCramerSupportXmtx(mtx: []f32, cols: usize, srcCol: usize, dstCol: usize, ret: []f32) !bool {
    const l: usize = mtx.len;
    const rows: usize = l / cols;
    if (rows >= cols) {
        prntNlStr("!! Warning getCramerSupportMtx requires the augmented matrix, mtx, should have one more column than rows. !!");
        return Error.AugmenedMatrixRequired;
    }

    if (srcCol != cols - 1) {
        prntNlStr("!! Warning getCramerSupportMtx requires augmented matrix, mtx, specify constants, cols - 1, to create the necessary support matrix. !!");
        return Error.AugmenedMatrixRequired;
    }

    if (dstCol < 0 or dstCol >= cols - 1) {
        prntNlStr("!! Warning getCramerSupportMtx requires return matrix, ret, should have column count, (cols - 1). !!");
        return Error.InvalidLengths;
    }

    var i: usize = 0;
    var j: usize = 0;
    const dstCols: usize = (cols - 1);
    const dstRows: usize = ret.len / dstCols;
    //std.debug.print("\nRows: {} Cols: {} DstRows: {} DstCols: {}", .{rows, cols, dstRows, dstCols});

    while (i < dstCols) : (i += 1) {
        j = 0;
        while (j < dstRows) : (j += 1) {
            if (i == dstCol) {
                //std.debug.print("\nSetting Idx: {} From Idx: {}", .{((j * dstCols) + i), ((j * cols) + srcCol)});
                ret[(j * dstCols) + i] = mtx[(j * cols) + srcCol];
            } else {
                //std.debug.print("\nSetting Idx: {} From Idx: {}", .{((j * dstCols) + i), ((j * cols) + i)});
                ret[(j * dstCols) + i] = mtx[(j * cols) + i];
            }
        }
    }

    return true;
}

test "XMTX: getCramerSupportMtx test" {
    //A = -1  2 -3  1
    //     2  0  1  0
    //     3 -4  4  2
    var A: [12]f32 = .{ -1, 2, -3, 1, 2, 0, 1, 0, 3, -4, 4, 2 };
    var A1: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var b: bool = false;
    const cols: usize = 4;
    const srcCol: usize = 3;
    var dstCol: usize = 0;
    var expA1: [9]f32 = .{ 1, 2, -3, 0, 0, 1, 2, -4, 4 };

    var start = try Instant.now();
    b = try getCramerSupportXmtx(&A, cols, srcCol, dstCol, &A1);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ngetCramerSupportMtx #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getCramerSupportMtx", elapsed1);

    std.debug.print("\nMatrix A{}:", .{dstCol});
    prntXmtxNl(&A1, 3);
    std.debug.print("\nMatrix Exp A{}:", .{dstCol});
    prntXmtxNl(&expA1, 3);
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&A1, &expA1));
    dstCol = 1;
    var expA2: [9]f32 = .{ -1, 1, -3, 2, 0, 1, 3, 2, 4 };

    start = try Instant.now();
    b = try getCramerSupportXmtx(&A, cols, srcCol, dstCol, &A1);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ngetCramerSupportMtx #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getCramerSupportMtx", elapsed1);

    std.debug.print("\nMatrix A{}:", .{dstCol});
    prntXmtxNl(&A1, 3);
    std.debug.print("\nMatrix Exp A{}:", .{dstCol});
    prntXmtxNl(&expA1, 3);
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&A1, &expA2));
    dstCol = 2;
    var expA3: [9]f32 = .{ -1, 2, 1, 2, 0, 0, 3, -4, 2 };

    start = try Instant.now();
    b = try getCramerSupportXmtx(&A, cols, srcCol, dstCol, &A1);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ngetCramerSupportMtx #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getCramerSupportMtx", elapsed1);

    std.debug.print("\nMatrix A{}:", .{dstCol});
    prntXmtxNl(&A1, 3);
    std.debug.print("\nMatrix Exp A{}:", .{dstCol});
    prntXmtxNl(&expA1, 3);
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&A1, &expA3));
    prntNl();
}

///Applies Cramer's rule to the given matrix, mtx, with mtxA represented the unaugmented version of matrix mtx, and mtxAi representing Cramer's rule's supporting matrix.
///
///  mtx = The augmented matrix to apply Cramer's rule to.
///
///  cols = The number of columns in the mtx matrix.
///
///  mtxA = The matrix representing the unaugmented version of mtx.
///
///  colsA = The number of columns in the mtxA matrix.
///
///  mtxAi = The supporting matrix used by Cramer's rule.
///
///  ret = An array of f32 values used to store the results of Cramer's rule.
///
///  returns = A Boolean value indicating if the operation was a success.
///
pub fn rslvCramersRule(mtx: []f32, cols: usize, mtxA: []f32, colsA: usize, mtxAi: []f32, ret: []f32) !bool {
    const l: usize = mtx.len;
    const rows: usize = l / cols;
    if (rows >= cols) {
        prntNlStr("\n!! Warning rslvCramersRule requires the augmented matrix, mtx, should have one more column than rows. !!");
        return Error.InvalidLengths;
    }

    const lA: usize = mtxA.len;
    const rowsA: usize = lA / colsA;
    if (rowsA != (cols - 1)) {
        prntNlStr("!! Warning getCramerSupportMtx requires augmented matrix, mtx, specify constants, cols - 1, to create the necessary support matrix. !!");
        return Error.AugmenedMatrixRequired;
    }

    var detAi: f32 = 0;
    var detA: f32 = 0;
    if (rowsA == 1) {
        detA = detXmtx1(@as(*[1]f32, @ptrCast(mtxA)));
    } else if (rowsA == 2) {
        detA = detXmtx2(@as(*[4]f32, @ptrCast(mtxA)));
    } else if (rowsA == 3) {
        detA = detXmtx3(@as(*[9]f32, @ptrCast(mtxA)));
    } else if (rowsA == 4) {
        detA = detXmtx4(@as(*[16]f32, @ptrCast(mtxA)));
    } else {
        prntNlStr("!! Warning rslvCramersRule this function only supports matrix dimensions of 1, 2, 3, 4 for argument mtxA. !!");
        return Error.UnsupportedMarixDimension;
    }

    //std.debug.print("\ndetA: {}", .{detA});

    var b: bool = false;
    var i: usize = 0;
    while (i < colsA) : (i += 1) {
        //create cramer support matrix Ai
        b = try getCramerSupportXmtx(mtx, cols, (cols - 1), i, mtxAi);
        if (!b) {
            prntNlStrArgs("!! Warning rslvCramersRule could not populate the Cramer's rule support matrix for column {}. !!", .{i});
            return Error.OperationFailed;
        }

        if (rowsA == 1) {
            detAi = detXmtx1(@as(*[1]f32, @ptrCast(mtxAi)));
        } else if (rowsA == 2) {
            detAi = detXmtx2(@as(*[4]f32, @ptrCast(mtxAi)));
        } else if (rowsA == 3) {
            detAi = detXmtx3(@as(*[9]f32, @ptrCast(mtxAi)));
        } else if (rowsA == 4) {
            detAi = detXmtx4(@as(*[16]f32, @ptrCast(mtxAi)));
        } else {
            prntNlStr("!! Warning rslvCramersRule This function only supports matrix dimensions of 1, 2, 3, 4 for Cramer's rule support matrix, mtxAi. !!");
            return Error.UnsupportedMarixDimension;
        }

        //std.debug.print("\nDetAI:{} DetA:{} RowsA:{}", .{detAi, detA, rowsA});
        ret[i] = detAi / detA;
    }

    return true;
}

test "XMTX: rslvCramersRule test" {
    //A = -1  2 -3  1
    //     2  0  1  0
    //     3 -4  4  2
    var mtx: [12]f32 = .{ -1, 2, -3, 1, 2, 0, 1, 0, 3, -4, 4, 2 };
    var A: [9]f32 = .{ -1, 2, -3, 2, 0, 1, 3, -4, 4 };
    var Ai: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var ret: [3]f32 = .{ 0, 0, 0 };
    var b: bool = false;
    const cols: usize = 4;
    const colsA: usize = 3;

    const start = try Instant.now();
    b = try rslvCramersRule(&mtx, cols, &A, colsA, &Ai, &ret);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nrslvCramersRule: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("rslvCramersRule", elapsed1);

    try std.testing.expectEqual(true, b);
    const expX: f32 = 4.0 / 5.0;
    const expY: f32 = -3.0 / 2.0;
    const expZ: f32 = -8.0 / 5.0;

    std.debug.print("\nCramer's Rule Return {} {} {}:", .{ expX, expY, expZ });
    prntXvecNl(&ret);
    prntNl();

    try std.testing.expectEqual(true, isEquF32(ret[0], expX, true));
    try std.testing.expectEqual(true, isEquF32(ret[1], expY, true));
    try std.testing.expectEqual(true, isEquF32(ret[2], expZ, true));
    prntNl();
}

///This function is responsible for changing the basis of the provided vector, vec, from the basis, basis, to chgBasis.
///
///  vec = The vector to convert from one basis to another.
///
///  basis = The matrix that constitutes the current basis for vector vec, and holds the resulting transformation matrix.
///
///  cols = The number of columns in the matrix, basis.
///
///  isStd = A Boolean value indicating if the basis is considered the standard basis.
///
///  chgBasis = The mtrix representing the basis to change to.
///
///  chgCols = The number of columns in the chgCols matrix.
///
///  chgIsStd = A Boolean value indicating if the change basis is considered the standard basis.
///
///  idtMtx = The identity associated with the basis and chgBasis matrices.
///
///  ret = The result, left-hand side, of the tranformation matrix calculation, should result in an identity matrix.
///
///  verbose = A Boolean indicating if the function is verbose.
///
///  returns = A Boolean value indicating if the operation was a success or not.
///
pub fn chgXvecBasis(vec: []f32, basis: []f32, cols: usize, isStd: bool, chgBasis: []f32, chgCols: usize, chgIsStd: bool, idtMtx: []f32, ret: []f32, nvec: []f32, verbose: bool) !bool {
    _ = idtMtx;
    _ = isStd;

    if (cols != chgCols) {
        prntNlStr("\n!! Warning expected column counts to match !!");
        return Error.InvalidLengths;
    }

    var b: bool = false;
    if (chgIsStd) {
        b = true;
        if (verbose) {
            prntNlStrArgs("chgVecBasis:change basis is identity so conversion matrix = basis {}", .{b});
        }
    } else {
        b = try getBasisCnvXmtx(basis, cols, chgBasis, chgCols, ret, verbose);
        if (verbose) {
            prntNlStrArgs("chgVecBasis:getBasisCnvMtx {}", .{b});
        }
    }

    if (verbose) {
        prntXmtxNl(basis, cols);
        prntNl();

        prntXmtxNl(ret, cols);
        prntNl();

        prntXmtxNl(vec, 1);
        prntNl();
    }

    if (!b) { // or equXmtx(ret, idtMtx) == false) {}
        prntNlStr("!! Warning matrix reduction failed !!");
        return Error.OperationFailed;
    }

    b = try tmsXmtx(basis, cols, vec, 1, nvec, 1);
    if (!b) {
        prntNlStr("!! Warning matrix multiplication failed !!");
        return Error.OperationFailed;
    }
    return true;
}

test "XMTX: chgXvecBasis test" {
    //B = {(1,0), (1,2)}
    //Bp = {(1,0), (0,1)}

    //B = | 1  1|
    //    | 0  2|
    //vec = [3 2]

    //Bp = | 1  0|
    //     | 0  1|
    //nvec = [? ?]

    var B: [4]f32 = .{ 1, 1, 0, 2 };
    var vec: [2]f32 = .{ 3, 2 };
    var exp: [2]f32 = .{ 5, 4 };
    var nvec: [2]f32 = .{ 0, 0 };

    const cols: usize = 2;
    var Bp: [4]f32 = .{ 1, 0, 0, 1 };
    const colsp: usize = 2;
    var idtMtx: [4]f32 = .{ 1, 0, 0, 1 };

    var b: bool = false;
    const vbose: bool = true;
    var ret: [4]f32 = .{ 0, 0, 0, 0 };

    const start = try Instant.now();
    b = try chgXvecBasis(&vec, &B, cols, false, &Bp, colsp, true, &idtMtx, &ret, &nvec, vbose);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nchgXvecBasis: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("chgXvecBasis", elapsed1);

    try std.testing.expectEqual(true, b);
    std.debug.print("\nchgXvecBasis ret {}", .{b});
    prntXmtxNl(&ret, cols);
    prntNl();

    std.debug.print("\nchgXvecBasis B", .{});
    prntXmtxNl(&B, cols);
    prntNl();

    std.debug.print("\nchgXvecBasis nvec", .{});
    prntXmtxNl(&nvec, 1);
    prntNl();

    try std.testing.expectEqual(true, equXmtx(&exp, &nvec));
    prntNl();
}

///Returns the basis conversion, tranformation, matrix for the given current basis and change of basis.
///
///  basis = The matrix that constitutes the current basis and holds the resulting transformation matrix.
///
///  cols = The number of columns in the basis matrix.
///
///  chgBasis = The matrix that constitues the new basis to find coordinates in.
///
///  ret = The result, left-hand side, of the tranformation matrix calculation, should result in an identity matrix.
///
///  verbose = A Boolean indicating if the function is verbose.
///
///  returns = A Boolean value indicating if this function was successful or not.
///
pub fn getBasisCnvXmtx(basis: []f32, cols: usize, chgBasis: []f32, chgCols: usize, ret: []f32, verbose: bool) !bool {
    if (cols != chgCols) {
        prntNlStr("!! Warning expected column counts to match !!");
        return Error.InvalidLengths;
    }

    //Reduce the provided matrix to reduced row escelon form using Gauss-Jordan Elimination and optionaly calculate the matrix inverse.
    //  mtx = The matrix to reduce.
    //  cols = The number of columns in the matrix.
    //  hasAug = A Boolean indicating if the matrix is an augmented matrix.
    //  ret = The matrix that holds the reduced matrix.
    //  hasIdtMtx = A Boolean indicating if the identity matrix has been provided.
    //  idtMtx = The identity matrix associated with the mtx matrix provided.
    //  dim = The number of matrix columns that must be zero for a zero row to exist.
    //  triagRdcOnly = A Boolean value indicating if the reduction operation should stop when the matrix is triangular.
    //  sclr = A pointer to a floating point variable that keeps track of the scalar multiplication performed against the matrix, mtx.
    //  returns = A Boolean value indicating  if the matrix was reduced successfuly.
    var b: bool = false;
    const hasAug: bool = false;
    const hasIdtMtx: bool = true;
    const dim: usize = cols;
    var sclr: f32 = 0.0;
    const triagRdcOnly: bool = false;
    b = try rdcXmtx(chgBasis, chgCols, hasAug, ret, hasIdtMtx, basis, dim, triagRdcOnly, &sclr);

    if (verbose) {
        prntNlStrArgs("getBasisCnvMtx:rdcXmtx {}", .{b});
        prntXmtxNl(ret, cols);
        prntNl();
    }

    if (!b) {
        prntNlStr("!! Warning could not reduce the transfomation matrix !!");
        return Error.OperationFailed;
    }

    return true;
}

test "XMTX: getBasisCnvXmtx test" {
    //B = {(1,0,0), (0, 1, 0), (0, 0, 1)};
    //B' = {(1,0,1), (0, -1, 2), (2, 3, -5)};

    //B = | 1  0  0|
    //    | 0  1  0|
    //    | 0  0  1|
    const vbose: bool = false;
    var ret: [9]f32 = std.mem.zeroes([9]f32); //.{};
    var B: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    const cols: usize = 3;

    //B' = | 1  0  2|
    //     | 0 -1  3|
    //     | 1  2 -5|
    var Bp: [9]f32 = .{ 1, 0, 2, 0, -1, 3, 1, 2, -5 };
    const colsp: usize = 3;
    var b: bool = false;

    const start = try Instant.now();
    b = try getBasisCnvXmtx(&B, cols, &Bp, colsp, &ret, vbose);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ngetBasisCnvXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getBasisCnvXmtx", elapsed1);

    try std.testing.expectEqual(true, b);
    std.debug.print("\ngetBasisCnvMtx test:getBasisCnvMtx ret {}", .{b});
    prntXmtxNl(&ret, cols);
    prntNl();

    std.debug.print("\ngetBasisCnvMtx test:getBasisCnvMtx B", .{});
    prntXmtxNl(&B, cols);
    prntNl();

    var exp: [9]f32 = .{ -1, 4, 2, 3, -7, -3, 1, -2, -1 };
    try std.testing.expectEqual(true, equXmtx(&exp, &B));
    prntNl();
}

///Returns a Boolean value indicating if the mtx argument is an orthogonal matrix. Uses matrix reduction to determine
/// if the matrix, mtx, is orthogonal.
///
///  mtx = The matrix to determine othogonality for.
///
///  cols = The number of columns in the matrix mtx.
///
///  ret = An empty matrix the same size as mtx used to hold return information from the reduce matrix call.
///
///  idtMtx = An identity matrix with the same dimensions as the mtx matrix.
///
///  trnMtx = An empty matrix the same size as mtx used to hold transpose matrix information.
///
///  return = A Boolean value indicating if mtx is orthogonal or not.
///
pub fn isOrthXmtx(mtx: []f32, cols: usize, ret: []f32, idtMtx: []f32, trnMtx: []f32) !bool {
    //get the inverse

    //Reduce the provided matrix to reduced row escelon form using Gauss-Jordan Elimination and optionaly calculate the matrix inverse.
    //  mtx = The matrix to reduce.
    //  cols = The number of columns in the matrix.
    //  hasAug = A Boolean indicating if the matrix is an augmented matrix.
    //  ret = The matrix that holds the reduced matrix.
    //  hasIdtMtx = A Boolean indicating if the identity matrix has been provided.
    //  idtMtx = The identity matrix associated with the mtx matrix provided.
    //  dim = The number of matrix columns that must be zero for a zero row to exist.
    //  triagRdcOnly = A Boolean value indicating if the reduction operation should stop when the matrix is triangular.
    //  sclr = A pointer to a floating point variable that keeps track of the scalar multiplication performed against the matrix, mtx.
    //  returns = A Boolean value indicating  if the matrix was reduced successfuly.
    var b: bool = false;
    const hasAug: bool = false;
    const hasIdtMtx: bool = true;
    const dim: usize = cols;
    var sclr: f32 = 0.0;
    const triagRdcOnly: bool = false;
    b = try rdcXmtx(mtx, cols, hasAug, ret, hasIdtMtx, idtMtx, dim, triagRdcOnly, &sclr);

    //get the transpose
    trnXmtx(mtx, cols, trnMtx);

    //check equality
    b = equXmtx(idtMtx, trnMtx);
    return b;
}

test "XMTX: isOrthXmtx test" {
    var mtx: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    const cols: usize = 3;
    var ret: [9]f32 = std.mem.zeroes([9]f32); //.{};
    var idtMtx: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var trnMtx: [9]f32 = std.mem.zeroes([9]f32); //.{};
    const expB: bool = true;
    var b: bool = false;

    const start = try Instant.now();
    b = try isOrthXmtx(&mtx, cols, &ret, &idtMtx, &trnMtx);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisOrthXmtx: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isOrthXmtx", elapsed1);

    try std.testing.expectEqual(expB, b);
    prntNl();
}

///Returns an enumeration value indicating the handedness of the 3 vector arguments assuming they form a basis,
///from Lengyel: Def: 4.1.2. If a handled error is encountered the function returns ERROR_ZERO or ERROR_INVALID_MATRIX
///depending on the type of error.
///
///  vecI = The i vector 3 of the basis.
///
///  vecJ = The j vector 3 of the basis.
///
///  vecK = The k vector 3 of the basis.
///
///  returns = A value from the BASIS_HAND enumeration, RIGHT, LEFT, ERROR_ZERO or ERROR_INVALID_MATRIX.
///
pub fn getBasisHndXvec3(vecI: *const [3]f32, vecJ: *const [3]f32, vecK: *const [3]f32) BASIS_HAND {
    //note these functions should have const parameters soon and this will need to be adjusted
    const iCrossJ: [3]f32 = crsPrdXvec3(vecI, vecJ);
    const iCrossJdotK: f32 = dotPrdXvec3(&iCrossJ, vecK);
    if (iCrossJdotK > 0) {
        return BASIS_HAND.RIGHT;
    } else if (iCrossJdotK < 0) {
        return BASIS_HAND.LEFT;
    } else {
        return BASIS_HAND.ERROR_ZERO;
    }
}

test "XMTX: getBasisHndXvec3 test" {
    var v1: [3]f32 = .{ 1, 0, 0 };
    var v2: [3]f32 = .{ 0, 1, 0 };
    var v3: [3]f32 = .{ 0, 0, 1 };
    var expResHnd: BASIS_HAND = BASIS_HAND.RIGHT;
    var resHnd: BASIS_HAND = BASIS_HAND.ERROR_ZERO;

    var start = try Instant.now();
    resHnd = getBasisHndXvec3(&v1, &v2, &v3);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ngetBasisHndXvec3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getBasisHndXvec3", elapsed1);

    try std.testing.expectEqual(expResHnd, resHnd);
    v1 = .{ 1, 0, 0 };
    v2 = .{ 0, 1, 0 };
    v3 = .{ 0, 0, -1 };
    expResHnd = BASIS_HAND.LEFT;

    start = try Instant.now();
    resHnd = getBasisHndXvec3(&v1, &v2, &v3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ngetBasisHndXvec3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getBasisHndXvec3", elapsed1);

    try std.testing.expectEqual(expResHnd, resHnd);
    v1 = .{ 0, 0, 0 };
    v2 = .{ 0, 0, 0 };
    v3 = .{ 0, 0, -1 };
    expResHnd = BASIS_HAND.ERROR_ZERO;

    start = try Instant.now();
    resHnd = getBasisHndXvec3(&v1, &v2, &v3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ngetBasisHndXvec3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getBasisHndXvec3", elapsed1);

    try std.testing.expectEqual(expResHnd, resHnd);
    prntNl();
}

///Returns an enumeration value indicating the handedness of the 3 vector arguments assuming they form a basis,
///from Lengyel: Def: 4.1.2. If a handled error is encountered the function returns ERROR_ZERO or ERROR_INVALID_MATRIX
///depending on the type of error.
///
///  vecI = The i vector 3 of the basis.
///
///  vecJ = The j vector 3 of the basis.
///
///  vecK = The k vector 3 of the basis.
///
///  ret = Stores a value from the BASIS_HAND enumeration, RIGHT, LEFT, ERROR_ZERO or ERROR_INVALID_MATRIX.
///
pub fn getBasisHndXvec3Ret(vecI: *const [3]f32, vecJ: *const [3]f32, vecK: *const [3]f32, ret: *BASIS_HAND) void {
    ret.* = getBasisHndXvec3(vecI, vecJ, vecK);
}

test "XMTX: getBasisHndXvec3Ret test" {
    var v1: [3]f32 = .{ 1, 0, 0 };
    var v2: [3]f32 = .{ 0, 1, 0 };
    var v3: [3]f32 = .{ 0, 0, 1 };
    var expResHnd: BASIS_HAND = BASIS_HAND.RIGHT;
    var resHnd: BASIS_HAND = BASIS_HAND.ERROR_ZERO;

    var start = try Instant.now();
    getBasisHndXvec3Ret(&v1, &v2, &v3, &resHnd);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ngetBasisHndXvec3Ret: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getBasisHndXvec3Ret", elapsed1);

    try std.testing.expectEqual(expResHnd, resHnd);
    v1 = .{ 1, 0, 0 };
    v2 = .{ 0, 1, 0 };
    v3 = .{ 0, 0, -1 };
    expResHnd = BASIS_HAND.LEFT;
    resHnd = BASIS_HAND.ERROR_ZERO;

    start = try Instant.now();
    getBasisHndXvec3Ret(&v1, &v2, &v3, &resHnd);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ngetBasisHndXvec3Ret #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getBasisHndXvec3Ret", elapsed1);

    try std.testing.expectEqual(expResHnd, resHnd);
    v1 = .{ 0, 0, 0 };
    v2 = .{ 0, 0, 0 };
    v3 = .{ 0, 0, -1 };
    expResHnd = BASIS_HAND.ERROR_ZERO;
    resHnd = BASIS_HAND.ERROR_ZERO;

    start = try Instant.now();
    getBasisHndXvec3Ret(&v1, &v2, &v3, &resHnd);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ngetBasisHndXvec3Ret #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getBasisHndXvec3Ret", elapsed1);

    try std.testing.expectEqual(expResHnd, resHnd);
    prntNl();
}

///Returns an enumeration value indicating the handedness of the 3 vector arguments assuming they form a basis,
///from Lengyel: Def: 4.1.2. If a handled error is encountered the function returns ERROR_ZERO or ERROR_INVALID_MATRIX
///depending on the type of error.
///
///  mtx = The matrix to use as the basis.
///
///  cols = The number of columns in mtx.
///
///  return = A value from the BASIS_HAND enmeration indicating the result.
///
pub fn getBasisHndXmtx3(mtx: *const [9]f32, cols: usize) BASIS_HAND {
    if (cols != 3 or mtx.len != 9) {
        return BASIS_HAND.ERROR_INVALID_MATRIX;
    }
    var vecI: [3]f32 = .{ 0, 0, 0 };
    cpyLessColRowXmtx3(mtx, &vecI, 0, cols, 0, 1, cols, cols);

    var vecJ: [3]f32 = .{ 0, 0, 0 };
    cpyLessColRowXmtx3(mtx, &vecJ, 0, cols, 1, 2, cols, cols);

    var vecK: [3]f32 = .{ 0, 0, 0 };
    cpyLessColRowXmtx3(mtx, &vecK, 0, cols, 2, 3, cols, cols);
    return getBasisHndXvec3(&vecI, &vecJ, &vecK);
}

test "XMTX: getBasisHndXmtx3 test" {
    var m1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var expResHnd: BASIS_HAND = BASIS_HAND.RIGHT;
    var resHnd: BASIS_HAND = BASIS_HAND.ERROR_ZERO;

    var start = try Instant.now();
    resHnd = getBasisHndXmtx3(&m1, 3);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ngetBasisHndXmtx3 #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getBasisHndXmtx3", elapsed1);

    try std.testing.expectEqual(expResHnd, resHnd);
    m1 = .{ 1, 0, 0, 0, 1, 0, 0, 0, -1 };
    expResHnd = BASIS_HAND.LEFT;
    resHnd = BASIS_HAND.ERROR_ZERO;

    start = try Instant.now();
    resHnd = getBasisHndXmtx3(&m1, 3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ngetBasisHndXmtx3 #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getBasisHndXmtx3", elapsed1);

    try std.testing.expectEqual(expResHnd, resHnd);
    m1 = .{ 0, 0, 0, 0, 0, 0, 0, 0, -1 };
    expResHnd = BASIS_HAND.ERROR_ZERO;
    resHnd = BASIS_HAND.LEFT;

    start = try Instant.now();
    resHnd = getBasisHndXmtx3(&m1, 3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ngetBasisHndXmtx3 #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getBasisHndXmtx3", elapsed1);

    try std.testing.expectEqual(expResHnd, resHnd);
    prntNl();
}

///Returns an enumeration value indicating the handedness of the 3 vector arguments assuming they form a basis,
///from Lengyel: Def: 4.1.2. If a handled error is encountered the function returns ERROR_ZERO or ERROR_INVALID_MATRIX
///depending on the type of error.
///
///  mtx = The matrix to use as the basis.
///
///  cols = The number of columns in mtx.
///
///  return = Stores a value from the BASIS_HAND enmeration indicating the result.
///
pub fn getBasisHndXmtx3Ret(mtx: *const [9]f32, cols: usize, ret: *BASIS_HAND) void {
    ret.* = getBasisHndXmtx3(mtx, cols);
}

test "XMTX: getBasisHndXmtx3Ret test" {
    var m1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var expResHnd: BASIS_HAND = BASIS_HAND.RIGHT;
    var resHnd: BASIS_HAND = BASIS_HAND.ERROR_ZERO;

    var start = try Instant.now();
    getBasisHndXmtx3Ret(&m1, 3, &resHnd);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ngetBasisHndXmtx3Ret #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getBasisHndXmtx3Ret", elapsed1);

    try std.testing.expectEqual(expResHnd, resHnd);
    m1 = .{ 1, 0, 0, 0, 1, 0, 0, 0, -1 };
    expResHnd = BASIS_HAND.LEFT;
    resHnd = BASIS_HAND.ERROR_ZERO;

    start = try Instant.now();
    getBasisHndXmtx3Ret(&m1, 3, &resHnd);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ngetBasisHndXmtx3Ret #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getBasisHndXmtx3Ret", elapsed1);

    try std.testing.expectEqual(expResHnd, resHnd);
    m1 = .{ 0, 0, 0, 0, 0, 0, 0, 0, -1 };
    expResHnd = BASIS_HAND.ERROR_ZERO;
    resHnd = BASIS_HAND.LEFT;

    start = try Instant.now();
    getBasisHndXmtx3Ret(&m1, 3, &resHnd);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ngetBasisHndXmtx3Ret #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("getBasisHndXmtx3Ret", elapsed1);

    try std.testing.expectEqual(expResHnd, resHnd);
    prntNl();
}

///Returns a Boolean indicating if a 3x3 matrix is right handed or not.
///
///  mtx = The matrix to check handedness for.
///
///  cols = The numbers of columns in mtx.
///
///  returns = A Boolean indicating if mtx is right handed.
///
pub fn isRghtHandedXmtx3(mtx: *const [9]f32, cols: usize) bool {
    if (getBasisHndXmtx3(mtx, cols) == BASIS_HAND.RIGHT) {
        return true;
    }
    return false;
}

test "XMTX: isRghtHandedXmtx3 test" {
    var m1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var expResHnd: bool = true;
    var resHnd: bool = false;

    var start = try Instant.now();
    resHnd = isRghtHandedXmtx3(&m1, 3);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisRghtHandedXmtx3 #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isRghtHandedXmtx3", elapsed1);

    try std.testing.expectEqual(expResHnd, resHnd);
    m1 = .{ 1, 0, 0, 0, 1, 0, 0, 0, -1 };
    expResHnd = false;
    resHnd = false;

    start = try Instant.now();
    resHnd = isRghtHandedXmtx3(&m1, 3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisRghtHandedXmtx3 #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isRghtHandedXmtx3", elapsed1);

    try std.testing.expectEqual(expResHnd, resHnd);
    m1 = .{ 0, 0, 0, 0, 0, 0, 0, 0, -1 };
    expResHnd = false;
    resHnd = false;

    start = try Instant.now();
    resHnd = isRghtHandedXmtx3(&m1, 3);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nisRghtHandedXmtx3 #3: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isRghtHandedXmtx3", elapsed1);

    try std.testing.expectEqual(expResHnd, resHnd);
    prntNl();
}

///Compares two function names, []const u8, and returns a number indicating their alphabetical relationship, -1 if a < b, 0 if a = b, 1 if a > b.
///
///  a = The left-hand function name to compare.
///
///  b = The right-hand function name to compare.
///
///  return = A number indicating the alphabetical relationship between the two strings.
///
pub fn cmpFncName(a: []const u8, b: []const u8) f32 {
    const alen = a.len;
    const blen = b.len;
    var l: usize = 0;

    if (alen < blen) {
        l = alen;
    } else {
        l = blen;
    }

    for (0..l) |i| {
        if (a[i] < b[i]) {
            return -1;
        } else if (a[i] > b[i]) {
            return 1;
        } else if (a[i] == b[i]) {
            continue;
        }
    }

    if (alen < blen) {
        return -1;
    } else if (alen > blen) {
        return 1;
    } else {
        return 0;
    }
}

test "XMTX: cmpFncName test" {
    const a: []const u8 = "testing";
    const b: []const u8 = "example";
    var res: f32 = -3.0;
    var exp: f32 = 1.0;

    var start = try Instant.now();
    res = cmpFncName(a, b);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ncmpFncName #1: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cmpFncName", elapsed1);

    exp = -1.0;
    start = try Instant.now();
    res = cmpFncName(b, a);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\ncmpFncName #2: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("cmpFncName", elapsed1);

    try std.testing.expectEqual(exp, res);
    prntNl();
}

///Calculates the magnitude of a vector in a general inner product space.
///
///  vec = The vector to calculate the inner product for.
///
///  prdct = A function that takes two vectors as an argument and returns a f32 value determining the general inner product space.
///
///  returns = An f32 value representing the magnitude of vector vec in the general inner product space.
///
pub fn magInrPrdctSpcXvec(vec: []f32, prdct: *const fn (l: []f32, r: []f32) f32) f32 {
    return std.math.sqrt(prdct(vec, vec));
}

test "XMTX: magInrPrdctSpcXvec test" {
    const prdct: *const fn (l: []f32, r: []f32) f32 = dotPrdXvec;
    var u: [3]f32 = .{ 1, 2, 3 };
    var v: [3]f32 = .{ 2, 4, 6 };
    const v1: f32 = magInrPrdctSpcXvec(&u, prdct);
    const v2: f32 = magInrPrdctSpcXvec(&v, prdct);
    const e1: f32 = 3.74165749e+00;
    const e2: f32 = 7.48331499e+00;

    prntNlStrArgs("Found val 1: {}", .{v1});
    prntNlStrArgs("Found exp 1: {}", .{e1});
    try std.testing.expectEqual(v1, e1);
    prntNl();

    prntNlStrArgs("Found val 2: {}", .{v2});
    prntNlStrArgs("Found exp 2: {}", .{e2});
    try std.testing.expectEqual(v2, e2);
    prntNl();
}

///Calculates the distance of two vectors in a general inner product spaces.
///
///  vecL = The left hand vector in the inner product distance calculation, vecL - vecR.
///
///  vecR = The right hand vector in the inner product distance calculation, vecL - vecR.
///
///  prdct = A function that takes two vectors as an argument and returns a f32 value determining the general inner product space.
///
///  returns = An f32 value representing the distance between two vectors in a general inner product space or -1.0 if there's an error.
///
pub fn dstInrPrdctSpcXvec(vecL: []f32, vecR: []f32, prdct: *const fn (l: []f32, r: []f32) f32, alloc: *const std.mem.Allocator) !f32 {
    if (vecL.len != vecR.len) {
        prntNlStr("dstInrPrdctSpcXvec: Error, vectors u and v must be the same size, returning -1.0.");
        return Error.InvalidLengths;
    }

    const t: []f32 = try alloc.*.alloc(f32, vecL.len);
    defer alloc.*.free(t);

    diff2Xvec(@constCast(t), vecL, vecR);
    return magInrPrdctSpcXvec(@constCast(t), prdct);
}

test "XMTX: dstInrPrdctSpcXvec test" {
    var alloc = std.testing.allocator;
    const prdct: *const fn (l: []f32, r: []f32) f32 = dotPrdXvec;
    var u: [3]f32 = .{ 1, 0, 0 };
    var v: [3]f32 = .{ 0, 0, 0 };
    var val: f32 = try dstInrPrdctSpcXvec(&u, &v, prdct, &alloc);
    var exp: f32 = 1.0;
    prntNlStrArgs("Found val 1: {}", .{val});
    prntNlStrArgs("Found exp 1: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();

    u = .{ 1, 1, 1 };
    v = .{ 1, 0, 0 };
    val = try dstInrPrdctSpcXvec(&u, &v, prdct, &alloc);
    exp = 1.41421353e+00;
    prntNlStrArgs("Found val 2: {}", .{val});
    prntNlStrArgs("Found exp 2: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();
}

///Calculates the distance of two vectors in a general inner product spaces.
///
///  vecL = The left hand vector in the inner product distance calculation, vecL - vecR.
///
///  vecR = The right hand vector in the inner product distance calculation, vecL - vecR.
///
///  prdct = A function that takes two vectors as an argument and returns a f32 value determining the general inner product space.
///
///  returns = An f32 value representing the distance between two vectors in a general inner product space or -1.0 if there's an error.
///
pub fn dstXvec(vecL: []f32, vecR: []f32, alloc: *const std.mem.Allocator) !f32 {
    if (vecL.len != vecR.len) {
        prntNlStr("dstXvec: Error, vectors u and v must be the same length.");
        return Error.InvalidLengths;
    }

    const t: []f32 = try alloc.*.alloc(f32, vecL.len);
    defer alloc.*.free(t);

    diff2Xvec(@constCast(t), vecL, vecR);
    return magXvec(@constCast(t));
}

test "XMTX: dstXvec test" {
    var alloc = std.testing.allocator;
    var u: [3]f32 = .{ 1, 0, 0 };
    var v: [3]f32 = .{ 0, 0, 0 };
    var val: f32 = try dstXvec(&u, &v, &alloc);
    var exp: f32 = 1.0;
    prntNlStrArgs("Found val 1: {}", .{val});
    prntNlStrArgs("Found exp 1: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();

    u = .{ 1, 1, 1 };
    v = .{ 1, 0, 0 };
    val = try dstXvec(&u, &v, &alloc);
    exp = 1.41421353e+00;
    prntNlStrArgs("Found val 2: {}", .{val});
    prntNlStrArgs("Found exp 2: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();
}

///Calculates the angle in radians between two vectors in a general inner product space.
///
///  vecL = The left-hand vector to use in the angle calculation.
///
///  vecR = The right-hand vector to use in the angle calculation.
///
///  prdct = A function that takes two vectors as an argument and returns a f32 value determining the general inner product space.
///
///  returns = The angle between vecL and vecR in radians.
///
pub fn aglBtwnInrPrdctSpcXvec(vecL: []f32, vecR: []f32, prdct: *const fn (l: []f32, r: []f32) f32) f32 {
    const dotP = prdct(vecL, vecR);
    const mag1 = magInrPrdctSpcXvec(vecL, prdct);
    const mag2 = magInrPrdctSpcXvec(vecR, prdct);
    const cosA: f32 = (dotP / (mag1 * mag2));
    const arcCosA = std.math.acos(cosA);
    prntNlStrArgs("aglBtwnInrPrdctSpcXvec: Found cosA: {} arcCosA: {}", .{ cosA, arcCosA });
    return arcCosA;
}

test "XMTX: aglBtwnInrPrdctSpcXvec test" {
    var v1: [3]f32 = .{ 1, 0, 0 };
    var v2: [3]f32 = .{ 0, 0, 1 };
    var angle: f32 = -1.0;

    const start = try Instant.now();
    angle = aglBtwnInrPrdctSpcXvec(&v1, &v2, dotPrdXvec);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\naglBtwnInrPrdctSpcXvec: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("aglBtwnInrPrdctSpcXvec", elapsed1);

    const exp: f32 = 1.57079625e+00;
    prntXvecNl(&v1);
    std.debug.print("\nAngle RAD: {}", .{angle});
    std.debug.print("\nAngle DEG: {}", .{rad2Deg(angle)});
    try std.testing.expectEqual(exp, angle);
    prntNl();
}

///Determines if the vectors u, v, and w are members of a general inner product space.
///
///  u = A vector used to determine if there's an inner product space.
///
///  v = A vector used to determine if there's an inner product space.
///
///  w = A vector used to determine if there's an inner product space.
///
///  c = A scalar used to determine if there's an inner product space.
///
///  prdct = A function that takes two vector arguments and returns an f32 value.
///
///  returns = A Boolean value indicating that prdct creates an inner product space.
///
pub fn isInrPrdctSpc(u: []f32, v: []f32, w: []f32, c: f32, prdct: *const fn (l: []f32, r: []f32) f32, alloc: *const std.mem.Allocator) !bool {
    if (u.len != v.len or u.len != w.len or v.len != w.len) {
        prntNlStrArgs("isInrPrdctSpc: test 0 u, v, and w must have the same number of components.", .{});
        return Error.InvalidLengths;
    }

    const v1: f32 = prdct(u, v);
    const v2: f32 = prdct(v, u);
    if (v1 != v2) {
        prntNlStrArgs("isInrPrdctSpc: test 1 (<u,v> = <v,u>) failed {} != {}", .{ v1, v2 });
        return false;
    }

    const t: []f32 = try alloc.*.alloc(f32, u.len);
    defer alloc.*.free(t);
    //cpyXvec(w, t);
    clrXvec(t);
    sum2Xvec(t, v, w);

    const v3: f32 = prdct(u, t);
    const v4: f32 = prdct(u, w);
    const v5: f32 = (v1 + v4);

    if (v3 != v5) {
        prntNlStrArgs("isInrPrdctSpc: test 2 (<u, v+w> = <u,v> + <u,w>) failed {} != {}", .{ v3, v5 });
        return false;
    }

    clrXvec(t);
    cpyXvec(u, t);
    mulXvec(t, c);
    const v6: f32 = (c * v1);
    const v7: f32 = prdct(t, v);

    if (v6 != v7) {
        prntNlStrArgs("isInrPrdctSpc: test 3 (c<u,v> = <cu,v>) failed {} != {}", .{ v6, v7 });
        //alloc.*.free(t);
        return false;
    }

    const v8: f32 = prdct(v, v);
    clrXvec(t);
    const v9: f32 = prdct(t, t);

    if (v8 < 0) {
        prntNlStrArgs("isInrPrdctSpc: test 4 (<v,v> > 0) failed {} != {}", .{ v8, 0 });
        //alloc.*.free(t);
        return false;
    }

    if (v9 != 0) {
        prntNlStrArgs("isInrPrdctSpc: test 5 (<v,v> = 0 iff v = 0) failed {} != {}", .{ v9, 0 });
        //alloc.*.free(t);
        return false;
    }

    //alloc.*.free(t);
    return true;
}

test "XMTX: isInrPrdctSpc test" {
    var u: [3]f32 = .{ 1, 0, 0 };
    var v: [3]f32 = .{ 0, 1, 0 };
    var w: [3]f32 = .{ 0, 0, 1 };
    const c: f32 = 3.0;
    var b: bool = false;
    const alloc = std.testing.allocator;
    const exp: bool = true;

    const start = try Instant.now();
    b = try isInrPrdctSpc(&u, &v, &w, c, dotPrdXvec, &alloc);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nisInrPrdctSpc: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("isInrPrdctSpc", elapsed1);

    prntNlStrArgs("isInrPrdctSpc: b: {}", .{b});
    prntNlStrArgs("isInrPrdctSpc: exp: {}", .{exp});
    try std.testing.expectEqual(exp, b);
    prntNl();
}

///Calculates the general inner product given vectors u, v, and a function pointer, prdct.
///
///  u = A vector used to calculate the inner product space.
///
///  v = A vector used to calculate the inner product space.
///
///  prdct = A function that takes two vector arguments and returns an f32 value.
///
///  returns = A f32 value that is the calculated general innter product space for the given vectors u and v.
///
pub fn inrPrdct(u: []f32, v: []f32, prdct: *const fn (l: []f32, r: []f32) f32) f32 {
    return prdct(u, v);
}

test "XMTX: inrPrdct test" {
    const u: [3]f32 = .{ 1, 0, 0 };
    const v: [3]f32 = .{ 0, 1, 0 };
    var b: f32 = -1000.00;
    const exp: f32 = 0;

    const start = try Instant.now();
    b = inrPrdct(@constCast(&u), @constCast(&v), dotPrdXvec);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ninrPrdct: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("inrPrdct", elapsed1);

    prntNlStrArgs("inrPrdct: b: {}", .{b});
    prntNlStrArgs("inrPrdct: exp: {}", .{exp});
    try std.testing.expectEqual(exp, b);
    prntNl();
}

///Tests the vector arguments, u and v, with the given general inner product space, prdct, to verify the Cauchy-Schwarz Inequality theorem.
///
///  u = A vector used in this theorem test.
///
///  v = A vector used in this theorem test.
///
///  prdct = A function that takes two vector arguments and returns an f32 value. This function defines the inner product space.
///
///  returns = A boolean value indicating if the Cauchy-Schwarz Inequality theorem holds for the given arguments or an error value indicating the error encountered.
///
pub fn tstCauchySchwarzIneq(u: []f32, v: []f32, prdct: *const fn (l: []f32, r: []f32) f32) bool {
    const v1: f32 = absF32(prdct(u, v));
    const v2: f32 = magInrPrdctSpcXvec(u, prdct);
    const v3: f32 = magInrPrdctSpcXvec(v, prdct);
    const v4: f32 = (v2 * v3);
    return (v1 < v4);
}

test "XMTX: tstCauchySchwarzIneq test" {
    var u: [2]f32 = .{ 5, 12 };
    var v: [2]f32 = .{ 3, 4 };
    const prdct: *const fn (l: []f32, r: []f32) f32 = dotPrdXvec;
    var exp: bool = true;
    var val: bool = tstCauchySchwarzIneq(&u, &v, prdct);
    prntNlStrArgs("Found val 1: {}", .{val});
    prntNlStrArgs("Found exp 1: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();

    u = [2]f32{ 5, 15 };
    v = [2]f32{ 3, 5 };
    exp = true;
    val = tstCauchySchwarzIneq(&u, &v, prdct);
    prntNlStrArgs("Found val 2: {}", .{val});
    prntNlStrArgs("Found exp 2: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();

    u = [2]f32{ 5, 15 };
    v = [2]f32{ 0, 0 };
    exp = false;
    val = tstCauchySchwarzIneq(&u, &v, prdct);
    prntNlStrArgs("Found val 3: {}", .{val});
    prntNlStrArgs("Found exp 3: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();
}

///Tests the vector arguments, u and v, with the given general inner product space, prdct, to verify the Triangle Inequality theorem.
///
///  u = A vector used in this theorem test.
///
///  v = A vector used in this theorem test.
///
///  prdct = A function that takes two vector arguments and returns an f32 value. This function defines the inner product space.
///
///  alloc = A memory allocator.
///
///  returns = A boolean value indicating if the Triangle Inequality theorem holds for the given arguments or an error value indicating the error encountered.
///            Returns false if any triangle side is length zero.
///
pub fn tstTriangleIneq(u: []f32, v: []f32, prdct: *const fn (l: []f32, r: []f32) f32, alloc: *const std.mem.Allocator) !bool {
    const t: []f32 = try alloc.*.alloc(f32, u.len);
    defer alloc.*.free(t);
    clrXvec(t);
    sum2Xvec(t, u, v);
    const v1: f32 = magInrPrdctSpcXvec(t, prdct);
    const v2: f32 = magInrPrdctSpcXvec(u, prdct);
    const v3: f32 = magInrPrdctSpcXvec(v, prdct);
    const v4: f32 = (v2 + v3);
    if (v1 == 0 or v2 == 0 or v3 == 0) {
        return false;
    } else {
        return (v1 <= v4);
    }
}

test "XMTX: tstTriangleIneq test" {
    var u: [2]f32 = .{ 5, 5 };
    var v: [2]f32 = .{ 3, 3 };
    const prdct: *const fn (l: []f32, r: []f32) f32 = dotPrdXvec;
    const alloc = std.testing.allocator;
    var exp: bool = true;
    var val: bool = try tstTriangleIneq(&u, &v, prdct, &alloc);
    prntNlStrArgs("Found val 1: {}", .{val});
    prntNlStrArgs("Found exp 1: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();

    u = [2]f32{ 1, 1 };
    v = [2]f32{ -1, -1 };
    exp = false;
    val = try tstTriangleIneq(&u, &v, prdct, &alloc);
    prntNlStrArgs("Found val 2: {}", .{val});
    prntNlStrArgs("Found exp 2: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();

    u = [2]f32{ 1, 1 };
    v = [2]f32{ 0, 0 };
    exp = false;
    val = try tstTriangleIneq(&u, &v, prdct, &alloc);
    prntNlStrArgs("Found val 3: {}", .{val});
    prntNlStrArgs("Found exp 3: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();
}

///Tests the vector arguments, u and v, with the given general inner product space, prdct, to verify the Pythagorean theorem.
///
///  u = A vector used in this theorem test.
///
///  v = A vector used in this theorem test.
///
///  prdct = A function that takes two vector arguments and returns an f32 value. This function defines the inner product space.
///
///  alloc = A memory allocator.
///
///  returns = A boolean value indicating if the Pythagorean theorem holds for the given arguments or an error value indicating the error encountered.
///            Returns false if any triangle side is length zero.
///
pub fn tstPythagoreanTheorem(u: []f32, v: []f32, prdct: *const fn (l: []f32, r: []f32) f32, alloc: *const std.mem.Allocator) !bool {
    const t: []f32 = try alloc.*.alloc(f32, u.len);
    defer alloc.*.free(t);
    clrXvec(t);
    sum2Xvec(t, u, v);
    var v1: f32 = magInrPrdctSpcXvec(t, prdct);
    const v2: f32 = magInrPrdctSpcXvec(u, prdct);
    const v3: f32 = magInrPrdctSpcXvec(v, prdct);
    v1 = (v1 * v1);
    const v4: f32 = (v2 * v2) + (v3 * v3);
    if (v1 == 0 or v2 == 0 or v3 == 0) {
        return false;
    } else {
        return (isEquF32(v1, v4, false));
    }
}

test "XMTX: tstPythagoreanTheorem test" {
    var u: [2]f32 = .{ 0, 5 };
    var v: [2]f32 = .{ 3, 0 };
    const prdct: *const fn (l: []f32, r: []f32) f32 = dotPrdXvec;
    const alloc = std.testing.allocator;
    var exp: bool = true;
    var val: bool = try tstPythagoreanTheorem(&u, &v, prdct, &alloc);
    prntNlStrArgs("Found val 1: {}", .{val});
    prntNlStrArgs("Found exp 1: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();

    u = [2]f32{ 1, 0 };
    v = [2]f32{ 0, -1 };
    exp = true;
    val = try tstPythagoreanTheorem(&u, &v, prdct, &alloc);
    prntNlStrArgs("Found val 2: {}", .{val});
    prntNlStrArgs("Found exp 2: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();

    u = [2]f32{ 3, 3 };
    v = [2]f32{ 0, 0 };
    exp = false;
    val = try tstPythagoreanTheorem(&u, &v, prdct, &alloc);
    prntNlStrArgs("Found val 3: {}", .{val});
    prntNlStrArgs("Found exp 3: {}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();
}

///Projects vector P onto vector Q in a general inner product space. Alters argument vecQ.
///
///  vecP = The vector to project onto.
///
///  vecQ = The vector to be projected on.
///
///  prdct = A function that takes two vector arguments and returns an f32 value.
///
///  returns = The parameter of the vector P that projects onto Q, stored in vecQ.
///
pub fn projXvec_VecP_Onto_VecQ_InrPrdctSpc(vecP: []f32, vecQ: []f32, prdct: *const fn (l: []f32, r: []f32) f32) []f32 {
    //Larson, Edwards: Chapter 5: Definition of Orthogonal projection: pg 272
    const dotProd: f32 = prdct(vecP, vecQ);
    const mag: f32 = magXvec(vecQ);
    const magSqr: f32 = (mag * mag);
    const mul: f32 = (dotProd / magSqr);
    mulXvec(vecQ, mul);
    return vecQ;
}

test "XMTX: projXvec_VecP_Onto_VecQQ_InrPrdctSpc test" {
    //Q = 2 2 1
    //P = 1 -2 0
    var Q: [3]f32 = .{ 2, 2, 1 };
    var P: [3]f32 = .{ 1, -2, 0 };
    var exp: [3]f32 = .{ (-4.0 / 9.0), (-4.0 / 9.0), (-2.0 / 9.0) };

    const start = try Instant.now();
    const ret: []f32 = projXvec_VecP_Onto_VecQ_InrPrdctSpc(&P, &Q, dotPrdXvec);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nprojXvec_VecP_Onto_VecQ_InrPrdctSpc: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("projXvec_VecP_Onto_VecQ_InrPrdctSpc", elapsed1);

    clnXmtx(ret);
    clnXmtx(&exp);

    prntXvecNl(ret);
    prntNl();

    prntXvecNl(&exp);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(ret, &exp));
    prntNl();
}

///Projects vector P onto vector unit vector Q in a general inner product space. Alters argument vecQ.
///
///  vecP = The vector to project onto.
///
///  vecQ = The unit vector to be projected on.
///
///  returns = The parameters of the vector P that projects onto Q, stored in vecQ.
///
/// Note:
///
/// progUontoV = progPontoQ; if v is a unit vector then <v,v> = ||v||^2 = 1,
///
/// which simplifies the formula to <u,v>v
///
pub fn projUnitXvec_VecP_Onto_VecQ_InrPrdctSpc(vecP: []f32, vecQ: []f32, prdct: *const fn (l: []f32, r: []f32) f32) []f32 {
    //Larson, Edwards: Chapter 5: remark: pg 272
    const dotProd: f32 = prdct(vecP, vecQ);
    mulXvec(vecQ, dotProd);
    return vecQ;
}

test "XMTX: projUnitXvec_VecP_Onto_VecQ_InrPrdctSpc test" {
    //Q = 2 2 1
    //P = 1 -2 0
    var Q: [3]f32 = .{ 2, 2, 1 };
    nrmXvec(&Q);
    var P: [3]f32 = .{ 1, -2, 0 };
    nrmXvec(&P);

    var ret1: []f32 = undefined;
    var ret2: []f32 = undefined;
    const exp: [3]f32 = .{ -1.76676996e-02, -1.76676996e-02, -8.83384980e-03 };

    var start = try Instant.now();
    ret1 = projXvec_VecP_Onto_VecQ_InrPrdctSpc(&P, &Q, dotPrdXvec);
    var end = try Instant.now();
    var elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nprojXvec_VecP_Onto_VecQ_InrPrdctSpc: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("projXvec_VecP_Onto_VecQ_InrPrdctSpc", elapsed1);

    start = try Instant.now();
    ret2 = projUnitXvec_VecP_Onto_VecQ_InrPrdctSpc(&P, &Q, dotPrdXvec);
    end = try Instant.now();
    elapsed1 = @floatFromInt(end.since(start));
    std.debug.print("\nprojUnitXvec_VecP_Onto_VecQ_InrPrdctSpc: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("projUnitXvec_VecP_Onto_VecQ_InrPrdctSpc", elapsed1);

    clnXvec(ret1);
    clnXvec(ret2);
    clnXvec(@constCast(&exp));

    prntXvecNl(ret1);
    prntNl();
    prntXvecNl(ret2);
    prntNl();
    prntXvecNl(@constCast(&exp));
    prntNl();

    try std.testing.expectEqual(true, equXvecWrkr(ret1, ret2, false));
    try std.testing.expectEqual(true, equXvecWrkr(ret1, @constCast(&exp), false));
    try std.testing.expectEqual(true, equXvecWrkr(ret2, @constCast(&exp), false));
    prntNl();
}

///Returns the area of the triangle given with R2 coordinates and a default value of 1 for the 3rd, Z, component.
///The format of mtx is not checked. The caller is responsible for making sure the matrix is in the correct format.
///
///  Parameter mtx format:
///
///  |x1 y1 1|
///
///  |x2 y2 1|
///
///  |x3 y3 1|
///
///  mtx = The matrix of points to use for the calulation of the area of a triangle.
///
///  returns = A value indicating the area of the given triangle.
///
pub fn areaOfTriangle(mtx: *const [9]f32) f32 {
    return absF32((detXmtx3(@constCast(mtx)) * (1.0 / 2.0)));
}

test "XMTX: areaOfTriangle test" {
    const mtx: [9]f32 = .{ 1, 0, 1, 2, 2, 1, 4, 3, 1 };
    const exp: f32 = (3.0 / 2.0);
    var res: f32 = 0.0;

    const start = try Instant.now();
    res = areaOfTriangle(&mtx);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nareaOfTriangle: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("areaOfTriangle", elapsed1);

    try std.testing.expectEqual(exp, res);
    prntNl();
}

///Returns the volume of the tetrahedron given with R3 coordinates and a default value of 1 for the 4th, K, component.
///The format of mtx is not checked. The caller is responsible for making sure the matrix is in the correct format.
///
///  Parameter mtx format:
///
///  |x1 y1 z1 1|
///
///  |x2 y2 z2 1|
///
///  |x3 y3 y3 1|
///
///  |x4 y4 y4 1|
///
///  mtx = The matrix of points to use for the calulation of the volume of a tetrahedron.
///
///  returns = A value indicating the volume of the given tetrahedron.
///
pub fn volumeOfTetrahedron(mtx: *const [16]f32) f32 {
    return absF32((detXmtx4(@constCast(mtx)) * (1.0 / 6.0)));
}

test "XMTX: volumeOfTetrahedron test" {
    const mtx: [16]f32 = .{ 0, 4, 1, 1, 4, 0, 0, 1, 3, 5, 2, 1, 2, 2, 5, 1 };
    const exp: f32 = 12.0;
    var res: f32 = 0.0;

    const start = try Instant.now();
    res = volumeOfTetrahedron(&mtx);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\nvolumeOfTetrahedron: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("volumeOfTetrahedron", elapsed1);

    try std.testing.expectEqual(exp, res);
    prntNl();
}

///Returns the orthogonal set of vectors by storing their values in basis prime, the bsPrm function argument.
///
///  basis = A set of vectors that are considered a basis in the general inner product space defined by the argument prdct.
///
///  bsPrm = A matrix the same size as the basis argument used to hold the results of the algorithm.
///
///  cols = The number of columns in the matrix, basis, and subsequently bsPrm.
///
///  alloc = A constant pointer to a memory allocator.
///
///  prdct = A constant pointer to a function, with the given signature, that is used for the general inner product space.
///
///  returns = Usually void unless there was an error allocating memory for the algorithm.
///
pub fn gramSchmidtOthogonal(basis: []f32, bsPrm: []f32, cols: usize, alloc: *const std.mem.Allocator, prdct: *const fn (l: []f32, r: []f32) f32) !void {
    const rows: usize = basis.len / cols;
    var row: usize = 0;
    var sRow: usize = 0;
    var s: usize = 0;
    var e: usize = 0;
    var sM: usize = 0;
    var eM: usize = 0;
    var Vrow: []f32 = try alloc.*.alloc(f32, cols);
    var Wrow: []f32 = try alloc.*.alloc(f32, cols);
    var tmpM: []f32 = try alloc.*.alloc(f32, cols);
    var top: f32 = 0.0;
    var bot: f32 = 0.0;
    var i: usize = 0;

    if (VERBOSE) {
        prntNl();
        prntNlStr("Start Inner");
        prntNlStr("Basis:");
        prntXmtxNl(basis, cols);
    }

    while (row < rows) : (row += 1) {
        sM = (row * cols);
        eM = (sM + cols);

        //Vrow = basis[sM..eM];
        i = 0;
        while (i < cols) : (i += 1) {
            Vrow[i] = basis[sM + i];
        }

        //tmpM = Vrow[0..cols];
        i = 0;
        while (i < cols) : (i += 1) {
            tmpM[i] = Vrow[i];
        }

        sRow = 0;

        if (VERBOSE) {
            prntNl();
            prntNlStrArgs("Row: {} sM: {} eM: {}", .{ row, sM, eM });
            prntNlStr("Vrow:");
            prntXvecNl(Vrow);
            prntNlStr("tmpM:");
            prntXvecNl(tmpM);
        }

        const lt: i32 = @as(i32, @intCast(row)) - 1;
        while ((sRow <= lt and lt >= 0) or (sRow == 0 and lt == 0)) : (sRow += 1) {
            s = (sRow * cols);
            e = (s + cols);

            //Wrow = bsPrm[s..e];
            i = 0;
            while (i < cols) : (i += 1) {
                Wrow[i] = bsPrm[s + i];
            }

            top = prdct(Vrow, Wrow);
            bot = prdct(Wrow, Wrow);
            mulXvec(Wrow, (top / bot));
            diff1Xvec(tmpM, Wrow);

            if (VERBOSE) {
                prntNl();
                prntNlStrArgs("sRow: {} s: {} e: {}", .{ sRow, s, e });
                prntNlStr("Wrow:");
                prntXvecNl(Wrow);
                prntNlStrArgs("Top: {} Botttom: {} Top/Bottom: {}", .{ top, bot, (top / bot) });
                prntNlStr("tmpM:");
                prntXvecNl(tmpM);
            }
        }

        if (VERBOSE) {
            prntNl();
            prntNlStr("End Inner 1");
            prntNlStr("Vrow:");
            prntXvecNl(Vrow);
            prntNlStr("tmpM:");
            prntXvecNl(tmpM);
        }

        //copy the tmpM values into basis prime matrix, bsPrm
        //bsPrm[sM..eM] = tmpM[0..cols];
        //diff1Xvec(Vrow, tmpM);
        i = 0;
        while (i < cols) : (i += 1) {
            bsPrm[sM + i] = tmpM[i];
        }

        if (VERBOSE) {
            prntNl();
            prntNlStr("End Inner 2");
            prntNlStr("basis Prime:");
            prntXmtxNl(bsPrm, cols);
            //prntNlStr("tmpM:");
            //prntXvecNl(tmpM);
        }
    }

    defer alloc.*.free(Vrow);
    defer alloc.*.free(Wrow);
    defer alloc.*.free(tmpM);
}

test "XMTX: gramSchmidtOthogonal test" {
    //Larson, Edwards Chapter 5, Example 6 page 282 - 283 1st half
    var basis: [4]f32 = .{ 1, 1, 0, 1 };
    const cols: usize = 2;
    const alloc = std.testing.allocator;
    const prdct: *const fn (l: []f32, r: []f32) f32 = dotPrdXvec;
    var res: [4]f32 = .{ 0, 0, 0, 0 };
    const exp: [4]f32 = .{ 1, 1, -0.5, 0.5 };

    const start = try Instant.now();
    try gramSchmidtOthogonal(&basis, &res, cols, &alloc, prdct);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ngramSchmidtOthogonal: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("gramSchmidtOthogonal", elapsed1);

    prntNlStr("Basis:");
    prntXmtxNl(&basis, cols);
    prntNlStr("Res:");
    prntXmtxNl(&res, cols);
    prntNlStr("Exp:");
    prntXmtxNl(@constCast(&exp), cols);

    try std.testing.expectEqual(true, equXvecWrkr(&res, @constCast(&exp), false));
    prntNl();
}

///Returns the orthonormal set of vectors by storing their values in basis prime, the bsPrm function argument.
///
///  basis = A set of vectors that are considered a basis in the general inner product space defined by the argument prdct.
///
///  bsPrm = A matrix the same size as the basis argument used to hold the results of the algorithm.
///
///  cols = The number of columns in the matrix, basis, and subsequently bsPrm.
///
///  alloc = A constant pointer to a memory allocator.
///
///  prdct = A constant pointer to a function, with the given signature, that is used for the general inner product space.
///
///  returns = Usually void unless there was an error allocating memory for the algorithm.
///
pub fn gramSchmidtOthonormal(basis: []f32, bsPrm: []f32, cols: usize, alloc: *const std.mem.Allocator, prdct: *const fn (l: []f32, r: []f32) f32) !void {
    try gramSchmidtOthogonal(basis, bsPrm, cols, alloc, prdct);

    const rows: usize = basis.len / cols;
    var row: usize = 0;
    var sM: usize = 0;
    var eM: usize = 0;
    var tmpM: []f32 = try alloc.*.alloc(f32, cols);
    var i: usize = 0;

    while (row < rows) : (row += 1) {
        sM = (row * cols);
        eM = (sM + cols);

        i = 0;
        while (i < cols) : (i += 1) {
            tmpM[i] = bsPrm[sM + i];
        }

        nrmXvec(tmpM);

        i = 0;
        while (i < cols) : (i += 1) {
            bsPrm[sM + i] = tmpM[i];
        }
    }

    defer alloc.*.free(tmpM);
}

test "XMTX: gramSchmidtOthonormal test" {
    //Larson, Edwards Chapter 5, Example 6 page 282 - 283 2nd half
    var basis: [4]f32 = .{ 1, 1, 0, 1 };
    const cols: usize = 2;
    const alloc = std.testing.allocator;
    const prdct: *const fn (l: []f32, r: []f32) f32 = dotPrdXvec;
    var res: [4]f32 = .{ 0, 0, 0, 0 };
    const exp: [4]f32 = .{ std.math.sqrt(2.0) / 2.0, std.math.sqrt(2.0) / 2.0, -std.math.sqrt(2.0) / 2.0, std.math.sqrt(2.0) / 2.0 };

    const start = try Instant.now();
    try gramSchmidtOthonormal(&basis, &res, cols, &alloc, prdct);
    const end = try Instant.now();
    const elapsed1: f64 = @floatFromInt(end.since(start));
    std.debug.print("\ngramSchmidtOthonormal: Time elapsed is: {d:.3}ms, {d:.3}ns", .{ elapsed1 / time.ns_per_ms, elapsed1 });
    addExecTime("gramSchmidtOthonormal", elapsed1);

    prntNlStr("Basis:");
    prntXmtxNl(&basis, cols);
    prntNlStr("Res:");
    prntXmtxNl(&res, cols);
    prntNlStr("Exp:");
    prntXmtxNl(@constCast(&exp), cols);

    try std.testing.expectEqual(true, equXvecWrkr(&res, @constCast(&exp), false));
    prntNl();
}

///Returns the coordinate relative to the given basis and relative coordinate, coordRel, for the given inner product space.
///
///  basis: A matrix of vectors that form the given basis.
///
///  cols: The number of columns in the given matrix, basis.
///
///  coordRel: The coordinate relative to the basis.
///
///  alloc: Allocation used to create a temporary vector for the calculation.
///
///  prdct: The function to use to calculate the inner product space.
///
pub fn coordRelOrthBasisInrPrdctSpc(basis: []f32, cols: usize, coordRel: []f32, res: []f32, alloc: *const std.mem.Allocator, prdct: *const fn (l: []f32, r: []f32) f32) !void {
    if (cols != coordRel.len) {
        prntNlStr("!! coordRelOrthBasis: Error, the vector coordRel must have the same length as the matrix basis has cols !!");
        return Error.InvalidLengths;
    } else if (basis.len % cols != 0) {
        prntNlStr("!! coordRelOrthBasis: Error, the matrix basis length must be evenly divisible by the given number of cols !!");
        return Error.InvalidLengths;
    } else if (cols != res.len) {
        prntNlStr("!! coordRelOrthBasis: Error, the vector res must have the same length as the matrix basis has cols !!");
        return Error.InvalidLengths;
    }

    const rows: usize = basis.len / cols;
    var row: usize = 0;
    var s: usize = 0;
    //var e: usize = 0;
    var Vrow: []f32 = try alloc.*.alloc(f32, cols);
    var i: usize = 0;

    //clear out the res matrix
    //clrXvec(res);

    while (row < rows) : (row += 1) {
        s = (row * cols);
        //e = (s + cols);

        //Vrow = basis[s..e];
        i = 0;
        while (i < cols) : (i += 1) {
            Vrow[i] = basis[s + i];
        }

        res[row] = prdct(Vrow, coordRel);
        //mulXvec(Vrow, prdct(Vrow, coordRel));
        //sum1Xvec(res, Vrow);
    }

    defer alloc.*.free(Vrow);
}

test "XMTX: coordRelOrthBasisInrPrdctSpc test" {
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
    const cols: usize = 2;
    var res1: [2]f32 = .{ 0.0, 0.0 };
    var exp1: [2]f32 = .{ ((4.0 * std.math.sqrt(13.0)) / 13.0), ((7.0 * std.math.sqrt(13.0)) / 13.0) };

    try coordRelOrthBasisInrPrdctSpc(&B1, cols, &x1, &res1, &alloc, dotPrdXvec);

    prntNlStr("Problem 13:");
    prntNlStr("Basis B1:");
    prntXmtxNl(&B1, cols);
    prntNlStr("Vector x1:");
    prntXvecNl(&x1);
    prntNlStrArgs("Cols: {}", .{cols});
    prntNlStr("Res1:");
    prntXvecNl(&res1);
    prntNlStr("Exp1:");
    prntXvecNl(&exp1);

    try std.testing.expectEqual(true, equXvecWrkr(&exp1, &res1, false));
    prntNl();
}

///Returns the coordinate relative to the given orthogonal basis and relative coordinate, coordRel.
///
///  basis: A matrix of vectors that form the given basis.
///
///  cols: The number of columns in the given matrix, basis.
///
///  coordRel: The coordinate relative to the basis.
///
///  alloc: Allocation used to create a temporary vector for the calculation.
///
pub fn coordRelOrthBasis(basis: []f32, cols: usize, coordRel: []f32, res: []f32, alloc: *const std.mem.Allocator) !void {
    if (cols != coordRel.len) {
        prntNlStr("!! coordRelOrthBasis: Error, the vector coordRel must have the same length as the matrix basis has cols !!");
        return Error.InvalidLengths;
    } else if (basis.len % cols != 0) {
        prntNlStr("!! coordRelOrthBasis: Error, the matrix basis length must be evenly divisible by the given number of cols !!");
        return Error.InvalidLengths;
    } else if (cols != res.len) {
        prntNlStr("!! coordRelOrthBasis: Error, the vector res must have the same length as the matrix basis has cols !!");
        return Error.InvalidLengths;
    }

    const rows: usize = basis.len / cols;
    var row: usize = 0;
    var s: usize = 0;
    //var e: usize = 0;
    var Vrow: []f32 = try alloc.*.alloc(f32, cols);
    var i: usize = 0;

    //clear out the res matrix
    clrXvec(res);

    while (row < rows) : (row += 1) {
        s = (row * cols);
        //e = (s + cols);

        //Vrow = basis[s..e];
        i = 0;
        while (i < cols) : (i += 1) {
            Vrow[i] = basis[s + i];
        }

        res[row] = dotPrdXvec(Vrow, coordRel);
        //mulXvec(Vrow, dotPrdXvec(Vrow, coordRel));
        //sum1Xvec(res, Vrow);
    }

    defer alloc.*.free(Vrow);
}

test "XMTX: coordRelOrthBasis test" {
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
    const cols: usize = 2;
    var res1: [2]f32 = .{ 0.0, 0.0 };
    var exp1: [2]f32 = .{ ((4.0 * std.math.sqrt(13.0)) / 13.0), ((7.0 * std.math.sqrt(13.0)) / 13.0) };

    try coordRelOrthBasis(&B1, cols, &x1, &res1, &alloc);

    prntNlStr("Problem 13:");
    prntNlStr("Basis B1:");
    prntXmtxNl(&B1, cols);
    prntNlStr("Vector x1:");
    prntXvecNl(&x1);
    prntNlStrArgs("Cols: {}", .{cols});
    prntNlStr("Res1:");
    prntXvecNl(&res1);
    prntNlStr("Exp1:");
    prntXvecNl(&exp1);

    try std.testing.expectEqual(true, equXvecWrkr(&exp1, &res1, false));
    prntNl();
}

///Finds the projection of vector V, vecV, onto the Subspace S, defined by the basis matrix mtxBasis.
///
///  vecV: The vector that we're finding the projection for, onto a subspace.
///  
///  mtxBasis: A matrix of vectors that form the given basis.
///
///  colsBasis: The number of columns in the given matrix, basis.
///
///  res: A matrix used to store the result of the calculation.
///
///  alloc: Allocation used to create a temporary vector for the calculation.
///
pub fn projXvec_VecV_Onto_SubspaceS(vecV: []f32, mtxBasis: []f32, colsBasis: usize, res: []f32, alloc: *const std.mem.Allocator) !void {
    const lenBasis: usize = mtxBasis.len;
    const rowsBasis: usize = lenBasis / colsBasis;
    const mtxOrthoBasis: []f32 = try alloc.*.alloc(f32, lenBasis);
    try gramSchmidtOthonormal(mtxBasis, mtxOrthoBasis, colsBasis, alloc, dotPrdXvec);

    var val: f32 = 0.0;
    var cnt: usize = 0;
    const vecU: []f32 = try alloc.*.alloc(f32, colsBasis);

    while(cnt < rowsBasis) : (cnt += 1) {
        cpyLessColRowXmtx(mtxOrthoBasis, @constCast(vecU), 0, colsBasis, cnt, (cnt + 1), colsBasis, colsBasis);
        val = dotPrdXvec(vecV, vecU);
        mulXvec(vecU, val);
        sum1Xvec(res, vecU);
    }

    defer alloc.*.free(mtxOrthoBasis);
    defer alloc.*.free(vecU);
}

test "XMTX: projXvec_VecV_Onto_SubspaceS test" {
    //Find the projection of the vector v = [1, 1, 3] 
    //onto the subspace S of R^3 spanned by the vectors
    //w1 = [0, 3, 1] and w2 = [2, 0, 0]
    //which normalize to [0, 3/sqrt(10), 1/sqrt(10)], [1, 0, 0]

    const alloc = std.testing.allocator;
    var vecV: [3]f32 = .{1, 1, 3};
    var mtxBasis: [6]f32 = .{0, 3.0 / std.math.sqrt(10.0), 1.0 / std.math.sqrt(10.0), 1, 0, 0};
    const colsBasis: usize = 3;
    var res: [3]f32 = .{0, 0, 0};
    var exp: [3]f32 = .{1.0, (9.0 / 5.0), (3.0 / 5.0)};
    try projXvec_VecV_Onto_SubspaceS(&vecV, &mtxBasis, colsBasis, &res, &alloc);

    prntNlStr("Example 5:");
    prntNlStr("Basis:");
    prntXmtxNl(&mtxBasis, colsBasis);
    prntNlStr("Vector V:");
    prntXvecNl(&vecV);
    prntNlStr("Res:");
    prntXvecNl(&res);
    prntNlStr("Exp:");
    prntXvecNl(&exp);

    try std.testing.expectEqual(true, equXvecWrkr(&exp, &res, false));
    prntNl();    
}

///Finds the least squares solution given the matrix A, mtxA, and the vector B, vecB.
///
///  mtxA: The matrix A, mtxA, that is used as the set of data driving the least squares solution.
/// 
///  colA: The number of columns in the matrix A, mtxA.
/// 
///  vecB: The target vector to find the least squares solution that is closest to.
///  
///
///  res: A matrix used to store the result of the calculation.
///
///  alloc: Allocation used to create a temporary vector for the calculation.
///
pub fn leastSquaresSol(mtxA: []f32, colsA: usize, vecB: []f32, res: []f32, alloc: *const std.mem.Allocator) !bool {
    const lenA: usize = mtxA.len;
    const rowsA: usize = (lenA / colsA);

    const lenB: usize = vecB.len;
    const colsB: usize = 1;
    const rowsB: usize = (lenB / colsB);

    const mtxATrn: []f32 = try alloc.*.alloc(f32, lenA);
    const lenATrn: usize = mtxATrn.len;
    const colsATrn: usize = rowsA;
    const rowsATrn: usize = colsA;

    const mtxATrnA: []f32 = try alloc.*.alloc(f32, (rowsATrn * colsA));
    const lenATrnA: usize = mtxATrnA.len;
    const colsATrnA: usize = colsA;
    const rowsATrnA: usize = rowsATrn;

    const mtxATrnB: []f32 = try alloc.*.alloc(f32, (rowsATrn * colsB));
    const lenATrnB: usize = mtxATrnB.len;
    const colsATrnB: usize = colsB;
    const rowsATrnB: usize = rowsATrn;    

    var sclr: f32 = 0.0;
    var b: bool = false;

    const mtxAug: []f32 = try alloc.*.alloc(f32, (mtxATrnA.len + mtxATrnB.len));
    const lenAug: usize = mtxAug.len;
    const colsAug: usize = (colsATrnA + 1);
    const rowsAug: usize = (rowsATrnA);

    const idtMtx: []f32 = try alloc.*.alloc(f32, mtxATrnA.len);

    clrXmtx(mtxATrn);
    clrXmtx(mtxATrnA);
    clrXmtx(mtxATrnB);    
    clrXmtx(mtxAug);
    
    trnXmtxRect(@constCast(mtxA), colsA, @constCast(mtxATrn), colsATrn);
    if(VERBOSE or true) {
        prntNlStr("mtxATrn AAA:");
        prntNlStrArgs("LenATrn: {}, ColsATrn: {}, RowsAtrn: {}", .{lenATrn, colsATrn, rowsATrn});    
        prntNl(); 
        prntXmtx(@constCast(mtxATrn), colsATrn);        
    }

    b = try tmsXmtx(@constCast(mtxATrn), colsATrn, @constCast(mtxA), colsA, @constCast(mtxATrnA), colsATrnA);
    if(!b) {
        return Error.OperationFailed;
    }

    b = try tmsXmtx( @constCast(mtxATrn), colsATrn, @constCast(vecB), colsB, @constCast(mtxATrnB), colsATrnB);
    if(!b) {
        return Error.OperationFailed;        
    }

    cpyXmtxSqr(@constCast(mtxATrnA), colsATrnA, @constCast(mtxAug), colsAug, 0, colsATrnA, 0, rowsATrnA, 0, 0);
    cpyXmtxSqr(@constCast(mtxATrnB), colsATrnB, @constCast(mtxAug), colsAug, 0, colsATrnB, 0, rowsATrnB, colsATrnA, 0);

    defer alloc.*.free(mtxATrnA);
    defer alloc.*.free(mtxATrnB);
    defer alloc.*.free(mtxATrn);
    defer alloc.*.free(mtxAug);
    defer alloc.*.free(idtMtx);    

    b = try rdcXmtxInl(@constCast(mtxAug), colsAug, true, false, idtMtx, (colsAug - 1), false, &sclr);
    if(VERBOSE) {
        prntNlStr("MtxA:");
        prntNlStrArgs("LenA: {}, ColsA: {}, RowsA: {}", .{lenA, colsA, rowsA});
        prntNl(); 
        prntXmtx(@constCast(mtxA), colsA);

        prntNlStr("VecB:");
        prntNlStrArgs("LenB: {}, ColsB: {}, RowsB: {}", .{lenB, colsB, rowsB});    
        prntNl(); 
        prntXmtx(vecB, colsB);

        prntNlStr("mtxATrn:");
        prntNlStrArgs("LenATrn: {}, ColsATrn: {}, RowsAtrn: {}", .{lenATrn, colsATrn, rowsATrn});    
        prntNl(); 
        prntXmtx(@constCast(mtxATrn), colsATrn);

        prntNlStr("mtxATrnA:");
        prntNlStrArgs("LenATrnA: {}, ColsATrnA: {}, RowsAtrnA: {}", .{lenATrnA, colsATrnA, rowsATrnA});    
        prntNl(); 
        prntXmtx(@constCast(mtxATrnA), colsATrnA);

        prntNlStr("mtxATrnB:");
        prntNlStrArgs("LenATrnB: {}, ColsATrnB: {}, RowsAtrnB: {}", .{lenATrnB, colsATrnB, rowsATrnB});    
        prntNl(); 
        prntXmtx(@constCast(mtxATrnB), colsATrnB);

        prntNlStr("mtxAug:");
        prntNlStrArgs("LenAug: {}, ColsAug: {}, RowsAug: {}", .{lenAug, colsAug, rowsAug});    
        prntNl(); 
        prntXmtx(@constCast(mtxAug), colsAug);
    }

    //_ = rowsB;
    //_ = rowsAug;

    if(b) {
        cpyXmtxSqr(@constCast(mtxAug), colsAug, res, 1, colsATrnA, (colsATrnA + 1), 0, rowsATrnA, 0, 0);
        return true;
    } else {
        return Error.OperationFailed;
    }

    return false;
}

test "XMTX: leastSquaresSol test" {
    //Section 5.4 Example 7
    //Find the solution to the least squares problem Ax = b.
    //| 1  1 | | c0 | = | 0 |
    //| 1  2 | | c1 |   | 1 |
    //| 1  3 |          | 3 |

    const alloc = std.testing.allocator;
    var mtxA: [6]f32 = .{1, 1, 1, 2, 1, 3};
    var vecB: [3]f32 = .{0, 1, 3};
    const colsA: usize = 2;
    var res: [2]f32 = .{0, 0};
    var exp: [2]f32 = .{(-5.0 / 3.0), (3.0 / 2.0)};
    var b: bool = false;

    b = try leastSquaresSol(&mtxA, colsA, &vecB, &res, &alloc);

    prntNlStr("Example 7:");
    prntNlStr("mtxA:");
    prntNl();
    prntXmtxNl(&mtxA, colsA);
    prntNlStr("Vector B:");
    prntXvecNl(&vecB);
    prntNlStr("Res:");
    prntXvecNl(&res);
    prntNlStr("Exp:");
    prntXvecNl(&exp);
    prntNlStrArgs("LeastSquaresSol Result: {}", .{b});

    try std.testing.expectEqual(true, equXvecWrkr(&exp, &res, false));
    prntNl();
}

//Compile function execution summary
test "XMTX: sortExecTimeList process" {
    //collapse exec time entries into a string hashmap
    var t1: ?ExecTime = null;
    var map = std.StringHashMap(ExecTime).init(std.testing.allocator);
    var et1: ExecTime = undefined;
    for (0..execTimesCnt) |i| {
        et1 = execTimes[i];
        t1 = map.get(et1.fncName);
        if (t1 == null) {
            et1.fncCnt = 1;
            et1.fncAvg = et1.fncTime;
            try map.put(et1.fncName, et1);
        } else {
            t1.?.fncCnt += 1.0;
            t1.?.fncTime += et1.fncTime;
            t1.?.fncAvg = et1.fncTime / et1.fncCnt;
            try map.put(t1.?.fncName, t1.?);
        }
    }
    defer map.deinit();

    //copy string hash map into an array
    prntNl();
    var klen: usize = 0;
    var keys: [1024]ExecTime = undefined;
    var iterator = map.iterator();
    while (iterator.next()) |execTime| {
        keys[klen] = .{
            .fncName = execTime.key_ptr.*,
            .fncTime = execTime.value_ptr.*.fncTime,
            .fncCnt = execTime.value_ptr.*.fncCnt,
            .fncAvg = execTime.value_ptr.*.fncAvg,
        };
        klen += 1;
        //prntNlStrArgs("MapList: {s}: {} {d:.3}ms {d:.3}ns", .{ execTime.key_ptr.*, execTime.value_ptr.*.fncCnt, execTime.value_ptr.*.fncAvg / time.ns_per_ms, execTime.value_ptr.*.fncAvg });
    }
    prntNl();

    //Test string compare
    //const tt1: *const [6:0]u8 = "absF32";
    //const tt2: *const [6:0]u8 = "absF32";
    //const tt3: *const [9:0]u8 = "absF32Ref";
    //const aa1: f32 = cmpFncName(tt1, tt2);
    //const aa2: f32 = cmpFncName(tt1, tt3);
    //const aa3: f32 = cmpFncName(tt3, tt1);
    //prntNlStrArgs("Answer1: {}", .{aa1});
    //prntNlStrArgs("Answer2: {}", .{aa2});
    //prntNlStrArgs("Answer3: {}", .{aa3});

    //sot exec times
    const n: usize = klen;
    var temp: ExecTime = undefined;
    for (0..n) |i| {
        for (1..(n - i)) |j| {
            if (cmpFncName(keys[j - 1].fncName, keys[j].fncName) == 1) {
                //swap elements
                temp = keys[j - 1];
                keys[j - 1] = keys[j];
                keys[j] = temp;
            }
        }
    }

    //List exec times
    prntNlStr("Function Execution Times List:");
    for (0..klen) |i| {
        prntNlStrArgs("{s}:\tCount: {}\tAvg: {d:.3}ms {d:.3}ns", .{ keys[i].fncName, keys[i].fncCnt, keys[i].fncAvg / time.ns_per_ms, keys[i].fncAvg });
    }
    prntNl();

    //var sum: f64 = execTimes[0].fncTime;
    //var avg: f64 = 0.0;
    //var cnt: f64 = 1.0;
    //var crntFncName: []u8 = @constCast(execTimes[0].fncName);
    //prntNlStrArgs("{s}: {d:.3}ms {d:.3}ns", .{ execTimes[0].fncName, execTimes[0].fncTime / time.ns_per_ms, execTimes[0].fncTime });
    //for (1..execTimesCnt) |i| {
    //    const execTime = execTimes[i];
    //    if (eqlFncName(execTime.fncName, crntFncName)) {
    //        sum += execTime.fncTime;
    //        cnt += 1.0;
    //    } else {
    //        avg = sum / cnt;
    //        prntNlStrArgs("AVG: {s}: {} {d:.3}ms {d:.3}ns", .{ crntFncName, cnt, avg / time.ns_per_ms, avg });
    //
    //        crntFncName = @constCast(execTime.fncName);
    //        sum = execTime.fncTime;
    //        cnt = 1.0;
    //        avg = 0.0;
    //    }
    //    prntNlStrArgs("{s}: {d:.3}ms {d:.3}ns", .{ execTime.fncName, execTime.fncTime / time.ns_per_ms, execTime.fncTime });
    //}
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
//START THEOREMS
//Mathematics for 3D Game Programming and Computer Graphics - Eric Lengyel - 3rd Edition
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

test "XMTX: MF3D - Lengyel: Theorem 2.1 test" {
    prntNl();
    std.debug.print("Theorem 2.1:\n", .{});
    var P: [2]f32 = .{ 1, 3 }; //P
    var Q: [2]f32 = .{ 2, 6 }; //Q
    var R: [2]f32 = .{ 3, 9 }; //R
    var La: [2]f32 = .{ 0, 0 }; //P + Q
    var Ra: [2]f32 = .{ 0, 0 }; //Q + P
    sum2Xvec(&La, &P, &Q);
    sum2Xvec(&Ra, &Q, &P);
    std.debug.print("XMTX: P + Q = Q + P:\n", .{});
    try std.testing.expectEqual(true, equXvec(&La, &Ra));

    P = .{ 1, 3 }; //P
    Q = .{ 2, 6 }; //Q
    R = .{ 3, 9 }; //R
    La = .{ 0, 0 }; //(P + Q) + R
    Ra = .{ 0, 0 }; //P + (Q + R)
    sum2Xvec(&La, &P, &Q);
    sum1Xvec(&La, &R);
    sum2Xvec(&Ra, &Q, &R);
    sum1Xvec(&Ra, &P);
    std.debug.print("XMTX: (P + Q) + R = P + (Q + R):\n", .{});
    try std.testing.expectEqual(true, equXvec(&La, &Ra));

    const a: f32 = 3;
    const b: f32 = 5;
    const c: f32 = a * b;
    P = .{ 1, 3 }; //P
    Q = .{ 2, 6 }; //Q
    R = .{ 3, 9 }; //R
    La = .{ 1, 1 }; //a * (b * Q)
    Ra = .{ 1, 1 }; //(a * b) * Q
    mulXvec(&La, b);
    mulXvec(&La, a);
    mulXvec(&Ra, c);
    std.debug.print("XMTX: a * (b * Q) = (a * b) * Q:\n", .{});
    try std.testing.expectEqual(true, equXvec(&La, &Ra));

    P = .{ 1, 3 }; //P
    Q = .{ 2, 6 }; //Q
    R = .{ 3, 9 }; //R
    La = .{ 0, 0 }; //a * (P + Q)
    Ra = .{ 0, 0 }; //(a * P) + (a * Q)
    sum2Xvec(&La, &P, &Q);
    mulXvec(&La, a);
    sum1Xvec(&Ra, &P);
    mulXvec(&Ra, a);
    var t1: [2]f32 = .{ 0, 0 };
    sum1Xvec(&t1, &Q);
    mulXvec(&t1, a);
    sum1Xvec(&Ra, &t1);
    std.debug.print("a * (P + Q) = (a * P) + (a * Q):\n", .{});
    try std.testing.expectEqual(true, equXvec(&La, &Ra));

    P = .{ 1, 3 }; //P
    Q = .{ 2, 6 }; //Q
    R = .{ 3, 9 }; //R
    La = .{ 0, 0 }; //(a + b) * P
    Ra = .{ 0, 0 }; //(a * P) + (b * P)
    const d: f32 = (a + b);
    _ = d;
    std.debug.print("(a + b) * P = (a * P) + (b * P):\n", .{});
    try std.testing.expectEqual(true, equXvec(&La, &Ra));
}

test "XMTX: MF3D - Lengyel: Theorem 2.2 test" {
    prntNl();
    std.debug.print("Theorem 2.2:\n", .{});

    var P: [2]f32 = .{ 1, 3 }; //P
    var Q: [2]f32 = .{ 2, 6 }; //Q
    var La: f32 = 0; //||P||
    var Ra: f32 = 0; // 0
    La = magXvec(&P);
    Ra = 0;
    std.debug.print("||P|| >= 0:\n", .{});
    try std.testing.expectEqual(true, (La >= Ra));

    P = .{ 1, 3 }; //P
    La = magXvec(&P); //||P||
    Ra = 0;
    var b: bool = isZeroXvec(&P); //isZero
    std.debug.print("||P|| = 0 if isZerO(P) = true:\n", .{});
    try std.testing.expectEqual(true, (La > Ra));

    P = .{ 0, 0 }; //P
    La = magXvec(&P); //||P||
    b = isZeroXvec(&P); //isZero
    Ra = 0;
    try std.testing.expectEqual(true, (La == Ra));

    P = .{ 1, 3 };
    var a: f32 = 7; //a
    var c: f32 = absF32(a); //|a|
    Q = .{ 1, 3 };
    mulXvec(&Q, a); //a * P
    a = magXvec(&Q); //||a * P||
    La = a;
    a = magXvec(&P); //||P||
    Ra = c * a; //|a| * ||P||
    std.debug.print("||a * P|| = |a| * ||P||:\n", .{});
    try std.testing.expectEqual(true, (La == Ra));

    P = .{ 1, 3 };
    Q = .{ 3, 9 };
    var R: [2]f32 = .{ 1, 3 };
    sum1Xvec(&P, &Q);
    La = magXvec(&P); //||P + Q||
    a = magXvec(&R); //||P||
    c = magXvec(&Q); //||Q||
    Ra = (a + c); //||P|| + ||Q||
    std.debug.print("||P + Q|| <= ||P|| + ||Q||:\n", .{});
}

test "XMTX: MF3D - Lengyel: Theorem 2.4 test" {
    prntNl();
    std.debug.print("Theorem 2.4:\n", .{});
    var P: [2]f32 = .{ 1, 3 };
    var Q: [2]f32 = .{ 6, 9 };
    const dotP = dotPrdXvec(&P, &Q);
    const mag1 = magXvec(&P);
    const mag2 = magXvec(&Q);
    const cosA: f32 = (dotP / (mag1 * mag2));
    const arcCosA = std.math.acos(cosA);
    std.debug.print("dotProd(P, Q) = ||P|| * ||Q|| * cos(a):\n", .{});
    try std.testing.expectEqual(cosA, std.math.cos(arcCosA));
}

test "XMTX: MF3D - Lengyel: Theorem 2.5 test" {
    prntNl();
    std.debug.print("Theorem 2.5:\n", .{});
    var P: [2]f32 = .{ 1, 3 }; //P
    var Q: [2]f32 = .{ 6, 9 }; //Q
    var dotP1: f32 = dotPrdXvec(&P, &Q);
    var dotP2: f32 = dotPrdXvec(&Q, &P);
    std.debug.print("dotProd(P, Q) = dotProd(Q, P):\n", .{});
    try std.testing.expectEqual(dotP1, dotP2);

    const a: f32 = 21;
    P = .{ 1, 3 }; //P
    Q = .{ 6, 9 }; //Q
    var R: [2]f32 = .{ 1, 1 };
    cpyXvec(&P, &R); //R = P
    prntXvecNl(&R);
    mulXvec(&R, a); //a * P
    dotP1 = dotPrdXvec(&R, &Q); //dotProd(a*P, Q)
    prntXvecNl(&R);
    std.debug.print("dotProd(a*P, Q): {}\n", .{dotP1});
    std.debug.print("dotProd(P, Q):: {}\n", .{dotPrdXvec(&P, &Q)});
    dotP2 = (a * dotPrdXvec(&P, &Q));
    std.debug.print("a * dotProd(a*P, Q): {} == {}\n", .{ dotP2, dotP1 });
    std.debug.print("dotProd(a*P, Q) = a * dotProd(P, Q):\n", .{});
    try std.testing.expectEqual(true, (dotP1 == dotP2));

    P = .{ 1, 3 }; //P
    Q = .{ 6, 9 }; //Q
    R = .{ -1, -13 }; //R
    sum1Xvec(&R, &Q);
    dotP1 = dotPrdXvec(&P, &R);
    var S: [2]f32 = .{ -1, -13 };
    dotP2 = dotPrdXvec(&P, &Q) + dotPrdXvec(&P, &S);
    std.debug.print("dotProd(P, (Q + R)) = dotProd(P, Q) + dotProd(P, R):\n", .{});
    try std.testing.expectEqual(dotP1, dotP2);

    P = .{ 1, 3 }; //P
    Q = .{ 6, 9 }; //Q
    const mag1: f32 = magXvec(&P);
    dotP1 = dotPrdXvec(&P, &P);
    std.debug.print("dotProd(P, P) = ||P|| * ||P||:\n", .{});
    try std.testing.expectEqual(dotP1, (mag1 * mag1));

    P = .{ 1, 3 }; //P
    Q = .{ 6, 9 }; //Q
    dotP1 = dotPrdXvec(&P, &Q);
    const mag2: f32 = magXvec(&Q);
    std.debug.print("|dotProd(P, Q)| <= ||P|| * ||Q||:\n", .{});
    try std.testing.expectEqual((dotP1 <= (mag1 * mag2)), true);
}

test "XMTX: MF3D - Lengyel: Theorem 2.7 test" {
    prntNl();
    var P: [3]f32 = .{ 1, 3, 6 }; //P
    var Q: [3]f32 = .{ 6, 9, 12 }; //Q
    var R: [3]f32 = try crsPrdXvec(&P, &Q);
    const dp1: f32 = dotPrdXvec(&R, &Q);
    const dp2: f32 = dotPrdXvec(&R, &P);
    std.debug.print("(P X Q) . Q = (P X Q)  . P = 0:\n", .{});
    try std.testing.expectEqual(dp1, dp2);
    try std.testing.expectEqual(dp1, 0);
    try std.testing.expectEqual(dp2, 0);
}

test "XMTX: MF3D - Lengyel: Theorem 2.9 test" {
    prntNl();
    var v1: [3]f32 = .{ 1, 3, 6 }; //P
    var v2: [3]f32 = .{ 6, 9, 12 }; //Q
    var v3: [3]f32 = .{ 2, 4, 8 }; //R
    v3 = try crsPrdXvec(&v2, &v1);
    var v4: [3]f32 = try crsPrdXvec(&v1, &v2);
    mulXvec(&v4, -1);
    std.debug.print("Q X P = -(P X Q):\n", .{});
    try std.testing.expectEqual(true, equXvec(&v3, &v4));

    v1 = .{ 1, 3, 6 }; //P
    v2 = .{ 6, 9, 12 }; //Q
    v3 = .{ 2, 4, 8 }; //R
    const a: f32 = 3;
    _ = cpyXvec(&v1, &v3);
    mulXvec(&v3, a);
    v3 = try crsPrdXvec(&v3, &v2);
    v4 = try crsPrdXvec(&v1, &v2);
    mulXvec(&v4, a);
    std.debug.print("(aP) X Q = a(P X Q):\n", .{});
    try std.testing.expectEqual(true, equXvec(&v3, &v4));

    v1 = .{ 1, 3, 6 }; //P
    v2 = .{ 6, 9, 12 }; //Q
    v3 = .{ 2, 4, 8 }; //R
    cpyXvec(&v2, &v4);
    sum1Xvec(&v4, &v3);
    v4 = try crsPrdXvec(&v1, &v4);
    var v5: [3]f32 = try crsPrdXvec(&v1, &v2);
    var v6: [3]f32 = try crsPrdXvec(&v1, &v3);
    sum1Xvec(&v5, &v6);
    prntXvecNl(&v4);
    prntXvecNl(&v5);
    std.debug.print("P X (Q + R) = (P X Q) + (P X R):\n", .{});
    try std.testing.expectEqual(true, equXvec(&v4, &v5));

    v1 = .{ 1, 3, 6 }; //P
    v2 = .{ 6, 9, 12 }; //Q
    v3 = .{ 2, 4, 8 }; //R
    v4 = try crsPrdXvec(&v1, &v1);
    v5 = .{ 0, 0, 0 };
    std.debug.print("P X P = [0,0,0]:\n", .{});
    try std.testing.expectEqual(true, equXvec(&v4, &v5));

    v1 = .{ 1, 3, 6 }; //P
    v2 = .{ 6, 9, 12 }; //Q
    v3 = .{ 2, 4, 8 }; //R
    v4 = try crsPrdXvec(&v1, &v2);
    const c1 = dotPrdXvec(&v4, &v3);
    v5 = try crsPrdXvec(&v3, &v1);
    const c2 = dotPrdXvec(&v5, &v2);
    v6 = try crsPrdXvec(&v2, &v3);
    const c3 = dotPrdXvec(&v6, &v1);
    std.debug.print("(P X Q) . R = (R X P) . Q = (Q X R) . P:\n", .{});
    try std.testing.expectEqual(true, (c1 == c2 and c1 == c3));

    v4 = try crsPrdXvec(&v2, &v1);
    v4 = try crsPrdXvec(&v1, &v4);
    v5 = try crsPrdXvec(&v1, &v2);
    v5 = try crsPrdXvec(&v5, &v1);

    const d: f32 = (magXvec(&v1) * magXvec(&v1));
    cpyXvec(&v2, &v6);
    mulXvec(&v6, d);
    const e = dotPrdXvec(&v1, &v2);
    var v7: [3]f32 = .{ 1, 1, 1 };
    cpyXvec(&v1, &v7);
    mulXvec(&v7, e);
    diff1Xvec(&v6, &v7);
    std.debug.print("P X (Q X P) = P X Q X P = (||P|| * ||P||) * Q - (P . Q) * P:\n", .{});
    try std.testing.expectEqual(true, equXvec(&v4, &v5));
    try std.testing.expectEqual(true, equXvec(&v4, &v6));
    try std.testing.expectEqual(true, equXvec(&v5, &v6));
}

test "XMTX: MF3D - Lengyel: Theorem 2.14 test" {
    prntNl();
    std.debug.print("Theorem 2.14:\n", .{});
    var v1: [2]f32 = .{ 1, 0 };
    var v2: [2]f32 = .{ 0, 1 };
    const dotP: f32 = dotPrdXvec(&v1, &v2);
    const val: f32 = 0;
    std.debug.print("(P . Q) = 0 if P and Q are linearly independent:\n", .{});
    try std.testing.expectEqual(val, dotP);
}

test "XMTX: MF3D - Lengyel: Theorem 3.1 test" {
    prntNl();
    var F: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var G: [9]f32 = .{ 1, 0, 1, 0, 1, 0, 1, 0, 1 };
    var FG: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var GF: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    cpyXmtx(&F, &FG);
    cpyXmtx(&G, &GF);
    sum1Xvec(&FG, &G);
    sum1Xvec(&GF, &F);
    std.debug.print("XMTX: F + G = G + F:\n", .{});
    try std.testing.expectEqual(true, equXmtx(&FG, &GF));

    F = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    G = .{ 1, 0, 1, 0, 1, 0, 1, 0, 1 };
    var H: [9]f32 = .{ 3, 1, 6, 7, 0, 5, 0, 9, 10 };
    FG = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var GH: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var FGH: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    cpyXmtx(&F, &FG);
    sum1Xvec(&FG, &G);
    sum1Xvec(&FG, &H);
    cpyXmtx(&G, &GH);
    sum1Xvec(&GH, &H);
    cpyXmtx(&F, &FGH);
    sum1Xvec(&FGH, &GH);
    std.debug.print("XMTX: (F + G) + H = F + (G + H):\n", .{});
    try std.testing.expectEqual(true, equXmtx(&FG, &FGH));

    const a: f32 = 13;
    const b: f32 = 21;
    F = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    G = .{ 1, 0, 1, 0, 1, 0, 1, 0, 1 };
    H = .{ 3, 1, 6, 7, 0, 5, 0, 9, 10 };
    FG = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    GH = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    FGH = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    cpyXmtx(&F, &FG);
    mulXvec(&FG, b);
    mulXvec(&FG, a);
    cpyXmtx(&F, &GH);
    mulXvec(&GH, (a * b));
    std.debug.print("XMTX: (a * (b * F)) = ((a * b) * F):\n", .{});
    try std.testing.expectEqual(true, equXmtx(&FG, &GH));

    F = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    G = .{ 1, 0, 1, 0, 1, 0, 1, 0, 1 };
    H = .{ 3, 1, 6, 7, 0, 5, 0, 9, 10 };
    FG = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    GH = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    FGH = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    cpyXmtx(&F, &FG);
    sum1Xvec(&FG, &G);
    mulXvec(&FG, a);
    cpyXmtx(&F, &GH);
    mulXvec(&GH, a);
    cpyXmtx(&G, &FGH);
    mulXvec(&FGH, a);
    sum1Xvec(&GH, &FGH);
    std.debug.print("XMTX: (a * (F + G)) = ((a * F) + (a * G)):\n", .{});
    try std.testing.expectEqual(true, equXmtx(&FG, &GH));

    const c: f32 = (a + b);
    F = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    G = .{ 1, 0, 1, 0, 1, 0, 1, 0, 1 };
    H = .{ 3, 1, 6, 7, 0, 5, 0, 9, 10 };
    FG = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    GH = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    FGH = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    cpyXmtx(&F, &FG);
    mulXvec(&FG, c);
    cpyXmtx(&F, &GH);
    mulXvec(&GH, a);
    cpyXmtx(&F, &FGH);
    mulXvec(&FGH, b);
    sum1Xvec(&GH, &FGH);
    std.debug.print("XMTX: ((a + b) * F) = ((a * F) + (b * F)):\n", .{});
    try std.testing.expectEqual(true, equXmtx(&FG, &GH));
}

test "XMTX: MF3D - Lengyel: Theorem 3.2 test" {
    prntNl();
    const a: f32 = 13;
    var F: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var G: [9]f32 = .{ 1, 0, 1, 0, 1, 0, 1, 0, 1 };
    var FG: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var GF: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var GH: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var RES: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var b: bool = false;
    cpyXmtx(&F, &FG);
    mulXvec(&FG, a);
    b = try tmsXmtx(&FG, 3, &G, 3, &RES, 3);
    try std.testing.expectEqual(true, b);
    cpyXmtx(&RES, &FG);

    RES = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    cpyXmtx(&F, &GF);
    cpyXmtx(&G, &GH);
    b = try tmsXmtx(&GF, 3, &GH, 3, &RES, 3);
    try std.testing.expectEqual(true, b);
    cpyXmtx(&RES, &GF);
    mulXvec(&GF, a);
    std.debug.print("XMTX: ((a * F) * G) = (a * (F * G)):\n", .{});
    try std.testing.expectEqual(true, equXmtx(&FG, &GF));

    F = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    G = .{ 1, 0, 1, 0, 1, 0, 1, 0, 1 };
    var H: [9]f32 = .{ 1, 2, 1, 2, 1, 3, 1, 3, 1 };
    FG = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    GF = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    GH = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    RES = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    b = try tmsXmtx(&F, 3, &G, 3, &RES, 3);
    try std.testing.expectEqual(true, b);
    cpyXmtx(&RES, &FG);

    RES = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    b = try tmsXmtx(&FG, 3, &H, 3, &RES, 3);
    try std.testing.expectEqual(true, b);
    cpyXmtx(&RES, &FG); //LEFT

    RES = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    b = try tmsXmtx(&G, 3, &H, 3, &RES, 3);
    try std.testing.expectEqual(true, b);
    cpyXmtx(&RES, &GH);

    RES = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    b = try tmsXmtx(&F, 3, &GH, 3, &RES, 3);
    try std.testing.expectEqual(true, b);
    cpyXmtx(&RES, &GH); //RIGHT
    std.debug.print("XMTX: ((F * G) * H) = (F * (G * H)):\n", .{});
    try std.testing.expectEqual(true, equXmtx(&FG, &GH));

    F = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    G = .{ 1, 0, 1, 0, 1, 0, 1, 0, 1 };
    FG = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    GF = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    GH = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    H = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    RES = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    b = try tmsXmtx(&F, 3, &G, 3, &RES, 3);
    try std.testing.expectEqual(true, b);
    cpyXmtx(&RES, &GH);
    trnXmtx(&GH, 3, &H); //Left

    RES = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    trnXmtx(&G, 3, &GF);
    trnXmtx(&F, 3, &FG);
    b = try tmsXmtx(&GF, 3, &FG, 3, &RES, 3);
    try std.testing.expectEqual(true, b);
    cpyXmtx(&RES, &GF); //Right
    std.debug.print("XMTX: ((F * G) ^ T) = (G^T * F^T):\n", .{});
    try std.testing.expectEqual(true, equXmtx(&H, &GF));
}

test "XMTX: MF3D - Lengyel: Theorem 3.9 test" {
    prntNl();
    var m1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 0, 0, 0 };
    var idtM1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var b: bool = isZeroXmtx(&m1, 3);
    try std.testing.expectEqual(true, b);

    var sclr: f32 = 0.0;
    b = try rdcXmtxInl(&m1, 3, false, false, &idtM1, 3, false, &sclr);
    try std.testing.expectEqual(false, b);

    b = isIdtXmtx(&m1, 3);
    try std.testing.expectEqual(false, b);

    //true
    m1 = .{ 2, 3, 8, 6, 0, -3, -1, 3, 2 };

    idtM1 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    b = isZeroXmtx(&m1, 3);
    prntXmtxNl(&m1, 3);
    prntNl();
    try std.testing.expectEqual(false, b);

    sclr = 0.0;
    b = try rdcXmtxInl(&m1, 3, false, false, &idtM1, 3, false, &sclr);
    prntXmtxNl(&m1, 3);
    prntNl();
    try std.testing.expectEqual(true, b);

    b = isIdtXmtx(&m1, 3);
    prntXmtxNl(&m1, 3);
    prntNl();
    try std.testing.expectEqual(true, b);
}

test "XMTX: MF3D - Lengyel: Theorem 3.10 test" {
    prntNl();
    //M has inverse iff transp(M) has inverse
    var origM1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 0, 0, 0 };
    var m1: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 0, 0, 0 };
    var trnM1: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var idtM1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var b: bool = false;

    //First test
    trnXmtx(&m1, 3, &trnM1);
    cpyXmtx(&trnM1, &m1);
    idtM1 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    b = isZeroXmtx(&m1, 3);
    prntXmtxNl(&m1, 3);
    prntNl();
    try std.testing.expectEqual(true, b);

    var sclr: f32 = 0.0;
    b = try rdcXmtxInl(&m1, 3, false, false, &idtM1, 3, false, &sclr);
    prntXmtxNl(&m1, 3);
    prntNl();
    try std.testing.expectEqual(false, b);

    b = isIdtXmtx(&m1, 3);
    prntXmtxNl(&m1, 3);
    prntNl();
    try std.testing.expectEqual(false, b);

    //Second test
    cpyXmtx(&origM1, &m1);
    idtM1 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    b = isZeroXmtx(&m1, 3);
    prntXmtxNl(&m1, 3);
    prntNl();
    try std.testing.expectEqual(true, b);

    sclr = 0.0;
    b = try rdcXmtxInl(&m1, 3, false, false, &idtM1, 3, false, &sclr);
    prntXmtxNl(&m1, 3);
    prntNl();
    try std.testing.expectEqual(false, b);

    b = isIdtXmtx(&m1, 3);
    prntXmtxNl(&m1, 3);
    prntNl();
    try std.testing.expectEqual(false, b);

    m1 = .{ 2, 3, 8, 6, 0, -3, -1, 3, 2 };

    origM1 = .{ 2, 3, 8, 6, 0, -3, -1, 3, 2 };

    trnM1 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    idtM1 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    b = false;

    //First test
    trnXmtx(&m1, 3, &trnM1);
    cpyXmtx(&trnM1, &m1);
    idtM1 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    b = isZeroXmtx(&m1, 3);
    prntXmtxNl(&m1, 3);
    prntNl();
    try std.testing.expectEqual(false, b);

    sclr = 0.0;
    b = try rdcXmtxInl(&m1, 3, false, false, &idtM1, 3, false, &sclr);
    prntXmtxNl(&m1, 3);
    prntNl();
    try std.testing.expectEqual(true, b);

    b = isIdtXmtx(&m1, 3);
    prntXmtxNl(&m1, 3);
    prntNl();
    try std.testing.expectEqual(true, b);

    //Second test
    cpyXmtx(&origM1, &m1);
    idtM1 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    b = isZeroXmtx(&m1, 3);
    prntXmtxNl(&m1, 3);
    prntNl();
    try std.testing.expectEqual(false, b);

    sclr = 0.0;
    b = try rdcXmtxInl(&m1, 3, false, false, &idtM1, 3, false, &sclr);
    prntXmtxNl(&m1, 3);
    prntNl();
    try std.testing.expectEqual(true, b);

    b = isIdtXmtx(&m1, 3);
    prntXmtxNl(&m1, 3);
    prntNl();
    try std.testing.expectEqual(true, b);
}

test "XMTX: MF3D - Lengyel: Theorem 3.11 test" {
    prntNl();
    //if F and G are n X n invertible matrices then (F * G) is invertible and (F * G) ^ (-1) = G^(-1) * F^(-1)
    var origF: [9]f32 = .{ 2, 3, 8, 6, 0, -3, -1, 3, 2 };

    var F: [9]f32 = .{ 2, 3, 8, 6, 0, -3, -1, 3, 2 };

    var origG: [9]f32 = .{ 8, 4, 3, -5, 6, -2, 7, 9, -8 };

    var G: [9]f32 = .{ 8, 4, 3, -5, 6, -2, 7, 9, -8 };

    var FG: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var invFG: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var invF: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var invG: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var invGinvF: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var b: bool = false;

    b = try tmsXmtx(&F, 3, &G, 3, &FG, 3);
    std.debug.print("BBB FG\n", .{});
    prntXmtxNl(&FG, 3);
    prntNl();
    try std.testing.expectEqual(true, b);

    var sclr: f32 = 0.0;
    b = try rdcXmtxInl(&FG, 3, false, true, &invFG, 3, false, &sclr);
    std.debug.print("BBB Reduced FG\n", .{});
    prntXmtxNl(&FG, 3);
    prntNl();
    std.debug.print("BBB invFG\n", .{});
    prntXmtxNl(&invFG, 3);
    prntNl();
    try std.testing.expectEqual(true, b);

    sclr = 0.0;
    b = try rdcXmtxInl(&F, 3, false, true, &invF, 3, false, &sclr);
    std.debug.print("CCC Reduced F\n", .{});
    prntXmtxNl(&F, 3);
    prntNl();
    std.debug.print("CCC invF\n", .{});
    prntXmtxNl(&invF, 3);
    prntNl();
    cpyXmtx(&origF, &F);
    try std.testing.expectEqual(true, b);

    sclr = 0.0;
    b = try rdcXmtxInl(&G, 3, false, true, &invG, 3, false, &sclr);
    std.debug.print("DDD Reduced G\n", .{});
    prntXmtxNl(&G, 3);
    prntNl();
    std.debug.print("DDD invG\n", .{});
    prntXmtxNl(&invG, 3);
    prntNl();
    cpyXmtx(&origG, &G);
    try std.testing.expectEqual(true, b);

    clrXmtx(&invGinvF);
    b = try tmsXmtx(&invG, 3, &invF, 3, &invGinvF, 3);
    std.debug.print("EEE invG * invF\n", .{});
    prntXmtxNl(&invGinvF, 3);
    prntNl();
    try std.testing.expectEqual(true, b);

    clnXmtx(&invFG);
    clnXmtx(&invGinvF);

    //if F and G are n X n invertible matrices then (F * G) is invertible and (F * G) ^ (-1) = G^(-1) * F^(-1)
    std.debug.print("FIN invFG\n", .{});
    prntXmtxNl(&invFG, 3);
    prntNl();
    std.debug.print("FIN invGinvF\n", .{});
    prntXmtxNl(&invGinvF, 3);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&invFG, &invGinvF));
}

test "XMTX: MF3D - Lengyel: Theorem 3.15 test" {
    prntNl();
    //matrix M has an inverse if rows of M are linrarly independent
    //First Test
    var origM: [9]f32 = .{ 2, 3, 8, 6, 0, -3, -1, 3, 2 };

    var M: [9]f32 = .{ 2, 3, 8, 6, 0, -3, -1, 3, 2 };

    var origInvM: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var invM: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var b: bool = false;

    var sclr: f32 = 0.0;
    b = try rdcXmtxInl(&M, 3, false, true, &invM, 3, false, &sclr);
    std.debug.print("AAA Reduced M\n", .{});
    prntXmtxNl(&M, 3);
    prntNl();
    std.debug.print("AAA invM\n", .{});
    prntXmtxNl(&invM, 3);
    prntNl();
    try std.testing.expectEqual(true, b);

    b = isIdtXmtx(&M, 3);
    try std.testing.expectEqual(true, b);

    b = idnfXmtx(&M, 3.0, MTX_OPS.MTX_IS_LIN_INDP);
    try std.testing.expectEqual(true, b);

    //Second Test
    origM = .{ 2, 3, 8, 4, 6, 16, -1, 3, 2 };

    M = .{ 2, 3, 8, 4, 6, 16, -1, 3, 2 };

    origInvM = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    invM = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    b = false;

    sclr = 0.0;
    b = try rdcXmtxInl(&M, 3, false, true, &invM, 3, false, &sclr);
    std.debug.print("BBB Reduced M\n", .{});
    prntXmtxNl(&M, 3);
    prntNl();
    std.debug.print("BBB invM\n", .{});
    prntXmtxNl(&invM, 3);
    prntNl();
    try std.testing.expectEqual(false, b);

    b = isIdtXmtx(&M, 3);
    try std.testing.expectEqual(false, b);

    b = idnfXmtx(&M, 3.0, MTX_OPS.MTX_IS_LIN_INDP);
    try std.testing.expectEqual(false, b);
}

test "XMTX: MF3D - Lengyel: Theorem 3.17 test" {
    prntNl();
    const detRow: usize = 0;
    var m2: [4]f32 = .{ 5, 6, 8, 9 };
    const cols: usize = 2;
    const alloc: std.mem.Allocator = std.testing.allocator;
    var v: f32 = try detXmtx(&m2, cols, &alloc, detRow);
    var exp: f32 = -3.0;
    std.debug.print("Found detXmtx 2x2 value: {}\n", .{v});
    try std.testing.expectEqual(exp, v);

    m2 = .{ 5, 6, 8, 9 };
    var m3: [4]f32 = .{ 0, 0, 0, 0 };
    addSclMulXmtxRows(0, 1, 2, 2, &m2, &m3);
    v = try detXmtx(&m3, cols, &alloc, detRow);
    exp = -3.0;
    std.debug.print("Found detXmtx 2x2 after add scalar multiple copy of row: {}\n", .{v});
    try std.testing.expectEqual(exp, v);

    m2 = .{ 5, 6, 8, 9 };
    m3 = .{ 0, 0, 0, 0 };
    sclMulXmtxRows(0, 2, cols, &m2, &m3);
    v = try detXmtx(&m3, cols, &alloc, detRow);
    exp = (-3.0 * 2);
    std.debug.print("Found detXmtx 2x2 after add scalar multiple to row: {}\n", .{v});
    try std.testing.expectEqual(exp, v);

    m2 = .{ 5, 6, 8, 9 };
    m3 = .{ 0, 0, 0, 0 };
    altXmtxRows(0, 1, cols, &m2, &m3);
    v = try detXmtx(&m3, cols, &alloc, detRow);
    exp = 3.0;
    std.debug.print("Found detXmtx 2x2 after alternating two rows: {}\n", .{v});
    try std.testing.expectEqual(exp, v);
}

test "XMTX: MF3D - Lengyel: Corollary 3.18 test" {
    prntNl();
    const detRow: usize = 0;
    var m2: [4]f32 = .{ 5, 6, 5, 6 };
    const cols: usize = 2;
    const alloc: std.mem.Allocator = std.testing.allocator;
    const v: f32 = try detXmtx(&m2, cols, &alloc, detRow);
    const exp: f32 = 0.0;
    std.debug.print("Found detXmtx 2x2 value: {}\n", .{v});
    try std.testing.expectEqual(exp, v);
}

test "XMTX: MF3D - Lengyel: Theorem 3.19 test" {
    prntNl();
    //An n X n matrix M is invertible iff detM != 0
    const detRow: usize = 0;
    var m2: [4]f32 = .{ 5, 6, 5, 6 };
    var cols: usize = 2;
    const alloc: std.mem.Allocator = std.testing.allocator;
    var v: f32 = try detXmtx(&m2, cols, &alloc, detRow);
    var exp: f32 = 0.0;
    std.debug.print("AAA Found detXmtx 2x2 value: {}\n", .{v});
    try std.testing.expectEqual(exp, v);

    m2 = .{ 5, 6, 5, 6 };
    var m3: [4]f32 = .{ 1, 0, 0, 1 };
    var sclr: f32 = 0.0;
    var b: bool = try rdcXmtxInl(&m2, cols, false, true, &m3, 2, false, &sclr);
    try std.testing.expectEqual(false, b);

    b = isIdtXmtx(&m2, 2);
    try std.testing.expectEqual(false, b);

    m2 = .{ 5, 6, 8, 9 };
    cols = 2;
    v = try detXmtx(&m2, cols, &alloc, detRow);
    exp = -3.0;
    std.debug.print("BBB Found detXmtx 2x2 value: {}\n", .{v});
    try std.testing.expectEqual(exp, v);

    m2 = .{ 5, 6, 8, 9 };
    m3 = .{ 1, 0, 0, 1 };
    sclr = 0.0;
    b = try rdcXmtxInl(&m2, cols, false, true, &m3, 2, false, &sclr);
    try std.testing.expectEqual(true, b);

    prntNl();
    prntXmtxNl(&m2, 2);

    clnXmtx(&m2);
    prntXmtxNl(&m2, 2);

    b = isIdtXmtx(&m2, 2);
    try std.testing.expectEqual(true, b);
}

test "XMTX: MF3D - Lengyel: Theorem 3.20 test" {
    prntNl();
    //det(FG) = det(F)det(G)
    const detRow: usize = 0;
    const cols: usize = 3;
    var F: [9]f32 = .{ 1, 2, 3, 0, 1, 4, 5, 6, 0 };

    var G: [9]f32 = .{ 4, 6, 5, 1, 2, 3, 0, 4, 6 };

    var FG: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    const alloc: std.mem.Allocator = std.testing.allocator;
    const b: bool = try tmsXmtx(&F, cols, &G, cols, &FG, cols);
    try std.testing.expectEqual(true, b);

    const detFG: f32 = try detXmtx(&FG, cols, &alloc, detRow);
    const detF: f32 = try detXmtx(&F, cols, &alloc, detRow);
    const detG: f32 = try detXmtx(&G, cols, &alloc, detRow);
    try std.testing.expectEqual(detFG, (detF * detG));
}

test "XMTX: MF3D - Lengyel: Theorem 3.21 test" {
    prntNl();
    //Let F be an n X n matrix and let the entries of n X n matrix G be defined as
    //Gij = Cji(F) * (1 / detF)
    //where Cji(F) is the cofactor of (F^Tij) then G = F^-1
    //Gij = Cij(F^T) * (1 / detF)

    //Cij(H) = (-1)^(i + j)detH
    std.debug.print("Test 1:\n", .{});
    const alloc: std.mem.Allocator = std.testing.allocator;
    var F: [4]f32 = .{ 5, 6, 8, 9 };
    var invF: [4]f32 = .{ 0, 0, 0, 0 };
    var cols: usize = 2;

    std.debug.print("F Matrix:\n", .{});
    clnXmtx(&F);
    prntXmtxNl(&F, cols);

    const detF = try detXmtx(&F, cols, &alloc, 0);
    std.debug.print("detF = {}\n", .{detF});

    F = .{ 5, 6, 8, 9 };
    var idtF: [4]f32 = .{ 1, 0, 0, 1 };
    var sclr: f32 = 0.0;
    var b: bool = try rdcXmtxInl(&F, cols, false, true, &idtF, 2, false, &sclr);
    try std.testing.expectEqual(true, b);

    std.debug.print("Calculated inverse F (should be identity matrix):\n", .{});
    clnXmtx(&F);
    prntXmtxNl(&F, cols);
    try std.testing.expectEqual(true, isIdtXmtx(&F, cols));

    std.debug.print("Calculated inverse F (should be inverse matrix):\n", .{});
    clnXmtx(&idtF);
    prntXmtxNl(&idtF, cols);

    F = .{ 5, 6, 8, 9 };
    getInvFromDet2(&F, detF, &invF);

    std.debug.print("Generated inverse from detF (should be inverse matrix):\n", .{});
    clnXmtx(&invF);
    prntXmtxNl(&invF, cols);
    try std.testing.expectEqual(true, equXmtx(&idtF, &invF));

    std.debug.print("Test 2:\n", .{});
    var F2: [9]f32 = .{ 2, 3, 8, 6, 0, -3, -1, 3, 2 };
    var invF2: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    cols = 3;

    std.debug.print("F2 Matrix:\n", .{});
    clnXmtx(&F2);
    prntXmtxNl(&F2, cols);

    const detF2 = try detXmtx(&F2, cols, &alloc, 0);
    std.debug.print("detF2 = {}\n", .{detF2});

    F2 = .{ 2, 3, 8, 6, 0, -3, -1, 3, 2 };
    var idtF2: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    sclr = 0.0;
    b = try rdcXmtxInl(&F2, cols, false, true, &idtF2, 3, false, &sclr);
    try std.testing.expectEqual(true, b);

    std.debug.print("Calculated inverse F2 (should be identity matrix):\n", .{});
    clnXmtx(&F2);
    prntXmtxNl(&F2, cols);
    try std.testing.expectEqual(true, isIdtXmtx(&F2, cols));

    std.debug.print("Calculated inverse F2 (should be inverse matrix):\n", .{});
    clnXmtx(&idtF2);
    prntXmtxNl(&idtF2, cols);

    F2 = .{ 2, 3, 8, 6, 0, -3, -1, 3, 2 };
    getInvFromDet3(&F2, detF2, &invF2);

    std.debug.print("Generated inverse from detF2 (should be inverse matrix):\n", .{});
    clnXmtx(&invF2);
    prntXmtxNl(&invF2, cols);
    try std.testing.expectEqual(true, equXmtx(&idtF2, &invF2));
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//STOP THEOREMS
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

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
//START PROBLEMS
//Elementary Linear Algebra - Larson, Edwards- 4th Edition
//Chapters 1.0 - 5.1
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

test "XMTX: ELA - Larson, Edwards: 1.2 Example 3 test" {
    prntNl();
    var m1: [12]f32 = .{ 1, -2, 3, 9, -1, 3, 0, -4, 2, -5, 5, 17 };

    var retM1: [12]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    var idtM1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    const dim: usize = 3; //overrides the has augment column difference of 1 and controls the zero row check
    const cols: usize = 4;
    var b: bool = false; //holds the result of the operation

    const hasAug: bool = true; //toggles isAugmented flag for the reduction function
    const hasIdt: bool = true; //indicates an matrix was provided to calculate and hold the inverse of m1.
    const triag: bool = false; //A Boolean value indicating if the reduction operation should stop when the matrix is triangular.

    var sclr: f32 = 0.0;
    b = try rdcXmtx(&m1, cols, hasAug, &retM1, hasIdt, &idtM1, dim, triag, &sclr);
    try std.testing.expectEqual(true, b);

    std.debug.print("Matrix M1:\n", .{});
    prntXmtxNl(&m1, cols);
    prntNl();

    std.debug.print("Matrix Ret:\n", .{});
    prntXmtxNl(&retM1, cols);
    prntNl();

    std.debug.print("Matrix Inv:\n", .{});
    prntXmtxNl(&idtM1, dim);
    prntNl();

    cpyLessXmtx(&retM1, &idtM1, cols, dim);
    std.debug.print("Copy Matrix M1:\n", .{});
    clnXmtx(&idtM1);
    prntXmtxNl(&idtM1, dim);
    prntNl();
    try std.testing.expectEqual(true, isIdtXmtx(&idtM1, dim));
}

test "XMTX: ELA - Larson, Edwards: 1.2 Example 4 test" {
    prntNl();
    var m1: [20]f32 = .{ 0, 1, 1, -2, -3, 1, 2, -1, 0, 2, 2, 4, 1, -3, -2, 1, -4, -7, -1, -19 };

    var retM1: [20]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    var idtM1: [16]f32 = .{ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };

    const dim: usize = 4; //overrides the has augment column difference of 1 and controls the zero row check
    const cols: usize = 5;
    var b: bool = false;
    var sclr: f32 = 0.0;
    b = try rdcXmtx(&m1, cols, true, &retM1, true, &idtM1, dim, false, &sclr);
    try std.testing.expectEqual(true, b);

    std.debug.print("Matrix M1:\n", .{});
    prntXmtxNl(&m1, cols);
    prntNl();

    std.debug.print("Matrix Ret:\n", .{});
    prntXmtxNl(&retM1, cols);
    prntNl();

    std.debug.print("Matrix Inv:\n", .{});
    prntXmtxNl(&idtM1, dim);
    prntNl();

    cpyLessXmtx(&retM1, &idtM1, cols, dim);
    std.debug.print("Copy Matrix M1:\n", .{});
    prntXmtxNl(&idtM1, dim);
    clnXmtx(&idtM1);
    prntNl();
    try std.testing.expectEqual(true, isIdtXmtx(&idtM1, dim));
}

test "XMTX: ELA - Larson, Edwards: 1.2 Example 6 test" {
    prntNl();
    var m1: [16]f32 = .{ 1, -1, 2, 4, 1, 0, 1, 6, 2, -3, 5, 4, 3, 2, -1, 1 };

    var retM1: [16]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    var idtM1: [16]f32 = .{ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };

    const dim: usize = 3; //the number of columns corresponding to variables when square and augmented = false
    const cols: usize = 4;
    var b: bool = false;
    var sclr: f32 = 0.0;
    b = try rdcXmtx(&m1, cols, false, &retM1, true, &idtM1, dim, false, &sclr);
    try std.testing.expectEqual(false, b);

    std.debug.print("Matrix M1:\n", .{});
    prntXmtxNl(&m1, cols);
    prntNl();

    std.debug.print("Matrix Ret:\n", .{});
    prntXmtxNl(&retM1, cols);
    prntNl();

    std.debug.print("Matrix Inv:\n", .{});
    prntXmtxNl(&idtM1, cols);
    prntNl();

    try std.testing.expectEqual(false, isIdtXmtx(&m1, cols));
}

test "XMTX: ELA - Larson, Edwards: 1.2 Example 7 test" {
    prntNl();
    var m1: [12]f32 = .{ 1, -2, 3, 9, -1, 3, 0, -4, 2, -5, 5, 17 };

    var retM1: [12]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    var idtM1: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    const dim: usize = 3; //the number of columns corresponding to variables when square and augmented = false
    const cols: usize = 4;
    var b: bool = false;
    var sclr: f32 = 0.0;
    b = try rdcXmtx(&m1, cols, true, &retM1, true, &idtM1, dim, false, &sclr);
    try std.testing.expectEqual(true, b);

    std.debug.print("Matrix M1:\n", .{});
    prntXmtxNl(&m1, cols);
    prntNl();

    std.debug.print("Matrix Ret:\n", .{});
    prntXmtxNl(&retM1, cols);
    prntNl();

    std.debug.print("Matrix Inv:\n", .{});
    prntXmtxNl(&idtM1, dim);
    prntNl();

    cpyLessXmtx(&retM1, &idtM1, cols, dim);
    std.debug.print("Copy Matrix M1:\n", .{});
    prntXmtxNl(&idtM1, dim);
    clnXmtx(&idtM1);
    prntNl();
    try std.testing.expectEqual(true, isIdtXmtx(&idtM1, dim));
}

test "XMTX: ELA - Larson, Edwards: 1.2 Example 8 test" {
    prntNl();
    var m1: [8]f32 = .{ 2, 4, -2, 0, 3, 5, 0, 1 };

    var retM1: [8]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0 };

    var idtM1: [4]f32 = .{ 1, 0, 0, 1 };

    const dim: usize = 2; //overrides the has augment column difference of 1 and controls the zero row check
    const cols: usize = 4;
    var b: bool = false;
    var sclr: f32 = 0.0;
    b = try rdcXmtx(&m1, cols, true, &retM1, true, &idtM1, dim, false, &sclr);
    try std.testing.expectEqual(true, b);

    std.debug.print("Matrix M1:\n", .{});
    prntXmtxNl(&m1, cols);
    prntNl();

    std.debug.print("Matrix Ret:\n", .{});
    prntXmtxNl(&retM1, cols);
    prntNl();

    std.debug.print("Matrix Inv:\n", .{});
    prntXmtxNl(&idtM1, dim);
    prntNl();

    cpyLessXmtx(&retM1, &idtM1, cols, dim);
    std.debug.print("Copy Matrix M1:\n", .{});
    prntXmtxNl(&idtM1, dim);
    clnXmtx(&idtM1);
    prntNl();
    try std.testing.expectEqual(true, isIdtXmtx(&idtM1, dim));
}

test "XMTX: ELA - Larson, Edwards: 1.2 Example 9 test" {
    prntNl();
    var m1: [8]f32 = .{ 1, -1, 3, 0, 2, 1, 3, 0 };

    var retM1: [8]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0 };

    var idtM1: [4]f32 = .{ 1, 0, 0, 1 };

    const dim: usize = 2; //overrides the has augment column difference of 1 and controls the zero row check
    const cols: usize = 4;
    var b: bool = false;
    var sclr: f32 = 0.0;
    b = try rdcXmtx(&m1, cols, true, &retM1, true, &idtM1, dim, false, &sclr);
    try std.testing.expectEqual(true, b);

    std.debug.print("Matrix M1:\n", .{});
    prntXmtxNl(&m1, cols);
    prntNl();

    std.debug.print("Matrix Ret:\n", .{});
    prntXmtxNl(&retM1, cols);
    prntNl();

    std.debug.print("Matrix Inv:\n", .{});
    prntXmtxNl(&idtM1, dim);
    prntNl();

    cpyLessXmtx(&retM1, &idtM1, cols, dim);
    std.debug.print("Copy Matrix M1:\n", .{});
    prntXmtxNl(&idtM1, dim);
    clnXmtx(&idtM1);
    prntNl();
    try std.testing.expectEqual(true, isIdtXmtx(&idtM1, dim));
}

test "XMTX: ELA - Larson, Edwards: 1.2 Problem 19 test" {
    prntNl();
    var m1: [6]f32 = .{ 1, 2, 7, 2, 1, 8 };

    var retM1: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };

    var idtM1: [4]f32 = .{ 1, 0, 0, 1 };

    const dim: usize = 2; //overrides the has augment column difference of 1 and controls the zero row check
    const cols: usize = 3;
    var b: bool = false;
    var sclr: f32 = 0.0;
    b = try rdcXmtx(&m1, cols, true, &retM1, true, &idtM1, dim, false, &sclr);
    try std.testing.expectEqual(true, b);

    std.debug.print("Matrix M1:\n", .{});
    prntXmtxNl(&m1, cols);
    prntNl();

    std.debug.print("Matrix Ret:\n", .{});
    prntXmtxNl(&retM1, cols);
    prntNl();

    std.debug.print("Matrix Inv:\n", .{});
    prntXmtxNl(&idtM1, dim);
    prntNl();

    cpyLessXmtx(&retM1, &idtM1, cols, dim);
    std.debug.print("Copy Matrix M1:\n", .{});
    prntXmtxNl(&idtM1, dim);
    clnXmtx(&idtM1);
    prntNl();
    try std.testing.expectEqual(true, isIdtXmtx(&idtM1, dim));
}

test "XMTX: ELA - Larson, Edwards: 1.2 Problem 20 test" {
    prntNl();
    var m1: [6]f32 = .{ 2, 6, 16, -2, -6, -16 };

    var retM1: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };

    var idtM1: [4]f32 = .{ 1, 0, 0, 1 };

    const dim: usize = 2; //overrides the has augment column difference of 1 and controls the zero row check
    const cols: usize = 3;
    var b: bool = false;
    var sclr: f32 = 0.0;
    b = try rdcXmtx(&m1, cols, true, &retM1, true, &idtM1, dim, false, &sclr);
    try std.testing.expectEqual(true, b);

    std.debug.print("Matrix M1:\n", .{});
    prntXmtxNl(&m1, cols);
    prntNl();

    std.debug.print("Matrix Ret:\n", .{});
    prntXmtxNl(&retM1, cols);
    prntNl();

    std.debug.print("Matrix Inv:\n", .{});
    prntXmtxNl(&idtM1, dim);
    prntNl();

    cpyLessXmtx(&retM1, &idtM1, cols, dim);
    std.debug.print("Copy Matrix M1:\n", .{});
    prntXmtxNl(&idtM1, dim);
    clnXmtx(&idtM1);
    prntNl();
    try std.testing.expectEqual(false, isIdtXmtx(&idtM1, dim));
}

test "XMTX: ELA - Larson, Edwards: 2.1 Example 2 test" {
    prntNl();
    var m1: [4]f32 = .{ -1, 2, 0, 1 };
    var m2: [4]f32 = .{ 1, 3, -1, 2 };
    var m3: [4]f32 = .{ 0, 0, 0, 0 };
    sum2Xvec(&m3, &m1, &m2);
    var exp1: [4]f32 = .{ 0, 5, -1, 3 };
    try std.testing.expectEqual(true, equXmtx(&exp1, &m3));

    var m4: [6]f32 = .{ 0, 1, -2, 1, 2, 3 };
    var m5: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    var m6: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    sum2Xvec(&m6, &m4, &m5);
    var exp2: [6]f32 = .{ 0, 1, -2, 1, 2, 3 };
    try std.testing.expectEqual(true, equXmtx(&exp2, &m6));

    var m7: [3]f32 = .{ 1, -3, -2 };
    var m8: [3]f32 = .{ -1, 3, 2 };
    var m9: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&m9, &m7, &m8);
    var exp3: [3]f32 = .{ 0, 0, 0 };
    try std.testing.expectEqual(true, equXmtx(&exp3, &m9));
}

test "XMTX: ELA - Larson, Edwards: 2.1 Example 3 test" {
    prntNl();
    var m1: [9]f32 = .{ 1, 2, 4, -3, 0, -1, 2, 1, 2 };
    var m2: [9]f32 = .{ 2, 0, 0, 1, -4, 3, -1, 3, 2 };

    //3A
    var exp1: [9]f32 = .{ 3, 6, 12, -9, 0, -3, 6, 3, 6 };
    mulXvec(&m1, 3);
    try std.testing.expectEqual(true, equXmtx(&exp1, &m1));

    //-B
    exp1 = .{ -2, 0, 0, -1, 4, -3, 1, -3, -2 };
    clnXmtx(&m2);
    mulXvec(&m2, -1);
    try std.testing.expectEqual(true, equXmtx(&exp1, &m2));

    //3A - B
    m1 = .{ 1, 2, 4, -3, 0, -1, 2, 1, 2 };
    m2 = .{ 2, 0, 0, 1, -4, 3, -1, 3, 2 };
    mulXvec(&m1, 3);
    diff1Xvec(&m1, &m2);
    exp1 = .{ 1, 6, 12, -10, 4, -6, 7, 0, 4 };
    clnXmtx(&m1);
    try std.testing.expectEqual(true, equXmtx(&exp1, &m1));
}

test "XMTX: ELA - Larson, Edwards: 2.1 Example 4, E5 test" {
    prntNl();
    //A
    var m1: [6]f32 = .{ -1, 3, 4, -2, 5, 0 };
    var m2: [4]f32 = .{ -3, 2, -4, 1 };
    var ret1: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    var exp1: [6]f32 = .{ -9, 1, -4, 6, -15, 10 };
    var b: bool = false;
    b = try tmsXmtx(&m1, 2, &m2, 2, &ret1, 2);
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&exp1, &ret1));

    //B
    var m3: [4]f32 = .{ 3, 4, -2, 5 };
    var m4: [4]f32 = .{ 1, 0, 0, 1 };
    var ret2: [4]f32 = .{ 0, 0, 0, 0 };
    var exp2: [4]f32 = .{ 3, 4, -2, 5 };
    b = try tmsXmtx(&m3, 2, &m4, 2, &ret2, 2);
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&exp2, &ret2));

    //C
    m3 = .{ 1, 2, 1, 1 };
    m4 = .{ -1, 2, 1, -1 };
    ret2 = .{ 0, 0, 0, 0 };
    exp2 = .{ 1, 0, 0, 1 };
    b = try tmsXmtx(&m3, 2, &m4, 2, &ret2, 2);
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&exp2, &ret2));

    //D
    var m5: [3]f32 = .{ 1, -2, -3 };
    var m6: [3]f32 = .{ 2, -1, 1 };
    var ret3: [1]f32 = .{0};
    var exp3: [1]f32 = .{1};
    std.debug.print("\nQuestion D:\n", .{});
    b = try tmsXmtx(&m5, 3, &m6, 1, &ret3, 1);

    prntXmtxNl(&m5, 3);
    prntNl();

    prntXmtxNl(&m6, 1);
    prntNl();

    prntXmtxNl(&ret3, 1);
    prntNl();

    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&exp3, &ret3));

    //E
    m5 = .{ 1, -2, -3 };
    m6 = .{ 2, -1, 1 };
    var ret4: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var exp4: [9]f32 = .{ 2, -4, -6, -1, 2, 3, 1, -2, -3 };
    std.debug.print("\nQuestion E:\n", .{});
    b = try tmsXmtx(&m6, 1, &m5, 3, &ret4, 3);

    prntXmtxNl(&m5, 3);
    prntNl();

    prntXmtxNl(&m6, 1);
    prntNl();

    prntXmtxNl(&ret4, 3);
    prntNl();

    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&exp4, &ret4));
}

test "XMTX: ELA - Larson, Edwards: 2.1 Example 6 test" {
    prntNl();
    var m1: [8]f32 = .{ 1, -2, 1, 0, 2, 3, -2, 0 };
    var ret: [8]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0 };
    var idtM1: [4]f32 = .{ 1, 0, 0, 1 };
    var m2: [4]f32 = .{ 0, 0, 0, 0 };
    var b: bool = false;
    var sclr: f32 = 0.0;
    b = try rdcXmtx(&m1, 4, true, &ret, false, &idtM1, 2, false, &sclr);

    prntXmtxNl(&m1, 4);
    prntNl();

    prntXmtxNl(&ret, 4);
    prntNl();

    prntXmtxNl(&idtM1, 2);
    prntNl();

    cpyLessXmtx(&ret, &m2, 4, 2);

    prntXmtxNl(&m2, 2);
    prntNl();

    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&m2, &idtM1));
}

test "XMTX: ELA - Larson, Edwards: 2.1 Problem 1 test" {
    prntNl();
    var a: [4]f32 = .{ 1, -1, 2, -1 };
    var b: [4]f32 = .{ 2, -1, -1, 8 };
    var aPb: [4]f32 = .{ 0, 0, 0, 0 };
    var exp: [4]f32 = .{ 3, -2, 1, 7 };

    //Find A + B, A - B, 2A, and 2A - B
    //A + B
    sum2Xvec(&aPb, &a, &b);
    try std.testing.expectEqual(true, equXmtx(&aPb, &exp));

    //A - B
    a = .{ 1, -1, 2, -1 };
    b = .{ 2, -1, -1, 8 };
    aPb = .{ 0, 0, 0, 0 };
    exp = .{ -1, 0, 3, -9 };
    diff2Xvec(&aPb, &a, &b);

    prntXmtxNl(&aPb, 2);
    prntNl();

    prntXmtxNl(&exp, 2);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&aPb, &exp));

    //2A
    a = .{ 1, -1, 2, -1 };
    b = .{ 2, -1, -1, 8 };
    aPb = .{ 0, 0, 0, 0 };
    exp = .{ 2, -2, 4, -2 };
    cpyXvec(&a, &aPb);
    mulXvec(&aPb, 2);
    try std.testing.expectEqual(true, equXmtx(&aPb, &exp));

    //2A - B
    exp = .{ 0, -1, 5, -10 };
    diff1Xvec(&aPb, &b);
    try std.testing.expectEqual(true, equXmtx(&aPb, &exp));
}

test "XMTX: ELA - Larson, Edwards: 2.1 Problem 3 test" {
    prntNl();
    var a: [6]f32 = .{ 6, -1, 2, 4, -3, 5 };
    var b: [6]f32 = .{ 1, 4, -1, 5, 1, 10 };
    var aPb: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    var exp: [6]f32 = .{ 7, 3, 1, 9, -2, 15 };

    //Find A + B, A - B, 2A, and 2A - B
    //A + B
    sum2Xvec(&aPb, &a, &b);
    try std.testing.expectEqual(true, equXmtx(&aPb, &exp));

    //A - B
    a = .{ 6, -1, 2, 4, -3, 5 };
    b = .{ 1, 4, -1, 5, 1, 10 };
    aPb = .{ 0, 0, 0, 0, 0, 0 };
    exp = .{ 5, -5, 3, -1, -4, -5 };
    diff2Xvec(&aPb, &a, &b);

    prntXmtxNl(&aPb, 2);
    prntNl();

    prntXmtxNl(&exp, 2);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&aPb, &exp));

    //2A
    a = .{ 6, -1, 2, 4, -3, 5 };
    b = .{ 1, 4, -1, 5, 1, 10 };
    aPb = .{ 0, 0, 0, 0, 0, 0 };
    exp = .{ 12, -2, 4, 8, -6, 10 };
    cpyXvec(&a, &aPb);
    mulXvec(&aPb, 2);
    try std.testing.expectEqual(true, equXmtx(&aPb, &exp));

    //2A - B
    exp = .{ 11, -6, 5, 3, -7, 0 };
    diff1Xvec(&aPb, &b);
    try std.testing.expectEqual(true, equXmtx(&aPb, &exp));
}

test "XMTX: ELA - Larson, Edwards: 2.1 Problem 9 test" {
    prntNl();
    //Find AB and BA
    //A = 1, 2, 4, 2
    //B = 2, -1, -1, 8

    var a: [4]f32 = .{ 1, 2, 4, 2 };
    var b: [4]f32 = .{ 2, -1, -1, 8 };
    var aTb: [4]f32 = .{ 0, 0, 0, 0 };
    var exp: [4]f32 = .{ 0, 15, 6, 12 };
    var res: bool = false;
    res = try tmsXmtx(&a, 2, &b, 2, &aTb, 2);
    try std.testing.expectEqual(true, res);
    try std.testing.expectEqual(true, equXmtx(&aTb, &exp));

    aTb = .{ 0, 0, 0, 0 };
    exp = .{ -2, 2, 31, 14 };
    res = false;
    res = try tmsXmtx(&b, 2, &a, 2, &aTb, 2);
    try std.testing.expectEqual(true, res);
    try std.testing.expectEqual(true, equXmtx(&aTb, &exp));
}

test "XMTX: ELA - Larson, Edwards: 2.2 Example 1 test" {
    prntNl();
    var m0: [3]f32 = .{ 0, 0, 0 };
    var m1: [3]f32 = .{ 1, 2, -3 };
    var m2: [3]f32 = .{ -1, -1, 2 };
    var m3: [3]f32 = .{ 0, 1, 4 };
    var m4: [3]f32 = .{ 2, -3, -2 };
    sum2Xvec(&m0, &m1, &m2);
    sum2Xvec(&m0, &m3, &m4);
    var exp: [3]f32 = .{ 2, -1, 1 };

    prntXmtxNl(&m0, 3);
    prntNl();

    prntXmtxNl(&exp, 3);
    prntNl();

    try std.testing.expectEqual(true, equXmtx(&m0, &exp));
}

test "XMTX: ELA - Larson, Edwards: 2.2 Example 2 test" {
    prntNl();
    //Find X = (1/3)(B - A)
    var a: [4]f32 = .{ 1, -2, 0, 3 };
    var b: [4]f32 = .{ -3, 4, 2, 1 };
    var exp: [4]f32 = .{ (-4.0 / 3.0), 2, (2.0 / 3.0), (-2.0 / 3.0) };
    diff1Xvec(&b, &a);
    mulXvec(&b, (1.0 / 3.0));

    prntXmtxNl(&b, 2);
    prntNl();

    prntXmtxNl(&exp, 2);
    prntNl();

    try std.testing.expectEqual(true, equXmtx(&b, &exp));
}

test "XMTX: ELA - Larson, Edwards: 2.2 Example 3 test" {
    prntNl();
    //(AB)C = A(BC)
    var a: [4]f32 = .{ 1, -2, 2, -1 };
    var b: [6]f32 = .{ 1, 0, 2, 3, -2, 1 };
    var c: [6]f32 = .{ -1, 0, 3, 1, 2, 4 };
    var ret1: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    var ret2: [4]f32 = .{ 0, 0, 0, 0 };
    var exp: [4]f32 = .{ 17, 4, 13, 14 };
    var res: bool = false;

    res = try tmsXmtx(&a, 2, &b, 3, &ret1, 3);
    try std.testing.expectEqual(true, res);

    prntXmtxNl(&ret1, 2);
    prntNl();

    res = try tmsXmtx(&ret1, 3, &c, 2, &ret2, 2);
    try std.testing.expectEqual(true, res);

    prntXmtxNl(&ret2, 2);
    prntNl();

    prntXmtxNl(&exp, 2);
    prntNl();

    try std.testing.expectEqual(true, equXmtx(&ret2, &exp));
}

test "XMTX: ELA - Larson, Edwards: 2.2 Example 4 test" {
    prntNl();
    //AB != BA
    var a: [4]f32 = .{ 1, 3, 2, -1 };
    var b: [4]f32 = .{ 2, -1, 0, 2 };
    var ret1: [4]f32 = .{ 0, 0, 0, 0 };
    var ret2: [4]f32 = .{ 0, 0, 0, 0 };
    var exp1: [4]f32 = .{ 2, 5, 4, -4 };
    var exp2: [4]f32 = .{ 0, 7, 4, -2 };
    var res: bool = false;

    res = try tmsXmtx(&a, 2, &b, 2, &ret1, 2);
    try std.testing.expectEqual(true, res);

    prntXmtxNl(&ret1, 2);
    prntNl();

    res = try tmsXmtx(&b, 2, &a, 2, &ret2, 2);
    try std.testing.expectEqual(true, res);

    prntXmtxNl(&ret2, 2);
    prntNl();

    try std.testing.expectEqual(true, equXmtx(&ret1, &exp1));
    try std.testing.expectEqual(true, equXmtx(&ret2, &exp2));
    try std.testing.expectEqual(false, equXmtx(&exp2, &exp1));
    try std.testing.expectEqual(false, equXmtx(&ret2, &ret1));
}

test "XMTX: ELA - Larson, Edwards: 2.2 Example 5 test" {
    prntNl();
    //Show that AC = BC yet A != B
    //A = 1 3 0 1
    //B = 2 4 2 3
    //C = 1 -2 -1 2

    var a: [4]f32 = .{ 1, 3, 0, 1 };
    var b: [4]f32 = .{ 2, 4, 2, 3 };
    var c: [4]f32 = .{ 1, -2, -1, 2 };
    var ret1: [4]f32 = .{ 0, 0, 0, 0 };
    var ret2: [4]f32 = .{ 0, 0, 0, 0 };
    var exp: [4]f32 = .{ -2, 4, -1, 2 };
    var res: bool = false;

    res = try tmsXmtx(&a, 2, &c, 2, &ret1, 2);
    try std.testing.expectEqual(true, res);

    prntXmtxNl(&ret1, 2);
    prntNl();

    res = try tmsXmtx(&b, 2, &c, 2, &ret2, 2);
    try std.testing.expectEqual(true, res);

    prntXmtxNl(&ret2, 2);
    prntNl();

    try std.testing.expectEqual(true, equXmtx(&ret1, &exp));
    try std.testing.expectEqual(true, equXmtx(&ret2, &exp));
    try std.testing.expectEqual(false, equXmtx(&a, &b));
}

test "XMTX: ELA - Larson, Edwards: 2.2 Example 6 test" {
    prntNl();
    //A * I = A
    //A
    var a: [6]f32 = .{ 3, -2, 4, 0, -1, 1 };
    var b: [4]f32 = .{ 1, 0, 0, 1 };
    var ret1: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    var exp1: [6]f32 = .{ 3, -2, 4, 0, -1, 1 };
    var res: bool = false;

    res = try tmsXmtx(&a, 2, &b, 2, &ret1, 2);
    try std.testing.expectEqual(true, res);

    prntXmtxNl(&ret1, 2);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&ret1, &exp1));
    try std.testing.expectEqual(true, equXmtx(&a, &exp1));

    //B
    //C = 1 0 0    -2    -2
    //    0 1 0  *  1 =   1
    //    0 0 1     4     4
    var c: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var d: [3]f32 = .{ -2, 1, 4 };
    var ret2: [3]f32 = .{ 0, 0, 0 };
    var exp2: [3]f32 = .{ -2, 1, 4 };

    res = try tmsXmtx(&c, 3, &d, 1, &ret2, 1);
    try std.testing.expectEqual(true, res);

    prntXmtxNl(&ret2, 3);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&ret2, &exp2));
    try std.testing.expectEqual(true, equXmtx(&d, &exp2));
}

test "XMTX: ELA - Larson, Edwards: 2.2 Example 7 test" {
    prntNl();
    //Find A^3 = 2 -1 = -4 -1
    //           3  0    3 -1
    var a: [4]f32 = .{ 2, -1, 3, 0 };
    var ret1: [4]f32 = .{ 0, 0, 0, 0 };
    var ret2: [4]f32 = .{ 0, 0, 0, 0 };
    var exp: [4]f32 = .{ -4, -1, 3, -6 };
    var res: bool = false;

    res = try tmsXmtx(&a, 2, &a, 2, &ret1, 2);
    try std.testing.expectEqual(true, res);

    prntXmtxNl(&a, 2);
    prntNl();

    prntXmtxNl(&ret1, 2);
    prntNl();

    res = try tmsXmtx(&ret1, 2, &a, 2, &ret2, 2);
    try std.testing.expectEqual(true, res);

    prntXmtxNl(&a, 2);
    prntNl();

    prntXmtxNl(&ret2, 2);
    prntNl();

    try std.testing.expectEqual(true, equXmtx(&ret2, &exp));
}

test "XMTX: ELA - Larson, Edwards: 2.2 Example 8 test" {
    prntNl();
    //Find the transpose of each matrix
    //A = 2 => 2 8
    //    8
    //
    //B = 1  2  3 => 1  4  7
    //    4  5  6    2  5  8
    //    7  8  9    3  6  9
    //
    //C = 1  2  0 => 1  2  0
    //    2  1  0    2  1  0
    //    0  0  1    0  0  1
    //
    //D = 0  1 => 0  2  1
    //    2  4    1  4 -1
    //    1 -1

    var m1: [2]f32 = .{ 2, 8 };
    var r1: [2]f32 = .{ 0, 0 };
    var e1: [2]f32 = .{ 2, 8 };

    var m2: [9]f32 = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    var r2: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var e2: [9]f32 = .{ 1, 4, 7, 2, 5, 8, 3, 6, 9 };

    var m3: [9]f32 = .{ 1, 2, 0, 2, 1, 0, 0, 0, 1 };
    var r3: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var e3: [9]f32 = .{ 1, 2, 0, 2, 1, 0, 0, 0, 1 };

    var m4: [6]f32 = .{ 0, 1, 2, 4, 1, -1 };
    var r4: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    var e4: [6]f32 = .{ 0, 2, 1, 1, 4, -1 };

    trnXmtx(&m1, 1, &r1);
    trnXmtx(&m2, 3, &r2);
    trnXmtxRect(&m3, 3, &r3, 3);
    trnXmtxRect(&m4, 2, &r4, 3);

    prntXmtxNl(&m4, 2);
    prntNl();

    prntXmtxNl(&r4, 3);
    prntNl();

    try std.testing.expectEqual(true, equXmtx(&r1, &e1));
    try std.testing.expectEqual(true, equXmtx(&r2, &e2));
    try std.testing.expectEqual(true, equXmtx(&r3, &e3));
    try std.testing.expectEqual(true, equXmtx(&r4, &e4));
}

test "XMTX: ELA - Larson, Edwards: 2.2 Example 9 test" {
    prntNl();
    //Show that (A * B)^T = B^T * A^T
    //A = 2  1 -1
    //   -1  0  3
    //    0 -2  1
    //
    //B = 3  1
    //    2 -1
    //    3  0

    var a: [9]f32 = .{ 2, 1, -2, -1, 0, 3, 0, -2, 1 };
    var b: [6]f32 = .{ 3, 1, 2, -1, 3, 0 };

    var ab: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    var expAb: [6]f32 = .{ 2, 1, 6, -1, -1, 2 };

    var aT: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var expAt: [9]f32 = .{ 2, -1, 0, 1, 0, -2, -2, 3, 1 };

    var bT: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    var expBt: [6]f32 = .{ 3, 2, 3, 1, -1, 0 };

    var abT: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    var expAbt: [6]f32 = .{ 2, 6, -1, 1, -1, 2 };

    var btAt: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    var expBtAt: [6]f32 = .{ 2, 6, -1, 1, -1, 2 };
    var res: bool = false;
    res = try tmsXmtx(&a, 3, &b, 2, &ab, 2);
    try std.testing.expectEqual(true, res);

    prntXmtxNl(&a, 3);
    prntNl();

    prntXmtxNl(&b, 2);
    prntNl();

    prntXmtxNl(&ab, 2);
    prntNl();

    prntXmtxNl(&expAb, 2);
    prntNl();
    try std.testing.expectEqual(true, res);
    try std.testing.expectEqual(true, equXmtx(&ab, &expAb));

    trnXmtx(&a, 3, &aT);
    prntXmtxNl(&aT, 3);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&aT, &expAt));

    trnXmtxRect(&b, 2, &bT, 3);
    prntXmtxNl(&bT, 3);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&bT, &expBt));

    trnXmtxRect(&ab, 2, &abT, 3);
    try std.testing.expectEqual(true, equXmtx(&abT, &expAbt));

    res = try tmsXmtx(&bT, 3, &aT, 3, &btAt, 3);
    prntXmtxNl(&btAt, 3);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&btAt, &expBtAt));
}

test "XMTX: ELA - Larson, Edwards: 2.2 Example 10 test" {
    prntNl();
    //Find AA^T
    //Show AA^T = (AA^T)^T
    //A = 1  3
    //    0 -2
    //   -2 -1
    var a: [6]f32 = .{ 1, 3, 0, -2, -2, -1 };
    var aT: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };
    var expAt: [6]f32 = .{ 1, 0, -2, 3, -2, -1 };
    var aAt: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var expAat: [9]f32 = .{ 10, -6, -5, -6, 4, 2, -5, 2, 5 };
    var res: bool = false;
    trnXmtxRect(&a, 2, &aT, 3);
    try std.testing.expectEqual(true, equXmtx(&aT, &expAt));

    res = try tmsXmtx(&a, 2, &aT, 3, &aAt, 3);
    try std.testing.expectEqual(true, res);
    try std.testing.expectEqual(true, equXmtx(&aAt, &expAat));
}

test "XMTX: ELA - Larson, Edwards: 2.2 Problem 1, 3, 5 test" {
    prntNl();
    //A = 1  2
    //    3  4
    //
    //B = 0  1
    //   -1  2
    //
    //a = 3
    //b = -4
    //
    //1. aA + bB
    //3. ab(B)
    //5. (a - b)(A - B)

    var A: [4]f32 = .{ 1, 2, 3, 4 };
    var B: [4]f32 = .{ 0, 1, -1, 2 };
    var exp: [4]f32 = .{ 0, 0, 0, 0 };
    const a: f32 = 3;
    const b: f32 = -4;
    const ab: f32 = a * b;

    //1. aA + bB
    var aApBb: [4]f32 = .{ 0, 0, 0, 0 };
    exp = .{ 3, 2, 13, 4 };
    var aA: [4]f32 = .{ 0, 0, 0, 0 };
    cpyXmtx(&A, &aA);
    mulXvec(&aA, a);

    var bB: [4]f32 = .{ 0, 0, 0, 0 };
    cpyXmtx(&B, &bB);
    mulXvec(&bB, b);
    sum2Xvec(&aApBb, &aA, &bB);
    try std.testing.expectEqual(true, equXmtx(&aApBb, &exp));

    //3. ab(B)
    var abB: [4]f32 = .{ 0, 0, 0, 0 };
    B = .{ 0, 1, -1, 2 };
    exp = .{ 0, -12, 12, -24 };
    cpyXmtx(&B, &abB);
    mulXvec(&abB, ab);
    try std.testing.expectEqual(true, equXmtx(&abB, &exp));

    //5. (a - b)(A - B)
    A = .{ 1, 2, 3, 4 };
    B = .{ 0, 1, -1, 2 };
    const amb: f32 = (a - b);
    var aMb: [4]f32 = .{ 0, 0, 0, 0 };
    exp = .{ 7, 7, 28, 14 };
    diff2Xvec(&aMb, &A, &B);
    mulXvec(&aMb, amb);
    try std.testing.expectEqual(true, equXmtx(&aMb, &exp));
}

test "XMTX: ELA - Larson, Edwards: 2.2 Problem 19, 21 test" {
    prntNl();
    //A = 1  2
    //    0 -1
    //
    //I = 1  0
    //    0  1
    //
    //19. A^2
    var A: [4]f32 = .{ 1, 2, 0, -1 };
    var I: [4]f32 = .{ 1, 0, 0, 1 };
    var AA: [4]f32 = .{ 0, 0, 0, 0 };
    var b: bool = false;
    cpyXmtx(&A, &AA);
    b = try tmsXmtx(&A, 2, &A, 2, &AA, 2);
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&AA, &I));

    var iPa: [4]f32 = .{ 0, 0, 0, 0 };
    sum2Xvec(&iPa, &I, &A);

    var aIpa: [4]f32 = .{ 0, 0, 0, 0 };
    var exp: [4]f32 = .{ 2, 2, 0, 0 };
    b = try tmsXmtx(&I, 2, &iPa, 2, &aIpa, 2);
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&exp, &aIpa));
}

test "XMTX: ELA - Larson, Edwards: 2.3 Example 1 test" {
    prntNl();
    //A = -1  2
    //    -1  1
    //
    //B =  1 -2
    //     1 -1
    //Show that AB = BA = I
    var A: [4]f32 = .{ -1, 2, -1, 1 };
    var B: [4]f32 = .{ 1, -2, 1, -1 };
    var AB: [4]f32 = .{ 0, 0, 0, 0 };
    var BA: [4]f32 = .{ 0, 0, 0, 0 };
    var I: [4]f32 = .{ 1, 0, 0, 1 };
    var b: bool = false;

    b = try tmsXmtx(&A, 2, &B, 2, &AB, 2);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&B, 2, &A, 2, &BA, 2);
    try std.testing.expectEqual(true, b);

    try std.testing.expectEqual(true, equXmtx(&AB, &I));
    try std.testing.expectEqual(true, equXmtx(&BA, &I));
}

test "XMTX: ELA - Larson, Edwards: 2.3 Example 2 test" {
    prntNl();
    //Find the inverse of A and verify it is correct.
    //A = 1  4
    //   -1  3
    //
    //A^-1 = -3 -4
    //        1  1
    //
    var A: [4]f32 = .{ 1, 4, -1, -3 };
    var I: [4]f32 = .{ 1, 0, 0, 1 };
    var B: [4]f32 = .{ 0, 0, 0, 0 };
    var exp: [4]f32 = .{ -3, -4, 1, 1 };
    var res: bool = false;

    var sclr: f32 = 0.0;
    res = try rdcXmtx(&A, 2, false, &B, true, &I, 2, false, &sclr);
    try std.testing.expectEqual(true, res);

    prntXmtxNl(&B, 2);
    prntNl();

    prntXmtxNl(&I, 2);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&exp, &I));
}

test "XMTX: ELA - Larson, Edwards: 2.3 Example 3 test" {
    prntNl();
    //A = 1 -1  0
    //    1  0 -1
    //   -6  2  3
    //
    //Find A^-1
    var A: [9]f32 = .{ 1, -1, 0, 1, 0, -1, -6, 2, 3 };
    var I: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var R: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var exp: [9]f32 = .{ -2, -3, -1, -3, -3, -1, -2, -4, -1 };
    var b: bool = false;

    var sclr: f32 = 0.0;
    b = try rdcXmtx(&A, 3, false, &R, true, &I, 3, false, &sclr);

    prntXmtxNl(&A, 3);
    prntNl();

    prntXmtxNl(&I, 3);
    prntNl();
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&exp, &I));
}

test "XMTX: ELA - Larson, Edwards: 2.3 Example 4 test" {
    prntNl();
    //A = 1  2  0
    //    3 -1  2
    //   -2  3  -2
    //
    //Find A^-1
    var A: [9]f32 = .{ 1, 2, 0, 3, -1, 2, -2, 3, -2 };
    var I: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var R: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var b: bool = false;

    var sclr: f32 = 0.0;
    b = try rdcXmtx(&A, 3, false, &R, true, &I, 3, false, &sclr);

    prntXmtxNl(&A, 3);
    prntNl();

    prntXmtxNl(&I, 3);
    prntNl();
    try std.testing.expectEqual(false, b);
}

test "XMTX: ELA - Larson, Edwards: 2.3 Example 5 test" {
    prntNl();
    //A = 3 -1 -2 2
    //B = 3 -1 -6 2
    //Find A^-1 and B^-1
    var A: [4]f32 = .{ 3, -1, -2, 2 };
    var aI: [4]f32 = .{ 0, 0, 0, 0 };
    var B: [4]f32 = .{ 3, -1, -6, 2 };
    var bI: [4]f32 = .{ 0, 0, 0, 0 };
    var I1: [4]f32 = .{ 1, 0, 0, 1 };
    var I2: [4]f32 = .{ 1, 0, 0, 1 };
    var exp: [4]f32 = .{ (1.0 / 2.0), (1.0 / 4.0), (1.0 / 2.0), (3.0 / 4.0) };
    var b: bool = false;

    var sclr: f32 = 0.0;
    b = try rdcXmtx(&A, 2, false, &aI, true, &I1, 2, false, &sclr);
    try std.testing.expectEqual(true, b);

    sclr = 0.0;
    b = try rdcXmtx(&B, 2, false, &bI, true, &I2, 2, false, &sclr);
    try std.testing.expectEqual(false, b);

    prntXmtxNl(&A, 2);
    prntNl();

    prntXmtxNl(&I1, 2);
    prntNl();

    clnXmtx(&I1);
    try std.testing.expectEqual(true, equXmtx(&exp, &I1));
}

test "XMTX: ELA - Larson, Edwards: 2.3 Example 6 test" {
    prntNl();
    //A = 1 1 2 4
    //Find A^-2 find (A^2)^-1 and find (A^-1)(A^-1)
    var A: [4]f32 = .{ 1, 1, 2, 4 };
    var aA: [4]f32 = .{ 0, 0, 0, 0 };
    var I: [4]f32 = .{ 1, 0, 0, 1 };
    var invA: [4]f32 = .{ 2.0, -5.0e-01, -1.0, 5.0e-01 };
    var invAa: [4]f32 = .{ 0, 0, 0, 0 };
    var expAa: [4]f32 = .{ 3, 5, 10, 18 };
    var b: bool = false;

    var expInvA: [4]f32 = .{ 2, (-1.0 / 2.0), -1, (1.0 / 2.0) };

    var expInvAa: [4]f32 = .{ (9.0 / 2.0), (-5.0 / 4.0), (-5.0 / 2.0), (3.0 / 4.0) };

    b = try tmsXmtx(&A, 2, &A, 2, &aA, 2);
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&expAa, &aA));

    var sclr: f32 = 0.0;
    b = try rdcXmtx(&A, 2, false, &invA, true, &I, 2, false, &sclr);
    try std.testing.expectEqual(true, b);

    clnXmtx(&I);
    clnXmtx(&invA);
    prntXmtxNl(&I, 2);
    prntNl();
    prntXmtxNl(&invA, 2);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&expInvA, &I));

    b = try tmsXmtx(&I, 2, &I, 2, &invAa, 2);
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&expInvAa, &invAa));

    A = .{ 1, 1, 2, 4 };
    b = try tmsXmtx(&A, 2, &A, 2, &aA, 2);
    I = .{ 1, 0, 0, 1 };
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&expAa, &aA));

    sclr = 0.0;
    b = try rdcXmtx(&aA, 2, false, &invAa, true, &I, 2, false, &sclr);
    try std.testing.expectEqual(true, b);
    cpyXmtx(&I, &invAa);
    try std.testing.expectEqual(true, equXmtx(&expInvAa, &invAa));
}

test "XMTX: ELA - Larson, Edwards: 2.3 Example 7 test" {
    prntNl();
    //A = 1  3  3
    //    1  4  3
    //    1  3  4
    //
    //B = 1  2  3
    //    1  3  3
    //    2  4  3
    var A: [9]f32 = .{ 1, 3, 3, 1, 4, 3, 1, 3, 4 };
    var B: [9]f32 = .{ 1, 2, 3, 1, 3, 3, 2, 4, 3 };
    var AB: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var invA: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var expInvA: [9]f32 = .{ 7, -3, -3, -1, 1, 0, -1, 0, 1 };
    var invB: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var expInvB: [9]f32 = .{ 1, -2, 1, -1, 1, 0, (2.0 / 3.0), 0, (-1.0 / 3.0) };
    var invAB: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var expInvAB: [9]f32 = .{ 8, -5, -2, -8, 4, 3, 5, -2, (-7.0 / 3.0) };
    var tmp: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var b: bool = false;
    const cols: usize = 3;

    b = try tmsXmtx(&A, cols, &B, cols, &AB, cols); //Set AB
    try std.testing.expectEqual(true, b);

    var sclr: f32 = 0.0;
    b = try rdcXmtx(&A, cols, false, &tmp, true, &invA, cols, false, &sclr); //Set invA
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&expInvA, &invA));

    sclr = 0.0;
    b = try rdcXmtx(&B, cols, false, &tmp, true, &invB, cols, false, &sclr); //Set invB
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&expInvB, &invB));

    sclr = 0.0;
    b = try rdcXmtx(&AB, cols, false, &tmp, true, &invAB, cols, false, &sclr); //Set invAB
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&expInvAB, &invAB));
}

test "XMTX: ELA - Larson, Edwards: 2.3 Example 8 test" {
    prntNl();
    //A = 2 3 1
    //    3 3 1
    //    2 4 1
    //solve for a: -1  1 -2
    //          b:  4  8  5
    //          c:  0  0  0
    var A: [9]f32 = .{ 2, 3, 1, 3, 3, 1, 2, 4, 1 };
    var invA: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var idtA: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var tmp: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    var sA: [3]f32 = .{ -1, 1, -2 };
    var sB: [3]f32 = .{ 4, 8, 5 };
    var sC: [3]f32 = .{ 0, 0, 0 };
    var sD: [3]f32 = .{ 0, 0, 0 };
    var b: bool = false;
    const cols: usize = 3;

    var sclr: f32 = 0.0;
    b = try rdcXmtx(&A, cols, false, &tmp, true, &invA, cols, false, &sclr); //Set invA
    try std.testing.expectEqual(true, b);

    prntNl();
    std.debug.print("XMTX: A\n", .{});
    prntXmtxNl(&A, cols);

    prntNl();
    std.debug.print("XMTX: TMP\n", .{});
    prntXmtxNl(&tmp, cols);

    prntNl();
    std.debug.print("XMTX: invA\n", .{});
    prntXmtxNl(&invA, cols);

    try std.testing.expectEqual(true, equXmtx(&idtA, &tmp));

    //ANSWERS
    var exp: [3]f32 = .{ 2, -1, -2 };
    b = try tmsXmtx(&invA, cols, &sA, 1, &sD, 1);
    try std.testing.expectEqual(true, b);
    prntNl();
    std.debug.print("XMTX: a:\n", .{});
    prntXmtxNl(&sD, 1);
    try std.testing.expectEqual(true, equXmtx(&exp, &sD));

    exp = .{ 4, 1, -7 };
    b = try tmsXmtx(&invA, cols, &sB, 1, &sD, 1);
    try std.testing.expectEqual(true, b);
    prntNl();
    std.debug.print("XMTX: b:\n", .{});
    prntXmtxNl(&sD, 1);
    try std.testing.expectEqual(true, equXmtx(&exp, &sD));

    exp = .{ 0, 0, 0 };
    b = try tmsXmtx(&invA, cols, &sC, 1, &sD, 1);
    try std.testing.expectEqual(true, b);
    prntNl();
    std.debug.print("XMTX: c:\n", .{});
    prntXmtxNl(&sD, 1);
    try std.testing.expectEqual(true, equXmtx(&exp, &sD));
}

test "XMTX: ELA - Larson, Edwards: 2.3 Problem 1, 3 test" {
    prntNl();
    //Show that B is the inverse of A
    //1:
    //A = 1 2   B = -2   1
    //    3 4       3/2 -1/2
    //3:
    //A = -2   2   3  B = 1/3 * -4  -5   3
    //     1  -1   0            -4  -8  -3
    //     0   1   4             1   2   0
    var b: bool = false;
    var A2: [4]f32 = .{ 1, 2, 3, 4 };
    var B2: [4]f32 = .{ -2, 1, (3.0 / 2.0), (-1.0 / 2.0) };
    var tmp2: [4]f32 = .{ 0, 0, 0, 0 };
    var idt2: [4]f32 = .{ 1, 0, 0, 1 };
    const cols2: usize = 2;

    var sclr: f32 = 0.0;
    b = try rdcXmtx(&A2, cols2, false, &tmp2, true, &idt2, cols2, false, &sclr);
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&B2, &idt2));

    var A3: [9]f32 = .{ -2, 2, 3, 1, -1, 0, 0, 1, 4 };
    var B3: [9]f32 = .{ -4, -5, 3, -4, -8, 3, 1, 2, 0 };
    var tmp3: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var idt3: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    const m: f32 = (1.0 / 3.0);
    const cols3: usize = 3;

    mulXvec(&B3, m);
    sclr = 0.0;
    b = try rdcXmtx(&A3, cols3, false, &tmp3, true, &idt3, cols3, false, &sclr);
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&B3, &idt3));
}

test "XMTX: ELA - Larson, Edwards: 2.3 Problem 5, 7, 9 test" {
    prntNl();
    //Find the inverse if it exists
    //5:
    //A = 1  2     invA = 7 -2
    //    3  7           -3  1
    //7:
    //B = -7  33   invB = -19 -33
    //     4 -19          -4  -7
    //9:
    //C = 1  1  1  invC =  1  1 -1
    //    3  5  4         -3  2 -1
    //    3  6  5          3 -3  2
    var A: [4]f32 = .{ 1, 2, 3, 7 };
    var invA: [4]f32 = .{ 7, -2, -3, 1 };
    var B: [4]f32 = .{ -7, 33, 4, -19 };
    var invB: [4]f32 = .{ -19, -33, -4, -7 };
    var tmp2: [4]f32 = .{ 0, 0, 0, 0 };
    var idt2: [4]f32 = .{ 1, 0, 0, 1 };
    var C: [9]f32 = .{ 1, 1, 1, 3, 5, 4, 3, 6, 5 };
    var invC: [9]f32 = .{ 1, 1, -1, -3, 2, -1, 3, -3, 2 };
    var tmp3: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var idt3: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var b: bool = false;

    var sclr: f32 = 0.0;
    b = try rdcXmtx(&A, 2, false, &tmp2, true, &idt2, 2, false, &sclr);
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&invA, &idt2));

    prntNl();
    std.debug.print("XMTX: A:\n", .{});
    prntXmtxNl(&idt2, 2);

    idt2 = .{ 1, 0, 0, 1 };
    sclr = 0.0;
    b = try rdcXmtx(&B, 2, false, &tmp2, true, &idt2, 2, false, &sclr);
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&invB, &idt2));

    prntNl();
    std.debug.print("XMTX: B:\n", .{});
    prntXmtxNl(&idt2, 2);

    sclr = 0.0;
    b = try rdcXmtx(&C, 3, false, &tmp3, true, &idt3, 3, false, &sclr);
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, equXmtx(&invC, &idt3));

    prntNl();
    std.debug.print("XMTX: C:\n", .{});
    prntXmtxNl(&idt3, 3);
}

test "XMTX: ELA - Larson, Edwards: 2.4 Example 3 test" {
    prntNl();
    //A = 0  1  3  5
    //    1 -3  0  2
    //    2 -6  2  0
    var A: [12]f32 = .{ 0, 1, 3, 5, 1, -3, 0, 2, 2, -6, 2, 0 };
    var tmp3: [12]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var idt3: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var b: bool = false;

    var sclr: f32 = 0.0;
    b = try rdcXmtx(&A, 4, true, &tmp3, true, &idt3, 3, false, &sclr);
    try std.testing.expectEqual(true, b);

    prntNl();
    std.debug.print("XMTX: A:\n", .{});
    prntXmtxNl(&A, 3);

    prntNl();
    std.debug.print("XMTX: TMP3:\n", .{});
    prntXmtxNl(&tmp3, 3);

    prntNl();
    std.debug.print("XMTX: IDT3:\n", .{});
    prntXmtxNl(&idt3, 3);
}

test "XMTX: ELA - Larson, Edwards: 2.5 Example 5, 6 test" {
    prntNl();
    //Matrix Cryptography
    //0 = _     6 = F       12 = L      18 = R      24 = X
    //1 = A     7 = G       13 = M      19 = S      25 = Y
    //2 = B     8 = H       14 = N      20 = T      26 = Z
    //3 = C     9 = I       15 = O      21 = U
    //4 = D     10 = J      16 = P      22 = V
    //5 = E     11 = K      17 = Q      23 = W

    //Clear text matrices
    //[13  5  5] [20  0 13] [5  0 13] [15 14  4] [1 25  0]
    //  M  E  E    T  _  M   E  _  M    O  N  D   A  Y  _

    //A = 1  -2   2
    //   -1   1   3
    //    1  -1  -4
    var A: [9]f32 = .{ 1, -2, 2, -1, 1, 3, 1, -1, -4 };
    var tmp3: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var idt3: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var b: bool = false;
    var w1: [3]f32 = .{ 13, 5, 5 };
    var w2: [3]f32 = .{ 20, 0, 13 };
    var w3: [3]f32 = .{ 5, 0, 13 };
    var w4: [3]f32 = .{ 15, 14, 4 };
    var w5: [3]f32 = .{ 1, 25, 0 };

    var eW1: [3]f32 = .{ 0, 0, 0 };
    var eW2: [3]f32 = .{ 0, 0, 0 };
    var eW3: [3]f32 = .{ 0, 0, 0 };
    var eW4: [3]f32 = .{ 0, 0, 0 };
    var eW5: [3]f32 = .{ 0, 0, 0 };

    var cW1: [3]f32 = .{ 0, 0, 0 };
    var cW2: [3]f32 = .{ 0, 0, 0 };
    var cW3: [3]f32 = .{ 0, 0, 0 };
    var cW4: [3]f32 = .{ 0, 0, 0 };
    var cW5: [3]f32 = .{ 0, 0, 0 };

    //Get inverse matrix
    var sclr: f32 = 0.0;
    b = try rdcXmtx(&A, 3, false, &tmp3, true, &idt3, 3, false, &sclr);
    try std.testing.expectEqual(true, b);

    prntNl();
    std.debug.print("XMTX: A:\n", .{});
    prntXmtxNl(&A, 3);

    prntNl();
    std.debug.print("XMTX: TMP3:\n", .{});
    prntXmtxNl(&tmp3, 3);

    prntNl();
    std.debug.print("XMTX: IDT3:\n", .{});
    prntXmtxNl(&idt3, 3);

    A = .{ 1, -2, 2, -1, 1, 3, 1, -1, -4 };

    //Encode word 1
    b = try tmsXmtx(&w1, 3, &A, 3, &eW1, 3);
    try std.testing.expectEqual(true, b);

    prntNl();
    std.debug.print("XMTX: W1:\n", .{});
    prntXmtxNl(&w1, 3);

    prntNl();
    std.debug.print("XMTX: ENC_W1:\n", .{});
    prntXmtxNl(&eW1, 3);

    //Encode word 2
    b = try tmsXmtx(&w2, 3, &A, 3, &eW2, 3);
    try std.testing.expectEqual(true, b);

    prntNl();
    std.debug.print("XMTX: W2:\n", .{});
    prntXmtxNl(&w2, 3);

    prntNl();
    std.debug.print("XMTX: ENC_W2:\n", .{});
    prntXmtxNl(&eW2, 3);

    //Encode word 3
    b = try tmsXmtx(&w3, 3, &A, 3, &eW3, 3);
    try std.testing.expectEqual(true, b);

    prntNl();
    std.debug.print("XMTX: W3:\n", .{});
    prntXmtxNl(&w3, 3);

    prntNl();
    std.debug.print("XMTX: ENC_W3:\n", .{});
    prntXmtxNl(&eW3, 3);

    //Encode word 4
    b = try tmsXmtx(&w4, 3, &A, 3, &eW4, 3);
    try std.testing.expectEqual(true, b);

    prntNl();
    std.debug.print("XMTX: W4:\n", .{});
    prntXmtxNl(&w4, 3);

    prntNl();
    std.debug.print("XMTX: ENC_W4:\n", .{});
    prntXmtxNl(&eW4, 3);

    //Encode word 5
    b = try tmsXmtx(&w5, 3, &A, 3, &eW5, 3);
    try std.testing.expectEqual(true, b);

    prntNl();
    std.debug.print("XMTX: W5:\n", .{});
    prntXmtxNl(&w5, 3);

    prntNl();
    std.debug.print("XMTX: ENC_W5:\n", .{});
    prntXmtxNl(&eW5, 3);

    //Decode word 1
    b = try tmsXmtx(&eW1, 3, &idt3, 3, &cW1, 3);
    try std.testing.expectEqual(true, b);

    prntNl();
    std.debug.print("XMTX: DEC_W1:\n", .{});
    prntXmtxNl(&cW1, 3);

    //Decode word 2
    b = try tmsXmtx(&eW2, 3, &idt3, 3, &cW2, 3);
    try std.testing.expectEqual(true, b);

    prntNl();
    std.debug.print("XMTX: DEC_W2:\n", .{});
    prntXmtxNl(&cW2, 3);

    //Decode word 3
    b = try tmsXmtx(&eW3, 3, &idt3, 3, &cW3, 3);
    try std.testing.expectEqual(true, b);

    prntNl();
    std.debug.print("XMTX: DEC_W3:\n", .{});
    prntXmtxNl(&cW3, 3);

    //Decode word 4
    b = try tmsXmtx(&eW4, 3, &idt3, 3, &cW4, 3);
    try std.testing.expectEqual(true, b);

    prntNl();
    std.debug.print("XMTX: DEC_W4:\n", .{});
    prntXmtxNl(&cW4, 3);

    //Decode word 5
    b = try tmsXmtx(&eW5, 3, &idt3, 3, &cW5, 3);
    try std.testing.expectEqual(true, b);

    prntNl();
    std.debug.print("XMTX: DEC_W5:\n", .{});
    prntXmtxNl(&cW5, 3);

    try std.testing.expectEqual(true, equXmtx(&w1, &cW1));
    try std.testing.expectEqual(true, equXmtx(&w2, &cW2));
    try std.testing.expectEqual(true, equXmtx(&w3, &cW3));
    try std.testing.expectEqual(true, equXmtx(&w4, &cW4));
    try std.testing.expectEqual(true, equXmtx(&w5, &cW5));
}

test "XMTX: ELA - Larson, Edwards: 2.5 Example 15, 17 test" {
    prntNl();
    //0 = _     6 = F       12 = L      18 = R      24 = X
    //1 = A     7 = G       13 = M      19 = S      25 = Y
    //2 = B     8 = H       14 = N      20 = T      26 = Z
    //3 = C     9 = I       15 = O      21 = U
    //4 = D     10 = J      16 = P      22 = V
    //5 = E     11 = K      17 = Q      23 = W

    //15
    //SEL       L_C      ONS        OLI       DAT      ED_
    //[19 5 12] [12 0 3] [15 14 19] [15 12 9] [4 1 20] [5 4 0]
    //A = 1  -1   0
    //    1   0  -1
    //   -6   2   3
    var A: [9]f32 = .{ 1, -2, 2, -1, 1, 3, 1, -1, -4 };
    var tmpA: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var idtA: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    var aW1: [3]f32 = .{ 19, 5, 12 };
    var aW2: [3]f32 = .{ 12, 0, 3 };
    var aW3: [3]f32 = .{ 15, 14, 19 };
    var aW4: [3]f32 = .{ 15, 12, 9 };
    var aW5: [3]f32 = .{ 4, 1, 20 };
    var aW6: [3]f32 = .{ 5, 4, 0 };

    var aEw1: [3]f32 = .{ 0, 0, 0 };
    var aEw2: [3]f32 = .{ 0, 0, 0 };
    var aEw3: [3]f32 = .{ 0, 0, 0 };
    var aEw4: [3]f32 = .{ 0, 0, 0 };
    var aEw5: [3]f32 = .{ 0, 0, 0 };
    var aEw6: [3]f32 = .{ 0, 0, 0 };

    var aCw1: [3]f32 = .{ 0, 0, 0 };
    var aCw2: [3]f32 = .{ 0, 0, 0 };
    var aCw3: [3]f32 = .{ 0, 0, 0 };
    var aCw4: [3]f32 = .{ 0, 0, 0 };
    var aCw5: [3]f32 = .{ 0, 0, 0 };
    var aCw6: [3]f32 = .{ 0, 0, 0 };

    var b: bool = false;
    var sclr: f32 = 0.0;
    b = try rdcXmtx(&A, 3, false, &tmpA, true, &idtA, 3, false, &sclr);
    try std.testing.expectEqual(true, b);

    //Encode
    b = try tmsXmtx(&aW1, 3, &A, 3, &aEw1, 3);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&aW2, 3, &A, 3, &aEw2, 3);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&aW3, 3, &A, 3, &aEw3, 3);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&aW4, 3, &A, 3, &aEw4, 3);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&aW5, 3, &A, 3, &aEw5, 3);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&aW6, 3, &A, 3, &aEw6, 3);
    try std.testing.expectEqual(true, b);

    //Decode
    b = try tmsXmtx(&aEw1, 3, &idtA, 3, &aCw1, 3);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&aEw2, 3, &idtA, 3, &aCw2, 3);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&aEw3, 3, &idtA, 3, &aCw3, 3);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&aEw4, 3, &idtA, 3, &aCw4, 3);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&aEw5, 3, &idtA, 3, &aCw5, 3);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&aEw6, 3, &idtA, 3, &aCw6, 3);
    try std.testing.expectEqual(true, b);

    try std.testing.expectEqual(true, equXmtx(&aW1, &aCw1));
    try std.testing.expectEqual(true, equXmtx(&aW2, &aCw2));
    try std.testing.expectEqual(true, equXmtx(&aW3, &aCw3));
    try std.testing.expectEqual(true, equXmtx(&aW4, &aCw4));
    try std.testing.expectEqual(true, equXmtx(&aW5, &aCw5));
    try std.testing.expectEqual(true, equXmtx(&aW6, &aCw6));

    //17
    //CO     ME     _H    OM      E_    SO      ON
    //[3 15] [13 5] [0 8] [15 13] [5 0] [19 15] [15 14]
    //B = 1   2
    //    3   5
    var B: [4]f32 = .{ 1, 2, 3, 5 };
    var tmpB: [4]f32 = .{ 0, 0, 0, 0 };
    var idtB: [4]f32 = .{ 1, 0, 0, 1 };

    var bW1: [2]f32 = .{ 3, 15 };
    var bW2: [2]f32 = .{ 13, 5 };
    var bW3: [2]f32 = .{ 0, 8 };
    var bW4: [2]f32 = .{ 15, 13 };
    var bW5: [2]f32 = .{ 5, 0 };
    var bW6: [2]f32 = .{ 19, 15 };
    var bW7: [2]f32 = .{ 15, 14 };

    var bEw1: [2]f32 = .{ 0, 0 };
    var bEw2: [2]f32 = .{ 0, 0 };
    var bEw3: [2]f32 = .{ 0, 0 };
    var bEw4: [2]f32 = .{ 0, 0 };
    var bEw5: [2]f32 = .{ 0, 0 };
    var bEw6: [2]f32 = .{ 0, 0 };
    var bEw7: [2]f32 = .{ 0, 0 };

    var bCw1: [2]f32 = .{ 0, 0 };
    var bCw2: [2]f32 = .{ 0, 0 };
    var bCw3: [2]f32 = .{ 0, 0 };
    var bCw4: [2]f32 = .{ 0, 0 };
    var bCw5: [2]f32 = .{ 0, 0 };
    var bCw6: [2]f32 = .{ 0, 0 };
    var bCw7: [2]f32 = .{ 0, 0 };

    //Encode
    sclr = 0.0;
    b = try rdcXmtx(&B, 2, false, &tmpB, true, &idtB, 2, false, &sclr);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&bW1, 2, &B, 2, &bEw1, 2);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&bW2, 2, &B, 2, &bEw2, 2);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&bW3, 2, &B, 2, &bEw3, 2);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&bW4, 2, &B, 2, &bEw4, 2);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&bW5, 2, &B, 2, &bEw5, 2);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&bW6, 2, &B, 2, &bEw6, 2);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&bW7, 2, &B, 2, &bEw7, 2);
    try std.testing.expectEqual(true, b);

    //Decode
    b = try tmsXmtx(&bEw1, 2, &idtB, 2, &bCw1, 2);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&bEw2, 2, &idtB, 2, &bCw2, 2);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&bEw3, 2, &idtB, 2, &bCw3, 2);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&bEw4, 2, &idtB, 2, &bCw4, 2);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&bEw5, 2, &idtB, 2, &bCw5, 2);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&bEw6, 2, &idtB, 2, &bCw6, 2);
    try std.testing.expectEqual(true, b);

    b = try tmsXmtx(&bEw7, 2, &idtB, 2, &bCw7, 2);
    try std.testing.expectEqual(true, b);

    try std.testing.expectEqual(true, equXmtx(&bW1, &bCw1));
    try std.testing.expectEqual(true, equXmtx(&bW2, &bCw2));
    try std.testing.expectEqual(true, equXmtx(&bW3, &bCw3));
    try std.testing.expectEqual(true, equXmtx(&bW4, &bCw4));
    try std.testing.expectEqual(true, equXmtx(&bW5, &bCw5));
    try std.testing.expectEqual(true, equXmtx(&bW6, &bCw6));
    try std.testing.expectEqual(true, equXmtx(&bW7, &bCw7));
}

test "XMTX: ELA - Larson, Edwards: 3.1 Example 1, 2, 3, 4, 5 test" {
    prntNl();
    //Example 1
    var A: [4]f32 = .{ 2, -3, 1, 2 };
    var detA: f32 = detXmtx2(&A);
    try std.testing.expectEqual(true, isEquF32(detA, 7.0, true));

    A = .{ 2, 1, 4, 2 };
    detA = detXmtx2(&A);
    try std.testing.expectEqual(true, isEquF32(detA, 0.0, true));

    A = .{ 0, 3, 2, 4 };
    detA = detXmtx2(&A);
    try std.testing.expectEqual(true, isEquF32(detA, -6.0, true));

    //Example 2
    var B: [9]f32 = .{ 0, 2, 1, 3, -1, 2, 4, 0, 1 };
    var expB: [9]f32 = .{ -1, 5, 4, -2, -4, 8, 5, 3, -6 };
    var cofSign: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var cofB: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var retB: [4]f32 = .{ 0, 0, 0, 0 };
    var b: bool = false;
    var row: usize = 0;
    var col: usize = 0;
    const rows: usize = 3;
    const cols: usize = 3;

    while (row < rows) : (row += 1) {
        col = 0;
        while (col < cols) : (col += 1) {
            prntNl();
            std.debug.print("Row: {} Col: {}\n", .{ row, col });

            cofSign[(row * cols) + col] = cofXmtxSign(row, col, true);

            b = cofXmtx(&B, cols, row, col, &retB, 2);
            try std.testing.expectEqual(true, b);

            std.debug.print("XMTX: B:\n", .{});
            prntXmtxNl(&B, cols);

            std.debug.print("XMTX: retB:\n", .{});
            prntXmtxNl(&retB, 2);

            cofB[(row * cols) + col] = detXmtx2(&retB) * cofSign[(row * cols) + col];
            std.debug.print("Sign: {} Val: {}\n", .{ cofSign[(row * cols) + col], cofB[(row * cols) + col] });
        }
    }

    try std.testing.expectEqual(true, equXmtx(&cofB, &expB));

    //Example 3
    const alloc: std.mem.Allocator = std.testing.allocator;
    var C: [9]f32 = .{ 0, 2, 1, 3, -1, 2, 4, 0, 1 };
    const detC: f32 = try detXmtx(&C, cols, &alloc, 0);
    try std.testing.expectEqual(true, isEquF32(detC, 14.0, true));

    //Example 4
    var D: [16]f32 = .{ 1, -2, 3, 0, -1, 1, 0, 2, 0, 2, 0, 3, 3, 4, 0, -2 };
    const detD: f32 = try detXmtx(&D, 4, &alloc, 0);
    try std.testing.expectEqual(true, isEquF32(detD, 39.0, true));

    //Example 5
    var E: [9]f32 = .{ 0, 2, 1, 3, -1, 2, 4, -4, 1 };
    const detE: f32 = try detXmtx(&E, cols, &alloc, 0);
    try std.testing.expectEqual(true, isEquF32(detE, 2.0, true));
}

test "XMTX: ELA - Larson, Edwards: 3.1 Example 6 test" {
    prntNl();
    var A: [16]f32 = .{ 2, 0, 0, 0, 4, -2, 0, 0, -5, 6, 1, 0, 1, 5, 3, 3 };
    var B: [25]f32 = .{ -1, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, -2 };
    var detA: f32 = 0;
    var detB: f32 = 0;
    detA = detDiagXmtx(&A, 4);
    detB = detDiagXmtx(&B, 5);
    std.debug.print("detA: {} detB: {}\n", .{ detA, detB });
    try std.testing.expectEqual(true, isEquF32(detA, -12.0, true));
    try std.testing.expectEqual(true, isEquF32(detB, 48.0, true));
}

test "XMTX: ELA - Larson, Edwards: 3.1 Problem 1, 3, 5 test" {
    prntNl();
    //Find the determinant of the given matrix
    //1
    //A = 1
    var A: [1]f32 = .{1};
    const detA: f32 = detXmtx1(&A);
    const expA: f32 = 1;
    try std.testing.expectEqual(true, isEquF32(detA, expA, true));

    //2
    //B = 2,1 ,3,4
    var B: [4]f32 = .{ 2, 1, 3, 4 };
    const detB: f32 = detXmtx2(&B);
    const expB: f32 = 5;
    try std.testing.expectEqual(true, isEquF32(detB, expB, true));

    //3
    //C = 5,2 ,-6,3
    var C: [4]f32 = .{ 5, 2, -6, 3 };
    const detC: f32 = detXmtx2(&C);
    const expC: f32 = 27;
    try std.testing.expectEqual(true, isEquF32(detC, expC, true));
}

test "XMTX: ELA - Larson, Edwards: 3.1 Problem 13, 15 test" {
    prntNl();
    //Find the minors and cofactors of the given matrix.
    //13
    //A = 1,2 ,3,4
    //Example 2
    var B: [4]f32 = .{ 1, 2, 3, 4 };
    var expB: [4]f32 = .{ 4, -3, -2, 1 };
    var cofSignB: [4]f32 = .{ 0, 0, 0, 0 };
    var cofB: [4]f32 = .{ 0, 0, 0, 0 };
    var retB: [1]f32 = .{0};
    var b: bool = false;
    var row: usize = 0;
    var col: usize = 0;
    var rows: usize = 2;
    var cols: usize = 2;

    while (row < rows) : (row += 1) {
        col = 0;
        while (col < cols) : (col += 1) {
            prntNl();
            std.debug.print("B: Row: {} Col: {}\n", .{ row, col });

            cofSignB[(row * cols) + col] = cofXmtxSign(row, col, true);

            b = cofXmtx(&B, cols, row, col, &retB, 1);
            try std.testing.expectEqual(true, b);

            //std.debug.print("XMTX: B:\n", .{});
            //prntXmtx(&B, 3);

            std.debug.print("XMTX: retB:\n", .{});
            prntXmtxNl(&retB, 2);

            cofB[(row * cols) + col] = detXmtx1(&retB) * cofSignB[(row * cols) + col];
            std.debug.print("SignB: {} Val: {}\n", .{ cofSignB[(row * cols) + col], cofB[(row * cols) + col] });
        }
    }

    try std.testing.expectEqual(true, equXmtx(&cofB, &expB));

    //15
    //B = -3,2,1 ,4,5,6 ,2,-3,1
    //Example 2
    var C: [9]f32 = .{ -3, 2, 1, 4, 5, 6, 2, -3, 1 };
    var expC: [9]f32 = .{ 23, 8, -22, -5, -5, -5, 7, 22, -23 };
    var cofSignC: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var cofC: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var retC: [4]f32 = .{ 0, 0, 0, 0 };
    row = 0;
    col = 0;
    rows = 3;
    cols = 3;

    while (row < rows) : (row += 1) {
        col = 0;
        while (col < cols) : (col += 1) {
            prntNl();
            std.debug.print("C: Row: {} Col: {}\n", .{ row, col });

            cofSignC[(row * cols) + col] = cofXmtxSign(row, col, true);

            b = cofXmtx(&C, cols, row, col, &retC, 2);
            try std.testing.expectEqual(true, b);

            //std.debug.print("XMTX: C:\n", .{});
            //prntXmtx(&C, 3);

            std.debug.print("XMTX: retC:\n", .{});
            prntXmtxNl(&retC, 2);

            cofC[(row * cols) + col] = detXmtx2(&retC) * cofSignC[(row * cols) + col];
            std.debug.print("SignC: {} Val: {}\n", .{ cofSignC[(row * cols) + col], cofC[(row * cols) + col] });
        }
    }

    try std.testing.expectEqual(true, equXmtx(&cofC, &expC));
}

test "XMTX: ELA - Larson, Edwards: 3.2 Example 1, 2, 3 test" {
    prntNl();
    //Show that the elementary row operation is true.
    //Example 1
    //A: Row interchange
    //detA = |2, -3, 1, 4| = 11
    //detB = |1, 4, 2, -3| = -11
    var A: [4]f32 = .{ 2, -3, 1, 4 };
    var detA: f32 = 0;
    detA = detXmtx2(&A);
    var B: [4]f32 = .{ 1, 4, 2, -3 };
    var detB: f32 = 0;
    altXmtxRows(0, 1, 2, &A, &B);
    detB = detXmtx2(&B);
    var exp: f32 = 11;
    try std.testing.expectEqual(true, isEquF32(exp, detA, true));
    try std.testing.expectEqual(true, isEquF32(-1.0 * exp, detB, true));

    //B: Add scalar (-2) to row
    //detA = |1, -3, 2, -4| = 2
    //detB = |1, -3, 0, 2| = 2
    exp = 2;
    var sclr: f32 = -2;
    var C: [4]f32 = .{ 1, -3, 2, -4 };
    var detC: f32 = 0;
    detC = detXmtx2(&C);
    var D: [4]f32 = .{ 1, -3, 0, 2 };
    addSclMulXmtxRows(0, 1, sclr, 2, &C, &D);
    var detD: f32 = 0;
    detD = detXmtx2(&D);
    try std.testing.expectEqual(true, isEquF32(exp, detC, true));
    try std.testing.expectEqual(true, isEquF32(exp, detD, true));

    //C: Multiply first row of A by 1/2
    //detA = |2, -8, -2, 9| = 2
    //detB = |1, -4, -2, 9| = 1
    sclr = 0.5;
    var E: [4]f32 = .{ 2, -8, -2, 9 };
    var detE: f32 = 0;
    detE = detXmtx2(&E);
    var F: [4]f32 = .{ 1, -4, -2, 9 };
    sclMulXmtxRows(0, sclr, 2, &E, &F);
    var detF: f32 = 0;
    detF = detXmtx2(&F);
    exp = 2.0;
    try std.testing.expectEqual(true, isEquF32(exp, detE, true));
    try std.testing.expectEqual(true, isEquF32(exp * sclr, detF, true));

    //Example 2
    //Find the determinant of, reduce then use diagonal determinant function
    //A = 2,-3,10 ,1,2,-2 ,0,1,-3
    var G: [9]f32 = .{ 2, -3, 10, 1, 2, -2, 0, 1, -3 };
    var retG: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var idtG: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var b: bool = false;
    std.debug.print("AAA START DetG:\n", .{});
    sclr = 0.0;
    b = try rdcXmtx(&G, 3, false, &retG, true, &idtG, 3, true, &sclr);
    try std.testing.expectEqual(true, b);
    var detG: f32 = detTriangXmtx(&retG, 3);
    prntXmtxNl(&retG, 3);
    std.debug.print("AAA STOP DetG: {}\n", .{detG * sclr});
    exp = -7.0; //CHECK THIS why is it 1.0 and not -7.0
    try std.testing.expectEqual(true, isEquF32(exp, (detG * sclr), true));

    const alloc: std.mem.Allocator = std.testing.allocator;
    detG = try detXmtx(&G, 3, &alloc, 0);
    std.debug.print("BBB DetG: {}\n", .{detG});
    try std.testing.expectEqual(true, isEquF32(exp, detG, true));

    //Example 3
    //Find the determinant of, reduce then use diagonal determinant function
    //A = -1,2,2 ,3,-6,4 ,5,-10,-3
    var H: [9]f32 = .{ -1, 2, 2, 3, -6, 4, 5, -10, -3 };
    var retH: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var idtH: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    sclr = 0.0;
    b = try rdcXmtx(&H, 3, false, &retH, true, &idtH, 3, true, &sclr);
    var detH: f32 = 0;

    if (b) {
        try std.testing.expectEqual(true, b);
        detH = detTriangXmtx(&retH, 3);
    } else {
        detH = 0;
    }

    exp = 0.0;
    try std.testing.expectEqual(true, isEquF32(exp, detH, true));
}

test "XMTX: ELA - Larson, Edwards: 3.2 Problem 15, 17, 19, 21 test" {
    prntNl();
    //Evaluate the determinant
    //15 = 1, 7,-3 ,1,3,1 ,4,8,1 = 28
    //17 = 2,-1,-1 ,1,3,2 ,1,1,3 = 17
    //19 = 4,3,-2 ,5,4,1 ,-2,3,4 = -60
    //21 = 5,-8,0 ,9,7,4 ,-8,7,1 = 223
    const cols: usize = 3;
    const alloc: std.mem.Allocator = std.testing.allocator;
    var A: [9]f32 = .{ 1, 7, -3, 1, 3, 1, 4, 8, 1 };
    var detA: f32 = 0;
    detA = try detXmtx(&A, cols, &alloc, 0);
    const expA: f32 = 28;
    try std.testing.expectEqual(true, isEquF32(expA, detA, true));

    var B: [9]f32 = .{ 2, -1, -1, 1, 3, 2, 1, 1, 3 };
    var detB: f32 = 0;
    detB = try detXmtx(&B, cols, &alloc, 0);
    const expB: f32 = 17;
    try std.testing.expectEqual(true, isEquF32(expB, detB, true));

    var C: [9]f32 = .{ 4, 3, -2, 5, 4, 1, -2, 3, 4 };
    var detC: f32 = 0;
    detC = try detXmtx(&C, cols, &alloc, 0);
    const expC: f32 = -60;
    try std.testing.expectEqual(true, isEquF32(expC, detC, true));

    var D: [9]f32 = .{ 5, -8, 0, 9, 7, 4, -8, 7, 1 };
    var detD: f32 = 0;
    detD = try detXmtx(&D, cols, &alloc, 0);
    const expD: f32 = 223;
    try std.testing.expectEqual(true, isEquF32(expD, detD, true));
}

test "XMTX: ELA - Larson, Edwards: 3.3 Example 1, 2, 3, 4 test" {
    prntNl();
    const alloc: std.mem.Allocator = std.testing.allocator;
    const cols: f32 = 3;

    //Find |A|, |B|, and |AB|
    //1
    //A = 1,-2,2 ,0,3,2 ,1,0,1
    //B = 2,0,1 ,0,-1,-2 ,3,1,-2
    var A: [9]f32 = .{ 1, -2, 2, 0, 3, 2, 1, 0, 1 };
    var B: [9]f32 = .{ 2, 0, 1, 0, -1, -2, 3, 1, -2 };
    var retA: [9]f32 = std.mem.zeroes([9]f32); //.{};
    var detA: f32 = 0.0;
    var detRet: f32 = 0.0;
    const expA: f32 = -7.0;
    var detB: f32 = 0.0;
    const expB: f32 = 11.0;

    detA = try detXmtx(&A, cols, &alloc, 0);
    try std.testing.expectEqual(true, isEquF32(expA, detA, true));
    detB = try detXmtx(&B, cols, &alloc, 0);
    try std.testing.expectEqual(true, isEquF32(expB, detB, true));
    const b: bool = try tmsXmtx(&A, cols, &B, cols, &retA, cols);
    try std.testing.expectEqual(true, b);
    detRet = try detXmtx(&retA, cols, &alloc, 0);
    try std.testing.expectEqual(true, isEquF32(expB * expA, detRet, true));

    //Find the determinant
    //2
    //A = 10,-20,40 ,30,0,50 ,-20,-30,10
    var C: [9]f32 = .{ 10, -20, 40, 30, 0, 50, -20, -30, 10 };
    var detC: f32 = 0.0;
    const expC: f32 = 5000;
    detC = try detXmtx(&C, cols, &alloc, 0);
    try std.testing.expectEqual(true, isEquF32(expC, detC, true));

    //Which has an inverse?
    //3
    //A = 0,2,-1 ,3,-2,1 ,3,2,-1
    //B = 0,2,-1 ,3,-2,1 ,3,2,1
    var D: [9]f32 = .{ 0, 2, -1, 3, -2, 1, 3, 2, -1 };
    var E: [9]f32 = .{ 0, 2, -1, 3, -2, 1, 3, 2, 1 };
    var detD: f32 = 0.0;
    var detE: f32 = 0.0;
    const expD: f32 = 0.0;
    const expE: f32 = -12.0;
    detD = try detXmtx(&D, cols, &alloc, 0);
    detE = try detXmtx(&E, cols, &alloc, 0);
    std.debug.print("DetD: {} DetE: {}\n", .{ detD, detE });
    try std.testing.expectEqual(true, isEquF32(expD, detD, true));
    try std.testing.expectEqual(true, isEquF32(expE, detE, true));
}

test "XMTX: ELA - Larson, Edwards: 3.3 Problem 1, 3, 5, 7 test" {
    prntNl();

    //Find |A|, |B|, AB, |AB| then verify that |A||B| = |AB|
    //1
    //A = -2,1 ,4,-2
    //B = 1,1 ,0,-1
    const alloc: std.mem.Allocator = std.testing.allocator;
    var b: bool = false;
    var A: [4]f32 = .{ -2, 1, 4, -2 };
    var B: [4]f32 = .{ 1, 1, 0, -1 };
    var AB: [4]f32 = .{ 0, 0, 0, 0 };
    var detA: f32 = 0.0;
    var detB: f32 = 0.0;
    const expA: f32 = 0.0;
    const expB: f32 = -1.0;
    var detAB: f32 = 0.0;
    const expAB: f32 = 0.0;

    detA = try detXmtx(&A, 2, &alloc, 0);
    try std.testing.expectEqual(true, isEquF32(expA, detA, true));
    detB = try detXmtx(&B, 2, &alloc, 0);
    try std.testing.expectEqual(true, isEquF32(expB, detB, true));
    b = try tmsXmtx(&A, 2, &B, 2, &AB, 2);
    try std.testing.expectEqual(true, b);
    detAB = try detXmtx(&AB, 2, &alloc, 0);
    try std.testing.expectEqual(true, isEquF32(expAB, detAB, true));

    //3
    //C = -1,2,1 ,1,0,1 ,0,1,0
    //D = -1,0,0 ,0,2,0 ,0,0,3
    var C: [9]f32 = .{ -1, 2, 1, 1, 0, 1, 0, 1, 0 };
    var D: [9]f32 = .{ -1, 0, 0, 0, 2, 0, 0, 0, 3 };
    var CD: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var detC: f32 = 0.0;
    var detD: f32 = 0.0;
    var detCD: f32 = 0.0;
    const expC: f32 = 2.0;
    const expD: f32 = -6.0;
    const expCD: f32 = -12.0;

    detC = try detXmtx(&C, 3, &alloc, 0);
    try std.testing.expectEqual(true, isEquF32(expC, detC, true));
    detD = try detXmtx(&D, 3, &alloc, 0);
    try std.testing.expectEqual(true, isEquF32(expD, detD, true));
    b = try tmsXmtx(&C, 3, &D, 3, &CD, 3);
    try std.testing.expectEqual(true, b);
    detCD = try detXmtx(&CD, 3, &alloc, 0);
    try std.testing.expectEqual(true, isEquF32(expCD, detCD, true));

    //Given that |cA| = c^n * |A| for an n X n matrix
    //5
    //A = 4,2 ,6,-8
    var sclr: f32 = 2.0;
    var sclrN: f32 = 4.0;
    var E: [4]f32 = .{ 4 / sclr, 2 / sclr, 6 / sclr, -8 / sclr };
    var detE: f32 = 0.0;
    const expE: f32 = -44.0;

    detE = try detXmtx(&E, 2, &alloc, 0);
    try std.testing.expectEqual(true, isEquF32(expE, detE * sclrN, true));

    //7
    //A = -3,6,9 ,6,9,12 ,9,12,15
    sclr = 3.0;
    sclrN = 27.0;
    var F: [9]f32 = .{ -3 / sclr, 6 / sclr, 9 / sclr, 6 / sclr, 9 / sclr, 12 / sclr, 9 / sclr, 12 / sclr, 15 / sclr };
    var detF: f32 = 0.0;
    const expF: f32 = 54.0;

    detF = try detXmtx(&F, 3, &alloc, 0);
    try std.testing.expectEqual(true, isEquF32(expF, detF * sclrN, true));
}

test "XMTX: ELA - Larson, Edwards: 3.4 Example 1, 3, 5 test" {
    prntNl();

    //Example 1: pg. 138
    //Let A = |1,4 ,2,3|   x1 = |1, 1| x2 = |2, -1|
    //Verify that l1 is an eigenvalue of A corresponding to x1 and that l2 = -1 is an eigenvalue of A corresponding to x2.
    var mtx: [4]f32 = .{ 1, 4, 2, 3 };
    var evs2: EigVal2 = .{};
    var b: bool = false;
    const expEv1: f32 = 5.0;
    const expEv2: f32 = -1.0;

    b = fndEigVal2(&mtx, &evs2);
    try std.testing.expectEqual(true, b);
    std.debug.print("Found Eigen Values:\n", .{});
    prntXvecNl(&evs2.eignVals);
    try std.testing.expectEqual(true, isEquF32(expEv1, evs2.eignVals[0], true));
    try std.testing.expectEqual(true, isEquF32(expEv2, evs2.eignVals[1], true));

    //Example 3: pg. 140
    //Find the adjoint of the matrix A = |-1,3,2 ,0,-2,1 ,1,0,-2|
    var A: [9]f32 = .{ -1, 3, 2, 0, -2, 1, 1, 0, -2 };
    var adjA: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var expAdjA: [9]f32 = .{ 4, 6, 7, 1, 0, 1, 2, 3, 2 };

    adjA = adjXmtx3(&A);
    std.debug.print("Matrix A:\n", .{});
    prntXmtxNl(&A, 3);
    std.debug.print("Matrix adjA:\n", .{});
    prntXmtxNl(&adjA, 3);
    std.debug.print("Matrix Expected adjA:\n", .{});
    prntXmtxNl(&expAdjA, 3);
    try std.testing.expectEqual(true, equXmtx(&adjA, &expAdjA));

    //Example 5: pg. 144
    //Use Cramer's rule to solve the following system of equations: B = |4,-2,10 ,3,-5,11|
    var B: [6]f32 = .{ 4, -2, 10, 3, -5, 11 };
    const cols: usize = 3;
    var BA: [4]f32 = .{ 4, -2, 3, -5 };
    const colsA: usize = 2;
    var BAi: [4]f32 = .{ 0, 0, 0, 0 };
    var Bret: [2]f32 = .{ 0, 0 };
    const expX1: f32 = 2;
    const expX2: f32 = -1;

    b = try rslvCramersRule(&B, cols, &BA, colsA, &BAi, &Bret);
    try std.testing.expectEqual(true, b);
    std.debug.print("Found Cramer's Rule Return Values:\n", .{});
    prntXvecNl(&Bret);
    try std.testing.expectEqual(true, isEquF32(expX1, Bret[0], true));
    try std.testing.expectEqual(true, isEquF32(expX2, Bret[1], true));
}

test "XMTX: ELA - Larson, Edwards: 3.4 Theorem 10 test" {
    prntNl();

    //If A is an n * n matrix, then
    //A^-1 = ((1 / detA) * adjA)

    var A: [9]f32 = .{ -1, 3, 2, 0, -2, 1, 1, 0, -2 };
    var adjA: [9]f32 = std.mem.zeroes([9]f32); //.{};
    var detA: f32 = 0.0;
    var expA: [9]f32 = .{ 4, 6, 7, 1, 0, 1, 2, 3, 2 };
    var invA: [9]f32 = std.mem.zeroes([9]f32); //.{};
    var expInvA: [9]f32 = .{ (4.0 / 3.0), 2, (7.0 / 3.0), (1.0 / 3.0), 0, (1.0 / 3.0), (2.0 / 3.0), 1, (2.0 / 3.0) };

    adjA = adjXmtx3(&A);
    detA = detXmtx3(&A);
    cpyXmtx(&adjA, &invA);
    mulXvec(&invA, (1.0 / detA));
    try std.testing.expectEqual(true, equXmtx(&adjA, &expA));

    std.debug.print("Adjoint Matrix Test:\n", .{});
    prntXmtxNl(&invA, 3);
    try std.testing.expectEqual(true, equXmtx(&invA, &expInvA));
}

test "XMTX: ELA - Larson, Edwards: 3.4 Theorem 11 test" {
    prntNl();
    //A = -1  2 -3  1
    //     2  0  1  0
    //     3 -4  4  2

    var mtx: [12]f32 = .{ -1, 2, -3, 1, 2, 0, 1, 0, 3, -4, 4, 2 };
    var A: [9]f32 = .{ -1, 2, -3, 2, 0, 1, 3, -4, 4 };
    var Ai: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var ret: [3]f32 = .{ 0, 0, 0 };
    var b: bool = false;
    const cols: usize = 4;
    const colsA: usize = 3;

    b = try rslvCramersRule(&mtx, cols, &A, colsA, &Ai, &ret);
    try std.testing.expectEqual(true, b);

    const expX: f32 = 4.0 / 5.0;
    const expY: f32 = -3.0 / 2.0;
    const expZ: f32 = -8.0 / 5.0;

    std.debug.print("Cramer's Rule Return {} {} {}:\n", .{ expX, expY, expZ });
    prntXvecNl(&ret);
    prntNl();
    try std.testing.expectEqual(true, isEquF32(ret[0], expX, true));
    try std.testing.expectEqual(true, isEquF32(ret[1], expY, true));
    try std.testing.expectEqual(true, isEquF32(ret[2], expZ, true));
}

test "XMTX: ELA - Larson, Edwards: 3.4 Problem 1, 3, 5, 7 test" {
    prntNl();

    //Problem 1: pg. 150
    //Verify that l1 is an eigenvalue of A and that xi is a corresponding eigenvector:
    //A = |1,2 ,0,-3|   l1 = 1, xi = |1,0|
    var A1: [4]f32 = .{ 1, 2, 0, -3 };
    const l1: f32 = 1;
    var xi1: [2]f32 = .{ 1, 0 };
    var evs1: EigVal2 = .{};
    var b: bool = false;

    b = fndEigVal2(&A1, &evs1);
    try std.testing.expectEqual(true, b);
    try std.testing.expectEqual(true, isEquF32(l1, evs1.eignVals[0], true));

    //A1 * xi1 = l1 * xi1
    //|1  2| * |1  0| = |1  0| => 1 * |1  0| = |1  0|
    //|0  3|
    var ret11: [2]f32 = .{ 0, 0 };
    b = try tmsXmtx(&A1, 2, &xi1, 1, &ret11, 1);
    var ret12: [2]f32 = .{ 1, 0 };
    mulXvec(&ret12, evs1.eignVals[0]);
    try std.testing.expectEqual(true, equXmtx(&ret11, &ret12));

    //Problem 3: pg. 150
    //Verify that l1 is an eigenvalue of A and that xi is a corresponding eigenvector:
    //A = |1,1,1 ,0,1,0 ,1,1,1|   l1 = 2, xi = |1,0,1|
    var A2: [9]f32 = .{ 1, 1, 1, 0, 1, 0, 1, 1, 1 };
    const l2: f32 = 2;
    var xi2: [3]f32 = .{ 1, 0, 1 };
    var evs2: EigVal3 = .{};
    var factors: [1000]f32 = std.mem.zeroes([1000]f32); //.{};
    b = false;

    b = fndEigVal3(&A2, &evs2, &factors);
    try std.testing.expectEqual(true, b);
    std.debug.print("Problem 3: pg. 150 Eigen Values:\n", .{});
    prntXvecNl(&evs2.eignVals);
    try std.testing.expectEqual(true, isEquF32(l2, evs2.eignVals[1], true));

    //A2 * xi2 = l2 * xi2
    //|1  1  1| * |1  0  1| = |2  0  2| => 2 * |1  0  1| = |2  0  2|
    //|0  1  0|
    //|1  1  1|
    var ret21: [3]f32 = .{ 0, 0, 0 };
    b = try tmsXmtx(&A2, 3, &xi2, 1, &ret21, 1);
    var ret22: [3]f32 = .{ 1, 0, 1 };
    mulXvec(&ret22, evs2.eignVals[1]);
    try std.testing.expectEqual(true, equXmtx(&ret21, &ret22));

    //Problem 5: pg. 150
    //Find the a) characteristic equation b) eigenvalues, and c) eigenvectors of the given matrix:
    //A = |4,-5 ,2,-3|
    var A3: [4]f32 = .{ 4, -5, 2, -3 };
    const l3a: f32 = 2.0;
    const l3b: f32 = -1.0;
    var xi3a: [2]f32 = .{ 5, 2 };
    var xi3b: [2]f32 = .{ 1, 1 };
    var evs3: EigVal2 = .{};
    b = false;

    b = fndEigVal2(&A3, &evs3);
    try std.testing.expectEqual(true, b);
    std.debug.print("Problem 5: pg. 150 a) Characteristic Equation:\n", .{});
    prntPolyExp(&evs3.polyExp);
    std.debug.print("Problem 5: pg. 150 b) Eigen Values:\n", .{});
    prntXvecNl(&evs3.eignVals);
    try std.testing.expectEqual(true, isEquF32(l3a, evs3.eignVals[0], true));
    try std.testing.expectEqual(true, isEquF32(l3b, evs3.eignVals[1], true));

    //|4  -5| * |x1  = 2 * |x1
    //|2  -3|    x2|        x2|
    //2x1 - 5x2 = 0 => |5 -2|
    //2x1 - 5x2 = 0 => |5 -2|
    var ret31: [2]f32 = .{ 0, 0 };
    b = try tmsXmtx(&A3, 2, &xi3a, 1, &ret31, 1);
    var ret32: [2]f32 = .{ 5, 2 };
    mulXvec(&ret32, evs3.eignVals[0]);
    try std.testing.expectEqual(true, equXmtx(&ret31, &ret32));

    //|4  -5| * |x1  = -1 * |x1
    //|2  -3|    x2|        x2|
    //5x1 - 5x2 = 0 => |1  1|
    //2x1 - 2x2 = 0 => |1  1|
    ret31 = .{ 0, 0 };
    b = try tmsXmtx(&A3, 2, &xi3b, 1, &ret31, 1);
    ret32 = .{ 1, 1 };
    mulXvec(&ret32, evs3.eignVals[1]);
    try std.testing.expectEqual(true, equXmtx(&ret31, &ret32));

    //Problem 7: pg. 150
    //Find the a) characteristic equation b) eigenvalues, and c) eigenvectors of the given matrix:
    //A = |2,1 ,3,0|
    var A4: [4]f32 = .{ 2, 1, 3, 0 };
    const l4a: f32 = 3.0;
    const l4b: f32 = -1.0;
    var xi4a: [2]f32 = .{ 1, 1 };
    var xi4b: [2]f32 = .{ -1, 3 };
    var evs4: EigVal2 = .{};
    b = false;

    b = fndEigVal2(&A4, &evs4);
    try std.testing.expectEqual(true, b);
    std.debug.print("Problem 7: pg. 150 a) Characteristic Equation:\n", .{});
    prntPolyExp(&evs4.polyExp);
    std.debug.print("Problem 7: pg. 150 b) Eigen Values:\n", .{});
    prntXvecNl(&evs4.eignVals);
    try std.testing.expectEqual(true, isEquF32(l4a, evs4.eignVals[0], true));
    try std.testing.expectEqual(true, isEquF32(l4b, evs4.eignVals[1], true));

    //|2  1| * |x1  = 3 * |x1
    //|3  0|    x2|        x2|
    //2x1 + 1x2 = 3x1 => 2x1 - 3x1 + 1x2 = 0 => -1x1 + 1x2 = 0 => -1 * (1x1 - 1x2) = 0 => |1  1|
    //3x1 + 0x2 = 3x2 => 3x1 + 0x2 - 3x2 = 0 => 3x1 - 3x2 = 0 => 3 * (1x1 - 1x2) = 0 => |1  1|
    var ret41: [2]f32 = .{ 0, 0 };
    b = try tmsXmtx(&A4, 2, &xi4a, 1, &ret41, 1);
    var ret42: [2]f32 = .{ 1, 1 };
    mulXvec(&ret42, evs4.eignVals[0]);
    try std.testing.expectEqual(true, equXmtx(&ret41, &ret42));

    //|2  1| * |x1  = -1 * |x1
    //|3  0|    x2|        x2|
    //2x1 + 1x2 = -1x1 => 3x1 + 1x2 = 0 => |-1  3|
    //3x1 + 0x2 = -1x2 => 3x1 + 1x2 = 0 => |-1  3|
    ret41 = .{ 0, 0 };
    b = try tmsXmtx(&A4, 2, &xi4b, 1, &ret41, 1);
    ret42 = .{ -1, 3 };
    mulXvec(&ret42, evs4.eignVals[1]);
    try std.testing.expectEqual(true, equXmtx(&ret41, &ret42));
}

test "XMTX: ELA - Larson, Edwards: 3.4 Problem 15, 17, 29, 31 test" {
    prntNl();

    //Problem 15: pg. 151
    //Find the adjoint of the given matrix A then use the adjoint to find the determinant of A if possible.
    //A = |1,0,0 ,0,2,6 ,0,-4,-12|
    var A: [9]f32 = .{ 1, 0, 0, 0, 2, 6, 0, -4, -12 };
    var adjA: [9]f32 = std.mem.zeroes([9]f32); //.{};
    adjA = adjXmtx3(&A);
    var detA: f32 = detXmtx3(&A);
    var expInvA: [9]f32 = std.mem.zeroes([9]f32); //.{};
    var ODetA: f32 = 0;
    std.debug.print("Problem 15: pg. 151: detA:{} \n", .{detA});
    if (detA != 0) {
        std.debug.print("Possible\n", .{});
        ODetA = 1.0 / detA;
        mulXvec(&adjA, ODetA);
        try std.testing.expectEqual(true, equXmtx(&adjA, &expInvA));
    } else {
        std.debug.print("NOT Possible\n", .{});
    }

    //Problem 17: pg. 151
    //Find the adjoint of the given matrix A then use the adjoint to find the determinant of A if possible.
    //A = |-3,-5,-7 ,2,4,3 ,0,1,-1|
    A = .{ -3, -5, -7, 2, 4, 3, 0, 1, -1 };
    adjA = std.mem.zeroes([9]f32); //.{};
    adjA = adjXmtx3(&A);
    detA = detXmtx3(&A);
    expInvA = .{ 7.0 / 3.0, 4, -13.0 / 3.0, -2.0 / 3.0, -1, 5.0 / 3.0, -2.0 / 3.0, -1, 2.0 / 3.0 };
    ODetA = 0;
    std.debug.print("Problem 17: pg. 151: detA:{}\n", .{detA});
    if (detA != 0) {
        std.debug.print("Possible\n", .{});
        ODetA = 1.0 / detA;
        mulXvec(&adjA, ODetA);
        try std.testing.expectEqual(true, equXmtx(&adjA, &expInvA));
    } else {
        std.debug.print("NOT Possible\n", .{});
    }

    //Problem 29: pg. 151
    //Use Cramer's rule to solve the given system of linear equations if possible.
    // 1x1 +  2x2 =  5
    //-1x1 +  0x2 =  1
    var B: [6]f32 = .{ 1, 2, 5, -1, 1, 1 };
    const colsB: usize = 3;
    var BA: [4]f32 = .{ 1, 2, -1, 1 };
    const colsBA: usize = 2;
    var BAi: [4]f32 = std.mem.zeroes([4]f32); //.{};
    var retB: [2]f32 = .{ 0, 0 };
    var b: bool = false;
    const expB: [2]f32 = .{ 1, 2 };
    b = try rslvCramersRule(&B, colsB, &BA, colsBA, &BAi, &retB);
    std.debug.print("Problem 29: pg. 151 Values:\n", .{});
    prntXvecNl(&retB);
    try std.testing.expectEqual(true, isEquF32(expB[0], retB[0], true));
    try std.testing.expectEqual(true, isEquF32(expB[1], retB[1], true));

    //Problem 31: pg. 151
    //Use Cramer's rule to solve the given system of linear equations if possible.
    // 3x1 +  4x2 = -2
    // 5x1 +  3x2 =  4
    var C: [6]f32 = .{ 3, 4, -2, 5, 3, 4 };
    const colsC: usize = 3;
    var CA: [4]f32 = .{ 3, 4, 5, 3 };
    const colsCA: usize = 2;
    var CAi: [4]f32 = std.mem.zeroes([4]f32); //.{};
    var retC: [2]f32 = .{ 0, 0 };
    b = false;
    const expC: [2]f32 = .{ 2, -2 };
    b = try rslvCramersRule(&C, colsC, &CA, colsCA, &CAi, &retC);
    std.debug.print("Problem 31: pg. 151 Values:\n", .{});
    prntXvecNl(&retC);
    try std.testing.expectEqual(true, isEquF32(expC[0], retC[0], true));
    try std.testing.expectEqual(true, isEquF32(expC[1], retC[1], true));
}

test "XMTX: ELA - Larson, Edwards: 4.1 Theorem 1 test" {
    prntNl();

    //Properties of vector addition and scalar multiplication
    //Let u, v, and w be vectors in the plane, and let c and d be scalars.

    //ADDITION
    //1: u + v is a vector in the plane, closure under addition
    var u: [3]f32 = .{ 1, 1, 0 };
    var negU: [3]f32 = .{ -1, -1, 0 };
    var v: [3]f32 = .{ 2, 2, 0 };
    var uPv: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&uPv, &u, &v);
    try std.testing.expectEqual(true, isEquF32(0.0, uPv[2], true));

    //2: u + v = v + u, commutative property of addition
    var vPu: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&vPu, &v, &u);
    try std.testing.expectEqual(true, equXmtx(&uPv, &vPu));

    //3: (u + v) + w = u + (v + w), associative property of addition
    var w: [3]f32 = .{ 3, 3, 0 };
    var vPw: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&vPw, &v, &w);
    var uPvPPw: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&uPvPPw, &uPv, &w);
    var uPPvPw: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&uPPvPw, &u, &vPw);
    try std.testing.expectEqual(true, equXmtx(&uPvPPw, &uPPvPw));

    //4: u + 0 = u, additive identity
    var zer: [3]f32 = .{ 0, 0, 0 };
    var uPzer: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&uPzer, &u, &zer);
    try std.testing.expectEqual(true, equXmtx(&uPzer, &u));

    //5: u + (-u) = 0, additive inverse
    var uPnegU: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&uPnegU, &u, &negU);
    try std.testing.expectEqual(true, equXmtx(&zer, &uPnegU));

    //MULTIPLICATION
    //6: c * u, is a vector on the plane, closure under scalar multiplication
    const c: f32 = 3.0;
    var cTu: [3]f32 = .{ 0, 0, 0 };
    cpyXmtx(&u, &cTu);
    mulXvec(&cTu, c);
    try std.testing.expectEqual(true, isEquF32(0.0, cTu[2], true));

    //7: c * (u + v) = c * u + c * v, distributive property
    var cTuPv: [3]f32 = .{ 0, 0, 0 };
    cpyXmtx(&uPv, &cTuPv);
    mulXvec(&cTuPv, c);
    var cTv: [3]f32 = .{ 0, 0, 0 };
    cpyXmtx(&v, &cTv);
    mulXvec(&cTv, c);
    var cTuPcTv: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&cTuPcTv, &cTu, &cTv);
    try std.testing.expectEqual(true, equXmtx(&cTuPcTv, &cTuPv));

    //8: (c + d) * u = c * u + d * u, distributive property
    const d: f32 = 2.0;
    const cPd: f32 = c + d;
    var cPdTu: [3]f32 = .{ 0, 0, 0 };
    cpyXmtx(&u, &cPdTu);
    mulXvec(&cPdTu, cPd);
    var dTu: [3]f32 = .{ 0, 0, 0 };
    cpyXmtx(&u, &dTu);
    mulXvec(&dTu, d);
    var dTuPcTu: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&dTuPcTu, &cTu, &dTu);
    try std.testing.expectEqual(true, equXmtx(&cPdTu, &dTuPcTu));

    //9: c * (d * u) - (c * d) * u, scalar associative property
    var cTTdTu: [3]f32 = .{ 0, 0, 0 };
    cpyXmtx(&dTu, &cTTdTu);
    mulXvec(&cTTdTu, c);
    const cTd: f32 = c * d;
    var cTdTu: [3]f32 = .{ 0, 0, 0 };
    cpyXmtx(&u, &cTdTu);
    mulXvec(&cTdTu, cTd);
    try std.testing.expectEqual(true, equXmtx(&cTdTu, &cTTdTu));

    //10: 1 * u = u, scalar multiplicative identity
    const one: f32 = 1.0;
    var oneTu: [3]f32 = .{ 0, 0, 0 };
    cpyXmtx(&u, &oneTu);
    mulXvec(&oneTu, one);
    try std.testing.expectEqual(true, equXmtx(&u, &oneTu));
}

test "XMTX: ELA - Larson, Edwards: 4.1 Theorem 2 test" {
    prntNl();

    //Properties of vector addition and scalar multiplication
    //Let u, v, and w be vectors in Rn, and let c and d be scalars.

    //ADDITION
    //1: u + v is a vector in the plane, closure under addition
    var u: [3]f32 = .{ 1, 1, 0 };
    var negU: [3]f32 = .{ -1, -1, 0 };
    var v: [3]f32 = .{ 2, 2, 0 };
    var uPv: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&uPv, &u, &v);
    try std.testing.expectEqual(true, isEquF32(0.0, uPv[2], true));

    //2: u + v = v + u, commutative property of addition
    var vPu: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&vPu, &v, &u);
    try std.testing.expectEqual(true, equXmtx(&uPv, &vPu));

    //3: (u + v) + w = u + (v + w), associative property of addition
    var w: [3]f32 = .{ 3, 3, 0 };
    var vPw: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&vPw, &v, &w);
    var uPvPPw: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&uPvPPw, &uPv, &w);
    var uPPvPw: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&uPPvPw, &u, &vPw);
    try std.testing.expectEqual(true, equXmtx(&uPvPPw, &uPPvPw));

    //4: u + 0 = u, additive identity
    var zer: [3]f32 = .{ 0, 0, 0 };
    var uPzer: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&uPzer, &u, &zer);
    try std.testing.expectEqual(true, equXmtx(&uPzer, &u));

    //5: u + (-u) = 0, additive inverse
    var uPnegU: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&uPnegU, &u, &negU);
    try std.testing.expectEqual(true, equXmtx(&zer, &uPnegU));

    //MULTIPLICATION
    //6: c * u, is a vector on the plane, closure under scalar multiplication
    const c: f32 = 3.0;
    var cTu: [3]f32 = .{ 0, 0, 0 };
    cpyXmtx(&u, &cTu);
    mulXvec(&cTu, c);
    try std.testing.expectEqual(true, isEquF32(0.0, cTu[2], true));

    //7: c * (u + v) = c * u + c * v, distributive property
    var cTuPv: [3]f32 = .{ 0, 0, 0 };
    cpyXmtx(&uPv, &cTuPv);
    mulXvec(&cTuPv, c);
    var cTv: [3]f32 = .{ 0, 0, 0 };
    cpyXmtx(&v, &cTv);
    mulXvec(&cTv, c);
    var cTuPcTv: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&cTuPcTv, &cTu, &cTv);
    try std.testing.expectEqual(true, equXmtx(&cTuPcTv, &cTuPv));

    //8: (c + d) * u = c * u + d * u, distributive property
    const d: f32 = 2.0;
    const cPd: f32 = c + d;
    var cPdTu: [3]f32 = .{ 0, 0, 0 };
    cpyXmtx(&u, &cPdTu);
    mulXvec(&cPdTu, cPd);
    var dTu: [3]f32 = .{ 0, 0, 0 };
    cpyXmtx(&u, &dTu);
    mulXvec(&dTu, d);
    var dTuPcTu: [3]f32 = .{ 0, 0, 0 };
    sum2Xvec(&dTuPcTu, &cTu, &dTu);
    try std.testing.expectEqual(true, equXmtx(&cPdTu, &dTuPcTu));

    //9: c * (d * u) - (c * d) * u, scalar associative property
    var cTTdTu: [3]f32 = .{ 0, 0, 0 };
    cpyXmtx(&dTu, &cTTdTu);
    mulXvec(&cTTdTu, c);
    const cTd: f32 = c * d;
    var cTdTu: [3]f32 = .{ 0, 0, 0 };
    cpyXmtx(&u, &cTdTu);
    mulXvec(&cTdTu, cTd);
    try std.testing.expectEqual(true, equXmtx(&cTdTu, &cTTdTu));

    //10: 1 * u = u, scalar multiplicative identity
    const one: f32 = 1.0;
    var oneTu: [3]f32 = .{ 0, 0, 0 };
    cpyXmtx(&u, &oneTu);
    mulXvec(&oneTu, one);
    try std.testing.expectEqual(true, equXmtx(&u, &oneTu));
}

test "XMTX: ELA - Larson, Edwards: 4.1 Theorem 3 test" {
    prntNl();

    //Properties of additive identity and additive inverse
    //Let v be a vector in Rn, and let c be a scalar. Then the following properties are true.
    //1. The additive identity is unique. If v + u = v then u = 0.
    var a: [2]f32 = .{ 1, 2 };
    var b: [2]f32 = .{ 0, 0 };
    var c: [2]f32 = .{ 0, 0 };
    var zero: [2]f32 = .{ 0, 0 };
    var expC: [2]f32 = .{ 1, 2 };
    var res: bool = false;

    sum2Xvec(&c, &a, &b);
    res = equXvec(&expC, &c);
    try std.testing.expectEqual(true, res);
    res = equXvec(&zero, &b);
    try std.testing.expectEqual(true, res);

    //2. The additive inverse of v is unique. If v + u = 0, then u = -v.
    b = .{ -1, -2 };
    c = .{ 0, 0 };
    expC = .{ 0, 0 };
    res = false;

    sum2Xvec(&c, &a, &b);
    res = equXvec(&expC, &c);
    try std.testing.expectEqual(true, res);
    var d: [2]f32 = .{ 0, 0 };
    cpyXvec(&a, &d);
    mulXvec(&d, -1);
    res = equXvec(&d, &b);
    try std.testing.expectEqual(true, res);

    //3. 0v = 0
    const sclZer: f32 = 0.0;
    var e: [2]f32 = .{ 3, 4 };
    mulXvec(&e, sclZer);
    res = equXvec(&zero, &e);
    try std.testing.expectEqual(true, res);

    //4. c0 = 0
    const sclr: f32 = 3.14;
    var zero2: [2]f32 = .{ 0, 0 };
    mulXvec(&zero2, sclr);
    res = equXvec(&zero2, &zero);
    try std.testing.expectEqual(true, res);

    //5. If cv = 0, then c = 0 or v = 0
    const sclZ: f32 = 7.77;
    var z: [2]f32 = .{ 7, 77 };
    mulXvec(&z, sclZer);
    mulXvec(&zero, sclZ);
    res = equXvec(&zero, &z);
    try std.testing.expectEqual(true, res);

    //6. -(-v) = v
    var V: [2]f32 = .{ 10, 20 };
    var negV: [2]f32 = .{ -10, -20 };
    const sclNeg1: f32 = -1;
    mulXvec(&V, sclNeg1);
    res = equXvec(&V, &negV);
    try std.testing.expectEqual(true, res);
}

test "XMTX: ELA - Larson, Edwards: 4.2 Example 3 test" {
    prntNl();

    //Show that the set of all 2x3 matrices with the operations of matrix addition and scalar
    //multiplication is a vector space.

    var A: [6]f32 = .{ 1, 2, 3, 4, 5, 6 };
    var B: [6]f32 = .{ 2, 4, 6, 8, 10, 12 };
    var C: [6]f32 = .{ 0, 0, 0, 0, 0, 0 };

    //Closure under vector addition
    sum2Xvec(&C, &B, &A);
    try std.testing.expectEqual(true, isEquF32(C.len, A.len, true));

    //Closure under vector multiplication
    mulXvec(&A, 2);
    try std.testing.expectEqual(true, equXmtx(&A, &B));
}

test "XMTX: ELA - Larson, Edwards: 4.3 Theorem 4 test" {
    prntNl();

    //Properties of scalar multiplication
    //Let v be any element of a vector space V, let c be any scalar. Then the following properties are true:

    var A: [2]f32 = .{ 5, 10 };
    const c: f32 = 3.14;
    var ZERO: [2]f32 = .{ 0, 0 };
    const sclZero: f32 = 0.0;

    //1. 0 * v = zero, where v and zero are vectors and 0 is a scalar
    mulXvec(&A, sclZero);
    try std.testing.expectEqual(true, equXmtx(&A, &ZERO));

    //2. c * zero = zero, where c is a scalar and zero is a vector
    A = .{ 0, 0 };
    mulXvec(&A, c);
    try std.testing.expectEqual(true, equXmtx(&A, &ZERO));

    //3. If c * v = 0, then c = 0 or v = 0
    A = .{ 0, 0 };
    mulXvec(&A, c);
    try std.testing.expectEqual(true, equXmtx(&A, &ZERO));

    //4. -1 * v = -v
    const nOne: f32 = -1;
    A = .{ 5, 10 };
    var nA: [2]f32 = .{ -5, -10 };
    mulXvec(&A, nOne);
    try std.testing.expectEqual(true, equXmtx(&nA, &A));
}

test "XMTX: ELA - Larson, Edwards: 4.3 Theorem 5 test" {
    prntNl();

    //If W is a nonempty subset of a vector space V, then W is a subspace of V
    //if and only if the following closure conditions hold.

    //1. If u and v are in W, then u + v is in W.
    //Let W be the vector space of all singular 1x2 matrices where the second value is zero.
    var A: [2]f32 = .{ 1, 0 };
    var B: [2]f32 = .{ 5, 0 };
    sum1Xvec(&A, &B);
    try std.testing.expectEqual(true, isEquF32(A[1], B[1], true));

    //2. If u is in W and c is any scalar, then c * u is in W
    //Let W be the vector space of all singular 1x2 matrices where the second value is zero.
    var U: [2]f32 = .{ 10, 0 };
    const sclr: f32 = 3.14;
    mulXvec(&U, sclr);
    try std.testing.expectEqual(true, isEquF32(U[1], B[1], true));
}

test "XMTX: ELA - Larson, Edwards: 4.3 Example 2 test" {
    prntNl();

    //Example 2: pg. 179
    //Let W be the set of all 2x2 symmetic matrices. Show that W is a subspace of the vector space M2x2,
    //with the standard operations of matrix addition and scalar multilication.
    var A: [4]f32 = .{ 1, 0, 0, 1 };
    var B: [4]f32 = .{ 5, 0, 0, 5 };
    var C: [4]f32 = .{ 0, 0, 0, 0 };
    var b: bool = false;

    sum2Xvec(&C, &A, &B);
    b = isSymXmtx(&A, 2);
    try std.testing.expectEqual(true, b);
    b = isSymXmtx(&B, 2);
    try std.testing.expectEqual(true, b);
    b = isSymXmtx(&C, 2);
    try std.testing.expectEqual(true, b);

    const sclr: f32 = 4.456;
    mulXvec(&C, sclr);
    b = isSymXmtx(&C, 2);
    try std.testing.expectEqual(true, b);
}

test "XMTX: ELA - Larson, Edwards: 4.3 Theorem 6 test" {
    prntNl();

    //The intersection of two subspaces is a subspace
    //If V and W are both subspaces of a vector space U, then the intersection of V and W
    //(denoted by V int W) is also a subspace of U.

    //Subspace 1
    //1. If u and v are in W, then u + v is in W.
    //Let W be the vector space of all singular 1x2 matrices where the second value is zero.
    var A: [2]f32 = .{ 1, 0 };
    var B: [2]f32 = .{ 5, 0 };
    sum1Xvec(&A, &B);
    try std.testing.expectEqual(true, isEquF32(A[1], B[1], true));

    //2. If u is in W and c is any scalar, then c * u is in W
    //Let W be the vector space of all singular 1x2 matrices where the second value is zero.
    var U: [2]f32 = .{ 10, 0 };
    const sclr: f32 = 3.14;
    mulXvec(&U, sclr);
    try std.testing.expectEqual(true, isEquF32(U[1], B[1], true));

    //Subspace 2
    //Let the second subspace be the null space
    //V int W = null space
    var C: [2]f32 = .{ 0, 0 };
    var E: [2]f32 = .{ 0, 0 };
    const D: f32 = 1.234;

    sum1Xvec(&C, &E);
    try std.testing.expectEqual(true, isZeroXvec(&C));
    mulXvec(&C, D);
    try std.testing.expectEqual(true, isZeroXvec(&C));
}

test "XMTX: ELA - Larson, Edwards: 4.4 Example 8, 9, 10 test" {
    prntNl();

    //Example 8: pg. 191
    //Determine whether the following set of vectors in R3
    //is linearly dependent or independent.
    //S = {
    //       (1, 2, 3)    v1
    //       (0, 1, 2)    v2
    //       (-2, 0, 1)   v3
    //    }
    //c1*v1 + c2*v2 + c3*v3 = 0
    //1c1 + 0c2 - 2c3 = 0
    //2c1 + 1c2 + 0c3 = 0
    //3c1 + 2c2 + 1c3 = 0

    var mtx: [12]f32 = .{ 1, 0, -2, 0, 2, 1, 0, 0, 3, 2, 1, 0 };
    const cols: usize = 4;
    var hasAug: bool = true;
    var ret: [12]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const hasIdtMtx: bool = false;
    var idtMtx: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    const dim: usize = 3;
    const triagRdcOnly: bool = false;
    var sclr: f32 = 0.0;
    var b: bool = false;

    //Reduce the provided matrix to reduced row escelon form using Gauss-Jordan Elimination and optionaly calculate the matrix inverse.
    //  mtx = The matrix to reduce.
    //  cols = The number of columns in the matrix.
    //  hasAug = A Boolean indicating if the matrix is an augmented matrix.
    //  ret = The matrix that holds the reduced matrix.
    //  hasIdtMtx = A Boolean indicating if the identity matrix has been provided.
    //  idtMtx = The identity matrix associated with the mtx matrix provided.
    //  dim = The number of matrix columns that must be zero for a zero row to exist.
    //  triagRdcOnly = A Boolean value indicating if the reduction operation should stop when the matrix is triangular.
    //  sclr = A pointer to a floating point variable that keeps track of the scalar multiplication performed against the matrix, mtx.
    //  returns = A Boolean value indicating  if the matrix was reduced successfuly.
    b = try rdcXmtx(&mtx, cols, hasAug, &ret, hasIdtMtx, &idtMtx, dim, triagRdcOnly, &sclr);
    try std.testing.expectEqual(true, b);
    std.debug.print("RdcMtx Result Example 8: {}\n", .{b});
    prntXmtxNl(&ret, 4);
    prntNl();

    //Example 9: pg. 192
    //Determine whether the following set of vectors in P2 is linearly independent
    //or linearly dependent.
    //S = {
    //       (1 + 1x - 2x^2)  v1
    //       (2 + 5x - 1x^2)  v2
    //       (0 + 1x + 1x^2)  v3
    //    }
    //c1*v1 + c2*v2 + c3*v3 = 0
    //
    //c1*1 + c1*1x - c1*2x^2 = 0
    //c2*2 + c2*5x - c2*1x^2 = 0
    //c3*0 + c3*1x + c3*1x^2 = 0
    //| 1   2   0   0|
    //| 1   5   1   0|
    //|-2  -1   1   0|
    mtx = .{ 1, 2, 0, 0, 1, 5, 1, 0, -2, -1, 1, 0 };
    ret = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    b = try rdcXmtx(&mtx, cols, hasAug, &ret, hasIdtMtx, &idtMtx, dim, triagRdcOnly, &sclr);
    try std.testing.expectEqual(true, b);
    std.debug.print("RdcMtx Result Example 9: {}\n", .{b});
    prntXmtxNl(&ret, 4);
    prntNl();

    //Example 10: pg. 193
    //Determine whether the following set of vectors in 2x2 matrices is linearly independent
    //or linearly dependent.
    //S = {
    //       |2,1 ,0,1|  v1
    //       |3,0 ,2,1|  v2
    //       |1,0 ,2,0|  v3
    //    }
    //c1*v1 + c2*v2 + c3*v3 = 0
    //
    //2c1 + 3c2 + 1c3 = 0
    //1c1 + 0c2 + 0c3 = 0
    //0c1 + 2c2 + 2c3 = 0
    //1c1 + 1c2 + 0c3 = 0
    //| 2   3   1   0|
    //| 1   0   0   0|
    //| 0   2   2   0|
    //| 1   1   0   0|

    var mtxB: [16]f32 = .{ 2, 3, 1, 0, 1, 0, 0, 0, 0, 2, 2, 0, 1, 1, 0, 0 };
    hasAug = false;
    var retB: [16]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    b = try rdcXmtx(&mtxB, cols, hasAug, &retB, hasIdtMtx, &idtMtx, dim, triagRdcOnly, &sclr);
    try std.testing.expectEqual(false, b);
    std.debug.print("RdcMtx Result Example 10: {}\n", .{b});
    prntXmtxNl(&ret, 4);
    prntNl();
}

test "XMTX: ELA - Larson, Edwards: 4.7 Example 2, 3, 4, 5 test" {
    prntNl();

    //Example 2: pg 225
    //Finding a coordinate matrix relative to a standard basis.
    //B = {(1,0), (1,2)}
    //Bp = {(1,0), (0,1)}
    //v = [3 2]

    var B: [4]f32 = .{ 1, 1, 0, 2 };
    var vec: [2]f32 = .{ 3, 2 };
    var exp: [2]f32 = .{ 5, 4 };
    var nvec: [2]f32 = .{ 0, 0 };
    var cols: usize = 2;
    var Bp: [4]f32 = .{ 1, 0, 0, 1 };
    var colsp: usize = 2;
    var idtMtx: [4]f32 = .{ 1, 0, 0, 1 };
    var b: bool = false;
    var vbose: bool = true;
    var ret: [4]f32 = .{ 0, 0, 0, 0 };

    b = try chgXvecBasis(&vec, &B, cols, false, &Bp, colsp, true, &idtMtx, &ret, &nvec, vbose);
    try std.testing.expectEqual(true, b);
    std.debug.print("4.7 Example 2 ret {}\n", .{b});
    prntXmtxNl(&ret, cols);
    prntNl();
    std.debug.print("4.7 Example 2 B\n", .{});
    prntXmtxNl(&B, cols);
    prntNl();
    std.debug.print("4.7 Example 2 nvec\n", .{});
    prntXmtxNl(&nvec, 1);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&exp, &nvec));

    //Example 3: pg 226
    //Finding a coordinate matrix relative to a nonstandard basis
    //vec = [ 1  2 -1]
    //Bp = | 1  0  2|
    //     | 0 -1  3|
    //     | 1  2 -5|

    //B = | 1  0  0|
    //    | 0  1  0|
    //    | 0  0  1|

    //Bp -> B = [B | Bp] => [I | Bp]
    //Transformation matrix is Bp

    var B3: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var vec3: [3]f32 = .{ 1, 2, -1 };
    var exp3: [3]f32 = .{ 5, -8, -2 };
    var nvec3: [3]f32 = .{ 0, 0, 0 };
    var cols3: usize = 3;
    var Bp3: [9]f32 = .{ 1, 0, 2, 0, -1, 3, 1, 2, -5 };
    var colsp3: usize = 3;
    var idtMtx3: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    b = false;
    vbose = false;
    var ret3: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    b = try chgXvecBasis(&vec3, &B3, cols3, true, &Bp3, colsp3, false, &idtMtx3, &ret3, &nvec3, vbose);
    try std.testing.expectEqual(true, b);
    std.debug.print("4.7 Example 3 ret3 {}\n", .{b});
    prntXmtxNl(&ret3, cols3);
    prntNl();
    std.debug.print("4.7 Example 3 B3\n", .{});
    prntXmtxNl(&B3, cols3);
    prntNl();
    std.debug.print("4.7 Example 3 nvec3\n", .{});
    prntXmtxNl(&nvec3, 1);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&exp3, &nvec3));

    //Example 4: pg 231
    //Finding a transition matrix
    //Find the transition matrix from B to Bp for the following bases in R^3.
    //B = {(1,0,0), (0,1,0), (0,0,1)}
    //Bp = {(1,0,1), (0,-1,2), (2,3,-5)}
    //vec = [1, 1, 1]

    var expB3: [9]f32 = .{ -1, 4, 2, 3, -7, -3, 1, -2, -1 };
    B3 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    vec3 = .{ 1, 1, 1 };
    exp3 = .{ 0, 0, 0 };
    nvec3 = .{ 0, 0, 0 };
    cols3 = 3;
    Bp3 = .{ 1, 0, 2, 0, -1, 3, 1, 2, -5 };
    colsp3 = 3;
    idtMtx3 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    b = false;
    vbose = false;
    ret3 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    b = try getBasisCnvXmtx(&B3, cols3, &Bp3, colsp3, &ret3, vbose);
    try std.testing.expectEqual(true, b);
    std.debug.print("4.7 Example 4 ret3 {}\n", .{b});
    prntXmtxNl(&ret3, cols3);
    prntNl();
    std.debug.print("4.7 Example 4 B3\n", .{});
    prntXmtxNl(&B3, cols3);
    prntNl();
    std.debug.print("4.7 Example 4 nvec\n", .{});
    prntXmtxNl(&nvec3, 1);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&B3, &expB3));

    //Example 5: pg 232
    //Finding a transition matrix
    //Find the transition matrix from B to Bp for the following bases for R2.
    //B = {(-3,2),(4,-2)}
    //Bp = {(-1,2),(2,-2)}

    var expB: [4]f32 = .{ -1, 2, -2, 3 };
    B = .{ -3, 4, 2, -2 };
    vec = .{ 1, 1 };
    exp = .{ 0, 0 };
    nvec = .{ 0, 0 };
    cols = 2;
    Bp = .{ -1, 2, 2, -2 };
    colsp = 2;
    idtMtx = .{ 1, 0, 0, 1 };
    b = false;
    vbose = true;
    ret = .{ 0, 0, 0, 0 };

    b = try getBasisCnvXmtx(&B, cols, &Bp, colsp, &ret, vbose);
    try std.testing.expectEqual(true, b);
    std.debug.print("4.7 Example 5 ret5 {}\n", .{b});
    prntXmtxNl(&ret, cols);
    prntNl();
    std.debug.print("4.7 Example 5 B\n", .{});
    prntXmtxNl(&B, cols);
    prntNl();
    std.debug.print("4.7 Example 5 nvec\n", .{});
    prntXmtxNl(&nvec, 1);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&B, &expB));
}

test "XMTX: ELA - Larson, Edwards: 4.7 Problem 1, 3, 5 test" {
    prntNl();

    var b: bool = false;
    const vbose: bool = false;

    //Problem 1: pg 234
    //B = {(2, -1), (0, 1)}
    //vec = [4 1]
    var B2: [4]f32 = .{ 2, 0, -1, 1 };
    const cols2: usize = 2;
    var Bp2: [4]f32 = .{ 1, 0, 0, 1 };
    var vec2: [2]f32 = .{ 4, 1 };
    var nvec2: [2]f32 = .{ 0, 0 };
    var ret2: [4]f32 = .{ 0, 0, 0, 0 };
    var exp2: [2]f32 = .{ 8, -3 };
    var idtMtx2: [4]f32 = .{ 1, 0, 0, 1 };

    b = try chgXvecBasis(&vec2, &B2, cols2, false, &Bp2, cols2, true, &idtMtx2, &ret2, &nvec2, vbose);
    try std.testing.expectEqual(true, b);
    std.debug.print("4.7 Problem 1 ret2 {}\n", .{b});
    prntXmtxNl(&ret2, cols2);
    prntNl();
    std.debug.print("4.7 Problem 1 B2\n", .{});
    prntXmtxNl(&B2, cols2);
    prntNl();
    std.debug.print("4.7 Problem 1 nvec2\n", .{});
    prntXmtxNl(&nvec2, 1);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&exp2, &nvec2));

    //Problem 3: pg 234
    //B = {(1,0,1), (1,1,0), (0,1,1)}
    //vec = [2 3 1]
    var B3: [9]f32 = .{ 1, 1, 0, 0, 1, 1, 1, 0, 1 };
    const cols3: usize = 3;
    var Bp3: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    var vec3: [3]f32 = .{ 2, 3, 1 };
    var nvec3: [3]f32 = .{ 0, 0, 0 };
    var ret3: [9]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var exp3: [3]f32 = .{ 5, 4, 3 };
    var idtMtx3: [9]f32 = .{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    b = try chgXvecBasis(&vec3, &B3, cols3, false, &Bp3, cols3, true, &idtMtx3, &ret3, &nvec3, vbose);
    try std.testing.expectEqual(true, b);
    std.debug.print("4.7 Problem 3 ret3 {}\n", .{b});
    prntXmtxNl(&ret3, cols3);
    prntNl();
    std.debug.print("4.7 Problem 3 B3\n", .{});
    prntXmtxNl(&B3, cols3);
    prntNl();
    std.debug.print("4.7 Problem 3 nvec3\n", .{});
    prntXmtxNl(&nvec3, 1);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&exp3, &nvec3));

    //Problem 5: pg 234
    //B = {(0,0,0,1), (0,0,1,1), (0,1,1,1), (1,1,1,1)}
    //vec = [1 -2 3 -1]
    var B4: [16]f32 = .{ 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1 };
    const cols4: usize = 4;
    var Bp4: [16]f32 = .{ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };
    var vec4: [4]f32 = .{ 1, -2, 3, -1 };
    var nvec4: [4]f32 = .{ 0, 0, 0, 0 };
    var ret4: [16]f32 = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var exp4: [4]f32 = .{ -1, 2, 0, 1 };
    var idtMtx4: [16]f32 = .{ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };

    b = try chgXvecBasis(&vec4, &B4, cols4, false, &Bp4, cols4, true, &idtMtx4, &ret4, &nvec4, vbose);
    try std.testing.expectEqual(true, b);
    std.debug.print("4.7 Problem 5 ret4 {}\n", .{b});
    prntXmtxNl(&ret4, cols4);
    prntNl();
    std.debug.print("4.7 Problem 5 B4\n", .{});
    prntXmtxNl(&B4, cols4);
    prntNl();
    std.debug.print("4.7 Problem 5 nvec4\n", .{});
    prntXmtxNl(&nvec4, 1);
    prntNl();
    try std.testing.expectEqual(true, equXmtx(&exp4, &nvec4));
}

test "XMTX: ELA - Larson, Edwards: 5.1 Example 1, 2, 3, 4 test" {
    prntNl();

    //Example 1.a: pg 252
    //The length of a vector in Rn
    //Find the length of v = [0, -2, 1, 4, -2]
    var vec: [5]f32 = .{ 0, -2, 1, 4, -2 };
    var mag: f32 = magXvec(&vec);
    var exp: f32 = 5.0;
    try std.testing.expectEqual(true, isEquF32(mag, exp, true));

    //Example 2: pg 253
    //Find the unit vector in the direction of v = (3, -1, 2), and verify that this vector has length 1.
    var vec3: [3]f32 = .{ 3, -1, 2 };
    mag = magXvec(&vec3);
    std.debug.print("5.1 example 2 mag 1: {}\n", .{mag});
    //divXvec(&vec3, mag);
    mulXvec(&vec3, 1.0 / mag);
    mag = magXvec(&vec3);
    std.debug.print("5.1 example 2 mag 2: {}\n", .{mag});
    exp = 1.0;
    try std.testing.expectEqual(true, isEquF32(mag, exp, false));

    //Example 3: pg 254
    //Finding the distance between two vectors
    //The distance between u = (0,2,2) and v = (2,0,1)
    exp = 3;
    var vec3a: [3]f32 = .{ 0, 2, 2 };
    var vec3b: [3]f32 = .{ 2, 0, 1 };
    diff1Xvec(&vec3a, &vec3b);
    mag = magXvec(&vec3a);
    try std.testing.expectEqual(true, isEquF32(mag, exp, true));

    //Example 4: pg 255
    //Finding the dot product of two vectors
    //The dot product of u = (1,2,0,-3) and v = (3, -2, 4, 2)
    var v1: [4]f32 = .{ 1, 2, 0, -3 };
    var v2: [4]f32 = .{ 3, -2, 4, 2 };
    const dprd: f32 = dotPrdXvec(&v1, &v2);
    exp = -7;
    try std.testing.expectEqual(true, isEquF32(dprd, exp, true));
}

test "XMTX: ELA - Larson, Edwards: 5.1 Example Theorem 5.4 test" {
    prntNl();

    //The Cauchy-Shwarz Inequality
    //If u and v are vectors in Rn, then
    // |u * v| <= ||u|| * ||v||
    //where |u * v|, denotes the absolute value of the dot product of u and v
    var u: [3]f32 = .{ 1, -1, 3 };
    var v: [3]f32 = .{ 2, 0, -1 };
    const mag1: f32 = magXvec(&u);
    const mag2: f32 = magXvec(&v);
    const dprd: f32 = dotPrdXvec(&u, &v);
    var b: bool = false;
    if (absF32(dprd) > mag1 * mag2) {
        b = true;
    } else {
        b = false;
    }
    try std.testing.expectEqual(false, b);
}

test "XMTX: ELA - Larson, Edwards: 5.1 Example Theorem 5.5 test" {
    prntNl();

    //The Triangle Inequality
    //If u and v are vectors in Rn, then ||u + v|| <= ||u|| + ||v||
    var u: [3]f32 = .{ 1, -1, 3 };
    var v: [3]f32 = .{ 2, 0, -1 };
    var uPv: [3]f32 = .{ 0, 0, 0 };
    const mag1: f32 = magXvec(&u);
    const mag2: f32 = magXvec(&v);
    const mag3: f32 = magXvec(&uPv);
    var b: bool = false;
    if (mag3 > mag1 + mag2) {
        b = true;
    } else {
        b = false;
    }
    try std.testing.expectEqual(false, b);
}

test "XMTX: ELA - Larson, Edwards: 5.1 Problem 1, 3, 5 test" {
    //Chapter 5: Section 5.1: Problem 1,3,5: pg 275
    prntNl();

    var mag: f32 = 0;
    var exp: f32 = 0;

    //Find the length of the following vectors:

    //Problem 1: pg 262
    //v = [4, 3]
    var v2: [2]f32 = .{ 4, 3 };
    mag = magXvec(&v2);
    exp = 5.0;
    std.debug.print("5.1 Problem 1 mag: {}\n", .{mag});
    try std.testing.expectEqual(true, isEquF32(mag, exp, true));

    //Problem 3: pg 262
    //v = [1, 2, 2]
    var v3: [3]f32 = .{ 1, 2, 2 };
    mag = magXvec(&v3);
    exp = 3.0;
    std.debug.print("5.1 Problem 3 mag: {}\n", .{mag});
    try std.testing.expectEqual(true, isEquF32(mag, exp, true));

    //Problem 5: pg 262
    //v = [4, 0, -3, 5]
    var v4: [4]f32 = .{ 4, 0, -3, 5 };
    mag = magXvec(&v4);
    exp = 7.0710;
    std.debug.print("5.1 Problem 5 mag: {}\n", .{mag});
    try std.testing.expectEqual(true, isEquF32(mag, exp, false));
}

test "XMTX: ELA - Larson, Edwards: 5.1 Problem 7, 9 test" {
    //Chapter 5: Section 5.1: Problem 7, 9: pg 262
    //Find: (a) ||u||, (b) ||v||, (c) ||u + v||

    //7: u=(0, 4, 3), v=(1,-2,1)
    //Exp: (a) 5, (b) sqrt(6.0), (c) sqrt(21.0)
    const u: [3]f32 = .{ 0, 4, 3 };
    const v: [3]f32 = .{ 1, -2, 1 };
    const t: [3]f32 = .{ 0, 0, 0 };
    var val: f32 = 0;
    var exp: f32 = 0;

    //(a)
    val = magXvec(@constCast(&u));
    exp = 5;
    prntNlStrArgs("7a Found val: {any}", .{val});
    prntNlStrArgs("7a Found exp: {any}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();

    //(b)
    val = magXvec(@constCast(&v));
    exp = std.math.sqrt(6.0);
    prntNlStrArgs("7b Found val: {any}", .{val});
    prntNlStrArgs("7b Found exp: {any}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();

    //(c)
    sum2Xvec(@constCast(&t), @constCast(&u), @constCast(&v));
    val = magXvec(@constCast(&t));
    exp = std.math.sqrt(21.0);
    prntNlStrArgs("7c Found val: {any}", .{val});
    prntNlStrArgs("7c Found exp: {any}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();

    //9: u=(0, 4, 3), v=(1,-2,1)
    //Exp: (a) sqrt(6.0), (b) sqrt(11.0), (c) sqrt(13.0)
    const lu: [4]f32 = .{ 0, 1, -1, 2 };
    const lv: [4]f32 = .{ 1, 1, 3, 0 };
    const lt: [4]f32 = .{ 0, 0, 0, 0 };

    //(a)
    val = magXvec(@constCast(&lu));
    exp = std.math.sqrt(6.0);
    prntNlStrArgs("9a Found val: {any}", .{val});
    prntNlStrArgs("9a Found exp: {any}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();

    //(b)
    val = magXvec(@constCast(&lv));
    exp = std.math.sqrt(11.0);
    prntNlStrArgs("9b Found val: {any}", .{val});
    prntNlStrArgs("9b Found exp: {any}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();

    //(c)
    sum2Xvec(@constCast(&lt), @constCast(&lu), @constCast(&lv));
    val = magXvec(@constCast(&lt));
    exp = std.math.sqrt(13.0);
    prntNlStrArgs("9c Found val: {any}", .{val});
    prntNlStrArgs("9c Found exp: {any}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.1 Problem 11, 13 test" {
    //Chapter 5: Section 5.1: Problem 11, 13: pg 262
    //Find: (a) a unit vector in the direction of u
    //      (b) a unit vector in the direction opposite of u

    //11: u=(-5,12)
    //Exp: (a) ((-5.0 / 13.0), (12.0 / 13.0))
    //     (b) ((5.0 / 13.0), (-12.0 / 13.0))
    //(a)
    var u: [2]f32 = .{ -5.0, 12.0 };
    var exp: [2]f32 = .{ (-5.0 / 13.0), (12.0 / 13.0) };
    nrmXvec(&u);
    var val: [2]f32 = u;
    prntNlStrArgs("11a Found val: {any}", .{val});
    prntNlStrArgs("11a Found exp: {any}", .{exp});
    try std.testing.expectEqual(true, equXvecWrkr(&val, &exp, false));
    prntNl();

    //(b)
    u = .{ -5.0, 12.0 };
    var t: [2]f32 = .{ 0, 0 };
    diff1Xvec(&t, &u);
    exp = .{ (5.0 / 13.0), (-12.0 / 13.0) };
    nrmXvec(&t);
    val = t;
    prntNlStrArgs("11b Found val: {any}", .{val});
    prntNlStrArgs("11b Found exp: {any}", .{exp});
    try std.testing.expectEqual(true, equXvecWrkr(&val, &exp, false));
    prntNl();

    //13: u=(3,2,-5)
    //Exp: (a) ((3.0 / sqrt(38.0)), (2.0 / sqrt(38.0)), (-5.0 / sqrt(38.0)))
    //     (b) ((-3.0 / sqrt(38.0)), (-2.0 / sqrt(38.0)), (5.0 / sqrt(38.0)))
    var lu: [3]f32 = .{ 3, 2, -5 };
    var lexp: [3]f32 = .{ (3.0 / std.math.sqrt(38.0)), (2.0 / std.math.sqrt(38.0)), (-5.0 / std.math.sqrt(38.0)) };
    nrmXvec(&lu);
    var lval: [3]f32 = lu;
    prntNlStrArgs("13a Found val: {any}", .{lval});
    prntNlStrArgs("13a Found exp: {any}", .{lexp});
    try std.testing.expectEqual(true, equXvecWrkr(&lval, &lexp, false));
    prntNl();

    //(b)
    lu = .{ 3, 2, -5 };
    var lt: [3]f32 = .{ 0, 0, 0 };
    diff1Xvec(&lt, &lu);
    lexp = .{ (-3.0 / std.math.sqrt(38.0)), (-2.0 / std.math.sqrt(38.0)), (5.0 / std.math.sqrt(38.0)) };
    nrmXvec(&lt);
    lval = lt;
    prntNlStrArgs("13b Found val: {any}", .{lval});
    prntNlStrArgs("13b Found exp: {any}", .{lexp});
    try std.testing.expectEqual(true, equXvecWrkr(&lval, &lexp, false));
    prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.1 Problem 15 test" {
    //Chapter 5: Section 5.1: Problem 15: pg 262
    //For what values of c is ||c(1,2,3)|| = 1
    //15:
    const u = [3]f32{ 1, 2, 3 };
    const su: []f32 = @constCast(&u);
    const val = 1.0 / magXvec(@constCast(su));
    const exp = (1.0 / std.math.sqrt(14.0));

    prntNlStrArgs("15 Found val: {any}", .{val});
    prntNlStrArgs("15 Found exp: {any}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.1 Problem 17, 19 test" {
    //Chapter 5: Section 5.1: Problem 17, 19: pg 262

    //17:
    const u: [2]f32 = .{ 2.0 * std.math.sqrt(2.0), 2.0 * std.math.sqrt(2.0) };
    var val: f32 = magXvec(@constCast(&u));
    var exp: f32 = 4.0;

    prntNlStrArgs("17 Found val: {any}", .{val});
    prntNlStrArgs("17 Found exp: {any}", .{exp});
    try std.testing.expectEqual(true, isEquF32(val, exp, false));
    prntNl();

    //19:
    const v: [3]f32 = .{ 1, std.math.sqrt(3.0), 0 };
    val = magXvec(@constCast(&v));
    exp = 2.0;

    prntNlStrArgs("19 Found val: {any}", .{val});
    prntNlStrArgs("19 Found exp: {any}", .{exp});
    try std.testing.expectEqual(true, isEquF32(val, exp, false));
    prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.1 Problem 21 test" {
    //Chapter 5: Section 5.1: Problem 21: pg 262
    //Given (8,8,6) find v = same direction, 1/2 len, v = opposite direction, 1/4 len
    //21:
    //(a)
    const u: [3]f32 = .{ 8, 8, 6 };
    const l: f32 = magXvec(@constCast(&u));
    var val: f32 = l / 2.0;
    var v: [3]f32 = .{ 4, 4, 3 };
    var exp: f32 = magXvec(@constCast(&v));

    prntNlStrArgs("21a Found val: {any}", .{val});
    prntNlStrArgs("21a Found exp: {any}", .{exp});
    try std.testing.expectEqual(true, isEquF32(val, exp, false));
    prntNl();

    //(b)
    val = l / 4.0;
    v = .{ -2.0, -2.0, (-3.0 / 2.0) };
    exp = magXvec(@constCast(&v));

    prntNlStrArgs("21b Found val: {any}", .{val});
    prntNlStrArgs("21b Found exp: {any}", .{exp});
    try std.testing.expectEqual(true, isEquF32(val, exp, false));
    prntNl();
}

test "XMTX: ELA - Larson, Edwards: 5.1 Problem 23, 25 test" {
    //Chapter 5: Section 5.1: Problem 23, 25: pg 262
    //23:
    const u = [2]f32{ 1, -1 };
    const v = [2]f32{ -1, 1 };
    const alloc = std.testing.allocator;
    var val: f32 = try dstXvec(@constCast(&u), @constCast(&v), &alloc);
    var exp: f32 = (2.0 * std.math.sqrt(2.0));
    prntNlStrArgs("23 Found val: {any}", .{val});
    prntNlStrArgs("23 Found exp: {any}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();

    //25:
    const lu = [3]f32{ 1, 1, 2 };
    const lv = [3]f32{ -1, 3, 0 };
    val = try dstXvec(@constCast(&lu), @constCast(&lv), &alloc);
    exp = (2.0 * std.math.sqrt(3.0));
    prntNlStrArgs("25 Found val: {any}", .{val});
    prntNlStrArgs("25 Found exp: {any}", .{exp});
    try std.testing.expectEqual(val, exp);
    prntNl();
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//STOP PROBLEMS
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
