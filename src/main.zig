const std = @import("std");
const plot = @import("plotting.zig");

const POINTS: comptime_int = 2000;

const INIT_TIME: comptime_float = 0;
const MAX_TIME: comptime_float = 15;

const td = (MAX_TIME - INIT_TIME) / @as(comptime_float, POINTS);

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const alloc = gpa.allocator();
    defer std.debug.assert(gpa.deinit() == .ok);
    // var alloc = std.heap.c_allocator;

    // try projectile_motion_with_drag_improper(alloc);
    // try projectile_motion(alloc);
    try damped_harmonic_oscillator(alloc);
    // try oscillations(alloc);
}

// fn euler_free_fall(alloc: std.mem.Allocator) !void {
//     const g = 9.81;

//     var vs = [_]f64{0.0} ** POINTS;
//     var ts = [_]f64{0.0} ** POINTS;
//     ts[0] = INIT_TIME;

//     for (1..POINTS) |i| {
//         vs[i] = vs[i - 1] - td * g;
//         ts[i] = ts[i - 1] + td;
//     }

//     try createPlot(alloc, "./euler.png", &ts, &vs, "Euler's Method");
// }

fn projectile_motion(alloc: std.mem.Allocator) !void {
    const g: f64 = 9.8; // Gravitational acceleration
    const B_m: f64 = 4e-3; // Drag coefficient divided by mass

    // For storing points where particle goes below surface
    var wo_point_break: usize = undefined;
    var point_break: usize = undefined;

    // For storing x and y velocities. Initialised to zero
    // vys_w_drag and vxs_w_drag -> with drag
    // vys_wo_drg and vxs_wo_drag -> without drag
    var vys_w_drag = [_]f64{0.0} ** POINTS;
    vys_w_drag[0] = 20;
    var vxs_w_drag = [_]f64{0.0} ** POINTS;
    vxs_w_drag[0] = 5;

    var vys_wo_drag = [_]f64{0.0} ** POINTS;
    vys_wo_drag[0] = 20;
    var vxs_wo_drag = [_]f64{0.0} ** POINTS;
    vxs_wo_drag[0] = 5;

    // For stroing magnitude of velocities
    var v: f64 = undefined;

    // For storing x and y positions. Initialised to zero
    // ys_w_drag and xs_w_drag -> with drag
    // ys_wo_drg and xs_wo_drag -> without drag
    var xs_w_drag = [_]f64{0.0} ** POINTS;
    var ys_w_drag = [_]f64{0.0} ** POINTS;

    var xs_wo_drag = [_]f64{0.0} ** POINTS;
    var ys_wo_drag = [_]f64{0.0} ** POINTS;

    // Euler's method with drag
    for (1..POINTS) |i| {
        // Compute velocity magnitude
        v = std.math.sqrt(
            std.math.pow(
                f64,
                vys_w_drag[i - 1],
                2,
            ) + std.math.pow(
                f64,
                vxs_w_drag[i - 1],
                2,
            ),
        );

        vys_w_drag[i] = vys_w_drag[i - 1] - (B_m * vys_w_drag[i - 1] * v + g) * td;
        vxs_w_drag[i] = vxs_w_drag[i - 1] - B_m * vxs_w_drag[i - 1] * v * td;

        xs_w_drag[i] = xs_w_drag[i - 1] + vxs_w_drag[i - 1] * td;
        ys_w_drag[i] = ys_w_drag[i - 1] + vys_w_drag[i - 1] * td;

        if (ys_w_drag[i] < 0) {
            point_break = i;
            break;
        }
    }

    // Euler's method without drag
    for (1..POINTS) |i| {
        vys_wo_drag[i] = vys_wo_drag[i - 1] - g * td;
        vxs_wo_drag[i] = vxs_wo_drag[i - 1];

        xs_wo_drag[i] = xs_wo_drag[i - 1] + vxs_wo_drag[i - 1] * td;
        ys_wo_drag[i] = ys_wo_drag[i - 1] + vys_wo_drag[i - 1] * td;

        if (ys_wo_drag[i] < 0) {
            wo_point_break = i;
            break;
        }
    }

    // Plotting all data on single axis
    var plt = plot.Plot(f64).init(alloc);
    defer plt.deinit();

    plt.setXLabel("x");
    plt.setYLabel("y");
    plt.setTitle("Projectile Motion with Drag vs without Drag");
    plt.showGrid();

    try plt.addPlot(xs_wo_drag[0..wo_point_break], ys_wo_drag[0..wo_point_break], "Without Drag");
    try plt.addPlot(xs_w_drag[0..point_break], ys_w_drag[0..point_break], "With Drag");

    try plt.saveFig("./proj.png", .PNG);
}

fn damped_harmonic_oscillator(alloc: std.mem.Allocator) !void {
    var ts = [_]f64{0.0} ** POINTS;
    var omegas = [_]f64{0.0} ** POINTS;
    var thetas = [_]f64{0.0} ** POINTS;
    thetas[0] = std.math.degreesToRadians(f64, 20.0);

    const a2: comptime_float = 0.5;
    const a1: comptime_float = 1 - a2;
    const q11: comptime_float = 1 / (2 * a2);
    const p1: comptime_float = 1 / (2 * a2);

    var k1: f64 = 0.0;
    var k2: f64 = 0.0;
    var l1: f64 = 0.0;
    var l2: f64 = 0.0;
    for (1..POINTS) |i| {
        k1 = td * f(thetas[i - 1], omegas[i - 1], ts[i - 1] + p1 * td);
        l1 = td * G(thetas[i - 1], omegas[i - 1], ts[i - 1] + p1 * td);

        k2 = td * f(thetas[i - 1] + q11 * k1, omegas[i - 1] + q11 * l1, ts[i - 1] + p1 * td);
        l2 = td * G(thetas[i - 1] + q11 * k1, omegas[i - 1] + q11 * l1, ts[i - 1] + p1 * td);

        thetas[i] = thetas[i - 1] + (a1 * k1 + a2 * k2);
        omegas[i] = omegas[i - 1] + (a1 * l1 + a2 * l2);
        ts[i] = ts[i - 1] + td;
    }

    var plt = plot.Plot(f64).init(alloc);
    defer plt.deinit();

    try plt.addPlot(&ts, &omegas, null);
    plt.setXLabel("Time");
    plt.setYLabel("Thetas");
    plt.showGrid();

    try plt.saveFig("pendulum_with_drag.png", .PNG);
}

fn f(theta: f64, omega: f64, t: f64) f64 {
    _ = t;
    _ = theta;

    return omega;
}

pub fn G(theta: f64, omega: f64, t: f64) f64 {
    const Fd: comptime_float = 0.2;
    const q: comptime_float = 1.0;
    const Omega_d: comptime_float = 10.0;
    const l: comptime_float = 1.0;
    const g: comptime_float = 9.8;

    return (-(g / l) * theta) - (q * omega) + (Fd * std.math.sin(Omega_d * t));
}
