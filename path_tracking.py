"""
Path tracking for a car-like mobile robot: Pure Pursuit vs. Stanley vs. LQR.

Kinematic bicycle model (0.5 m wheelbase, +/-30 deg steering, 0.1 s steering lag) follows a
smooth 2-D path built from a cubic spline. Each controller is scored on cross-track error,
heading error and steering effort, at two speeds and with noisy localization.

Run:  python path_tracking.py           (prints the comparison, saves tracking.png)
"""
import argparse
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import solve_discrete_are

L, MAX_STEER, TAU, DT = 0.5, np.deg2rad(30), 0.1, 0.02

def wrap(a): return (a + np.pi) % (2 * np.pi) - np.pi

def make_path():
    wp = np.array([[0, 0], [6, 0], [10, 3], [12, 8], [9, 12], [4, 12], [1, 16], [5, 20], [12, 20], [16, 16]], float)
    t = np.r_[0, np.cumsum(np.linalg.norm(np.diff(wp, axis=0), axis=1))]
    sx, sy = CubicSpline(t, wp[:, 0]), CubicSpline(t, wp[:, 1])
    tt = np.linspace(0, t[-1], 4000)
    x, y = sx(tt), sy(tt)
    dx, dy, ddx, ddy = sx(tt, 1), sy(tt, 1), sx(tt, 2), sy(tt, 2)
    yaw = np.arctan2(dy, dx)
    kappa = (dx * ddy - dy * ddx) / np.power(dx * dx + dy * dy, 1.5)
    s = np.r_[0, np.cumsum(np.hypot(np.diff(x), np.diff(y)))]
    return dict(x=x, y=y, yaw=yaw, k=kappa, s=s)

def nearest(path, px, py, last):
    lo, hi = max(0, last - 50), min(len(path["x"]), last + 300)
    d = np.hypot(path["x"][lo:hi] - px, path["y"][lo:hi] - py)
    i = lo + int(np.argmin(d))
    # signed cross-track error: positive when the robot is left of the path
    dxp, dyp = px - path["x"][i], py - path["y"][i]
    e = -np.sin(path["yaw"][i]) * dxp + np.cos(path["yaw"][i]) * dyp
    return i, e

# ---------------------------------------------------------------- controllers
def pure_pursuit(st, path, i, e, v, p):
    ld = p.get("k_ld", 0.6) * v + p.get("ld0", 0.4)
    j = i
    while j < len(path["s"]) - 1 and path["s"][j] - path["s"][i] < ld: j += 1
    rx, ry = path["x"][j] - st[0], path["y"][j] - st[1]          # rear-axle reference point
    alpha = wrap(np.arctan2(ry, rx) - st[2])
    return np.arctan2(2 * L * np.sin(alpha), ld)

def stanley(st, path, i, e, v, p):
    fx, fy = st[0] + L * np.cos(st[2]), st[1] + L * np.sin(st[2])  # front axle
    jf, ef = nearest(path, fx, fy, i)
    return wrap(path["yaw"][jf] - st[2]) + np.arctan2(-p.get("k", 2.0) * ef, v + 0.3)

_lqr_cache = {}
def lqr(st, path, i, e, v, p):
    key = round(v, 2)
    if key not in _lqr_cache:                                     # lateral error model: [e, yaw error]
        A = np.array([[1, v * DT], [0, 1]]); B = np.array([[0], [v * DT / L]])
        Q, R = np.diag([4.0, 1.0]), np.array([[2.0]])
        P = solve_discrete_are(A, B, Q, R)
        _lqr_cache[key] = np.linalg.solve(R + B.T @ P @ B, B.T @ P @ A)
    K = _lqr_cache[key]
    ey = wrap(st[2] - path["yaw"][i])
    ff = np.arctan(L * path["k"][i])                              # curvature feed-forward
    return ff - float((K @ np.array([e, ey]))[0])

# ---------------------------------------------------------------- simulation
def run(ctrl, path, v, noise, rng, p=None):
    p = p or {}
    st = np.array([path["x"][0], path["y"][0] + 0.5, path["yaw"][0] + 0.1])   # start 0.5 m off the path
    steer, i = 0.0, 0
    cte, hdg, dsteer, xs, ys = [], [], [], [], []
    for _ in range(int(3 * path["s"][-1] / v / DT)):
        meas = st + (rng.normal(0, [0.03, 0.03, np.deg2rad(1)]) if noise else 0)
        i, e = nearest(path, meas[0], meas[1], i)
        cmd = np.clip(ctrl(meas, path, i, e, v, p), -MAX_STEER, MAX_STEER)
        new = steer + (cmd - steer) * DT / TAU                    # actuator lag
        dsteer.append(abs(new - steer) / DT); steer = new
        st = st + DT * np.array([v * np.cos(st[2]), v * np.sin(st[2]), v / L * np.tan(steer)])
        st[2] = wrap(st[2])
        ti, te = nearest(path, st[0], st[1], i)                   # score against the TRUE pose
        cte.append(te); hdg.append(wrap(st[2] - path["yaw"][ti])); xs.append(st[0]); ys.append(st[1])
        if ti >= len(path["x"]) - 5: break
    cte = np.array(cte); skip = int(2.0 / DT)                     # ignore the first 2 s (initial offset)
    far = np.nonzero(np.abs(cte) > 0.05)[0]                      # settling: last time error was above 5 cm
    gaps = np.nonzero(np.diff(far) > 1)[0]
    settle = (far[gaps[0]] if len(gaps) else far[-1]) * DT if len(far) else 0.0
    return dict(settle=settle, rms=np.sqrt(np.mean(cte[skip:] ** 2)), max=np.abs(cte[skip:]).max(),
                hdg=np.rad2deg(np.sqrt(np.mean(np.square(hdg[skip:])))), effort=np.rad2deg(np.mean(dsteer)),
                done=ti >= len(path["x"]) - 5, x=xs, y=ys, cte=cte)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--seed", type=int, default=0); a = ap.parse_args()
    path = make_path()
    ctrls = {"Pure Pursuit": pure_pursuit, "Stanley": stanley, "LQR + feed-forward": lqr}
    print(f"Path length {path['s'][-1]:.1f} m, max curvature {np.abs(path['k']).max():.2f} 1/m, "
          f"wheelbase {L} m, steering limit 30 deg, 0.1 s steering lag\n")
    print(f"{'scenario':<20}{'controller':<20}{'settle':>8}{'RMS CTE':>9}{'max CTE':>9}{'RMS heading':>13}{'steer rate':>12}")
    results = {}
    for v, noise in ((1.0, False), (2.0, False), (3.0, False), (2.0, True)):
        tag = f"{v:.0f} m/s" + (", noisy pose" if noise else "")
        for name, f in ctrls.items():
            r = run(f, path, v, noise, np.random.default_rng(a.seed))
            results[(tag, name)] = r
            flag = "" if r["done"] else "  (did not finish)"
            print(f"{tag:<20}{name:<20}{r['settle']:6.2f} s{r['rms'] * 100:7.1f} cm{r['max'] * 100:7.1f} cm{r['hdg']:10.2f} deg"
                  f"{r['effort']:8.1f} deg/s{flag}")
        print()
    best = min(ctrls, key=lambda n: results[("2 m/s", n)]["rms"])
    print(f"Lowest RMS cross-track error at 2 m/s: {best}")
    plot(path, results)

def plot(path, results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.5))
    a1.plot(path["x"], path["y"], "k-", lw=3, alpha=.3, label="reference path")
    cols = {"Pure Pursuit": "tab:red", "Stanley": "tab:blue", "LQR + feed-forward": "tab:green"}
    for (tag, name), r in results.items():
        if tag == "2 m/s":
            a1.plot(r["x"], r["y"], color=cols[name], lw=1.3, label=name)
            a2.plot(np.arange(len(r["cte"])) * DT, r["cte"] * 100, color=cols[name], lw=1.2, label=name)
    a1.set_aspect("equal"); a1.legend(fontsize=8); a1.set_title("Tracking at 2 m/s (starts 0.5 m off the path)")
    a2.set_xlabel("time (s)"); a2.set_ylabel("cross-track error (cm)"); a2.set_ylim(-25, 55)
    a2.axhline(0, color="k", lw=.5); a2.legend(fontsize=8); a2.set_title("Cross-track error at 2 m/s")
    fig.tight_layout(); fig.savefig("tracking.png", dpi=120)
    print("saved tracking.png")

if __name__ == "__main__":
    main()
