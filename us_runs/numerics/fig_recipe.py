"""Figure: what makes the transport runs fast (orthogonal TCV-X21 grid, 8 processes).

usage (from us_runs/): python numerics/fig_recipe.py <struct json> <out stem>
(the json is written by numerics/recipe_analyse.py for the runs in t23_recipe/struct)

Wall-clock time to steady state (within 1 % of the run's own final state) when ONE
ingredient of the recipe is changed. Each line starts at the recipe's time; runs
stopped after 1 h end in an arrow. Right-hand column: how much the change alters
the steady state (largest volume-weighted difference of n_i, p_i, p_e, n_n).
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paperstyle as ps  # noqa: E402

ps.setup()
ps.mpl.rcParams["text.latex.preamble"] += r"\usepackage{xcolor}"
plt = ps.plt
LIMIT_S = 3600.0

runs = {r["name"]: r for r in json.load(open(sys.argv[1]))}
if "bjacobi_ilu_clean" in runs:      # rerun on a quiet machine (the first shared it)
    runs["bjacobi_ilu"] = runs["bjacobi_ilu_clean"]
ref = runs["recipe"]
t_ref = ref["wall_1pct_s"]
# run-to-run spread of the wall time: identical runs (repeat of the recipe; a run whose
# change does not alter a single step) differ by this fraction
same = [runs[k]["wall_1pct_s"] for k in ("recipe_rep", "recipe_v18", "coulomb_orig") if k in runs
        and runs[k]["wall_1pct_s"] is not None]
SPREAD = max([abs(t / t_ref - 1) for t in same] + [0.05])

# (run name, label); a row with name None is a group header
ROWS = [(None, r"Jacobian"),
        ("jac_every_step", r"rebuilt every time step"),
        ("jac_every_newton", r"rebuilt every Newton iteration"),
        ("jac_persists", r"kept by PETSc across steps instead"),
        ("fd_operator", r"used as the operator (not matrix-free)"),
        (None, r"Linear solver"),
        ("bjacobi_ilu", r"block-Jacobi ILU(0) instead of LU"),
        ("ref_blas", r"reference BLAS instead of OpenBLAS"),
        (None, r"Time-step control"),
        ("threshold_dt", r"thresholds instead of PID control"),
        ("predictor_const", r"no extrapolated first guess"),
        ("no_ew", r"fixed linear tolerance (no Eisenstat--Walker)"),
        ("atol_1e-3", r"absolute tolerance $10^{-3}$ instead of $10^{-6}$"),
        (None, r"Slope limiter (van Albada in the recipe)"),
        ("lim_MinMod", r"MinMod"),
        ("lim_MC", r"MC (Hermes-3 default)"),
        ("lim_WENO3", r"WENO3"),
        ("lim_Superbee", r"Superbee"),
        ("lim_Upwind", r"first-order upwind"),
        (None, r"Model details"),
        ("coulomb_orig", r"Coulomb logarithm constant below 0.1\,eV"),
        ("no_limiter_lag", r"neutral flux limiter not lagged"),
        (None, r"All at once: Hermes-3's 2D tokamak example"),
        ("example_settings", r"its solver settings"),
        ("example_settings_mc", r"its solver settings and the MC limiter")]
NOTES = {"coulomb_orig": r"identical run: no cell cools below 0.1\,eV here"}
EXTRA = {"no_ew": lambda r: r"; within 0.1\,\%% only after %.0f\,s" % r["wall_01pct_s"]}
BOX = dict(boxstyle="square,pad=0.08", fc="white", ec="none")
present = [r for r in ROWS if r[0] is None or r[0] in runs]
ROWS = [r for i, r in enumerate(present)
        if r[0] is not None or (i + 1 < len(present) and present[i + 1][0] is not None)]

nrows = len(ROWS)
fig = plt.figure(figsize=(ps.WIDTH, 0.205 * nrows + 0.75))
ax = fig.add_axes([0.37, 0.75 / (0.205 * nrows + 0.75), 0.40, 1 - 0.95 / (0.205 * nrows + 0.75)])

y = 0.0
ys, labels = [], []
for name, label in ROWS:
    if name is None:
        y += 0.25
        ax.text(-0.015, y, r"\textbf{%s}" % label, transform=ax.get_yaxis_transform(), ha="right",
                va="center", fontsize=8, color=ps.INK)
        y += 1.0
        continue
    r = runs[name]
    ys.append(y)
    labels.append(label)
    if r["wall_1pct_s"] is not None:
        t = r["wall_1pct_s"]
        color = ps.MUTED if abs(t / t_ref - 1) <= SPREAD else (ps.ORANGE if t > t_ref else ps.BLUE)
        ax.plot([t_ref, t], [y, y], color=color, lw=1.5, solid_capstyle="butt", zorder=2)
        ax.plot([t], [y], marker="o", ms=5.5, color=color, mec="white", mew=0.8, zorder=3)
        fac = r"$\times%.1f$" % (t / t_ref) if t >= t_ref else r"$\times%.2f$" % (t / t_ref)
        if abs(t / t_ref - 1) <= SPREAD:
            fac = r"no change (within the run-to-run spread)"
        if name in EXTRA:
            fac += EXTRA[name](r)
        ax.text(max(t, t_ref) * 1.15, y, NOTES.get(name, fac), va="center", ha="left", fontsize=7.2,
                color=ps.INK if (name not in NOTES and abs(t / t_ref - 1) > SPREAD) else ps.SECONDARY,
                bbox=BOX, zorder=4)
    elif r.get("exit_code") == 124 or r.get("wall_run_s") is None:
        # stopped at the time limit
        ax.plot([t_ref, LIMIT_S], [y, y], color=ps.ORANGE, lw=1.5, solid_capstyle="butt", zorder=2)
        ax.annotate("", xy=(LIMIT_S * 1.6, y), xytext=(LIMIT_S * 0.98, y),
                    arrowprops=dict(arrowstyle="-|>", color=ps.ORANGE, lw=1.5, mutation_scale=8), zorder=3)
        ax.text(80, y - 0.05, r"no steady state in 1\,h (reached %.2g\,ms)" % r["t_last_ms"],
                va="center", ha="left", fontsize=7.2, color=ps.INK,
                bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none"))
    else:
        # the solver gave up
        t = r["wall_run_s"]
        ax.plot([t_ref, t], [y, y], color=ps.ORANGE, lw=1.5, solid_capstyle="butt", zorder=2)
        ax.plot([t], [y], marker="X", ms=6.5, color=ps.ORANGE, mec="white", mew=0.6, zorder=3)
        ax.text(t * 1.18, y, r"aborted after %.0f\,min at %.2g\,ms" % (t / 60, r["t_last_ms"]),
                va="center", ha="left", fontsize=7.2, color=ps.INK, bbox=BOX, zorder=4)
    # effect on the steady state, in its own column
    d = r.get("diff_ref_max")
    if r["complete"] and d is not None:
        txt = r"%.1f\,\%%" % (100 * d) if d >= 1e-3 else r"$<0.1$\,\%"
        ax.text(1.03, y, txt, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=7.2,
                color=ps.INK if d >= 1e-3 else ps.MUTED)
    y += 1.0
ax.text(1.03, -0.85, r"change of" "\n" r"steady state", transform=ax.get_yaxis_transform(),
        ha="left", va="bottom", fontsize=7.2, color=ps.SECONDARY, linespacing=1.0)
ax.axvline(t_ref, color=ps.SECONDARY, lw=0.8, zorder=1)
ax.axvspan(t_ref * (1 - SPREAD), t_ref * (1 + SPREAD), color=ps.RULE, alpha=0.45, lw=0, zorder=0)
ax.text(t_ref * 0.93, -0.85, r"recipe: %.0f\,s" % t_ref, ha="right", va="bottom", fontsize=7.2,
        color=ps.SECONDARY)
ax.set_yticks(ys)
ax.set_yticklabels(labels, fontsize=7.6)
ax.tick_params(axis="y", length=0, pad=4)
ax.set_ylim(y - 0.5, -1.0)
ax.set_xscale("log")
ax.set_xlim(12, LIMIT_S * 1.6)
TICKS = [(15, r"15\,s"), (30, r"30\,s"), (60, r"1\,min"), (120, r"2\,min"), (300, r"5\,min"),
         (900, r"15\,min"), (3600, r"1\,h")]
ax.set_xticks([t for t, _ in TICKS])
ax.set_xticklabels([l for _, l in TICKS])
ax.xaxis.set_minor_locator(ps.mpl.ticker.NullLocator())
ax.set_xlabel(r"wall-clock time to steady state (8 processes)")
ax.grid(axis="x", which="major", color=ps.RULE, lw=0.5, zorder=0)
for s in ("left", "right", "top"):
    ax.spines[s].set_visible(False)
ax.tick_params(axis="x", which="both", top=False)
ps.save(fig, sys.argv[2])
print("wrote", sys.argv[2])
