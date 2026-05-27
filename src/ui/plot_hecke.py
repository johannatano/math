import sys
import warnings
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
warnings.filterwarnings("ignore", category=UserWarning, module="bokeh")

from sympy import primerange
from bokeh.models import (ColumnDataSource, RangeTool, Spinner, Button, Range1d,
                          HoverTool, TapTool, DataTable, TableColumn, NumberFormatter, StringFormatter, Div, TextInput)
from bokeh.layouts import column, row
from bokeh.plotting import figure, curdoc
from nt.modular_forms import CuspForm, HeckeOperator
from nt.common import fmt_invariants
from tqdm import tqdm

STEP = 1000
state = {"pmin": 2, "pmax": 5000}

curdoc().theme = "dark_minimal"

N_input   = Spinner(title="Level N", value=500, low=1,   high=10000, step=1)
k_input   = Spinner(title="Weight k", value=2,  low=0,   high=10000, step=1)
btn_left  = Button(label="◀  -1000", width=110, button_type="default")
btn_right = Button(label="+1000  ▶", width=110, button_type="default")
dk_input  = TextInput(title="Filter D_K (blank = all)", value="", width=180)

source        = ColumnDataSource(dict(q=[], val=[], ref=[], val_raw=[], ref_raw=[]))
source_curves = ColumnDataSource(dict(q=[], nc=[], nss=[], np=[]))
source_traces = ColumnDataSource(dict(q=[], t=[], color=[]))
source_fpi    = ColumnDataSource(dict(q=[], fpi=[]))
source_detail = ColumnDataSource(dict(prop=[], num=[]))
source_torsion = ColumnDataSource(dict(label=[], f=[], h=[], hw=[], n_pts=[], value=[], inv=[], torsion_inv=[]))
source_dk_hist = ColumnDataSource(dict(dk=[], count=[], left=[], right=[]))

# global store for D_K per trace — updated in recompute, used to refresh histogram on pan
_all_trace_qs: list = []
_all_trace_dks: list = []

detail_range   = Range1d(start=state["pmin"], end=min(state["pmin"] + 500, state["pmax"]))
overview_range = Range1d(start=state["pmin"], end=state["pmax"])

p_detail = figure(
    height=350, width=900,
    x_axis_label="q (prime)", y_axis_label="a_q / q^((weight-1)/2)",
    tools="xpan,xwheel_zoom,reset",
    x_range=detail_range,
)
p_detail.scatter("q", "val", source=source, legend_label="Computed", color="steelblue", size=6)
# p_detail.scatter("q", "ref", source=source, legend_label="Sage ref", color="orange", size=4, marker="triangle")
p_detail.legend.location = "top_left"
p_detail.add_tools(HoverTool(tooltips=[
    ("q",            "@q{0}"),
    ("computed",     "@val{0.0000}  (@val_raw{0})"),
    ("sage ref",     "@ref{0.0000}  (@ref_raw{0})"),
]))
p_detail.add_tools(TapTool())

p_curves = figure(
    height=200, width=900,
    x_axis_label="q (prime)", y_axis_label="# curves",
    tools="xpan,xwheel_zoom,reset",
    x_range=detail_range,
)
p_curves.scatter("q", "nc",  source=source_curves, legend_label="#E",    color="steelblue", size=5)
p_curves.scatter("q", "nss", source=source_curves, legend_label="#E SS", color="tomato",    size=5, marker="triangle")
# p_curves.scatter("q", "np", source=source_curves, legend_label="#Points", color="gold", size=4, marker="square")
p_curves.legend.location = "top_left"
p_curves.add_tools(HoverTool(tooltips=[("q", "@q{0}"), ("#E", "@nc{0}"), ("#E SS", "@nss{0}"), ("#Points", "@np{0}")]))

p_traces = figure(
    height=400, width=900,
    x_axis_label="q (prime)", y_axis_label="trace t",
    tools="xpan,ypan,xwheel_zoom,ywheel_zoom,box_zoom,reset",
    x_range=detail_range,
)
p_traces.scatter("q", "t",   source=source_traces, color="color", size=4, alpha=0.7, legend_label="t")
p_traces.scatter("q", "fpi", source=source_fpi,    color="gold",  size=4, alpha=0.7, legend_label="f_pi")
p_traces.legend.location = "top_left"
p_traces.add_tools(HoverTool(tooltips=[("q", "@q{0}"), ("t", "@t{0}"), ("f_pi|N", "@color"), ("f_pi", "@fpi{0}")]))

p_overview = figure(
    height=130, width=900,
    x_axis_label="q", y_axis_type=None,
    tools="", toolbar_location=None,
    x_range=overview_range,
)
p_overview.line("q", "val", source=source, color="steelblue")

p_dk_hist = figure(
    height=300, width=900,
    x_axis_label="|D_K|", y_axis_label="count",
    tools="xpan,xwheel_zoom,reset",
)
p_dk_hist.scatter("dk", "count", source=source_dk_hist, color="steelblue", size=8)
p_dk_hist.add_tools(HoverTool(tooltips=[("D_K", "@dk{0}"), ("count", "@count{0}")]))
p_dk_hist.add_tools(TapTool())

range_tool = RangeTool(x_range=detail_range, start_gesture="pan")
range_tool.overlay.fill_color = "steelblue"
range_tool.overlay.fill_alpha = 0.2
p_overview.add_tools(range_tool)

hecke = HeckeOperator(CuspForm(int(N_input.value), int(k_input.value) + 2))


def recompute():
    pmin, pmax = state["pmin"], state["pmax"]
    N, k_off = int(N_input.value), int(k_input.value)
    weight = k_off + 2
    hecke.update_target(CuspForm(N, weight))
    dk_filter = int(dk_input.value) if dk_input.value.strip() else None
    global _all_trace_qs, _all_trace_dks
    qs, vals, refs, vals_raw, refs_raw, ncs, nsss, npts = [], [], [], [], [], [], [], []
    tq_qs, tq_ts, tq_colors = [], [], []
    tq_dks = []
    fpi_qs, fpis = [], []
    primes = list(primerange(pmin, pmax + 1))
    for prime in tqdm(primes, desc=f"N={N} w={weight} [{pmin},{pmax}]", unit="p", ncols=60):
        r = hecke.trace(prime, 1)
        norm = r.q ** ((weight - 1) / 2) if weight > 1 else 1
        # check if this q has any trace matching the D_K filter
        if dk_filter is not None:
            has_match = any(
                (ic := hecke.get_isogeny_class(r.q, t)) is not None and ic.field.D == dk_filter
                for t in (r.traces or [])
            )
            if not has_match:
                continue
        qs.append(r.q)
        vals.append(int(r.val) / norm)
        refs.append(int(r.val - r.eis_term) / norm)
        vals_raw.append(int(r.val))
        refs_raw.append(int(r.reference_val))
        ncs.append(r.num_curves)
        nsss.append(r.num_ss_curves)
        npts.append(r.num_points)
        if r.traces:
            for t in r.traces:
                ic = hecke.get_isogeny_class(r.q, t)
                if dk_filter is not None and (ic is None or ic.field.D != dk_filter):
                    continue
                if ic is not None and ic.field.D == -3:
                    color = "cyan"
                elif ic is not None and ic.field.D == -4:
                    color = "magenta"
                elif ic is not None and ic.f_pi % N == 0:
                    color = "limegreen"
                else:
                    color = "white"
                tq_qs.append(r.q)
                tq_ts.append(t)
                tq_colors.append(color)
                tq_dks.append(ic.field.D if ic is not None else 0)
                if ic is not None:
                    fpi_qs.append(r.q); fpis.append(ic.f_pi)
    if not qs:
        for src in (source, source_curves, source_traces, source_fpi):
            src.data = {k: [] for k in src.data}
        return

    source.data = dict(
        q=np.array(qs,        dtype=np.float64),
        val=np.array(vals,     dtype=np.float64),
        ref=np.array(refs,     dtype=np.float64),
        val_raw=np.array(vals_raw, dtype=np.float64),
        ref_raw=np.array(refs_raw, dtype=np.float64),
    )
    source_curves.data = dict(
        q=np.array(qs,    dtype=np.float64),
        nc=np.array(ncs,   dtype=np.float64),
        nss=np.array(nsss, dtype=np.float64),
        np=np.array(npts,  dtype=np.float64),
    )
    # save selection before replacing data (Bokeh clears it on .data assignment)
    prev = source_traces.selected.indices
    prev_q = int(source_traces.data["q"][prev[0]]) if prev and prev[0] < len(source_traces.data["q"]) else None
    prev_t = int(source_traces.data["t"][prev[0]]) if prev and prev[0] < len(source_traces.data["t"]) else None

    source_traces.data = dict(
        q=np.array(tq_qs, dtype=np.float64),
        t=np.array(tq_ts,  dtype=np.float64),
        color=tq_colors,
    )
    source_fpi.data = dict(q=np.array(fpi_qs, dtype=np.float64), fpi=np.array(fpis, dtype=np.float64))
    _all_trace_qs[:] = tq_qs
    _all_trace_dks[:] = tq_dks

    if prev_q is not None:
        for i, (q2, t2) in enumerate(zip(tq_qs, tq_ts)):
            if q2 == prev_q and t2 == prev_t:
                source_traces.selected.indices = [i]
                break

    overview_range.start = pmin
    overview_range.end   = pmax
    detail_range.start = max(pmin, detail_range.start)
    detail_range.end   = min(pmax, detail_range.end)
    if detail_range.start >= detail_range.end:
        detail_range.start = pmin
        detail_range.end   = min(pmin + 500, pmax)
    update_dk_hist()


def update_dk_hist():
    lo, hi = detail_range.start, detail_range.end
    from collections import Counter
    counts = Counter(
        dk for q, dk in zip(_all_trace_qs, _all_trace_dks)
        if lo <= q <= hi and dk != 0
    )
    if not counts:
        source_dk_hist.data = dict(dk=[], count=[], left=[], right=[])
        return
    dks    = sorted(counts, key=abs)
    source_dk_hist.data = dict(
        dk=[abs(d) for d in dks],
        count=[counts[d] for d in dks],
        left=[0] * len(dks),
        right=[counts[d] for d in dks],
    )


def shift(direction: int):
    delta = direction * STEP
    new_pmin = max(2, state["pmin"] + delta)
    actual_delta = new_pmin - state["pmin"]
    state["pmin"] = new_pmin
    state["pmax"] += actual_delta
    detail_range.start += actual_delta
    detail_range.end   += actual_delta
    recompute()


detail_table = DataTable(
    source=source_detail,
    columns=[
        TableColumn(field="prop", title="Property", formatter=StringFormatter()),
        TableColumn(field="num",  title="Value"),
    ],
    width=280, height=300, index_position=None,
)

torsion_table = DataTable(
    source=source_torsion,
    columns=[
        TableColumn(field="label",       title="class",      formatter=StringFormatter()),
        TableColumn(field="f",           title="f",          formatter=StringFormatter()),
        TableColumn(field="inv",         title="E inv",      formatter=StringFormatter()),
        TableColumn(field="torsion_inv", title="N-torsion",  formatter=StringFormatter()),
        TableColumn(field="h",           title="h(f)",       formatter=StringFormatter()),
        TableColumn(field="hw",          title="hw(f)",      formatter=StringFormatter()),
        TableColumn(field="n_pts",       title="#pts(N)",    formatter=StringFormatter()),
        TableColumn(field="value",       title="weighted",   formatter=StringFormatter()),
    ],
    width=700, height=400, index_position=None,
)


def on_trace_selected(attr, old, new):
    if not new:
        return
    indices = new
    idx = indices[0]
    q   = int(source_traces.data["q"][idx])
    t   = int(source_traces.data["t"][idx])
    ic  = hecke.get_isogeny_class(q, t)
    if ic is None:
        return
    rows = [
        ("q", int(q)),
        ("t", int(t)),
        ("D_pi", int(ic.D_pi)),
        ("D_K", int(ic.field.D)),
        ("f_pi", int(ic.f_pi)),
        ("n_pts", int(ic.n_pts)),
        ("size H(D)", int(ic.size)),
        ("ordinary", int(ic.ordinary)),
        ("is_quat", int(ic.is_quaternion)),
        ("unit_index w", int(ic.field._maximal_unit_index)),
    ]
    N = int(N_input.value)
    tor = ic.get_torsion(N)
    labels, fs, hs, hws, npts, vals, invs, torsion_invs = [], [], [], [], [], [], [], []
    header = f"t={t}  D_K={ic.field.D}  f_π={ic.f_pi}"
    for c in tor.conductor_levels:
        labels.append(header); header = ""
        fs.append(str(c.f))
        hs.append(str(ic.field.h(c.f)))
        hws.append(str(ic.field.hw(c.f)))
        npts.append(str(c.n_pts_exact_order))
        vals.append(str(c.value))
        invs.append(fmt_invariants(c.inv, color=False))
        torsion_invs.append(fmt_invariants(c.torsion_inv, color=False))
    source_torsion.data = dict(label=labels, f=fs, h=hs, hw=hws, n_pts=npts, value=vals, inv=invs, torsion_inv=torsion_invs)

    # clear first to force length change → Bokeh re-renders the table
    source_detail.data = dict(prop=[], num=[])
    source_detail.data = dict(prop=[r[0] for r in rows], num=[r[1] for r in rows])


def on_detail_selected(attr, old, new):
    if not new:
        return
    idx = new[0]
    q = int(source.data["q"][idx])
    N = int(N_input.value)
    # collect all isogeny classes for this q that have N-torsion
    all_fs, all_hs, all_hws, all_npts, all_vals, all_invs, all_torsion_invs = [], [], [], [], [], [], []
    all_rows = [("q", q)]
    all_labels, all_fs, all_hs, all_hws, all_npts, all_vals, all_invs, all_torsion_invs = [], [], [], [], [], [], [], []
    for t_val, ic in (hecke._curves_cache[q].isogeny_classes if q in hecke._curves_cache else []):
        if ic.n_pts % N != 0:
            continue
        tor = ic.get_torsion(N)
        all_rows += [("t", int(t_val)), ("D_K", int(ic.field.D)), ("f_pi", int(ic.f_pi))]
        header = f"t={t_val}  D_K={ic.field.D}  f_π={ic.f_pi}"
        for c in tor.conductor_levels:
            all_labels.append(header); header = ""
            all_fs.append(str(c.f))
            all_hs.append(str(ic.field.h(c.f)))
            all_hws.append(str(ic.field.hw(c.f)))
            all_npts.append(str(c.n_pts_exact_order))
            all_vals.append(str(c.value))
            all_invs.append(fmt_invariants(c.inv, color=False))
            all_torsion_invs.append(fmt_invariants(c.torsion_inv, color=False))
    source_torsion.data = dict(label=all_labels, f=all_fs, h=all_hs, hw=all_hws, n_pts=all_npts, value=all_vals, inv=all_invs, torsion_inv=all_torsion_invs)
    source_detail.data = dict(prop=[], num=[])
    source_detail.data = dict(prop=[r[0] for r in all_rows], num=[r[1] for r in all_rows])

source.selected.on_change("indices", on_detail_selected)
source_traces.selected.on_change("indices", on_trace_selected)
p_traces.add_tools(TapTool())

N_input.on_change("value",  lambda a, o, n: recompute())
k_input.on_change("value",  lambda a, o, n: recompute())
dk_input.on_change("value", lambda a, o, n: recompute())
detail_range.on_change("start", lambda a, o, n: update_dk_hist())
detail_range.on_change("end",   lambda a, o, n: update_dk_hist())

def on_dk_selected(attr, old, new):
    if not new:
        return
    idx = new[0]
    dk_abs = int(source_dk_hist.data["dk"][idx])
    dk_input.value = str(-dk_abs)
    source_dk_hist.selected.indices = []

source_dk_hist.selected.on_change("indices", on_dk_selected)
btn_left.on_click( lambda _: shift(-1))
btn_right.on_click(lambda _: shift(+1))

recompute()

plots = column(
    row(N_input, k_input, dk_input, btn_left, btn_right),
    p_detail,
    p_curves,
    p_traces,
    p_overview,
)
curdoc().add_root(column(
    row(plots, p_dk_hist),
    row(detail_table, torsion_table),
))
