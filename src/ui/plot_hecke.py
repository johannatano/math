import sys
import warnings
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
warnings.filterwarnings("ignore", category=UserWarning, module="bokeh")

from sympy import primerange
from bokeh.models import (ColumnDataSource, RangeTool, Spinner, Button, Range1d,
                          HoverTool, TapTool, DataTable, TableColumn, NumberFormatter, StringFormatter, Div)
from bokeh.layouts import column, row
from bokeh.plotting import figure, curdoc
from nt.modular_forms import CuspForm, HeckeOperator
from nt.common import fmt_invariants
from tqdm import tqdm

STEP = 1000
state = {"pmin": 2, "pmax": 50000}

curdoc().theme = "dark_minimal"

N_input   = Spinner(title="Level N", value=500, low=1,   high=10000, step=1)
k_input   = Spinner(title="Weight k", value=2,  low=0,   high=10000, step=1)
btn_left  = Button(label="◀  -1000", width=110, button_type="default")
btn_right = Button(label="+1000  ▶", width=110, button_type="default")

source        = ColumnDataSource(dict(q=[], val=[], ref=[], val_raw=[], ref_raw=[]))
source_curves = ColumnDataSource(dict(q=[], nc=[], nss=[], np=[]))
source_traces = ColumnDataSource(dict(q=[], t=[], color=[]))
source_fpi    = ColumnDataSource(dict(q=[], fpi=[]))
source_detail = ColumnDataSource(dict(prop=[], num=[]))
source_torsion = ColumnDataSource(dict(f=[], h=[], hw=[], n_pts=[], value=[], inv=[], torsion_inv=[]))

detail_range   = Range1d(start=state["pmin"], end=min(state["pmin"] + 500, state["pmax"]))
overview_range = Range1d(start=state["pmin"], end=state["pmax"])

p_detail = figure(
    height=350, width=900,
    x_axis_label="q (prime)", y_axis_label="a_q / q^((weight-1)/2)",
    tools="xpan,xwheel_zoom,reset",
    x_range=detail_range,
)
p_detail.scatter("q", "val", source=source, legend_label="Computed", color="steelblue", size=6)
p_detail.scatter("q", "ref", source=source, legend_label="Sage ref",  color="orange",   size=4, marker="triangle")
p_detail.legend.location = "top_left"
p_detail.add_tools(HoverTool(tooltips=[
    ("q",            "@q{0}"),
    ("computed",     "@val{0.0000}  (@val_raw{0})"),
    ("sage ref",     "@ref{0.0000}  (@ref_raw{0})"),
]))

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
    qs, vals, refs, vals_raw, refs_raw, ncs, nsss, npts = [], [], [], [], [], [], [], []
    tq_qs, tq_ts, tq_colors = [], [], []
    fpi_qs, fpis = [], []
    primes = list(primerange(pmin, pmax + 1))
    for prime in tqdm(primes, desc=f"N={N} w={weight} [{pmin},{pmax}]", unit="p", ncols=60):
        r = hecke.trace(prime, 1)
        norm = r.q ** ((weight - 1) / 2) if weight > 1 else 1
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
                if ic is not None:
                    fpi_qs.append(r.q); fpis.append(ic.f_pi)
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
    source_traces.data = dict(
        q=np.array(tq_qs, dtype=np.float64),
        t=np.array(tq_ts,  dtype=np.float64),
        color=tq_colors,
    )
    source_fpi.data = dict(q=np.array(fpi_qs, dtype=np.float64), fpi=np.array(fpis, dtype=np.float64))
    overview_range.start = pmin
    overview_range.end   = pmax
    detail_range.start = max(pmin, detail_range.start)
    detail_range.end   = min(pmax, detail_range.end)
    if detail_range.start >= detail_range.end:
        detail_range.start = pmin
        detail_range.end   = min(pmin + 500, pmax)


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
        TableColumn(field="f",           title="f"),
        TableColumn(field="inv",         title="E inv",      formatter=StringFormatter()),
        TableColumn(field="torsion_inv", title="N-torsion",  formatter=StringFormatter()),
        TableColumn(field="h",           title="h(f)"),
        TableColumn(field="hw",          title="hw(f)"),
        TableColumn(field="n_pts",       title="#pts(N)"),
        TableColumn(field="value",       title="weighted"),
    ],
    width=550, height=300, index_position=None,
)


def on_trace_selected(attr, old, new):
    indices = new if new else old
    if not indices:
        return
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
    fs, hs, hws, npts, vals, invs, torsion_invs = [], [], [], [], [], [], []
    for c in tor.conductor_levels:
        fs.append(str(c.f))
        hs.append(float(ic.field.h(c.f)))
        hws.append(float(ic.field.hw(c.f)))
        npts.append(int(c.n_pts_exact_order))
        vals.append(float(c.value))
        invs.append(fmt_invariants(c.inv, color=False))
        torsion_invs.append(fmt_invariants(c.torsion_inv, color=False))
    source_torsion.data = dict(f=fs, h=hs, hw=hws, n_pts=npts, value=vals, inv=invs, torsion_inv=torsion_invs)

    # clear first to force length change → Bokeh re-renders the table
    source_detail.data = dict(prop=[], num=[])
    source_detail.data = dict(prop=[r[0] for r in rows], num=[r[1] for r in rows])


source_traces.selected.on_change("indices", on_trace_selected)
p_traces.add_tools(TapTool())

N_input.on_change("value",  lambda a, o, n: recompute())
k_input.on_change("value",  lambda a, o, n: recompute())
btn_left.on_click( lambda _: shift(-1))
btn_right.on_click(lambda _: shift(+1))

recompute()

plots = column(
    row(N_input, k_input, btn_left, btn_right),
    p_detail,
    p_curves,
    p_traces,
    p_overview,
)
curdoc().add_root(row(plots, column(detail_table, torsion_table)))
