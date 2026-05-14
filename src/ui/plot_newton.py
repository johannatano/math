import sys
import warnings
import math
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
warnings.filterwarnings("ignore", category=UserWarning, module="bokeh")

from bokeh.models import ColumnDataSource, RangeTool, Spinner, Button, Range1d
from bokeh.layouts import column, row
from bokeh.plotting import figure, curdoc
from nt.modular_forms import CuspForm, HeckeOperator
from tqdm import tqdm
from nt.common import valuation
STEP = 1                          # k shifts by even steps (weights are even)
state = {"kmin": 0, "kmax": 50}

curdoc().theme = "dark_minimal"

p_input   = Spinner(title="Prime p", value=2, low=2,  high=10000, step=1)
n_input   = Spinner(title="Power n", value=8,  low=1,  high=20,    step=1)
N_input   = Spinner(title="Level N", value=1,  low=1,  high=10000, step=1)
t_input   = Spinner(title="Trace t", value=1,  low=-10000, high=10000, step=1)
btn_left  = Button(label="◀  k-2", width=110, button_type="default")
btn_right = Button(label="k+2  ▶", width=110, button_type="default")

source    = ColumnDataSource(dict(k=[], val=[], ref=[]))
source_hk = ColumnDataSource(dict(k=[], hk=[]))

detail_range   = Range1d(start=state["kmin"], end=state["kmax"])
overview_range = Range1d(start=state["kmin"], end=state["kmax"])

p_detail = figure(
    height=350, width=900,
    x_axis_label="weight k", y_axis_label="log10|a_{p^n}|",
    tools="xpan,xwheel_zoom,reset",
    x_range=detail_range,
)
p_detail.scatter("k", "val", source=source, legend_label="Computed", color="steelblue", size=8)
p_detail.scatter("k", "ref", source=source, legend_label="Sage ref",  color="orange",   size=5, marker="triangle")
p_detail.legend.location = "top_left"

p_overview = figure(
    height=130, width=900,
    x_axis_label="k", y_axis_type=None,
    tools="", toolbar_location=None,
    x_range=overview_range,
)
p_overview.line("k", "val", source=source, color="steelblue")

p_hk = figure(
    height=250, width=900,
    x_axis_label="weight k", y_axis_label="hk(t,k,q) / q^(k/2)",
    tools="xpan,xwheel_zoom,reset",
    x_range=detail_range,  # shared x-range with main plot
)
p_hk.scatter("k", "hk", source=source_hk, color="mediumpurple", size=6)

range_tool = RangeTool(x_range=detail_range, start_gesture="pan")
range_tool.overlay.fill_color = "steelblue"
range_tool.overlay.fill_alpha = 0.2
p_overview.add_tools(range_tool)

hecke = HeckeOperator(CuspForm(int(N_input.value), 2))


def recompute_hk():
    p = int(p_input.value)
    n = int(n_input.value)
    t = int(t_input.value)
    q = p ** n
    kmin, kmax = state["kmin"], state["kmax"]
    hks = []
    for k in range(kmin, kmax + 1):
        raw = sum(
            math.comb(k - j, j) * (-q) ** j * t ** (k - 2 * j)
            for j in range(k // 2 + 1)
        )
        norm = q ** (k / 2) if k > 0 else 1
        hks.append(raw / norm)
    source_hk.data = dict(
        k=np.array(range(kmin, kmax + 1), dtype=np.float64),
        hk=np.array(hks, dtype=np.float64),
    )


def recompute():
    p = int(p_input.value)
    n = int(n_input.value)
    N = int(N_input.value)
    kmin, kmax = state["kmin"], state["kmax"]

    ks, vals, refs = [], [], []
    for w in tqdm(range(kmin, kmax + 1), desc=f"p={p}^{n} N={N} k={kmin}..{kmax}", unit="k", ncols=60):
        hecke.update_target(CuspForm(N, w))
        r = hecke.trace(p, n)
        ks.append(w)
        v = int(r.val)
        vals.append(math.copysign(math.log10(abs(v)), v) if v != 0 else float("nan"))
        refs.append(valuation(int(r.reference_val), p) if r.reference_val != 0 else float("inf"))

    source.data = dict(
        k=np.array(ks,    dtype=np.float64),
        val=np.array(vals, dtype=np.float64),
        ref=np.array(refs, dtype=np.float64),
    )
    recompute_hk()
    overview_range.start = kmin
    overview_range.end   = kmax
    detail_range.start = max(kmin, detail_range.start)
    detail_range.end   = min(kmax, detail_range.end)
    if detail_range.start >= detail_range.end:
        detail_range.start = kmin
        detail_range.end   = kmax


def shift(direction: int):
    delta = direction * STEP
    new_kmin = max(2, state["kmin"] + delta)
    actual = new_kmin - state["kmin"]
    state["kmin"] = new_kmin
    state["kmax"] += actual
    detail_range.start += actual
    detail_range.end   += actual
    recompute()


p_input.on_change("value", lambda a, o, n: recompute())
n_input.on_change("value", lambda a, o, n: recompute())
N_input.on_change("value", lambda a, o, n: recompute())
t_input.on_change("value", lambda a, o, n: recompute_hk())
btn_left.on_click( lambda _: shift(-1))
btn_right.on_click(lambda _: shift(+1))

recompute()

curdoc().add_root(column(
    row(p_input, n_input, N_input, t_input, btn_left, btn_right),
    p_detail,
    p_hk,
    p_overview,
))
