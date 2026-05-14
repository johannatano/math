import sys
import warnings
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
warnings.filterwarnings("ignore", category=UserWarning, module="bokeh")

from sympy import primerange
from bokeh.models import ColumnDataSource, RangeTool, Spinner, Button, Range1d
from bokeh.layouts import column, row
from bokeh.plotting import figure, curdoc
from nt.modular_forms import CuspForm, HeckeOperator
from tqdm import tqdm

STEP = 1000
state = {"pmin": 2, "pmax": 5000}

curdoc().theme = "dark_minimal"

N_input   = Spinner(title="Level N", value=500, low=1,   high=10000, step=1)
k_input   = Spinner(title="Weight k", value=2,  low=0,   high=10000,  step=1)
btn_left  = Button(label="◀  -1000", width=110, button_type="default")
btn_right = Button(label="+1000  ▶", width=110, button_type="default")

source = ColumnDataSource(dict(q=[], val=[], ref=[]))

detail_range  = Range1d(start=state["pmin"], end=min(state["pmin"] + 500, state["pmax"]))
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
    weight = k_off + 2  # actual modular weight
    hecke.update_target(CuspForm(N, weight))
    qs, vals, refs = [], [], []
    primes = list(primerange(pmin, pmax + 1))
    for prime in tqdm(primes, desc=f"N={N} w={weight} [{pmin},{pmax}]", unit="p", ncols=60):
        r = hecke.trace(prime, 1)
        norm = r.q ** ((weight - 1) / 2) if weight > 1 else 1
        qs.append(r.q)
        vals.append(int(r.val) / norm)
        refs.append(int(r.reference_val) / norm)
    source.data = dict(
        q=np.array(qs, dtype=np.float64),
        val=np.array(vals, dtype=np.float64),
        ref=np.array(refs, dtype=np.float64),
    )
    overview_range.start = pmin
    overview_range.end   = pmax
    # clamp detail window to stay within the new data extent
    detail_range.start = max(pmin, detail_range.start)
    detail_range.end   = min(pmax, detail_range.end)
    if detail_range.start >= detail_range.end:
        detail_range.start = pmin
        detail_range.end   = min(pmin + 500, pmax)


def shift(direction: int):
    delta = direction * STEP
    new_pmin = max(2, state["pmin"] + delta)
    actual_delta = new_pmin - state["pmin"]  # real shift after clamping
    state["pmin"] = new_pmin
    state["pmax"] += actual_delta
    detail_range.start += actual_delta
    detail_range.end   += actual_delta
    recompute()


N_input.on_change("value",  lambda a, o, n: recompute())
k_input.on_change("value",  lambda a, o, n: recompute())
btn_left.on_click( lambda _: shift(-1))   # on_click passes click count as arg
btn_right.on_click(lambda _: shift(+1))

recompute()

curdoc().add_root(column(
    row(N_input, k_input, btn_left, btn_right),
    p_detail,
    p_overview,
))
