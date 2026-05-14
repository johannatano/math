import numpy as np

from bokeh.layouts import column
from bokeh.models import ColumnDataSource, RangeTool
from bokeh.plotting import figure, figure, show, curdoc
from bokeh.io import output_file, save
from bokeh.sampledata.stocks import AAPL

curdoc().theme = "dark_minimal"

dates = np.array(AAPL["date"], dtype=np.datetime64)
source = ColumnDataSource(data=dict(date=dates, close=AAPL["adj_close"]))

p = figure(
    height=300,
    width=800,
    tools="xpan,xwheel_zoom,reset",
    x_axis_type="datetime",
    x_axis_location="above",
    window_axis="x",
    x_range=(dates[1500], dates[2500]),
)

p.line("date", "close", source=source)
p.yaxis.axis_label = "Price"

select = figure(
    title="Drag the middle and edges of the selection box to change the range above",
    height=130,
    width=800,
    x_axis_type="datetime",
    y_axis_type=None,
    tools="",
    toolbar_location=None,
)
select.x_range.range_padding = 0
select.x_range.bounds = "auto"

range_tool = RangeTool(x_range=p.x_range, start_gesture="pan")
range_tool.overlay.fill_color = "navy"
range_tool.overlay.fill_alpha = 0.2

select.line("date", "close", source=source)
select.ygrid.grid_line_color = None
select.add_tools(range_tool)

output_file("plot.html")
save(column(p, select))
print("Saved: plot.html")
