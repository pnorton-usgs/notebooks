# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import numpy as np
from bokeh.plotting import *
import pandas as pd
N = 80
output_notebook()

# %%
from bokeh import plotting
colors = [
      "#1f77b4",
      "#ff7f0e", 
      "#2ca02c", "#98df8a",
      "#d62728", "#ff9896",
      "#9467bd", "#c5b0d5",
      "#8c564b", "#c49c94",
      "#e377c2", "#f7b6d2",
      "#7f7f7f", "#ffbb78",
      "#bcbd22", "#dbdb8d",
      "#17becf", "#9edae5"
    ]
def make_color_cycle():
    ic = iter(colors)
    while True:
        try:
            yield ic.next()
        except StopIteration as e:
            ic = iter(colors)


def multiplot(df, colsets, **kwargs):
    x_range = [None]

    x=df.index
    columns = dict([[c, df[c]] for c in df.columns])
    columns['x'] = df.index
    datasource = plotting.ColumnDataSource(columns)
    def default_plot(c__, title=False):
        defaults = dict(f=plotting.line)
        defaults.update(kwargs)
        if type(c__) == type({}):
            defaults.update(c__)
            colname = defaults.pop('col')
        else:
            colname = c__
        plotf = defaults.pop('f')
        l = plotf('x', colname, source=datasource,
                x_range=x_range[0], color=color_cycle.next(), 
                legend=colname, title=title or colname, **defaults)
        x_range[0] = l.x_range
        return l
    rows = []

    color_cycle = make_color_cycle()
    for col_or_colset in colsets:
        if type(col_or_colset) == type([]):
            plotting.hold(True)
            l = default_plot(col_or_colset[0])
            for c in col_or_colset[1:]:
                default_plot(c)
            plotting.hold(False)
        else:
            l = default_plot(col_or_colset)

        rows.append([l])
        plotting.figure()
    plotting.gridplot(rows)
    plotting.show()


# %%

x = np.linspace(0, 4*np.pi, N)

df = pd.DataFrame( index=x)
df['y'] = np.sin(x)
df['z'] = np.sin(2*x) + 1.75 * np.cos(x)
df['y2'] = np.sin(1.7*x)
df['w'] = np.tan(x)
multiplot(df, [dict(col='z', f=plotting.scatter), ['y', 'z'], 'w'], plot_width=900, plot_height=200, 
          f=plotting.scatter, tools="save, select, box_zoom, resize, crosshair")

# %%
multiplot(df, ['z', ['y', 'z'], 'w'], plot_width=900, plot_height=200, tools="save, box_zoom, resize, crosshair")

# %%
multiplot(df, ['y','y2'], plot_width=700, plot_height=200, tools="pan,box_zoom,previewsave,resize,crosshair")

# %%
#this plot causes an error, I'm not sure why.  Running it completely borks js on the page
multiplot(df, [['y','y2']], plot_width=700, plot_height=200, tools="pan,zoom,preview,resize,crosshair")

# %%
