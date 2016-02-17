from os.path import dirname, join

import pandas as pd
import numpy as np
import math
import pickle
from bokeh.plotting import figure, show, ColumnDataSource
from bokeh.models import CustomJS, VBox, HBox, Select
from bokeh.io import show, vform, curdoc, hplot, vplot


def expression_by_age(source_circle, source_line):
    tools=""

    factors = list(source_circle.to_df()['age'].unique())

    p = figure(title="Gene Expression by Age", x_range=factors, plot_height=500, plot_width=800, tools=tools)

    p.circle('age','rpkm', fill_alpha=0.2, size=5, source=source_circle)
    p.xaxis.major_label_orientation = math.pi/4

    p.line('age', 'rpkm', source=source_line)

    # change just some things about the x-grid
    p.xgrid.grid_line_color = None

    # change just some things about the y-grid
    p.ygrid.grid_line_alpha = 0.5
    p.ygrid.grid_line_dash = [6, 4]

    p.background_fill_color = '#FFFFFF'

    return p


def expression_by_structure(source):
    tools=""
    factors = list(source.to_df()['structure'].unique())

    r = figure(title="Gene Expression by Brain Region", x_range=factors, plot_height=700, plot_width=800, tools=tools)
    r.rect('structure',
           'median',
           width=0.6,
           height='iqr', source=source)

    #r.circle(x='structure', y='rpkm', size=0.1, source=source)
    r.segment('structure','lower','structure','upper', line_width=2, source=source)
    r.rect('structure', 'lower', 0.2, 0.01, source=source)
    r.rect('structure', 'upper', 0.2, 0.01, source=source)
    r.rect('structure', 'median', 0.6, 0.005, line_color="black", source=source)

    r.xaxis.major_label_orientation = math.pi/4

    return r

def update_plot(attrname, old, new):
    gene = str(select.value)
    src1 = ColumnDataSource(structure_stats[gene])
    src2 = ColumnDataSource(sources[gene])
    src3 = ColumnDataSource(line[gene])
    source1.data.update(src1.data)
    source2.data.update(src2.data)
    source3.data.update(src3.data)

genes = ['WDR4', 'RPN1', 'EPHA5', 'WDR7', 'ROBO3', 'MARCH7', 'ZNF665', 'CSMD2', 'FYN', 'SUV420H1', 'ETFB', 'KYNU']

with open('app/data/structure_stats.p', 'rb') as fp1:
    structure_stats = pickle.load(fp1)

with open('app/data/sources.p', 'rb') as fp2:
    sources = pickle.load(fp2)

with open('app/data/line.p', 'rb') as fp3:
    line = pickle.load(fp3)

gene = genes[0]

source1 = ColumnDataSource(structure_stats[gene])
source2 = ColumnDataSource(sources[gene])
source3 = ColumnDataSource(line[gene])

structure_plot = expression_by_structure(source1)
age_plot = expression_by_age(source2, source3)

plot = hplot(structure_plot,age_plot)

select = Select(title="Gene:", value=genes[0], options=genes)
select.on_change('value', update_plot)

curdoc().add_root(vplot(select, plot))
