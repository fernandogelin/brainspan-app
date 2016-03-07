from os.path import dirname, join

import pandas as pd
import numpy as np
import math
import pickle
import bokeh.palettes
from bokeh.plotting import figure, show, ColumnDataSource
from bokeh.models import CustomJS, VBox, HBox, Select, HoverTool
from bokeh.models.widgets import MultiSelect
from bokeh.io import show, vform, curdoc, hplot, vplot
from sklearn import preprocessing

def expression_by_age(source_circle, source_line):
    tools="save"
    factors = list(source_circle.to_df()['age'].unique())
    p = figure(title="Gene Expression by Age",
               background_fill_color="#EFE8E2",
               x_range=factors,
               plot_height=450,
               plot_width=650,
               tools=tools)
    p.circle('age','rpkm', fill_alpha=1, color="#F38630", size=4, source=source_circle)
    p.xaxis.major_label_orientation = math.pi/4
    p.line('age', 'rpkm_smooth', line_color="#3B8686", line_width=3, source=source_line)
    # change just some things about the x-grid
    p.xgrid.grid_line_color = None
    # change just some things about the y-grid
    p.ygrid.grid_line_alpha = 0.5
    p.ygrid.grid_line_dash = [6, 4]
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = "white"
    p.grid.grid_line_width = 2
    p.xaxis.major_label_text_font_size="7pt"
    return p

def expression_by_structure(source):
    tools="save"
    factors = list(source.to_df()['structure'].unique())
    r = figure(title="Gene Expression by Brain Region",
               background_fill_color="#EFE8E2",
               x_range=factors,
               plot_height=585,
               plot_width=650,
               tools=tools)
    r.rect('structure',
           'median',
            width=0.6,
            height='iqr',
            source=source)
    #r.circle(x='structure', y='rpkm', size=0.1, source=source)
    r.segment('structure','lower','structure','upper', line_width=2, source=source)
    r.rect('structure', 'lower', 0.2, 0.01, source=source)
    r.rect('structure', 'upper', 0.2, 0.01, source=source)
    r.rect('structure', 'median', 0.6, 0.005, line_color="black", source=source)
    r.xaxis.major_label_orientation = math.pi/4
    r.xgrid.grid_line_color = None
    r.ygrid.grid_line_color = "white"
    r.ygrid.grid_line_dash = [6, 4]
    r.grid.grid_line_width = 2
    r.xaxis.major_label_text_font_size="7pt"
    return r

def get_dataset(df, str_list):
    new_df = df[df['structure'].isin(list(str_list))]
    return ColumnDataSource(new_df)

def update_plot_structures(attrname, old, new):
    gene = select.value
    src = get_dataset(sources[gene], multi_select.value)
    source2.data.update(src.data)

def update_plot(attrname, old, new):
    gene = select.value
    src1 = ColumnDataSource(structure_stats[gene])
    src2 = ColumnDataSource(sources[gene])
    src3 = ColumnDataSource(line[gene])
    source1.data.update(src1.data)
    source2.data.update(src2.data)
    source3.data.update(src3.data)

genes = ['WDR4', 'RPN1', 'EPHA5', 'WDR7', 'ROBO3', 'MARCH7', 'ZNF665', 'CSMD2', 'FYN', 'SUV420H1', 'ETFB', 'KYNU']
#genes = ["ARHGAP19", "FMN2", "SCLY", "WDR4", "ZNF714", "ABCA3", "PRRC2B", "ZNF665", "ABR", "ADAM18", "ADAMDEC1", "AGRN", "AHNAK", "AHNAK2", "AIRE", "ALS2CR11", "AMBRA1", "ANKRD26", "AP5Z1", "APOPT1", "ARFGAP2", "ARHGEF18", "ARSD", "ASB17", "ASH2L", "ASMTL", "ASTE1", "ATF6B", "ATXN7L1", "AUTS2", "AXIN1", "BCLAF1", "BMP4", "BTBD11", "BTN2A1", "C11orf65", "C1orf177", "C1orf86", "C20orf26", "C2CD3", "C3orf30", "C5orf60", "C8orf33", "CAD", "CASC10", "CCDC158", "CCDC33", "CCDC85B", "CCNE2", "CCSER1", "CD97", "CDK15", "CEP68", "CGRRF1", "CHRD", "CHRM4", "CNBD1", "COL18A1", "COL6A6", "COX6C", "CPE", "CRELD2", "CRYBA2", "CSMD2", "CUBN", "DDX27", "DECR2", "DHRSX", "DLK1", "DNAH17", "DNMT3L", "DOCK4", "DSC3", "DSG2", "DSPP", "DYRK3", "DYRK4", "EHBP1L1", "EPG5", "EPHA5", "EPO", "EPRS", "ERMN", "ERRFI1", "ETFDH", "EZH2", "FAM208A", "FAM83H", "FBXO42", "FERMT2", "FHOD3", "FOXJ1", "FRAS1", "FREM3", "GAL3ST3", "GANAB", "GART", "GBA", "GCFC2", "GOLT1A", "GOT1L1", "GPR113", "GPR151", "GPRIN1", "GRM5", "GUCY1A3", "HCN4", "HEPHL1", "HIST1H3A", "IL4I1", "KIAA0586", "KIAA0825", "KIF22", "KLHDC2", "KRTAP10-11", "KRTAP10-7", "LAMA3", "LAMC3", "LAX1", "LETM2", "LINGO1", "LRRC26", "LRRC32", "LSG1", "LY75", "MAPK15", "MARCH7", "MBLAC1", "MCCC1", "MCF2L", "MFSD4", "MICAL3", "MMP1", "MRGPRX1", "MRPL28", "MS4A15", "MSLN", "MUC4", "MYBPHL", "MYO1C", "MYO7A", "NAALAD2", "NDUFA11", "NEO1", "NFATC4", "NID1", "NLRC4", "NOTUM", "NPAP1", "NPHP4", "NPHS1", "NUTM1", "OGDH", "OLAH", "OR10G7", "OR11H6", "OR4D9", "OR8D4", "OR9A4", "OR9I1", "OSBPL8", "OXTR", "PARM1", "PAX9", "PCSK5", "PDZD2", "PELP1", "PERM1", "PEX5L", "PIDD1", "PIK3C2B", "PITRM1", "PKHD1", "PKP3", "PLA2G4C", "PLEKHG4B", "POM121L12", "PPP1R15B", "PPP1R35", "PRR25", "PRR27", "PRR35", "PRSS23", "PRSS27", "PSMD14", "RAB11FIP1", "RFXAP", "RGS21", "RHBDF1", "RIMBP2", "RIN3", "ROBO3", "RPN1", "RPRD2", "RPS6KC1", "SCN9A", "SERPINA4", "SEZ6", "SFTPC", "SGPP2", "SHCBP1L", "SLC14A2", "SLC26A1", "SLC26A9", "SLC2A6", "SLC39A5", "SLC4A4", "SLIT3", "SMAD3", "SMG1", "SNX16", "SOX30", "SPATA18", "SPATS2", "SPIRE2", "SPOCK3", "SPTB", "SRD5A1", "SSH3", "SSPO", "SSRP1", "STBD1", "SYNE3", "SYNGR1", "SYT5", "TAGLN", "TAS1R2", "TCF7", "TECTA", "THAP3", "THOC1", "THPO", "TIGD4", "TLN1", "TMEM14E", "TMEM198", "TMEM249", "TMEM8B", "TNXB", "TRAF5", "TRAPPC10", "TTC16", "TTC21B", "TTN", "TUBGCP6", "TWISTNB", "UFSP1", "UGT2B15", "ULK1", "UMODL1", "UNC79", "USP35", "VWA3A", "VWDE", "WDR27", "WDR7", "WDR90", "WHSC1L1", "WRAP73", "XDH", "ZASP", "ZBTB41", "ZC3H7A", "ZDBF2", "ZFYVE9", "ZNF214", "ZNF423", "ZNF438", "ZNF7", "ZNF765", "ZNFX1", "ZSCAN31"]


with open('brainspan-app/data/structure_stats.p', 'rb') as fp1:
    structure_stats = pickle.load(fp1)

with open('brainspan-app/data/sources.p', 'rb') as fp2:
    sources = pickle.load(fp2)

with open('brainspan-app/data/line.p', 'rb') as fp3:
    line = pickle.load(fp3)

gene = genes[0]

source1 = ColumnDataSource(structure_stats[gene])
source2 = ColumnDataSource(sources[gene])
source3 = ColumnDataSource(line[gene])

structure_plot = expression_by_structure(source1)
age_plot = expression_by_age(source2, source3)

plot = hplot(structure_plot,age_plot)

multi_select = MultiSelect(title="Brain Region:",
                           options=list(sources[gene]['structure'].unique()))

multi_select.on_change('value', update_plot_structures)

select = Select(title="Gene:", value=genes[0], options=genes)
select.on_change('value', update_plot)



curdoc().add_root(vplot(select, plot, multi_select))
