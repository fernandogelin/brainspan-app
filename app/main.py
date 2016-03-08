from os.path import dirname, join

from math import pi
import pandas as pd
import numpy as np
import json
from pymongo import MongoClient
from scipy.signal import savgol_filter
from bokeh.plotting import figure, output_notebook, show, ColumnDataSource, vplot, output_server, hplot
from bokeh.models import CustomJS, VBox, HBox, Select, MultiSelect
from bokeh.io import output_file, show, vform, curdoc
from bokeh.charts import BoxPlot
from bokeh.palettes import RdPu9
from bokeh.models import HoverTool
from sklearn import preprocessing

categories=[u'8 pcw', u'9 pcw', u'12 pcw', u'13 pcw', u'16 pcw', u'17 pcw',
        u'19 pcw', u'21 pcw', u'24 pcw', u'25 pcw', u'26 pcw', u'35 pcw',
        u'37 pcw', u'4 mos', u'10 mos', u'1 yrs', u'2 yrs', u'3 yrs',
        u'4 yrs', u'8 yrs', u'11 yrs', u'13 yrs', u'15 yrs', u'18 yrs',
        u'19 yrs', u'21 yrs', u'23 yrs', u'30 yrs', u'36 yrs', u'37 yrs',
        u'40 yrs']
sorterIndex = dict(zip(categories,range(len(categories))))

def get_dataframes(gene, structures):
    query = {"gene":gene, "structure_name": {"$in": structures}}
    cursor = db.brainspan2.find(query)
    df = pd.DataFrame(list(cursor))
    df = df[['gene','age','structure_name','rpkm']]

    df_line = pd.pivot_table(df, values='rpkm', index='age', aggfunc=np.mean).to_frame().reindex(index=list(df['age'].unique()))
    df_line = df_line.reset_index()
    df_line['rank'] = df_line['age'].map(sorterIndex)
    df_line.sort_values(by='rank', ascending = True, inplace = True)
    df_line.drop('rank', 1, inplace = True)

    if len(structures)%2 == 0:
        window = len(structures) + 1
    else:
        window = len(structures)
    df_line['rpkm_smooth'] = savgol_filter(df_line['rpkm'], window, 3)

    return df, df_line

def get_boxplot_data(source_df):

    """ Takes a dataframe with expression for one gene, returns datasources for boxplot """
    structures = list(source_df.structure_name.unique())

    # find the quartiles and IQR for each category
    groups = source_df.groupby('structure_name')
    q1 = groups.quantile(q=0.25)
    q2 = groups.quantile(q=0.5)
    q3 = groups.quantile(q=0.75)
    iqr = q3 - q1
    upper = q3 + 1.5*iqr
    lower = q1 - 1.5*iqr

    # find the outliers for each category
    def outliers(group):
        structure = group.name
        return group[(group.rpkm > upper.loc[structure][0]) | (group.rpkm < lower.loc[structure][0])]['rpkm']

    out = groups.apply(outliers).dropna()

    # prepare outlier data for plotting, we need coordinates for every outlier.
    outx = []
    outy = []
    for structure in structures:
        # only add outliers if they exist
        if not out.loc[structure].empty:
            for value in out[structure]:
                outx.append(structure)
                outy.append(value)

    # if no outliers, shrink lengths of stems to be no longer than the minimums or maximums
    qmin = groups.quantile(q=0.00)
    qmax = groups.quantile(q=1.00)
    upper.rpkm = [min([x,y]) for (x,y) in zip(list(qmax.iloc[:,0]),upper.rpkm) ]
    lower.rpkm = [max([x,y]) for (x,y) in zip(list(qmin.iloc[:,0]),lower.rpkm) ]

    full_data_frame = pd.DataFrame({'q1': q1.rpkm,
                                    'q2': q2.rpkm,
                                    'q3': q3.rpkm,
                                    'upper': upper.rpkm,
                                    'lower': lower.rpkm,
                                    'top_box_1': (q3.rpkm+q2.rpkm)/2,
                                    'top_box_2': q3.rpkm-q2.rpkm,
                                    'bottom_box_1':(q2.rpkm+q1.rpkm)/2,
                                    'bottom_box_2': q2.rpkm-q1.rpkm})

    outliers = pd.DataFrame({'outx':outx, 'outy':outy})
    return ColumnDataSource(full_data_frame.reset_index()), ColumnDataSource(outliers)

def structure_boxplot(source, source_outliers):

    """ Takes datasource returns a boxplot for expression in brain regions """

    #structures = list(df.structure_name.unique())
    structures = source.data['structure_name']


    p = figure(tools="save",
               background_fill_color="#fafafa",
               title="",
               x_range=structures,
               plot_height=700,
               plot_width=1000)

    # stems
    p.segment(structures, 'upper', structures, 'q3', line_width=2, line_color='#3B8686', source=source)
    p.segment(structures, 'lower', structures, 'q1', line_width=2, line_color='#3B8686', source=source)

    # boxes
    p.rect(structures, 'top_box_1', 0.7, 'top_box_2',
        fill_color='#BEB790', line_width=2, line_color='#3B8686', source=source)
    p.rect(structures, 'bottom_box_1', 0.7, 'bottom_box_2',
        fill_color='#BEB790', line_width=2, line_color='#3B8686', source=source)

    # whiskers (almost-0 height rects simpler than segments)
    p.rect(structures, 'lower', 0.2, 0.01, line_color='#3B8686', source=source)
    p.rect(structures, 'upper', 0.2, 0.01, line_color='#3B8686', source=source)

    # outliers
    p.circle('outx', 'outy', size=6, color=RdPu9[0], alpha=0.3, source=source_outliers)

    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = "white"
    p.grid.grid_line_width = 2
    p.xaxis.major_label_text_font_size="12pt"
    p.xaxis.major_label_orientation = pi/4

    return p

def expression_by_age(source_circle, source_line):
    tools="save"
    factors = list(source_circle.to_df()['age'].unique())
    p = figure(title="Gene Expression by Age",
               background_fill_color="#fafafa",
               x_range=factors,
               plot_height=450,
               plot_width=1000,
               tools=tools)
    p.circle('age','rpkm', fill_alpha=1, color='#BEB790', size=4, source=source_circle)
    p.xaxis.major_label_orientation = pi/4
    p.line('age', 'rpkm_smooth', line_color='#3B8686', line_width=2, source=source_line)
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

def get_heatmap_data(genes):
    query = {"gene": {"$in": genes}}
    cursor = db.brainspan2.find(query)
    gene_hm_df = pd.DataFrame(list(cursor))

    heatmap_source = pd.DataFrame()
    for gene in genes:
        data = gene_hm_df[gene_hm_df['gene'] == gene]
        hm_data = pd.concat([data['gene'], data['structure_name'], data['age'], data['rpkm']], axis=1)

        rpkm = hm_data.rpkm.reshape(-1,1) #returns a numpy array
        min_max_scaler = preprocessing.MinMaxScaler()
        rpkm_scaled = min_max_scaler.fit_transform(rpkm)

        hm_data['rpkm_norm'] = rpkm_scaled
        heatmap_source = pd.concat([heatmap_source, hm_data])

    genes_hm_df = heatmap_source.groupby(['gene', 'age']).mean().reset_index()
    pivot_data_norm = genes_hm_df.pivot(index='gene', columns='age', values='rpkm_norm')
    pivot_data = genes_hm_df.pivot(index='gene', columns='age', values='rpkm')


    # Set up the data for plotting. We will need to have values for every
    # pair of age/gene names. Map the rate to a color and size.
    ages = list(gene_hm_df['age'].unique())
    colors = RdPu9
    col_rev = []
    for i in reversed(colors):
        col_rev.append(i)

    gene = []
    age = []
    color = []
    rpkm_norm = []
    rpkm = []
    size = []
    for a in ages:
        for g in genes:
            gene.append(g)
            age.append(a)
            rpkm_ = pivot_data_norm[a][g]
            size.append(rpkm_*30)
            rpkm.append(pivot_data[a][g])
            rpkm_norm.append(rpkm_)
            color.append(col_rev[min(int(rpkm_*10), 8)])

    source = ColumnDataSource(
        data=dict(gene=gene, age=age, color=color, rpkm=rpkm, rpkm_norm=rpkm_norm, size=size)
    )
    return source

def plot_heatmap(source, genes):
    # Generate the plot
    TOOLS = "hover,save"
    factors = list(source.to_df()['age'].unique())

    p = figure(title="Gene Expression by Age",
               x_range=factors, y_range=genes,
               x_axis_location="above",
               plot_width=1000,
               plot_height=400,
               tools=TOOLS,
               background_fill_color="#fafafa")

    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "7pt"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = pi/3

    p.circle("age", "gene", 1, 1, source=source,
           color="color", line_color=None, size="size")

    p.select_one(HoverTool).tooltips = [
        ('gene', '@gene'),
        ('age', '@age'),
        ('rpkm', '@rpkm'),
    ]
    return p

def update_plot(attrname, old, new):
    gene = str(select.value)
    structures = multi_select.value
    df1, df2 = get_dataframes(gene, structures)
    src1 = ColumnDataSource(df1)
    src2 = ColumnDataSource(df2)
    src3, src4 = get_boxplot_data(df1)

    source1.data.update(src1.data)
    source2.data.update(src2.data)
    source3.data.update(src3.data)
    source4.data.update(src4.data)

client = MongoClient("mongodb://localhost:27017")
db = client.brainspan

genes1 = ['WDR4', 'RPN1', 'EPHA5', 'WDR7', 'ROBO3', 'MARCH7', 'ZNF665', 'CSMD2', 'FYN', 'SUV420H1', 'ETFB', 'KYNU']
genes = ["ARHGAP19", "FMN2", "UBE2F", "WDR4", "ZNF714", "ABCA3", "PRRC2B", "ZNF665", "ABR", "ADAM18", "ADAMDEC1", "AGRN", "AHNAK", "AHNAK2", "AIRE", "ALS2CR11", "AMBRA1", "ANKRD26", "AP5Z1", "APOPT1", "ARFGAP2", "ARHGEF18", "ARSD", "ASB17", "ASH2L", "ASMTL", "ASTE1", "ATF6B", "ATXN7L1", "AUTS2", "AXIN1", "BCLAF1", "BMP4", "BTBD11", "BTN2A1", "C11orf65", "C1orf177", "C1orf86", "C20orf26", "C2CD3", "C3orf30", "C5orf60", "C8orf33", "CAD", "CASC10", "CCDC158", "CCDC33", "CCDC85B", "CCNE2", "CCSER1", "CD97", "CDK15", "CEP68", "CGRRF1", "CHRD", "CHRM4", "CNBD1", "COL18A1", "COL6A6", "COX6C", "CPE", "CRELD2", "CRYBA2", "CSMD2", "CUBN", "DDX27", "DECR2", "DHRSX", "DLK1", "DNAH17", "DNMT3L", "DOCK4", "DSC3", "DSG2", "DSPP", "DYRK3", "DYRK4", "EHBP1L1", "EPG5", "EPHA5", "EPO", "EPRS", "ERMN", "ERRFI1", "ETFDH", "EZH2", "FAM208A", "FAM83H", "FBXO42", "FERMT2", "FHOD3", "FOXJ1", "FRAS1", "FREM3", "GAL3ST3", "GANAB", "GART", "GBA", "GCFC2", "GOLT1A", "GOT1L1", "GPR113", "GPR151", "GPRIN1", "GRM5", "GUCY1A3", "HCN4", "HEPHL1", "HIST1H3A", "IL4I1", "KIAA0586", "KIAA0825", "KIF22", "KLHDC2", "KRTAP10-11", "KRTAP10-7", "LAMA3", "LAMC3", "LAX1", "LETM2", "LINGO1", "LRRC26", "LRRC32", "LSG1", "LY75", "MAPK15", "MARCH7", "MBLAC1", "MCCC1", "MCF2L", "MFSD4", "MICAL3", "MMP1", "MRGPRX1", "MRPL28", "MS4A15", "MSLN", "MUC4", "MYBPHL", "MYO1C", "MYO7A", "NAALAD2", "NDUFA11", "NEO1", "NFATC4", "NID1", "NLRC4", "NOTUM", "NPAP1", "NPHP4", "NPHS1", "NUTM1", "OGDH", "OLAH", "OR10G7", "OR11H6", "OR4D9", "OR8D4", "OR9A4", "OR9I1", "OSBPL8", "OXTR", "PARM1", "PAX9", "PCSK5", "PDZD2", "PELP1", "PERM1", "PEX5L", "PIDD1", "PIK3C2B", "PITRM1", "PKHD1", "PKP3", "PLA2G4C", "PLEKHG4B", "POM121L12", "PPP1R15B", "PPP1R35", "PRR25", "PRR27", "PRR35", "PRSS23", "PRSS27", "PSMD14", "RAB11FIP1", "RFXAP", "RGS21", "RHBDF1", "RIMBP2", "RIN3", "ROBO3", "RPN1", "RPRD2", "RPS6KC1", "SCN9A", "SERPINA4", "SEZ6", "SFTPC", "SGPP2", "SHCBP1L", "SLC14A2", "SLC26A1", "SLC26A9", "SLC2A6", "SLC39A5", "SLC4A4", "SLIT3", "SMAD3", "SMG1", "SNX16", "SOX30", "SPATA18", "SPATS2", "SPIRE2", "SPOCK3", "SPTB", "SRD5A1", "SSH3", "SSPO", "SSRP1", "STBD1", "SYNE3", "SYNGR1", "SYT5", "TAGLN", "TAS1R2", "TCF7", "TECTA", "THAP3", "THOC1", "THPO", "TIGD4", "TLN1", "TMEM14E", "TMEM198", "TMEM249", "TMEM8B", "TNXB", "TRAF5", "TRAPPC10", "TTC16", "TTC21B", "TTN", "TUBGCP6", "TWISTNB", "UFSP1", "UGT2B15", "ULK1", "UMODL1", "UNC79", "USP35", "VWA3A", "VWDE", "WDR27", "WDR7", "WDR90", "WHSC1L1", "WRAP73", "XDH", "ZASP", "ZBTB41", "ZC3H7A", "ZDBF2", "ZFYVE9", "ZNF214", "ZNF423", "ZNF438", "ZNF7", "ZNF765", "ZNFX1", "ZSCAN31"]
gene = genes[2]
structures = [u'occipital neocortex',
             u'primary motor-sensory cortex (samples)',
             u'amygdaloid complex',
             u'medial ganglionic eminence',
             u'posterior (caudal) superior temporal cortex (area 22c)',
             u'upper (rostral) rhombic lip',
             u'caudal ganglionic eminence',
             u'dorsal thalamus',
             u'anterior (rostral) cingulate (medial prefrontal) cortex',
             u'dorsolateral prefrontal cortex',
             u'orbital frontal cortex',
             u'lateral ganglionic eminence',
             u'inferolateral temporal cortex (area TEv, area 20)',
             u'hippocampus (hippocampal formation)',
             u'ventrolateral prefrontal cortex',
             u'parietal neocortex',
             u'temporal neocortex',
             u'primary auditory cortex (core)',
             u'primary visual cortex (striate cortex, area V1/17)',
             u'striatum',
             u'primary motor cortex (area M1, area 4)',
             u'posteroventral (inferior) parietal cortex',
             u'primary somatosensory cortex (area S1, areas 3,1,2)',
             u'cerebellum',
             u'cerebellar cortex',
             u'mediodorsal nucleus of thalamus']

df1, df2 = get_dataframes(gene, structures)
source1 = ColumnDataSource(df1)
source2 = ColumnDataSource(df2)
source3, source4 = get_boxplot_data(df1)

age_plot = expression_by_age(source1, source2)
structure_plot = structure_boxplot(source3, source4)

source5 = get_heatmap_data(genes1)
heatmap = plot_heatmap(source5, genes1)

plot = vplot(structure_plot,age_plot, heatmap)

multi_select = MultiSelect(title="Brain Regions:", value=structures,
                           options=structures)
multi_select.on_change('value', update_plot)

select = Select(title="Gene:", value=genes[2], options=genes)
select.on_change('value', update_plot)

curdoc().add_root(vplot(select, plot, multi_select))
