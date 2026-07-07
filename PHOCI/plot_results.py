#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from matplotlib.axes import Axes

from trackc.pl.bigwig import _make_multi_region_ax
from trackc.tl._getRegionsCmat import GenomeRegion

import matplotlib.colors as mcolors
from itertools import combinations
import itertools
from matplotlib.patches import Ellipse
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.patches import Shadow
from matplotlib import patches

from datetime import datetime

import pyBigWig

import numpy as np


# In[2]:


def gene_track(
    ax = None,
    bed12 = None,
    regions = None,
    track_type = "gene",
    show_label = True,
    pos_strand_gene_color = "#3366CC",
    neg_strand_gene_color = "#EECFA1",
    line = 1,
    gene_fontszie = 7,
    label = None,
    label_rotation = 0,
    label_fontsize = 16,
    ax_on = False,
):

    if isinstance(regions, list):
        line_GenomeRegions = pd.concat(
            [GenomeRegion(i).GenomeRegion2df() for i in regions]
        )
    else:
        line_GenomeRegions = GenomeRegion(regions).GenomeRegion2df()

    axs = _make_multi_region_ax(ax, line_GenomeRegions)
    line_GenomeRegions = line_GenomeRegions.reset_index()

    ax.set_ylabel(
        label,
        fontsize=label_fontsize,
        rotation=label_rotation,
        horizontalalignment="right",
        verticalalignment="center",
    )

    ax.set_xticks([])
    ax.set_xticklabels("")
    ax.set_yticks([])
    ax.set_yticklabels("")
    if ax_on == False:
        spines = ["top", "bottom", "left", "right"]
        for i in spines:
            ax.spines[i].set_color("none")

    if isinstance(bed12, str):
        bed12 = pd.read_table(bed12, sep="\t", header=None)

    bed12 = bed12.iloc[:, 0:12]
    bed12.columns = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "itemRgb",
        "blockCount",
        "blockSizes",
        "blockStarts",
    ]
    bed12["blockSizes"] = bed12["blockSizes"].str.rstrip(",")
    bed12["blockStarts"] = bed12["blockStarts"].str.rstrip(",")
    bed12["chrom"] = bed12["chrom"].astype(str)
    chrom_names = bed12["chrom"].unique()

    for ix, row in line_GenomeRegions.iterrows():
        if row["chrom"] not in chrom_names:
            raw_chr = row["chrom"]
            if row["chrom"].startswith("chr"):
                row["chrom"] = row["chrom"].lstrip("chr")
            else:
                row["chrom"] = "chr" + row["chrom"]
            if row["chrom"] not in chrom_names:
                print(f"{raw_chr} not in bigwig chroms!")
                return

        if track_type == "gene":
            _plot_gene(
                axs[ix],
                bed12,
                row["chrom"],
                row["fetch_start"],
                row["fetch_end"],
                needReverse=row["isReverse"],
                show_label=show_label,
                pos_strand_gene_color=pos_strand_gene_color,
                neg_strand_gene_color=neg_strand_gene_color,
                line=line,
                fontszie=gene_fontszie,
                ax_on=ax_on,
            )
        if track_type == "density":
            print("This gene type is developping")
        else:
            pass


def _plot_gene(
    ax,
    gene_bed,
    chrom,
    start,
    end,
    needReverse=False,
    show_label=True,
    pos_strand_gene_color="#3366CC",
    neg_strand_gene_color="#EECFA1",
    line=1,
    fontszie=5,
    ax_on=False,
):
    gene_bed = gene_bed[gene_bed["chrom"] == chrom]
    gene_bed_plot = gene_bed[
        ((gene_bed["start"] >= start) & (gene_bed["start"] <= end))
        | ((gene_bed["end"] >= start) & (gene_bed["end"] <= end))
    ]
    gene_bed_plot = gene_bed_plot.sort_values(by="end")
    # print(gene_bed_plot

    plot_gene_num = gene_bed_plot.shape[0]

    ii = 0
    head_length = (abs(end - start) / (line + 2)) / 5
    if line <= 3:
        head_length = (abs(end - start) / (line * 3)) / 10

    for i, row in gene_bed_plot.iterrows():
        # col = pos_strand_gene_color
        text_col = pos_strand_gene_color

        if row["strand"] == "-":
            # col = neg_strand_gene_color
            text_col = neg_strand_gene_color

        # text_col = col
        
        if row["strand"] == "+":
            plot_y = 0
            
        elif row["strand"] == "-":
            plot_y = 1
        
        #plot_y = ii % line

        ax.plot(
            (row["start"], row["end"]),
            (plot_y + 0.5, plot_y + 0.5),
            color="k",
            linewidth=1,
            solid_capstyle="butt",
        )
        starts = [int(x) for x in row["blockStarts"].split(",")]
        widths = [int(x) for x in row["blockSizes"].split(",")]

        ax.bar(
            x=starts,
            height=0.4,
            width=widths,
            bottom=plot_y + 0.3,
            edgecolor="k",
            linewidth=1,
            align="edge",
            color="k",
        )

        if row["start"] < start:
            row["start"] = start
        if row["end"] > end:
            row["end"] = end

        arrow_s = row["end"]
        dx = 0.3
        if row["strand"] == "-":
            arrow_s = row["start"]
            dx = -0.1
        ax.arrow(
            arrow_s,
            plot_y + 0.5,
            dx,
            0,
            overhang=0.5,
            width=0,
            head_width=0.28,
            head_length=head_length,
            length_includes_head=False,
            color=text_col,
            linewidth=0.5,
        )
        if isinstance(show_label, bool):
            if show_label == False:
                ii += 1
                continue
        if isinstance(show_label, str):
            if row["name"] != show_label:
                ii += 1
                continue
        if isinstance(show_label, list):
            if row["name"] not in show_label:
                ii += 1
                continue

        if (row["name"] in gene_bed_plot.iloc[-4:, :]["name"]) or (
            int(line / (ii + 1)) < 2
        ):
            ha = "center"
            genename = row["name"] + "  "
            xpos = (row["start"]+row["end"])/2
            if needReverse:
                ha = "center"
                genename = "  " + row["name"]
                # xpos = row['end']
            ypos = plot_y #+ 0.5
            if line == 1:
                xpos = row["start"] + abs(row["start"] - row["end"]) / 2
                ypos = 0.8
                ha = "center"
                
            if row["strand"] == "+":
                ax.text(
                xpos,
                ypos-0.2,
                genename + "  ",
                ha=ha,
                va="top",#"center",
                color=text_col,
                fontsize=fontszie,)

            elif row["strand"] == "-":
                ax.text(
                xpos,
                ypos+1.2,
                genename + "  ",
                ha=ha,
                va="bottom",#"center",
                color=text_col,
                fontsize=fontszie,)


        else:
            ha = "center"
            genename = "  " + row["name"]
            xpos = (row["start"]+row["end"])/2
            if needReverse:
                ha = "center"
                genename = row["name"] + "  "
                # xpos = row['start']

            ypos = plot_y #+ 0.5
            if line == 1:
                xpos = row["start"] + abs(row["start"] - row["end"]) / 2
                ypos = 0.8
                ha = "center"
                
            if row["strand"] == "+":
                ax.text(
                xpos,
                ypos-0.2,
                genename + "  ",
                ha=ha,
                va="top",#"center",
                color=text_col,
                fontsize=fontszie,)

            elif row["strand"] == "-":
                ax.text(
                xpos,
                ypos+1.2,
                genename + "  ",
                ha=ha,
                va="bottom",#"center",
                color=text_col,
                fontsize=fontszie,)

        ii += 1

    xlim_s = start
    xlim_e = end
    if needReverse == True:
        xlim_s = end
        xlim_e = start

    ax.set_xlim(xlim_s, xlim_e)
    ax.set_ylim(top=0, bottom=line)
    if plot_gene_num < line:
        ax.spines["bottom"].set_position(("data", plot_gene_num))

    if ax_on == False:
        spines = ["top", "bottom", "left", "right"]
        for i in spines:
            ax.spines[i].set_visible(False)
            # for i in ['left','top','right']:
            ax.spines[i].set_color("none")
            ax.spines[i].set_linewidth(0)
    ax.spines["bottom"].set_color("black")
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(bottom=True, top=False, left=False, right=False)
    ax.set_xticklabels("")
    ax.set_yticklabels("")


# In[ ]:


def pcolormesh_45deg(ax, matrix_c, matrix_c2, start_pos_vector, start, end):
    
    n = matrix_c.shape[0]
    # create rotation/scaling matrix
    t = np.array([[1, 0.5], [-1, 0.5]])
    # create coordinate matrix and transform it
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                            start_pos_vector)]), t)
    
    # this is to convert the indices into bp ranges
    x = matrix_a[:, 1].reshape(n , n )
    y = matrix_a[:, 0].reshape(n , n )
    
    cmap = plt.cm.Reds#coolwarm
    norm = mcolors.TwoSlopeNorm(vmin=0, vcenter=1.5, vmax=3) #vcenter=0.7,

    #plt.imshow(hic_data[:,:], cmap=cmap,norm=norm,  interpolation='nearest')
    
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), cmap=cmap, norm=norm)
    #im = ax.pcolormesh(x, y, np.flipud(matrix_c2), cmap=plt.cm.Greens, norm=norm)
    
    threshold = 0  # define your threshold
    for i in range(n):
        for j in range(n):
            if np.flipud(matrix_c2)[i, j] >= threshold:  # condition to highlight
                rect = Ellipse((x[i, j], y[i, j]), 1, 1.5, fill=False, edgecolor='yellow', linewidth=1)
                                
                ax.add_patch(rect)
    
    ax.set_xticks([])
    #ax.set_xticklabels([start, end], fontsize=25)
    #ax.tick_params(axis='x',which='major',direction='out',bottom=True,length=6, width = 1)
    ax.set_xlim(0,n-1)
    ax.set_ylim(0, n)
    
    ax.xaxis.set_ticks_position('bottom') 
    #ax.tick_params(length=6, width=2, direction='inout')
    
    # Hide y-axis
    # Hide y-axis and spines
    ax.get_yaxis().set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)   
    
    return im


def plot_hic_tri(matrix_data, matrix_data2, start, end, axs, cax):
    
    matrix_size = matrix_data.shape[0]
    
    triangle_matrix = np.triu(matrix_data)

    # 创建一个掩码，将下半部分设为 False
    mask = np.tri(matrix_size, k=-2)
    masked_matrix = np.ma.array(triangle_matrix, mask=mask)

    triangle_matrix2 = np.triu(matrix_data2)
    masked_matrix2 = np.ma.array(triangle_matrix2, mask=mask)
    
    
    # 调用 pcolormesh_45deg() 方法绘制旋转后的矩阵
    im = pcolormesh_45deg(axs, masked_matrix, masked_matrix2, np.arange(matrix_size), start, end)

    cbar = plt.colorbar(im, cax=cax, anchor = (-15, 0.5))


# In[ ]:


def load_bigwig_stats(file_path, chrom, start, end, bins=1000, stats_type="mean"):
    """加载 BigWig 数据并计算统计值"""
    with pyBigWig.open(file_path) as bw:
        return bw.stats(chrom, start, end, type=stats_type, nBins=bins)

def plot_track(ax, data, label, color, ylabel, xlim=(0, 1000)):
    """通用轨迹绘图"""
    ax.plot(data, label=label, linewidth=1, color=color)
    ax.set_xlim(*xlim)
    ax.set_ylim(ymin=0)
    ax.set_ylabel(ylabel, rotation=0, ha='right', va='center')
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.xaxis.set_visible(False)


# In[ ]:


def plot_apriori(target, graph_data, rs, es, start, end, start_bin, end_bin, count, file_paths, chrom, gene_index, gene_name, resample_support, resample_confidence, num_samples):
    
    atac = load_bigwig_stats(file_paths['atac'], chrom, start, end)
    h3k27ac = load_bigwig_stats(file_paths['h3k27ac'], chrom, start, end)
    h3k4me1 = load_bigwig_stats(file_paths['h3k4me1'], chrom, start, end)
    ctcf = load_bigwig_stats(file_paths['ctcf'], chrom, start, end)
    rad21 = load_bigwig_stats(file_paths['rad21'], chrom, start, end)

    bb = pyBigWig.open(file_paths['chrom_state'])
    c18 = bb.entries(chrom, start,end+1)

    bb = pyBigWig.open("annotation/encodeCcreCombined_hg38.bb")
    ccre = bb.entries(chrom, start,end+1)

    size = end_bin - start_bin +1

    bins = list(range(start_bin,end_bin+1))

    f = np.full((size, size), np.nan)

    # plot hic
    for i in range(0, graph_data.edge_index.shape[1]):
        k0 =  graph_data.edge_index[0,i]
        k1 =  graph_data.edge_index[1,i]

        if k0 > start_bin-1 and k0 < end_bin+1: 
            if k1 > start_bin-1 and k1 < end_bin+1:
                f[k0-start_bin, k1-start_bin] =  graph_data.edge_attr[i]


    f2 = np.full((size, size), np.nan)
    for rss in rs:
        for i in combinations(rss,2):
            f2[i[0]-start_bin,i[1]-start_bin] = f[i[0]-start_bin,i[1]-start_bin]
            f2[i[1]-start_bin,i[0]-start_bin] = f[i[1]-start_bin,i[0]-start_bin]            


    sorted_keys = sorted(count.keys(), key=lambda k: count[k], reverse=True)

    sorted_keys2 = sorted(count.keys(), key=lambda k: count[k], reverse=False)

    sorted_rs = [list(i) for i in sorted_keys]

    sorted_count = [count[i] for i in sorted_keys2]

    sn = start
    sstr = chrom+":"+"{:,}".format(sn)
    en = end
    estr = "{:,}".format(en)

    # 示例数据
    c18 = c18
    track1_positions = es
    data = sorted_rs

    barh_data = sorted_count

    # 创建图形和自定义布局
    fig = plt.figure(figsize=(14, 14))
    gs = gridspec.GridSpec(11, 2, width_ratios=[3, 1], height_ratios=[3.5, 0.6, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.5, 2.5])
    gs.update(hspace=0.30, wspace=0.15)  # 更新子图间距

    # 生成颜色映射
    cmap = plt.cm.Set2
    colors = []#[cmap(i / len(track1_positions)) for i in range(len(track1_positions))]
    for i in track1_positions:
        if i == target:
            colors.append("#d9333f")
        else:
            colors.append("#008899")


    ax0000 = plt.subplot(gs[0,0])
    inner_grid = gridspec.GridSpecFromSubplotSpec(1, 6, subplot_spec=gs[0, 1], wspace=0.0, hspace=0.0)
    cax = fig.add_subplot(inner_grid[0])
    plot_hic_tri(f, f2, sstr, estr, ax0000, cax)

    ax0001 = plt.subplot(gs[1,0])
    gene_track(ax=ax0001, bed12 = 'annotation/hg38.bed12', regions=[chrom+':'+str(sn)+'-'+str(en)], line=2)

    ax_ctcf = plt.subplot(gs[2, 0])
    plot_track(ax_ctcf, ctcf, label='CTCF', color="#674196", ylabel='CTCF')

    ax_rad21 = plt.subplot(gs[3, 0])
    plot_track(ax_rad21, rad21, label='RAD21', color="#e49e61", ylabel='RAD21')

    ax_atac = plt.subplot(gs[4, 0])
    plot_track(ax_atac, atac, label='ATAC', color="#c7b370", ylabel='ATAC')

    ax_h3k4me1 = plt.subplot(gs[5, 0])
    plot_track(ax_h3k4me1, h3k4me1, label='H3K4me1', color="#839b5c", ylabel='H3K4me1')

    ax_h3k27ac = plt.subplot(gs[6, 0])
    plot_track(ax_h3k27ac, h3k27ac, label='H3K27ac', color="#165e83", ylabel='H3K27ac')

    
    # 设定 colormaps
    cmap = plt.cm.get_cmap('tab20')
    cmap2 = plt.cm.get_cmap('tab20b')
    cmap3 = plt.cm.get_cmap('tab20c')

    # 颜色字典
    label_colors = {
        'Tx': cmap.colors[18], 
        'TxWk': cmap.colors[19], 
        'ReprPC': cmap.colors[2], 
        'ReprPCWk': cmap.colors[3], 
        'Quies': '#fbfaf5',
        'Het': cmap.colors[10],
        'ZNF/Rpts': cmap.colors[12],
        'EnhA1': cmap2.colors[0],
        'EnhA2': cmap2.colors[1],
        'EnhG1': cmap2.colors[2],
        'EnhG2': cmap2.colors[3],
        'EnhBiv': cmap3.colors[14],
        'EnhWk': cmap3.colors[15],
        'TssA': cmap2.colors[4],
        'TssBiv': cmap3.colors[8],
        'TssFlnk': cmap3.colors[9],
        'TssFlnkD': cmap3.colors[10],
        'TssFlnkU': cmap3.colors[11],
    }


    ax01 = plt.subplot(gs[7, 0])

    # 遍历数据并绘制方块
    for s, e, label in c18:
        label_key = label.split('\t')[0]  # 提取label的前半部分
        color = label_colors.get(label_key, 'black')  # 获取颜色，如果不存在则使用黑色
        rect = patches.Rectangle((s, 0), e - s, 1, color=color)
        ax01.add_patch(rect)

    # 设置轴的范围
    ax01.set_xlim(sn, en)
    ax01.set_ylim(0, 1)
    ax01.set_xticks([])
    ax01.set_yticks([])
    ax01.set_ylabel('ChromHMM\n 18-state',rotation=0, ha='right', va='center')
    ax01.spines['top'].set_visible(False)
    ax01.spines['right'].set_visible(False)
    ax01.spines['left'].set_visible(False)
    ax01.spines['bottom'].set_visible(False)


    ax10 = plt.subplot(gs[8, 0])

    for s, e, label in ccre:
        label_key = label.split('\t')[0]  # 提取label的前半部分

        a = label.split('\t')[7]

        if a == 'PLS':
            rect = patches.Rectangle((s, 0), e - s, 1, color='#c83c23')
            ax10.add_patch(rect)

    # 设置轴的范围
    ax10.set_xlim(sn, en)
    ax10.set_ylim(0, 1)
    ax10.set_xticks([])
    ax10.set_yticks([])
    ax10.set_ylabel('Promoter-Like\nSignature',rotation=0, ha='right', va='center')
    ax10.spines['top'].set_visible(False)
    ax10.spines['right'].set_visible(False)
    ax10.spines['left'].set_visible(False)
    ax10.spines['bottom'].set_visible(False)


    # Track1: 方块图
    ax1 = plt.subplot(gs[9, 0])

    ax1.add_patch(plt.Rectangle((start_bin, 0), end_bin, 1, facecolor='#c1e4e9', linewidth=1))
    ax1.add_patch(plt.Rectangle((gene_index[0], 0), gene_index[-1]-gene_index[0], 1, facecolor='#f6bfbc', linewidth=1))

    for i, position in enumerate(track1_positions):
        ax1.add_patch(plt.Rectangle((position, 0), 1, 1, facecolor=colors[i], edgecolor='black', linewidth=1))
        ax1.text(position + 0.5, -0.15, str(position), ha='center', va='top', fontsize=5)  # 在方块下方显示数字

    ax1.set_xlim(start_bin, end_bin+1)
    ax1.set_ylim(-0.5, 1.1)  # 确保显示完整
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.set_yticks([])
    ax1.set_xticks([start_bin, end_bin+1],[sstr, estr])  # 只保留首尾两个刻度
    #ax1.set_xticklabels()  # 设置刻度标签


    # 中间的方块集合图
    ax2 = plt.subplot(gs[10, 0])
    gap = 0.2  # 方块之间的间隙
    block_height = 1 - gap  # 方块的高度

    for i, sublist in enumerate(data):
        for j, pos in enumerate(sublist):
            # 调整方块位置，使其中心与y轴刻度对齐
            color_index = track1_positions.index(pos) if pos in track1_positions else 0
            ax2.add_patch(plt.Rectangle((pos, len(data) - i - 1 - 0.5 * block_height), 1, block_height, 
                                        facecolor=colors[color_index], edgecolor='black', linewidth=1))
    ax2.set_xlim(start_bin, end_bin+1)
    ax2.set_ylim(-0.5, len(data) - 0.5)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    #ax2.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)
    ax2.set_xticks([])
    ax2.set_xlabel(gene_name+"_"+str(target)+"\nApriori frequent sets (Support > "+str(resample_support)+", Confidence > "+str(resample_confidence)+", Samples = "+str(num_samples)+")", fontsize = 16)
    
    # Barh plot: 频数统计
    ax3 = plt.subplot(gs[10, 1], sharey=ax2)
    bars = ax3.barh(list(range(0,len(barh_data))), barh_data, facecolor='#f8b862')

    for bar in bars:
        shadow = Shadow(bar, ox=0.1, oy=-0.01, color='gray', alpha=0.5)
        ax3.add_patch(shadow)

    ax3.set_ylim(-0.5, len(data) - 0.5)  # 使y轴对齐
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    #ax3.xaxis.set_visible(False)
    ax3.set_xticks([])
    ax3.set_xlabel("Set Counts")
    ax3.yaxis.set_visible(False)

    for i, v in enumerate(barh_data):
        ax3.text(v + 0.5, i, str(v), color='black', verticalalignment='center')

        
    # 获取当前时间
    now = datetime.now()

    # 格式化为字符串
    time_str = now.strftime("%Y-%m-%d_%H-%M-%S")

    # 使用 tight_layout 确保布局正确
    plt.tight_layout()
    plt.savefig("apriori_plots/apriori_"+gene_name+"_"+str(target)+"_"+time_str+".png",dpi=300, bbox_inches='tight')

    plt.show()


# In[ ]:




