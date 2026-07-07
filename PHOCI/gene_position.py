#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import numpy as np
import pyBigWig


def get_gene_position(gene_name, cell_line):

    with open("annotation/gene_chr.pkl", "rb") as f:
        gene_chr = pickle.load(f)

    chrom = gene_chr[gene_name]

    gene_annotation = pd.read_pickle('annotation/'+chrom+'_bins_gene_id.pkl')

    gene_annotation['gene_name'] = gene_annotation['gene_name'].str.replace(',', '', 1)
    gene_annotation['gene_name'] = gene_annotation['gene_name'].str.replace(' ', '')

    gene_span = gene_annotation[gene_annotation['gene_name'] == gene_name].bin_id.tolist()

    with open("data/"+cell_line+"_hg38/input_graph/sliding_index","rb") as f:
        annotate = pickle.load(f)

    if chrom not in annotate.keys():
        print("Chr not in annotate file. Change annotate file.")

    anno = annotate[chrom]


    for i, ann in enumerate(anno):
        intersection = np.intersect1d(gene_span, ann)

        if len(intersection) >= 1:
            print(i, gene_span, ann[0], ann[-1])
            ann1 = ann
            i1 = i

            break


    sub_anno = gene_annotation[ann1[0]:ann1[-1]+1]
    sub_anno = sub_anno.reset_index()

    gene_index = []
    pls_index = []

    bb = pyBigWig.open("annotation/encodeCcreCombined_hg38.bb")

    for i in sub_anno.iterrows():
        if  gene_name in i[1].gene_name:
            gene_index.append(i[0])

            ccre = bb.entries(i[1].Chromosome, i[1].Start, i[1].End)

            if ccre == None:
                continue

            for s, e, label in ccre:
                label_key = label.split('\t')[0]  # 提取label的前半部分

                a = label.split('\t')[7]

                if a == 'PLS':
                    pls_index.append(i[0])



    if len(pls_index) == 0:
        print("Not found PLS for "+gene_name)

    bb.close()

    test_points = []
    for i in pls_index:
        test_points.append(i-1)
        test_points.append(i)
        test_points.append(i+1)

    test_points = list(set(test_points))

    
    return i1, chrom, test_points, sub_anno, gene_index



