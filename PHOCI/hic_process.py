#!/usr/bin/env python
# coding: utf-8

import numpy as np
import hicstraw
import glob
from config import config_test, test_cell_line


def hic_process():

    resolution = 5000
    
    hic_file = config_test["hic_file"]
    hic = hicstraw.HiCFile(hic_file)

    chrms = hic.getChromosomes()

    for chrom in chrms:
        if chrom.name == "All" or chrom.name == 'ALL' or chrom.name == 'all':
            continue
    
        matrix_object = hic.getMatrixZoomData(chrom.name, chrom.name, "observed", "KR", "BP", resolution)
        matrix = matrix_object.getRecords(0, chrom.length, 0, chrom.length)
    
        index = []
        attr = []
    
        for r in matrix:
            x = int(r.binX/resolution)
            y = int(r.binY/resolution)
            count = r.counts
    
            index.append([x,y])
            attr.append(count)
    
        index = np.array(index).T
        attr = np.array(attr)
    
        np.save(config_test["hic_dir_path"]+test_cell_line+"_chr_chr"+chrom.name+"_index_"+str(resolution), index)
        np.save(config_test["hic_dir_path"]+test_cell_line+"_chr_chr"+chrom.name+"_attr_"+str(resolution), attr)

