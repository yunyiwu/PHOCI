#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pyBigWig
import numpy as np

from config import config_test, test_cell_line


# In[ ]:


chroms = {"chr1":49788, 
          "chr2":48431, 
          "chr3":39633, 
          "chr4":38035, 
          "chr5":36271, 
          "chr6":34148, 
          "chr7":31867, 
          "chrX":31205, 
          "chr8":29015, 
          "chr9":27677,
          "chr11":27015,
          "chr10":26752, 
          "chr12":26652, 
          "chr13":22870, 
          "chr14":21376, 
          "chr15":20395, 
          "chr16":18045, 
          "chr17":16646, 
          "chr18":16052, 
          "chr20":12866, 
          "chr19":11721,
          "chr22":10160, 
          "chr21":9339}


# In[ ]:


def bigwig_process():

    bigwigs = config_test["bigwigs"]

    resolution = 5000

    for chrom in chroms.keys():
        signals = []
        for bigwig in bigwigs:
            bw = pyBigWig.open(config_test["bigwig_dir"]+bigwig+".bigWig")
        
            try:
                signal = bw.stats(chrom, 0, (chroms[chrom]+1)*resolution, type="mean", nBins=(chroms[chrom]+1))
                signals.append(signal)
            except:
                print("error in bigwig files:")
                print(chrom, bigwig)
                continue
        
            bw.close()

        signals = np.array(signals)

        np.save(config_test["feature_dir_path"]+test_cell_line+"_"+chrom+"_x", signals)


# In[ ]:





# In[ ]:




