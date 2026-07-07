#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from config import config_test
from hic_process import hic_process
from bigwig_process import bigwig_process
from data_process_for_prediction import data_process_all_chrom
import pickle

def run_data_process_for_prediction():
    # Step 1: Process hic data
    hic_process()

    # Step 2: Process bigwig data
    bigwig_process()

    # Step 3: Process data for prediction
    size = 1000  # Define the size for sliding window processing
    chr_index = {}

    for ch in config_test["test_chr_names"]:
        config_test["chr_names"] = [ch]

        try:
            sliding_data, sliding_index = data_process_all_chrom(config_test, size)
        except Exception as e:
            print(f"Error processing chromosome {ch}: {e}")
            continue

        # Save sliding data
        sliding_data_path = f"{config_test['input_graph_dir_path']}{ch}_sliding_data"
        with open(sliding_data_path, "wb") as f:
            pickle.dump(sliding_data, f)
        
        # Collect index data
        chr_index[ch] = sliding_index

    # Save the chromosome index data
    sliding_index_path = f"{config_test['input_graph_dir_path']}sliding_index"
    with open(sliding_index_path, "wb") as f:
        pickle.dump(chr_index, f)

def main():
    run_data_process_for_prediction()

if __name__ == "__main__":
    main()

