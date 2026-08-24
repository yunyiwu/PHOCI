#!/usr/bin/env python
# coding: utf-8

import pickle
from config import config_train
from hic_process import hic_process
from bigwig_process import bigwig_process
from porec_process import porec_process
from data_process_for_training import data_process_all_chrom


def run_data_process_for_training():
    # Step 1: Process hic data
    hic_process()

    # Step 2: Process bigwig data
    bigwig_process()

    # Step 3: Process porec data
    fc_ids = config_train.get("fc_ids", [])
    for fc_id in fc_ids:
        try:
            porec_process(file_id=fc_id)
        except Exception as e:
            print(f"Error processing Pore-C file {fc_id}: {e}")

    # Step 4: Process data for training
    size = 1000  # Define the size for sliding window processing
    chr_index = {}

    for ch in config_train["chr_names"]:
        config_train["chr_names"] = [ch]

        try:
            sliding_data, sliding_index = data_process_all_chrom(config_train, size)
        except Exception as e:
            print(f"Error processing chromosome {ch}: {e}")
            continue

        # Save sliding data
        sliding_data_path = f"{config_train['input_graph_dir_path']}{ch}_sliding_data"
        with open(sliding_data_path, "wb") as f:
            pickle.dump(sliding_data, f)
        
        # Collect index data
        chr_index[ch] = sliding_index

    # Save the chromosome index data
    sliding_index_path = f"{config_train['input_graph_dir_path']}sliding_index"
    with open(sliding_index_path, "wb") as f:
        pickle.dump(chr_index, f)


def main():
    run_data_process_for_training()


if __name__ == "__main__":
    main()
