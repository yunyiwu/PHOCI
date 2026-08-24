#!/usr/bin/env python
# coding: utf-8

import os
import pickle
from config import config_train, train_cell_line
from hic_process import hic_process
from bigwig_process import bigwig_process
from porec_process import porec_process
from data_process_for_training import data_process_all_chrom


def run_data_process_for_training():
    # Step 1: Process Hi-C data
    print("Step 1: Processing Hi-C data...")
    hic_process()

    # Step 2: Process BigWig feature data
    print("Step 2: Processing BigWig features...")
    bigwig_process()

    # Step 3: Process Pore-C hypergraph data for each FC ID in config_train
    print("Step 3: Processing Pore-C hypergraph data...")
    fc_ids = config_train.get("fc_ids", [])
    for fc_id in fc_ids:
        print(f"  -> Processing Pore-C for file_id: {fc_id}")
        try:
            porec_process(file_id=fc_id)
        except Exception as e:
            print(f"Error processing Pore-C file {fc_id}: {e}")

    # Step 4: Process sliding graph data for training
    print("Step 4: Generating sliding window graph datasets for training...")
    size = 1000  # Sliding window size

    # Ensure output graph directory exists
    input_graph_dir = config_train.get(
        "input_graph_dir_path", 
        f"/public/home/wuyy/{train_cell_line}_hg38/input_graph/"
    )
    os.makedirs(input_graph_dir, exist_ok=True)

    # Process each chromosome listed in training set
    for ch in config_train["chr_names"]:
        print(f"  -> Processing training data for chromosome: {ch}")
        # Temporarily isolate single chromosome for processing
        config_train_single = config_train.copy()
        config_train_single["chr_names"] = [ch]

        try:
            sliding_data = data_process_all_chrom(config_train_single, size)
        except Exception as e:
            print(f"Error processing training dataset for chromosome {ch}: {e}")
            continue

        # Save sliding data for current chromosome
        sliding_data_path = os.path.join(input_graph_dir, f"{ch}_train_sliding_data")
        with open(sliding_data_path, "wb") as f:
            pickle.dump(sliding_data, f)
            
        print(f"     Saved {len(sliding_data)} samples to {sliding_data_path}")

    print("All training data preprocessing completed successfully!")


def main():
    run_data_process_for_training()


if __name__ == "__main__":
    main()
