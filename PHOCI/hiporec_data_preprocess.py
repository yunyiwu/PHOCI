import os
import pickle
import pandas as pd
import pyranges as pr
from cooler.util import binnify
from config import config_train, train_cell_line


def porec_process(file_id: str):
    """
    Process Pore-C alignment file into hypergraphs using configuration from config_train.

    Parameters:
    -----------
    file_id : str
        File identifier (e.g., 'FC2').
    """
    resolution = 5000

    # 1. Load chromosome sizes
    chrom_sizes_file = "hg38_norm.chrom.sizes"
    chrom_dict = pd.read_csv(
        chrom_sizes_file,
        sep="\t",
        header=None,
        names=["chrom", "size"],
        index_col="chrom"
    )["size"]

    # 2. Get input alignment path from config_train
    # Path template: <porec_dir_path>/<cell_line>/<cell_line>_<file_id>_reads_alignment.csv.gz
    alignment_path = os.path.join(
        config_train["porec_dir_path"],
        f"{train_cell_line}_{file_id}_reads_alignment.csv.gz"
    )

    use_cols = ["read_name", "read_length", "read_start", "read_end", "strand", "chrom", "chrom_length", "start", "end", "MapQual"]
    alignment_df = pd.read_csv(alignment_path, compression="gzip", usecols=use_cols)

    # 3. Filter and clean alignment data
    alignment_df = alignment_df[alignment_df["MapQual"] != "MapQual"].copy()
    alignment_df["MapQual"] = alignment_df["MapQual"].astype(float)
    alignment_df = alignment_df[alignment_df["MapQual"] >= 10]

    int_cols = ["read_start", "read_end", "read_length", "start", "end"]
    alignment_df[int_cols] = alignment_df[int_cols].astype(int)

    # 4. Construct and ensure output directory exists based on config_train
    output_dir = os.path.join(
        config_train["porec_dir_path"],
        f"{file_id}_hyper"
    )
    os.makedirs(output_dir, exist_ok=True)

    # 5. Process chromosomes specified in config_train
    chrom_names = config_train.get("chr_names", ["chr10"])

    for chrom_name in chrom_names:
        if chrom_name == "chrY":
            continue

        alignment_df_sub = alignment_df[alignment_df["chrom"] == chrom_name]
        if alignment_df_sub.empty:
            continue

        # Bin the chromosome sequence
        chrom_dict_sub = chrom_dict.loc[[chrom_name]]
        bins_df = binnify(chrom_dict_sub, resolution)
        bins_df.index.name = "bin_id"

        bins_pr = pr.PyRanges(
            bins_df.reset_index().rename(columns={"start": "Start", "end": "End", "chrom": "Chromosome"}),
            int64=True
        )

        # Calculate fragment midpoints
        midpoint_pr = pr.PyRanges(
            alignment_df_sub[["read_name", "chrom", "start", "end"]]
            .assign(start=lambda x: ((x["start"] + x["end"]) * 0.5).round(0).astype(int))
            .eval("end = start + 1")
            .rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"}),
            int64=True
        )

        # Map midpoints to corresponding bin IDs
        midpoint_to_bin = midpoint_pr.join(bins_pr, how="left").df

        # Generate hypergraph via vectorized groupby
        grouped = midpoint_to_bin.groupby("read_name")["bin_id"].apply(tuple)
        hyper_graph = grouped[grouped.apply(len) > 1].tolist()

        # Save output pickle file to output_dir
        output_file_name = f"{train_cell_line}_{file_id}_chr_{chrom_name}_hypergraph_{resolution}"
        output_path = os.path.join(output_dir, output_file_name)

        with open(output_path, "wb") as f:
            pickle.dump(hyper_graph, f)

        print(f"Successfully processed {chrom_name}, saved to: {output_path}")


if __name__ == "__main__":
    # Example usage:
    porec_process(file_id="FC2")
