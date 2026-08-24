import glob
import os

def create_directories_if_not_exists(path):
    """
    Create directories if they don't exist.
    """
    if not os.path.exists(path):
        os.makedirs(path)

    return True

CELL_LINE_BIGWIGS = {
    "GM12878": [
        "ENCFF003DXG", "ENCFF340JIF", "ENCFF039JOT", "ENCFF564KBE", "ENCFF380LZI",
        "ENCFF683HCZ", "ENCFF599TRR", "ENCFF627OKN", "ENCFF479XIQ", "ENCFF601YET",
        "ENCFF931USZ", "ENCFF485CGE", "ENCFF200WHZ", "ENCFF571ZJJ", "ENCFF603BJO"
    ],
    "K562": [
        "ENCFF525ZRM", "ENCFF381NDD", "ENCFF928NWQ", "ENCFF761XBZ", "ENCFF440XMD",
        "ENCFF812HRW", "ENCFF937QUK", "ENCFF959YJV", "ENCFF605FAF", "ENCFF494WCA",
        "ENCFF544AVW", "ENCFF675GVW", "ENCFF124WLE", "ENCFF652NKM", "ENCFF754EAC"
    ],
    "A549": [
        "ENCFF242FAU", "ENCFF070DKP", "ENCFF702IOJ", "ENCFF160YWB", "ENCFF473XIC",
        "ENCFF142SPT", "ENCFF808VAQ", "ENCFF479HXK", "ENCFF417UUX", "ENCFF177CPK",
        "ENCFF375NRQ", "ENCFF109XKO", "ENCFF774RVE", "ENCFF498DXU", "ENCFF872SDF"
    ],
    "HepG2": [
        "ENCFF500VAH", "ENCFF022TZG", "ENCFF437XHN", "ENCFF576YVM", "ENCFF488DNL",
        "ENCFF754ROM", "ENCFF053ROV", "ENCFF057BKO", "ENCFF330AIV", "ENCFF253PND",
        "ENCFF655XBP", "ENCFF301SGJ", "ENCFF761IJZ", "ENCFF242MRW", "ENCFF664EJT"
    ],
    "H1": [
        "ENCFF493QWY", "ENCFF314KQD", "ENCFF345VHG", "ENCFF088MXE", "ENCFF488THD",
        "ENCFF183MHJ", "ENCFF084JKQ", "ENCFF860NVB", "ENCFF156JZY", "ENCFF296IBP",
        "ENCFF401PZS", "ENCFF648BTZ", "ENCFF933YTR", "ENCFF002NBT", "4DNFICPNO4M5"
    ]
}

config_train = {}
config_test = {}
chro = "chr0"

###########################set by users########################
train_cell_line = "GM12878"

config_train["train_test_cell_name"] = train_cell_line+"_"+test_cell_line
config_train["model_dir"] = "models/"+train_cell_line+"_"+test_cell_line+"_hypergcl"

config_train["hic_dir_path"] = "data/"+train_cell_line+"_hg38/hic_mcool/"
config_train["feature_dir_path"] = "data/"+train_cell_line+"_hg38/bigwig_features/"
config_train["input_graph_dir_path"] = "data/"+train_cell_line+"_hg38/input_graph/"
config_train["porec_dir_path"] = "raw_data/"+train_cell_line+"_hg38//hi_pore_c/"

config_train["bigwigs"] =  CELL_LINE_BIGWIGS[train_cell_line]

###########################set by users########################

config_train["dir_check"] = create_directories_if_not_exists(config_train["model_dir"])
config_train["dir_check"] = create_directories_if_not_exists(config_train["hic_dir_path"] )
config_train["dir_check"] = create_directories_if_not_exists(config_train["feature_dir_path"])
config_train["dir_check"] = create_directories_if_not_exists(config_train["input_graph_dir_path"])
config_train["dir_check"] = create_directories_if_not_exists(config_train["porec_dir_path"])

config_train["edge_file_name"] = train_cell_line+f"_chr_{chro}_index_5000.npy"
config_train["weight_file_name"] = train_cell_line+f"_chr_{chro}_attr_5000.npy"
config_train["feature_file_name"] = train_cell_line+f"_{chro}_x.npy"
config_train["hyperedge_file_name"] = train_cell_line+f"_{chro}_hypergraph_5000"

config_train["porec_dirs"] = glob.glob(config_train["porec_dir_path"]+"/*_hyper/")

config_train["learning_rate"] = 0.01
config_train["encoder_hidden_channel"] = 1024
config_train["encoder_out_channel"] = 128
config_train["num_layers"] = 5

config_train["chr_names"] = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                "chr10", "chr11", "chr12","chr16", "chr17", "chr18", "chr20", "chr21"]

config_train["test_chr_names"] = ["chr13", "chr14", "chr15", "chrX", "chr22"]

if train_cell_line == "GM12878":
    config_train["fc_ids"] = ["FC1", "FC2", "FC3", "FC4", "FC5", "FC6", "FC7"]

if train_cell_line == "K562":
    config_train["fc_ids"] = ["FC1", "FC2", "FC3", "FC4"]

###########################set by users########################
test_cell_line = "H1"

config_test["hic_file"] = "raw_data/H1_hg38/4DNFID162B9J.hic"
config_test["bigwig_dir"] = "raw_data/H1_hg38/bigwigs/"
config_test["bigwigs"] =  CELL_LINE_BIGWIGS[test_cell_line]

###########################set by users########################

config_test["hic_dir_path"] = "data/"+test_cell_line+"_hg38/hic_mcool/"
config_test["feature_dir_path"] = "data/"+test_cell_line+"_hg38/bigwig_features/"
config_test["input_graph_dir_path"] = "data/"+test_cell_line+"_hg38/input_graph/"

config_test["dir_check"] = create_directories_if_not_exists(config_test["feature_dir_path"])
config_test["dir_check"] = create_directories_if_not_exists(config_test["hic_dir_path"])
config_test["dir_check"] = create_directories_if_not_exists(config_test["input_graph_dir_path"])

config_test["edge_file_name"] = test_cell_line+f"_chr_{chro}_index_5000.npy"
config_test["weight_file_name"] = test_cell_line+f"_chr_{chro}_attr_5000.npy"
config_test["hyperedge_file_name"] = test_cell_line+f"_{chro}_hypergraph_5000"

config_test["feature_file_name"] = test_cell_line+f"_{chro}_x.npy"

config_test["chr_names"] = []

config_test["test_chr_names"] = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
                      "chr10", "chr11", "chr12","chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
                      "chr22", "chrX"]
