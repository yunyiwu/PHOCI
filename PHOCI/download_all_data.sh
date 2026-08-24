#!/bin/bash
# ==============================================================================
# Script: download_all_data.sh
# Description: Download Hi-C (.hic), Epigenomic bigWig tracks, and Hi-Pore-C 
#              alignment files, organizing them cleanly by cell line.
# ==============================================================================

set -e # Exit immediately if a command exits with a non-zero status

# Base output directory
BASE_DIR="./raw_data"

# Ensure base output directory exists
mkdir -p "${BASE_DIR}"

echo "=== Starting Full Dataset Download (Hi-C, Epigenomics & Hi-Pore-C) ==="

# ------------------------------------------------------------------------------
# 1. Epigenomic Files Arrays (Accession:Feature)
# ------------------------------------------------------------------------------
GM12878_EPI=(
  "ENCFF003DXG:H3K4me3" "ENCFF340JIF:H3K27ac" "ENCFF039JOT:H3K27me3"
  "ENCFF564KBE:H3K4me1" "ENCFF380LZI:H3K36me3" "ENCFF683HCZ:H3K9me3"
  "ENCFF599TRR:H3K9ac"  "ENCFF627OKN:H3K4me2"  "ENCFF479XIQ:H4K20me1"
  "ENCFF601YET:H2AFZ"   "ENCFF931USZ:H3K79me2" "ENCFF485CGE:CTCF"
  "ENCFF200WHZ:POLR2A"  "ENCFF571ZJJ:RAD21"    "ENCFF603BJO:ATAC"
)

K562_EPI=(
  "ENCFF525ZRM:H3K4me3" "ENCFF381NDD:H3K27ac" "ENCFF928NWQ:H3K27me3"
  "ENCFF761XBZ:H3K4me1" "ENCFF440XMD:H3K36me3" "ENCFF812HRW:H3K9me3"
  "ENCFF937QUK:H3K9ac"  "ENCFF959YJV:H3K4me2"  "ENCFF605FAF:H4K20me1"
  "ENCFF494WCA:H2AFZ"   "ENCFF544AVW:H3K79me2" "ENCFF675GVW:CTCF"
  "ENCFF124WLE:POLR2A"  "ENCFF652NKM:RAD21"    "ENCFF754EAC:ATAC"
)

A549_EPI=(
  "ENCFF242FAU:H3K4me3" "ENCFF070DKP:H3K27ac" "ENCFF702IOJ:H3K27me3"
  "ENCFF160YWB:H3K4me1" "ENCFF473XIC:H3K36me3" "ENCFF142SPT:H3K9me3"
  "ENCFF808VAQ:H3K9ac"  "ENCFF479HXK:H3K4me2"  "ENCFF417UUX:H4K20me1"
  "ENCFF177CPK:H2AFZ"   "ENCFF375NRQ:H3K79me2" "ENCFF109XKO:CTCF"
  "ENCFF774RVE:POLR2A"  "ENCFF498DXU:RAD21"    "ENCFF872SDF:ATAC"
)

HepG2_EPI=(
  "ENCFF500VAH:H3K4me3" "ENCFF022TZG:H3K27ac" "ENCFF437XHN:H3K27me3"
  "ENCFF576YVM:H3K4me1" "ENCFF488DNL:H3K36me3" "ENCFF754ROM:H3K9me3"
  "ENCFF053ROV:H3K9ac"  "ENCFF057BKO:H3K4me2"  "ENCFF330AIV:H4K20me1"
  "ENCFF253PND:H2AFZ"   "ENCFF655XBP:H3K79me2" "ENCFF301SGJ:CTCF"
  "ENCFF761IJZ:POLR2A"  "ENCFF242MRW:RAD21"    "ENCFF664EJT:ATAC"
)

H1_hESC_EPI=(
  "ENCFF493QWY:H3K4me3" "ENCFF314KQD:H3K27ac" "ENCFF345VHG:H3K27me3"
  "ENCFF088MXE:H3K4me1" "ENCFF488THD:H3K36me3" "ENCFF183MHJ:H3K9me3"
  "ENCFF084JKQ:H3K0ac"  "ENCFF860NVB:H3K4me2"  "ENCFF156JZY:H4K20me1"
  "ENCFF296IBP:H2AFZ"   "ENCFF401PZS:H3K79me2" "ENCFF648BTZ:CTCF"
  "ENCFF933YTR:POLR2A"  "ENCFF002NBT:RAD21"    "4DNFICPNO4M5:ATAC"
)

# ------------------------------------------------------------------------------
# 2. Hi-Pore-C GEO URLs Arrays
# ------------------------------------------------------------------------------
GM12878_POREC=(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6124010&format=file&file=GSM6124010%5FGM12878%5FFC1%5Freads%5Falignment%2Ecsv%2Egz"
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6124011&format=file&file=GSM6124011%5FGM12878%5FFC2%5Freads%5Falignment%2Ecsv%2Egz"
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6124012&format=file&file=GSM6124012%5FGM12878%5FFC3%5Freads%5Falignment%2Ecsv%2Egz"
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6124013&format=file&file=GSM6124013%5FGM12878%5FFC4%5Freads%5Falignment%2Ecsv%2Egz"
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6284586&format=file&file=GSM6284586%5FGM12878%5FFC5%5Freads%5Falignment%2Ecsv%2Egz"
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6284587&format=file&file=GSM6284587%5FGM12878%5FFC6%5Freads%5Falignment%2Ecsv%2Egz"
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6732990&format=file&file=GSM6732990%5FGM12878%5FFC7%5Freads%5Falignment%2Ecsv%2Egz"
)

K562_POREC=(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6284588&format=file&file=GSM6284588%5FK562%5FFC1%5Freads%5Falignment%2Ecsv%2Egz"
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6284589&format=file&file=GSM6284589%5FK562%5FFC2%5Freads%5Falignment%2Ecsv%2Egz"
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6284590&format=file&file=GSM6284590%5FK562%5FFC3%5Freads%5Falignment%2Ecsv%2Egz"
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6732991&format=file&file=GSM6732991%5FK562%5FFC4%5Freads%5Falignment%2Ecsv%2Egz"
)

# ------------------------------------------------------------------------------
# 3. Main Master Download Function
# ------------------------------------------------------------------------------
download_cell_line() {
  local cell_line=$1
  local hic_url=$2
  local epi_ref=$3
  local porec_ref=$4

  # De-reference array variables passed by name
  eval "local epi_files=("\${${epi_ref}[@]}")"
  eval "local porec_urls=("\${${porec_ref}[@]}")"

  local hic_dir="${BASE_DIR}/${cell_line}/hic"
  local epi_dir="${BASE_DIR}/${cell_line}/epigenomic"
  local porec_dir="${BASE_DIR}/${cell_line}/hi_pore_c"

  echo "========================================================"
  echo " Processing Cell Line: ${cell_line}"
  echo " Target Base: ${BASE_DIR}/${cell_line}"
  echo "========================================================"

  # [1/3] Download Hi-C File via direct URL
  if [[ -n "${hic_url}" ]]; then
    mkdir -p "${hic_dir}"
    echo "[1/3] Downloading Hi-C matrix for ${cell_line}..."
    local hic_filename=$(basename "${hic_url}")
    wget -c -q --show-progress -O "${hic_dir}/${hic_filename}" "${hic_url}"
  fi

  # [2/3] Download Epigenomic Tracks
  if [[ ${#epi_files[@]} -gt 0 ]]; then
    mkdir -p "${epi_dir}"
    echo "[2/2] Downloading Epigenomic tracks for ${cell_line}..."
    for item in "${epi_files[@]}"; do
      IFS=":" read -r acc feature <<< "${item}"
      out_file="${epi_dir}/${acc}.bigWig"
      
      if [[ "${acc}" == 4DN* ]]; then
        url="https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/aa08799c-f06e-450e-ab63-d8a2a6157b52/${acc}.bw"
      else
        url="https://www.encodeproject.org/files/${acc}/@@download/${acc}.bigWig"
      fi

      echo "  [+] Fetching ${feature} (${acc})..."
      wget -c -q --show-progress -O "${out_file}" "${url}"
    done
  fi

  # [3/3] Download Hi-Pore-C Files (with GEO redirect handling)
  if [[ ${#porec_urls[@]} -gt 0 ]]; then
    mkdir -p "${porec_dir}"
    echo "[3/3] Downloading Hi-Pore-C alignments for ${cell_line}..."
    for url in "${porec_urls[@]}"; do
      echo "  [+] Fetching GEO file..."
      # --trust-server-names automatically resolves URL encoding and 302 redirects from GEO
      wget -c -q --show-progress --trust-server-names -P "${porec_dir}" "${url}"
    done
  fi

  echo ""
}

# ------------------------------------------------------------------------------
# 4. Execute All Downloads
# ------------------------------------------------------------------------------
download_cell_line "GM12878"   "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/aa26f261-a88a-4cac-9118-ff8e90ab6f61/4DNFI9ZWZ5BS.hic"   "GM12878_EPI"   "GM12878_POREC"

download_cell_line "K562"   "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/93b8e020-b337-4685-81dc-40cc8a12b5e9/4DNFIXU2KPNQ.hic"   "K562_EPI"   "K562_POREC"

download_cell_line "A549"   "https://www.encodeproject.org/files/ENCFF689CUX/@@download/ENCFF689CUX.hic"   "A549_EPI"   ""

download_cell_line "HepG2"   "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/25104375-a588-46e6-a382-663cee6c332f/4DNFICSTCJQZ.hic"   "HepG2_EPI"   ""

download_cell_line "H1-hESC"   "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/bb3307fd-7162-477a-87c5-52f12d03befc/4DNFID162B9J.hic"   "H1_hESC_EPI"   ""

echo "=== All datasets downloaded successfully! ==="
