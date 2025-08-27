#!/bin/bash
#SBATCH -A your_account
#SBATCH -J Full_CnT_Pipeline
#SBATCH -p ct56
#SBATCH --ntasks=56      # Requesting more cores for the parallel counting step
#SBATCH --mem=128g       # Requesting more memory
#SBATCH -o %j.out
#SBATCH -e %j.err

# --- Script Safety Settings ---
set -e
set -u
set -o pipefail

################################################################################
#                                                                              #
#             PART 0: Global Parameters and Directory Structure Setup          #
#                                                                              #
################################################################################

echo "INFO: [PART 0] Setting up global parameters and directory structure..."

# === Global Parameters ===
# --- [1] Core Paths ---
proj_root=/work/your_account/CUTandTAG/TS250710005
pkg=/work/your_account/CUTandTAG/CUTandTag_pkg
filelist=$proj_root/filelist.txt

# --- [2] Directory Structure ---
# STAGE 1: Pre-processing outputs
raw_data_dir=$proj_root/00_raw_data
trimmed_fastq_dir=$proj_root/01_trimmed_fastq
bam_dir=$proj_root/02_alignment_bam
bw_dir=$proj_root/03_browser_tracks
peak_dir_base=$proj_root/04_peaks
seacr_dir=$peak_dir_base/SEACR_individual

# STAGE 2: Peak calling, merging, and counting outputs (from your new script)
macs3_individual_dir=$peak_dir_base/MACS3_individual
replicate_merged_dir=$peak_dir_base/MACS3_replicate_merged
mark_consensus_dir=$peak_dir_base/MACS3_mark_consensus
quant_dir=$proj_root/05_quantification_counts
log_dir=$proj_root/logs

# --- [3] Static & Software Files ---
gtf_file=$pkg/gencode.vM36.annotation.gtf
BLACKLIST_FILE="/work/j120885731/CUTandTAG/CUTandTag_pkg/mm39.excluderanges_2.bed"
export MACS3_PATH="$HOME/.local/bin/macs3" # As requested

# Export variables needed for parallel functions in STAGE 1
export proj_root pkg filelist raw_data_dir trimmed_fastq_dir bam_dir bw_dir \
        seacr_dir macs3_individual_dir log_dir gtf_file MACS3_PATH BLACKLIST_FILE

# Create the complete directory structure, including new sub-directories for BigWigs
mkdir -p "$trimmed_fastq_dir" "$bam_dir" "$log_dir" "$seacr_dir" \
         "$macs3_individual_dir" "$replicate_merged_dir" "$mark_consensus_dir" \
         "$quant_dir" "$log_dir" \
         "$bw_dir/raw_bw" "$bw_dir/cpm_normalized_bw" "$bw_dir/spike-in_normalized_bw"

# Clear previous summary file
> "$proj_root/QC_summary.txt"

################################################################################
#                                                                              #
#             STAGE 1: Sequence Pre-processing (Trimming to BigWig)            #
#                                                                              #
################################################################################

# === [STAGE 1] Environment Setup Function ===
# This function sets up the environment for trimming, alignment, and BAM processing.
init_env_stage1() {
  source /etc/profile.d/modules.sh
  module purge
  module load old-module
  module load biology/Trimmomatic/0.39
  module load biology/SAMTOOLS/1.18
  module load biology/bowtie2/2.4.2
  module load biology/BEDTOOLS/2.31.1
  module load biology/Picard/2.27.4
  export PATH=$PATH:/opt/ohpc/Taiwania3/pkg/biology/UCSC_Utilities/Utilities_20180515
}
export -f init_env_stage1

# === [STAGE 1] Trimming Function ===
run_trim() {
  init_env_stage1
  local sample_id="$1"; local i5="$2"; local i7="$3"
  local fq1_candidates=("$raw_data_dir/${sample_id}"*R1*.fastq.gz); local fq2_candidates=("$raw_data_dir/${sample_id}"*R2*.fastq.gz)
  
  # --- Robust File Check ---
  if [[ ! -f "${fq1_candidates[0]}" ]]; then echo "[ERROR] No R1 file found for pattern: ${raw_data_dir}/${sample_id}*R1*.fastq.gz" >&2; return 1; fi
  if [[ ${#fq1_candidates[@]} -ne 1 ]]; then echo "[ERROR] Found ${#fq1_candidates[@]} R1 files for ${sample_id}, expected 1." >&2; return 1; fi
  if [[ ! -f "${fq2_candidates[0]}" ]]; then echo "[ERROR] No R2 file found for pattern: ${raw_data_dir}/${sample_id}*R2*.fastq.gz" >&2; return 1; fi
  if [[ ${#fq2_candidates[@]} -ne 1 ]]; then echo "[ERROR] Found ${#fq2_candidates[@]} R2 files for ${sample_id}, expected 1." >&2; return 1; fi

  local fq1="${fq1_candidates[0]}"; local fq2="${fq2_candidates[0]}"
  local adaptor_file="$pkg/v2_Ad1_${i5}/v2_Ad2_${i7}.fa"
  if [[ ! -f "$adaptor_file" ]]; then adaptor_file="$pkg/default_adaptor.fa"; fi
  local trimmed_fq1="$trimmed_fastq_dir/${sample_id}_R1.paired.fastq.gz"; local trimmed_fq2="$trimmed_fastq_dir/${sample_id}_R2.paired.fastq.gz"
  local unpaired_fq1="$trimmed_fastq_dir/${sample_id}_R1.unpaired.fastq.gz"; local unpaired_fq2="$trimmed_fastq_dir/${sample_id}_R2.unpaired.fastq.gz"
  local log_file="$log_dir/${sample_id}.trimmomatic.log"
  echo "[INFO] STAGE 1: Starting Trimmomatic for sample: $sample_id"
  trimmomatic PE -threads 4 -phred33 "$fq1" "$fq2" "$trimmed_fq1" "$unpaired_fq1" "$trimmed_fq2" "$unpaired_fq2" \
    ILLUMINACLIP:"$adaptor_file":2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25 2> "$log_file"
}
export -f run_trim

# === [STAGE 1] Alignment Function ===
run_alignment() {
  init_env_stage1
  local sample_id=$1
  local fq1="$trimmed_fastq_dir/${sample_id}_R1.paired.fastq.gz"; local fq2="$trimmed_fastq_dir/${sample_id}_R2.paired.fastq.gz"
  local sam="$bam_dir/${sample_id}.sam"; local spike_sam="$bam_dir/${sample_id}_spikein.sam"
  local qc_temp_file="$log_dir/${sample_id}.qc.txt"
  if [[ ! -f "$fq1" || ! -f "$fq2" ]]; then echo "[WARN] Missing trimmed FASTQ for $sample_id. Skipping." >&2; return 1; fi
  > "$qc_temp_file"
  echo "[INFO] STAGE 1: Aligning ${sample_id} to mouse genome..."
  bowtie2 --local --very-sensitive-local --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 4 \
    -x "$pkg/index/mouse/mm39" -1 "$fq1" -2 "$fq2" -S "$sam" 2> "$log_dir/${sample_id}_bowtie2.log"
  local aln_rate=$(grep "overall alignment rate" "$log_dir/${sample_id}_bowtie2.log" | awk '{print $1}'); echo -e "${sample_id}\tmouse_alignment_rate\t${aln_rate}" >> "$qc_temp_file"
  echo "[INFO] STAGE 1: Aligning ${sample_id} to spike-in genome..."
  bowtie2 --local --very-sensitive-local --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 4 \
    --no-overlap --no-dovetail -x "$pkg/index/yeast/yeast" -1 "$fq1" -2 "$fq2" -S "$spike_sam" 2> "$log_dir/${sample_id}_spikein.log"
  local spike_rate=$(grep "overall alignment rate" "$log_dir/${sample_id}_spikein.log" | awk '{print $1}'); echo -e "${sample_id}\tspikein_alignment_rate\t${spike_rate}" >> "$qc_temp_file"
}
export -f run_alignment

# === [STAGE 1] Post-processing, BigWig, and Individual Peak Calling Function ===
run_postprocess_and_peak() {
  init_env_stage1
  local sample_id=$1
  
  # Define file paths
  local sam="$bam_dir/${sample_id}.sam"; local spike_sam="$bam_dir/${sample_id}_spikein.sam"
  local bam_sort="$bam_dir/${sample_id}.sorted.bam"; local bam_rmdup="$bam_dir/${sample_id}.sorted.rmDup.bam"
  local bam_final="$bam_dir/${sample_id}.final.bam" # Final, clean BAM
  local fragments_bed="$bam_dir/${sample_id}.fragments.bed"
  local qc_temp_file="$log_dir/${sample_id}.qc.txt"

  echo "[INFO] STAGE 1: Post-processing ${sample_id}..."
  samtools view -bS "$sam" | samtools sort -o "$bam_sort" -
  picard MarkDuplicates -I "$bam_sort" -O "$bam_rmdup" --REMOVE_DUPLICATES true --METRICS_FILE "$log_dir/${sample_id}_picard.rmDup.txt"
  samtools view -b -F 4 -f 0x2 -q 10 "$bam_rmdup" > "$bam_final"
  samtools index "$bam_final"
  
  samtools sort -n "$bam_final" | bedtools bamtobed -i - -bedpe | \
  awk 'BEGIN{OFS="\t"} $1==$4 && $6-$2 < 1000 {print $1, $2, $6}' | \
  sort -k1,1 -k2,2n -k3,3n > "$fragments_bed"

  # --- BigWig Generation & New folder classification ---
  
  # 1. Raw BigWig
  echo "[INFO] STAGE 1: Generating Raw BigWig for ${sample_id}..."
  local raw_bedgraph="$bam_dir/${sample_id}.raw.bedgraph"
  bedtools genomecov -bg -i "$fragments_bed" -g "$pkg/mm39.chrom.sizes" | sort -k1,1 -k2,2n > "$raw_bedgraph"
  bedGraphToBigWig "$raw_bedgraph" "$pkg/mm39.chrom.sizes" "$bw_dir/raw_bw/${sample_id}.raw.bw"

  # 2. Spike-in Normalized BigWig
  local spikein_reads=$(samtools view -c -F 0x04 "$spike_sam")
  local scale_factor=$(echo "10000 / (($spikein_reads / 2) + 1)" | bc -l) # Added +1 to avoid division by zero
  echo "[INFO] STAGE 1: Generating Spike-in normalized BigWig for ${sample_id} with scale factor ${scale_factor}..."
  local scaled_bedgraph="$bam_dir/${sample_id}.scaled.bedgraph"
  bedtools genomecov -bg -scale "$scale_factor" -i "$fragments_bed" -g "$pkg/mm39.chrom.sizes" | sort -k1,1 -k2,2n > "$scaled_bedgraph"
  bedGraphToBigWig "$scaled_bedgraph" "$pkg/mm39.chrom.sizes" "$bw_dir/spike-in_normalized_bw/${sample_id}.spikein.normalized.bw"

  # 3. CPM Normalized BigWig (NEW)
  local total_mapped_reads=$(samtools view -c -F 0x04 "$bam_final")
  local cpm_scale_factor=$(echo "1000000 / $total_mapped_reads" | bc -l)
  echo "[INFO] STAGE 1: Generating CPM normalized BigWig for ${sample_id} with scale factor ${cpm_scale_factor}..."
  local cpm_bedgraph="$bam_dir/${sample_id}.cpm.bedgraph"
  bedtools genomecov -bg -scale "$cpm_scale_factor" -i "$fragments_bed" -g "$pkg/mm39.chrom.sizes" | sort -k1,1 -k2,2n > "$cpm_bedgraph"
  bedGraphToBigWig "$cpm_bedgraph" "$pkg/mm39.chrom.sizes" "$bw_dir/cpm_normalized_bw/${sample_id}.cpm.normalized.bw"
  
  # --- Individual Peak Calling (SEACR & MACS3) ---
  echo "[INFO] STAGE 1: Running SEACR for ${sample_id}..."
  bash "$pkg/SEACR_1.3.sh" "$scaled_bedgraph" 0.01 non stringent "$seacr_dir/${sample_id}_seacr.peaks"

  echo "[INFO] STAGE 1: Running MACS3 and filtering for ${sample_id}..."
  local sample_macs3_outdir="${macs3_individual_dir}/${sample_id}"
  mkdir -p "$sample_macs3_outdir"
  "$MACS3_PATH" callpeak -t "$bam_final" -f BAMPE -g mm -q 0.01 --nomodel --shift 0 --extsize 200 \
    --outdir "$sample_macs3_outdir" -n "${sample_id}" --bdg > "${log_dir}/${sample_id}_macs3.log" 2>&1
  
  local macs3_raw_peak_file="${sample_macs3_outdir}/${sample_id}_peaks.narrowPeak"
  local final_filtered_peaks="${macs3_individual_dir}/${sample_id}.final_filtered_peaks.bed" # IMPORTANT output for STAGE 2
  bedtools intersect -v -a "$macs3_raw_peak_file" -b "$BLACKLIST_FILE" \
    | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"peak_"NR}' > "$final_filtered_peaks"

  # --- FrIP Score Calculation ---
  local total_reads=$(samtools view -c "$bam_final"); local reads_in_peaks=$(bedtools intersect -u -a "$bam_final" -b "$macs3_raw_peak_file" | wc -l)
  if [ "$total_reads" -gt 0 ]; then
    local frip_score=$(awk "BEGIN {print ($reads_in_peaks / $total_reads) * 100}")
    echo -e "${sample_id}\tFrIP_score_MACS3\t${frip_score}" >> "$qc_temp_file"
  fi
  
  # --- Clean up ---
  rm "$sam" "$spike_sam" "$bam_sort" "$bam_rmdup" "$scaled_bedgraph" "$raw_bedgraph" "$cpm_bedgraph"
  echo "[INFO] STAGE 1: Post-processing for ${sample_id} finished."
}
export -f run_postprocess_and_peak

# === [STAGE 1] Execution Logic ===
echo "========================================================"
echo "INFO: Starting STAGE 1: Per-Sample Pre-processing..."
echo "========================================================"
cat "$filelist" | parallel -j $SLURM_NTASKS --colsep '\s+' --joblog "$log_dir/parallel_trim.log" run_trim {1} {2} {3}
cut -f1 "$filelist" | parallel -j $SLURM_NTASKS --joblog "$log_dir/parallel_align.log" run_alignment {}
cut -f1 "$filelist" | parallel -j $SLURM_NTASKS --joblog "$log_dir/parallel_postprocess.log" run_postprocess_and_peak {}

echo "INFO: STAGE 1 Finished. Consolidating QC results..."
{ echo -e "SampleID\tMetric\tValue"; find "$log_dir" -name "*.qc.txt" -print0 | xargs -0 cat; } > "$proj_root/QC_summary.txt"

---

################################################################################
#                                                                              #
#      STAGE 2: Peak Merging, Consensus Calling & Quantification               #
#                                                                              #
################################################################################

echo "========================================================"
echo "INFO: Starting STAGE 2: Peak Merging and Quantification..."
echo "========================================================"

# === [STAGE 2] Environment Setup ===
echo "INFO: [STAGE 2] Setting up Conda environment..."
ml purge
ml load old-module
ml load biology/BEDTOOLS/2.31.1
ml load biology/Anaconda/Anaconda3
source activate featurecounts_env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate featurecounts_env
echo "INFO: [STAGE 2] Conda environment 'featurecounts_env' activated."

# === [STAGE 2] Part 2.1: Merge Biological Replicates ===
echo "INFO: [STAGE 2.1] Finding reproducible peaks from biological replicates..."
ALL_SAMPLES=($(find "${bam_dir}" -name "*.final.bam" -exec basename {} .final.bam \; | sort -u))
UNIQUE_BASE_SAMPLES=($(printf "%s\n" "${ALL_SAMPLES[@]}" | sed -E 's/_[0-9]+$//' | sort -u))

for base_name in "${UNIQUE_BASE_SAMPLES[@]}"; do
  replicate_files=($(find "$macs3_individual_dir" -type f -name "${base_name}_*.final_filtered_peaks.bed"))
  num_reps=${#replicate_files[@]}
  output_bed="${replicate_merged_dir}/${base_name}.reproducible_peaks.bed"

  echo "INFO: Processing group '${base_name}': found ${num_reps} replicate(s)."
  if [ "$num_reps" -gt 1 ]; then
    sorted_first_file=$(mktemp); sort -k1,1 -k2,2n "${replicate_files[0]}" > "$sorted_first_file"
    rest_files=("${replicate_files[@]:1}")
    bedtools intersect -a "$sorted_first_file" -b "${rest_files[@]}" -wa -sorted > "$output_bed"
    rm "$sorted_first_file"
    if [ ! -s "$output_bed" ]; then echo "  -> WARNING: No overlapping peaks found for ${base_name}."; fi
  elif [ "$num_reps" -eq 1 ]; then
    cp "${replicate_files[0]}" "$output_bed"
  else
    echo "  -> WARNING: No peak files found for group ${base_name}."
  fi
done

# === [STAGE 2] Part 2.2: Define Parallel Function for Counting ===
HISTONE_MARKS=($(printf "%s\n" "${ALL_SAMPLES[@]}" | cut -d'_' -f3 | sort -u))
echo "INFO: [STAGE 2.2] Found ${#HISTONE_MARKS[@]} unique marks to process in parallel: ${HISTONE_MARKS[*]}"

# Define how many histone marks to process in parallel.
NUM_PARALLEL_JOBS=5
THREADS_PER_JOB=$(( ${SLURM_NTASKS:-56} / NUM_PARALLEL_JOBS ))
if [ "$THREADS_PER_JOB" -lt 1 ]; then THREADS_PER_JOB=1; fi

process_histone_mark() {
  local histone_mark="$1"
  local log_file="${log_dir}/log_counting_${histone_mark}.log"
  {
    echo "INFO: Starting parallel job for MARK: ${histone_mark}"
    local MARK_CONSENSUS_DIR_SUB="${mark_consensus_dir}/${histone_mark}"
    mkdir -p "$MARK_CONSENSUS_DIR_SUB"
    local ALL_PEAKS_CAT_MARK="${MARK_CONSENSUS_DIR_SUB}/${histone_mark}_all_reproducible_peaks.bed"
    local MERGED_PEAKS_BED_MARK="${MARK_CONSENSUS_DIR_SUB}/${histone_mark}_consensus_peaks.bed"
    local FINAL_ANNOTATION_GTF_MARK="${MARK_CONSENSUS_DIR_SUB}/${histone_mark}_consensus_peaks.gtf"
    local FINAL_COUNTS_FILE="${quant_dir}/${histone_mark}_Raw_Counts.txt"

    echo "INFO: [1/4] Concatenating REPRODUCIBLE peaks for ${histone_mark}..."
    find "${replicate_merged_dir}" -type f -name "*_${histone_mark}_*.reproducible_peaks.bed" -print0 | xargs -0 --no-run-if-empty cat > "$ALL_PEAKS_CAT_MARK"
    if [ ! -s "$ALL_PEAKS_CAT_MARK" ]; then echo "WARNING: No reproducible peaks for ${histone_mark}. Exiting."; return; fi

    echo "INFO: [2/4] Creating consensus peaks for ${histone_mark}..."
    bedtools sort -i "$ALL_PEAKS_CAT_MARK" | bedtools merge -i stdin -c 4 -o count > "$MERGED_PEAKS_BED_MARK"

    echo "INFO: [3/4] Creating GTF for featureCounts for ${histone_mark}..."
    awk -v mark="$histone_mark" 'BEGIN{OFS="\t"} {peak_id=mark"_peak_"NR; print $1, "consensus_peak", "exon", $2+1, $3, ".", "+", ".", "gene_id \""peak_id"\"; transcript_id \""peak_id"\";"}' "$MERGED_PEAKS_BED_MARK" > "$FINAL_ANNOTATION_GTF_MARK"

    echo "INFO: [4/4] Running featureCounts for ${histone_mark}..."
    local MARK_BAM_FILES=()
    while IFS= read -r bam_path; do
      if [[ "$(basename "$bam_path")" == *_${histone_mark}_* ]]; then
        MARK_BAM_FILES+=("$bam_path")
      fi
    done < <(find "$bam_dir" -name "*.final.bam")
    
    if [ ${#MARK_BAM_FILES[@]} -eq 0 ]; then echo "WARNING: No BAMs for ${histone_mark}. Exiting."; return; fi
    
    featureCounts -a "$FINAL_ANNOTATION_GTF_MARK" -o "$FINAL_COUNTS_FILE" -F GTF \
      -T "$THREADS_PER_JOB" -p -g gene_id "${MARK_BAM_FILES[@]}"

    echo "INFO: Job for ${histone_mark} finished successfully."
  } > "$log_file" 2>&1
}
export -f process_histone_mark
export replicate_merged_dir mark_consensus_dir quant_dir THREADS_PER_JOB log_dir bam_dir

# === [STAGE 2] Part 2.3: Launch Parallel Counting Jobs ===
echo "INFO: [STAGE 2.3] Launching parallel counting jobs with GNU Parallel..."
parallel --jobs ${NUM_PARALLEL_JOBS} --joblog "${log_dir}/parallel_counting_runtimes.log" \
  --progress --eta process_histone_mark ::: "${HISTONE_MARKS[@]}"

---

## ????

```bash
# === Final Cleanup Function ===
clean_up() {
    echo "INFO: Cleaning up temporary files..."
    # Remove temporary .auc and .bed files from proj_root
    find "$proj_root" -maxdepth 1 -name "*.auc" -delete
    find "$proj_root" -maxdepth 1 -name "*.bed" -delete
    echo "INFO: Cleanup complete."
}

# --- Final Execution ---
echo ""
echo "========================================================"
echo "??? Full CUT&Tag Pipeline FINISHED! ???"
echo "Check ${log_dir}/parallel_counting_runtimes.log for job exit codes."
echo "Final count tables are in: ${quant_dir}"
echo "Next step: Use the count tables for DESeq2 analysis."
echo "========================================================"

# Run the cleanup function at the end
clean_up
