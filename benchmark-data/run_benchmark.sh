#!/usr/bin/env bash
set -euo pipefail

# Benchmark: rseqc-redux vs original RSeQC 5.0.4
# Uses chr22 subset of the official RSeQC test data

BAM="benchmark-data/chr22.bam"
BED="benchmark-data/chr22.bed"
CHROM="benchmark-data/hg19.chrom.sizes"
FULL_BED="benchmark-data/hg19_RefSeq.bed"

ORIG="/tmp/rseqc-original-venv/bin"
REDUX="uv run"

OUTDIR_ORIG="benchmark-data/output-original"
OUTDIR_REDUX="benchmark-data/output-redux"
PROFILE_DIR="benchmark-data/profiles"
RESULTS="benchmark-data/benchmark_results.tsv"

mkdir -p "$OUTDIR_ORIG" "$OUTDIR_REDUX" "$PROFILE_DIR"

echo -e "script\tversion\treal_sec\tuser_sec\tsys_sec\texit_code" > "$RESULTS"

run_timed() {
    local script="$1"
    local version="$2"
    local outfile="$3"
    shift 3
    local cmd=("$@")

    local timefile
    timefile=$(mktemp)

    echo "=== Running: $script ($version) ==="
    local exit_code=0

    # Use GNU/BSD time with portable format
    /usr/bin/time -p "${cmd[@]}" > "$outfile" 2> "$timefile" || exit_code=$?

    # Parse time output (last 3 lines from /usr/bin/time -p)
    local real_sec user_sec sys_sec
    real_sec=$(grep '^real' "$timefile" | awk '{print $2}' || echo "0")
    user_sec=$(grep '^user' "$timefile" | awk '{print $2}' || echo "0")
    sys_sec=$(grep '^sys' "$timefile" | awk '{print $2}' || echo "0")

    # Stderr output from the script goes to timefile too - save it
    cp "$timefile" "${outfile}.stderr"

    echo -e "${script}\t${version}\t${real_sec}\t${user_sec}\t${sys_sec}\t${exit_code}" >> "$RESULTS"
    rm -f "$timefile"
    echo "    -> ${real_sec}s (exit: ${exit_code})"
}

run_profiled() {
    local script="$1"
    local profile_out="$2"
    shift 2
    local cmd=("$@")

    echo "=== Profiling: $script ==="
    # Use cProfile to generate a .prof file for the redux version
    uv run python -m cProfile -o "$profile_out" "${cmd[@]}" > /dev/null 2>&1 || true
}

echo "============================================"
echo "  Benchmark: rseqc-redux vs RSeQC 5.0.4"
echo "  BAM: $BAM (chr22, ~430K reads)"
echo "  BED: $BED (1262 genes)"
echo "============================================"
echo ""

# ── 1. bam_stat ──
run_timed "bam_stat" "original" "$OUTDIR_ORIG/bam_stat.txt" \
    "$ORIG/bam_stat.py" -i "$BAM"
run_timed "bam_stat" "redux" "$OUTDIR_REDUX/bam_stat.txt" \
    $REDUX bam_stat -i "$BAM"

# ── 2. infer_experiment ──
run_timed "infer_experiment" "original" "$OUTDIR_ORIG/infer_experiment.txt" \
    "$ORIG/infer_experiment.py" -i "$BAM" -r "$BED"
run_timed "infer_experiment" "redux" "$OUTDIR_REDUX/infer_experiment.txt" \
    $REDUX infer_experiment -i "$BAM" -r "$BED"

# ── 3. read_distribution ──
run_timed "read_distribution" "original" "$OUTDIR_ORIG/read_distribution.txt" \
    "$ORIG/read_distribution.py" -i "$BAM" -r "$BED"
run_timed "read_distribution" "redux" "$OUTDIR_REDUX/read_distribution.txt" \
    $REDUX read_distribution -i "$BAM" -r "$BED"

# ── 4. read_duplication ──
run_timed "read_duplication" "original" "$OUTDIR_ORIG/read_duplication.txt" \
    "$ORIG/read_duplication.py" -i "$BAM" -o "$OUTDIR_ORIG/read_dup"
run_timed "read_duplication" "redux" "$OUTDIR_REDUX/read_duplication.txt" \
    $REDUX read_duplication -i "$BAM" -o "$OUTDIR_REDUX/read_dup"

# ── 5. read_GC ──
run_timed "read_GC" "original" "$OUTDIR_ORIG/read_GC.txt" \
    "$ORIG/read_GC.py" -i "$BAM" -o "$OUTDIR_ORIG/read_GC"
run_timed "read_GC" "redux" "$OUTDIR_REDUX/read_GC.txt" \
    $REDUX read_GC -i "$BAM" -o "$OUTDIR_REDUX/read_GC"

# ── 6. read_NVC ──
run_timed "read_NVC" "original" "$OUTDIR_ORIG/read_NVC.txt" \
    "$ORIG/read_NVC.py" -i "$BAM" -o "$OUTDIR_ORIG/read_NVC"
run_timed "read_NVC" "redux" "$OUTDIR_REDUX/read_NVC.txt" \
    $REDUX read_NVC -i "$BAM" -o "$OUTDIR_REDUX/read_NVC"

# ── 7. read_quality ──
run_timed "read_quality" "original" "$OUTDIR_ORIG/read_quality.txt" \
    "$ORIG/read_quality.py" -i "$BAM" -o "$OUTDIR_ORIG/read_quality"
run_timed "read_quality" "redux" "$OUTDIR_REDUX/read_quality.txt" \
    $REDUX read_quality -i "$BAM" -o "$OUTDIR_REDUX/read_quality"

# ── 8. clipping_profile ──
run_timed "clipping_profile" "original" "$OUTDIR_ORIG/clipping_profile.txt" \
    "$ORIG/clipping_profile.py" -i "$BAM" -o "$OUTDIR_ORIG/clipping" -s "PE"
run_timed "clipping_profile" "redux" "$OUTDIR_REDUX/clipping_profile.txt" \
    $REDUX clipping_profile -i "$BAM" -o "$OUTDIR_REDUX/clipping" -s "PE"

# ── 9. insertion_profile ──
run_timed "insertion_profile" "original" "$OUTDIR_ORIG/insertion_profile.txt" \
    "$ORIG/insertion_profile.py" -i "$BAM" -o "$OUTDIR_ORIG/insertion" -s "PE"
run_timed "insertion_profile" "redux" "$OUTDIR_REDUX/insertion_profile.txt" \
    $REDUX insertion_profile -i "$BAM" -o "$OUTDIR_REDUX/insertion" -s "PE"

# ── 10. deletion_profile ──
run_timed "deletion_profile" "original" "$OUTDIR_ORIG/deletion_profile.txt" \
    "$ORIG/deletion_profile.py" -i "$BAM" -o "$OUTDIR_ORIG/deletion"
run_timed "deletion_profile" "redux" "$OUTDIR_REDUX/deletion_profile.txt" \
    $REDUX deletion_profile -i "$BAM" -o "$OUTDIR_REDUX/deletion"

# ── 11. mismatch_profile ──
run_timed "mismatch_profile" "original" "$OUTDIR_ORIG/mismatch_profile.txt" \
    "$ORIG/mismatch_profile.py" -i "$BAM" -o "$OUTDIR_ORIG/mismatch" -l 51
run_timed "mismatch_profile" "redux" "$OUTDIR_REDUX/mismatch_profile.txt" \
    $REDUX mismatch_profile -i "$BAM" -o "$OUTDIR_REDUX/mismatch" -l 51

# ── 12. junction_annotation ──
run_timed "junction_annotation" "original" "$OUTDIR_ORIG/junction_annotation.txt" \
    "$ORIG/junction_annotation.py" -i "$BAM" -r "$BED" -o "$OUTDIR_ORIG/junc_annot"
run_timed "junction_annotation" "redux" "$OUTDIR_REDUX/junction_annotation.txt" \
    $REDUX junction_annotation -i "$BAM" -r "$BED" -o "$OUTDIR_REDUX/junc_annot"

# ── 13. junction_saturation ──
run_timed "junction_saturation" "original" "$OUTDIR_ORIG/junction_saturation.txt" \
    "$ORIG/junction_saturation.py" -i "$BAM" -r "$BED" -o "$OUTDIR_ORIG/junc_sat"
run_timed "junction_saturation" "redux" "$OUTDIR_REDUX/junction_saturation.txt" \
    $REDUX junction_saturation -i "$BAM" -r "$BED" -o "$OUTDIR_REDUX/junc_sat"

# ── 14. inner_distance ──
run_timed "inner_distance" "original" "$OUTDIR_ORIG/inner_distance.txt" \
    "$ORIG/inner_distance.py" -i "$BAM" -r "$BED" -o "$OUTDIR_ORIG/inner_dist"
run_timed "inner_distance" "redux" "$OUTDIR_REDUX/inner_distance.txt" \
    $REDUX inner_distance -i "$BAM" -r "$BED" -o "$OUTDIR_REDUX/inner_dist"

# ── 15. geneBody_coverage ──
run_timed "geneBody_coverage" "original" "$OUTDIR_ORIG/geneBody_coverage.txt" \
    "$ORIG/geneBody_coverage.py" -i "$BAM" -r "$BED" -o "$OUTDIR_ORIG/genebody"
run_timed "geneBody_coverage" "redux" "$OUTDIR_REDUX/geneBody_coverage.txt" \
    $REDUX geneBody_coverage -i "$BAM" -r "$BED" -o "$OUTDIR_REDUX/genebody"

# ── 16. tin ──
run_timed "tin" "original" "$OUTDIR_ORIG/tin.txt" \
    "$ORIG/tin.py" -i "$BAM" -r "$BED"
run_timed "tin" "redux" "$OUTDIR_REDUX/tin.txt" \
    $REDUX tin -i "$BAM" -r "$BED"

# ── 17. RPKM_saturation ──
run_timed "RPKM_saturation" "original" "$OUTDIR_ORIG/RPKM_saturation.txt" \
    "$ORIG/RPKM_saturation.py" -i "$BAM" -r "$BED" -o "$OUTDIR_ORIG/rpkm_sat"
run_timed "RPKM_saturation" "redux" "$OUTDIR_REDUX/RPKM_saturation.txt" \
    $REDUX RPKM_saturation -i "$BAM" -r "$BED" -o "$OUTDIR_REDUX/rpkm_sat"

# ── 18. FPKM_count ──
run_timed "FPKM_count" "original" "$OUTDIR_ORIG/FPKM_count.txt" \
    "$ORIG/FPKM_count.py" -i "$BAM" -r "$BED" -o "$OUTDIR_ORIG/fpkm_count"
run_timed "FPKM_count" "redux" "$OUTDIR_REDUX/FPKM_count.txt" \
    $REDUX FPKM_count -i "$BAM" -r "$BED" -o "$OUTDIR_REDUX/fpkm_count"

# ── 19. RNA_fragment_size ──
run_timed "RNA_fragment_size" "original" "$OUTDIR_ORIG/RNA_fragment_size.txt" \
    "$ORIG/RNA_fragment_size.py" -i "$BAM" -r "$BED"
run_timed "RNA_fragment_size" "redux" "$OUTDIR_REDUX/RNA_fragment_size.txt" \
    $REDUX RNA_fragment_size -i "$BAM" -r "$BED"

# ── 20. read_hexamer ──
# read_hexamer needs a refgene file (BED) — it reads sequences from BAM
run_timed "read_hexamer" "original" "$OUTDIR_ORIG/read_hexamer.txt" \
    "$ORIG/read_hexamer.py" -i "$BAM" -r "$BED"
run_timed "read_hexamer" "redux" "$OUTDIR_REDUX/read_hexamer.txt" \
    $REDUX read_hexamer -i "$BAM" -r "$BED"

# ── 21. split_bam ──
run_timed "split_bam" "original" "$OUTDIR_ORIG/split_bam.txt" \
    "$ORIG/split_bam.py" -i "$BAM" -r "$BED" -o "$OUTDIR_ORIG/split_bam"
run_timed "split_bam" "redux" "$OUTDIR_REDUX/split_bam.txt" \
    $REDUX split_bam -i "$BAM" -r "$BED" -o "$OUTDIR_REDUX/split_bam"

# ── 22. divide_bam ──
run_timed "divide_bam" "original" "$OUTDIR_ORIG/divide_bam.txt" \
    "$ORIG/divide_bam.py" -i "$BAM" -o "$OUTDIR_ORIG/divide" -n 2
run_timed "divide_bam" "redux" "$OUTDIR_REDUX/divide_bam.txt" \
    $REDUX divide_bam -i "$BAM" -o "$OUTDIR_REDUX/divide" -n 2

# ── 23. bam2wig ──
run_timed "bam2wig" "original" "$OUTDIR_ORIG/bam2wig.txt" \
    "$ORIG/bam2wig.py" -i "$BAM" -s "$CHROM" -o "$OUTDIR_ORIG/bam2wig" -t 1000000
run_timed "bam2wig" "redux" "$OUTDIR_REDUX/bam2wig.txt" \
    $REDUX bam2wig -i "$BAM" -s "$CHROM" -o "$OUTDIR_REDUX/bam2wig" -t 1000000

# ── 24. normalize_bigwig ──
# Needs a bigwig file — we'll create one from bam2wig output if available
# Skip if no bigwig output yet

# ── 25. bam2fq ──
run_timed "bam2fq" "original" "$OUTDIR_ORIG/bam2fq.txt" \
    "$ORIG/bam2fq.py" -i "$BAM" -o "$OUTDIR_ORIG/bam2fq"
run_timed "bam2fq" "redux" "$OUTDIR_REDUX/bam2fq.txt" \
    $REDUX bam2fq -i "$BAM" -o "$OUTDIR_REDUX/bam2fq"

echo ""
echo "============================================"
echo "  Now profiling redux scripts with cProfile"
echo "============================================"

# Profile the most important/slow scripts
for script_args in \
    "bam_stat:scripts/bam_stat.py:-i:$BAM" \
    "infer_experiment:scripts/infer_experiment.py:-i:$BAM:-r:$BED" \
    "read_distribution:scripts/read_distribution.py:-i:$BAM:-r:$BED" \
    "read_duplication:scripts/read_duplication.py:-i:$BAM:-o:$PROFILE_DIR/tmp_read_dup" \
    "read_GC:scripts/read_GC.py:-i:$BAM:-o:$PROFILE_DIR/tmp_read_GC" \
    "read_NVC:scripts/read_NVC.py:-i:$BAM:-o:$PROFILE_DIR/tmp_read_NVC" \
    "read_quality:scripts/read_quality.py:-i:$BAM:-o:$PROFILE_DIR/tmp_read_qual" \
    "clipping_profile:scripts/clipping_profile.py:-i:$BAM:-o:$PROFILE_DIR/tmp_clip:-s:PE" \
    "insertion_profile:scripts/insertion_profile.py:-i:$BAM:-o:$PROFILE_DIR/tmp_ins:-s:PE" \
    "deletion_profile:scripts/deletion_profile.py:-i:$BAM:-o:$PROFILE_DIR/tmp_del" \
    "mismatch_profile:scripts/mismatch_profile.py:-i:$BAM:-o:$PROFILE_DIR/tmp_mismatch:-l:51" \
    "junction_annotation:scripts/junction_annotation.py:-i:$BAM:-r:$BED:-o:$PROFILE_DIR/tmp_junc_annot" \
    "junction_saturation:scripts/junction_saturation.py:-i:$BAM:-r:$BED:-o:$PROFILE_DIR/tmp_junc_sat" \
    "inner_distance:scripts/inner_distance.py:-i:$BAM:-r:$BED:-o:$PROFILE_DIR/tmp_inner" \
    "geneBody_coverage:scripts/geneBody_coverage.py:-i:$BAM:-r:$BED:-o:$PROFILE_DIR/tmp_genebody" \
    "tin:scripts/tin.py:-i:$BAM:-r:$BED" \
    "RPKM_saturation:scripts/RPKM_saturation.py:-i:$BAM:-r:$BED:-o:$PROFILE_DIR/tmp_rpkm_sat" \
    "FPKM_count:scripts/FPKM_count.py:-i:$BAM:-r:$BED:-o:$PROFILE_DIR/tmp_fpkm" \
    "RNA_fragment_size:scripts/RNA_fragment_size.py:-i:$BAM:-r:$BED" \
    "split_bam:scripts/split_bam.py:-i:$BAM:-r:$BED:-o:$PROFILE_DIR/tmp_split" \
    "bam2fq:scripts/bam2fq.py:-i:$BAM:-o:$PROFILE_DIR/tmp_bam2fq" \
; do
    IFS=':' read -r name module args_str <<< "$script_args"
    # Convert colon-separated args to array
    IFS=':' read -ra args <<< "$args_str"

    echo "=== Profiling: $name ==="
    uv run python -m cProfile -o "$PROFILE_DIR/${name}.prof" "$module" "${args[@]}" > /dev/null 2>&1 || echo "    (profiling failed or script errored)"
done

echo ""
echo "============================================"
echo "  Benchmark complete! Results in:"
echo "    $RESULTS"
echo "    $PROFILE_DIR/*.prof"
echo "============================================"
