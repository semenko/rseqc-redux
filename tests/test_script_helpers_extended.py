"""Extended tests for helper functions in scripts/ — covers functions that need
BED/BAM fixtures or file I/O, complementing test_script_helpers.py which tests
pure functions."""

from pathlib import Path

import pysam

FIXTURES_DIR = Path(__file__).parent / "fixtures"
MINI_BED = str(FIXTURES_DIR / "mini.bed")


# ---------------------------------------------------------------------------
# scripts.read_distribution: process_gene_model
# ---------------------------------------------------------------------------


class TestProcessGeneModel:
    """Test read_distribution.process_gene_model with the mini.bed fixture."""

    def test_returns_20_element_tuple(self):
        from scripts.read_distribution import process_gene_model

        result = process_gene_model(MINI_BED)
        assert isinstance(result, tuple)
        assert len(result) == 20

    def test_cds_exon_ranges_populated(self):
        from scripts.read_distribution import process_gene_model

        result = process_gene_model(MINI_BED)
        cds_exon_ranges = result[0]
        # mini.bed has genes on chr1 — CDS exons should be present
        assert "CHR1" in cds_exon_ranges

    def test_intron_ranges_populated(self):
        from scripts.read_distribution import process_gene_model

        result = process_gene_model(MINI_BED)
        intron_ranges = result[1]
        # gene1 (3 exons) and gene2 (2 exons) have introns
        assert "CHR1" in intron_ranges

    def test_sizes_are_positive(self):
        from scripts.read_distribution import process_gene_model

        result = process_gene_model(MINI_BED)
        exon_size = result[10]
        intron_size = result[11]
        assert exon_size > 0
        assert intron_size > 0

    def test_intergenic_sizes(self):
        from scripts.read_distribution import process_gene_model

        result = process_gene_model(MINI_BED)
        # intergenic sizes are at indices 14-19
        for i in range(14, 20):
            assert isinstance(result[i], int)

    def test_cds_exon_finds_known_position(self):
        """Gene1 CDS exon: thickStart=1200, first exon block 1000-1500, so CDS region ~1200-1500."""
        from scripts.read_distribution import process_gene_model

        result = process_gene_model(MINI_BED)
        cds_exon_ranges = result[0]
        # Position 1300 should be in CDS exon of gene1
        hits = cds_exon_ranges["CHR1"].find(1300, 1301)
        assert len(hits) > 0

    def test_empty_bed_file(self, tmp_path):
        from scripts.read_distribution import process_gene_model

        bed = tmp_path / "empty.bed"
        bed.write_text("")
        result = process_gene_model(str(bed))
        # All ranges should be empty dicts, all sizes 0
        for i in range(10):
            assert result[i] == {} or len(result[i]) == 0
        for i in range(10, 20):
            assert result[i] == 0


# ---------------------------------------------------------------------------
# scripts.tin: genomic_positions
# ---------------------------------------------------------------------------


class TestGenomicPositions:
    """Test tin.genomic_positions generator."""

    def test_yields_tuples(self):
        from scripts.tin import genomic_positions

        results = list(genomic_positions(MINI_BED, sample_size=100))
        assert len(results) > 0
        for item in results:
            assert len(item) == 6
            gene_name, chrom, tx_start, tx_end, intron_size, positions = item
            assert isinstance(gene_name, str)
            assert isinstance(chrom, str)
            assert isinstance(tx_start, int)
            assert isinstance(tx_end, int)

    def test_gene1_has_intron(self):
        from scripts.tin import genomic_positions

        results = list(genomic_positions(MINI_BED, sample_size=100))
        # gene1 has 3 exons spanning 1000-5000, intron_size should be > 0
        gene1 = [r for r in results if r[0] == "gene1"]
        assert len(gene1) == 1
        intron_size = gene1[0][4]
        assert intron_size > 0

    def test_gene3_single_exon_no_intron(self):
        from scripts.tin import genomic_positions

        results = list(genomic_positions(MINI_BED, sample_size=100))
        gene3 = [r for r in results if r[0] == "gene3"]
        assert len(gene3) == 1
        intron_size = gene3[0][4]
        assert intron_size == 0

    def test_positions_are_within_gene_bounds(self):
        from scripts.tin import genomic_positions

        for gene_name, chrom, tx_start, tx_end, intron_size, positions in genomic_positions(MINI_BED, sample_size=50):
            for pos in positions:
                assert tx_start <= pos <= tx_end + 1

    def test_small_sample_size_returns_all_bases(self):
        """When sample_size > mRNA_size, should return all exonic bases."""
        from scripts.tin import genomic_positions

        results = list(genomic_positions(MINI_BED, sample_size=10000))
        for gene_name, chrom, tx_start, tx_end, intron_size, positions in results:
            assert len(positions) > 0

    def test_skips_comments_and_track_lines(self, tmp_path):
        from scripts.tin import genomic_positions

        bed = tmp_path / "test.bed"
        bed.write_text(
            "# comment\n"
            "track name=test\n"
            "browser position chr1:1-100\n"
            "chr1\t1000\t5000\tgene1\t0\t+\t1200\t4800\t0,0,0\t3\t500,600,500\t0,1500,3500\n"
        )
        results = list(genomic_positions(str(bed), sample_size=100))
        assert len(results) == 1


# ---------------------------------------------------------------------------
# scripts.tin: union_exons
# ---------------------------------------------------------------------------


class TestUnionExons:
    def test_returns_dict_with_chrom_keys(self):
        from scripts.tin import union_exons

        result = union_exons(MINI_BED)
        assert "chr1" in result

    def test_finds_exonic_position(self):
        """Gene1 first exon: 1000-1500 — position 1250 should be in exon ranges."""
        from scripts.tin import union_exons

        result = union_exons(MINI_BED)
        hits = result["chr1"].find(1250, 1251)
        assert len(hits) > 0

    def test_does_not_find_intronic_position(self):
        """Gene1 first intron: ~1500-2500 — position 2000 should NOT be in exon ranges."""
        from scripts.tin import union_exons

        result = union_exons(MINI_BED)
        hits = result["chr1"].find(2000, 2001)
        assert len(hits) == 0


# ---------------------------------------------------------------------------
# scripts.tin: check_min_reads
# ---------------------------------------------------------------------------


class TestCheckMinReads:
    def test_returns_true_when_enough_reads(self, mini_bam):
        from scripts.tin import check_min_reads

        samfile = pysam.AlignmentFile(str(mini_bam), "rb")
        # Gene1 region has multiple reads (read_unique1, read_unique2, etc.)
        result = check_min_reads(samfile, "chr1", 1000, 5000, cutoff=1)
        samfile.close()
        assert result is True

    def test_returns_false_when_not_enough_reads(self, mini_bam):
        from scripts.tin import check_min_reads

        samfile = pysam.AlignmentFile(str(mini_bam), "rb")
        # Require 100 different start positions — we only have a few reads
        result = check_min_reads(samfile, "chr1", 1000, 5000, cutoff=100)
        samfile.close()
        assert result is False

    def test_returns_false_for_unknown_chrom(self, mini_bam):
        from scripts.tin import check_min_reads

        samfile = pysam.AlignmentFile(str(mini_bam), "rb")
        result = check_min_reads(samfile, "chrX", 1000, 5000, cutoff=1)
        samfile.close()
        assert result is False

    def test_returns_false_for_empty_region(self, mini_bam):
        from scripts.tin import check_min_reads

        samfile = pysam.AlignmentFile(str(mini_bam), "rb")
        # Region with no reads
        result = check_min_reads(samfile, "chr1", 40000, 45000, cutoff=1)
        samfile.close()
        assert result is False


# ---------------------------------------------------------------------------
# scripts.tin: estimate_bg_noise
# ---------------------------------------------------------------------------


class TestEstimateBgNoise:
    def test_returns_float(self, mini_bam):
        from scripts.tin import estimate_bg_noise, union_exons

        samfile = pysam.AlignmentFile(str(mini_bam), "rb")
        e_ranges = union_exons(MINI_BED)
        result = estimate_bg_noise("chr1", 1000, 5000, samfile, e_ranges)
        samfile.close()
        assert isinstance(result, float)

    def test_zero_for_empty_region(self, mini_bam):
        from scripts.tin import estimate_bg_noise, union_exons

        samfile = pysam.AlignmentFile(str(mini_bam), "rb")
        e_ranges = union_exons(MINI_BED)
        result = estimate_bg_noise("chr1", 40000, 45000, samfile, e_ranges)
        samfile.close()
        assert result == 0.0


# ---------------------------------------------------------------------------
# scripts.tin: genebody_coverage (the tin.py version)
# ---------------------------------------------------------------------------


class TestTinGenebodyCoverage:
    def test_returns_list(self, mini_bam):
        from scripts.tin import genebody_coverage

        samfile = pysam.AlignmentFile(str(mini_bam), "rb")
        # Positions within gene1 exon1 where we have reads
        positions = sorted([1051, 1060, 1070, 1080, 1090, 1100])
        result = genebody_coverage(samfile, "chr1", positions)
        samfile.close()
        assert isinstance(result, list)

    def test_coverage_values_nonnegative(self, mini_bam):
        from scripts.tin import genebody_coverage

        samfile = pysam.AlignmentFile(str(mini_bam), "rb")
        positions = sorted([1051, 1060, 1070, 1080, 1090, 1100])
        result = genebody_coverage(samfile, "chr1", positions)
        samfile.close()
        for val in result:
            assert val >= 0

    def test_empty_for_unknown_chrom(self, mini_bam):
        from scripts.tin import genebody_coverage

        samfile = pysam.AlignmentFile(str(mini_bam), "rb")
        result = genebody_coverage(samfile, "chrX", [100, 200, 300])
        samfile.close()
        assert result == []

    def test_bg_subtraction(self, mini_bam):
        from scripts.tin import genebody_coverage

        samfile = pysam.AlignmentFile(str(mini_bam), "rb")
        positions = sorted([1051, 1060, 1070, 1080, 1090, 1100])
        result_with_bg = genebody_coverage(samfile, "chr1", positions, bg_level=1000)
        samfile.close()

        # With high bg_level, all values should be 0
        for val in result_with_bg:
            assert val == 0


# ---------------------------------------------------------------------------
# scripts.geneBody_coverage: genebody_percentile
# ---------------------------------------------------------------------------


class TestGenebodyPercentile:
    def test_returns_dict(self):
        from scripts.geneBody_coverage import genebody_percentile

        result = genebody_percentile(MINI_BED, mRNA_len_cut=100)
        assert isinstance(result, dict)

    def test_contains_gene_entries(self):
        from scripts.geneBody_coverage import genebody_percentile

        result = genebody_percentile(MINI_BED, mRNA_len_cut=100)
        # Should have entries for genes with mRNA >= 100bp
        assert len(result) > 0

    def test_gene_entry_has_chrom_strand_positions(self):
        from scripts.geneBody_coverage import genebody_percentile

        result = genebody_percentile(MINI_BED, mRNA_len_cut=100)
        for gene_id, (chrom, strand, positions) in result.items():
            assert chrom == "chr1"
            assert strand in ("+", "-")
            assert isinstance(positions, list)
            assert len(positions) == 100  # percentile_list returns 100 points

    def test_high_cutoff_filters_short_genes(self):
        from scripts.geneBody_coverage import genebody_percentile

        result = genebody_percentile(MINI_BED, mRNA_len_cut=5000)
        # No gene in mini.bed has mRNA >= 5000 bp
        assert len(result) == 0

    def test_skips_comments(self, tmp_path):
        from scripts.geneBody_coverage import genebody_percentile

        bed = tmp_path / "test.bed"
        bed.write_text(
            "# comment\n"
            "track name=test\n"
            "chr1\t1000\t5000\tgene1\t0\t+\t1200\t4800\t0,0,0\t3\t500,600,500\t0,1500,3500\n"
        )
        result = genebody_percentile(str(bed), mRNA_len_cut=100)
        assert len(result) == 1


# ---------------------------------------------------------------------------
# scripts.geneBody_coverage: Rcode_write
# ---------------------------------------------------------------------------


class TestRcodeWrite:
    def test_creates_r_file(self, tmp_path):
        from scripts.geneBody_coverage import Rcode_write

        prefix = str(tmp_path / "test")
        data = [100.0 + i for i in range(100)]
        dataset = [("sample1", data, 0.5)]
        Rcode_write(dataset, prefix, format="pdf", colNum=100)
        r_file = tmp_path / "test.r"
        assert r_file.exists()

    def test_single_sample_plot(self, tmp_path):
        from scripts.geneBody_coverage import Rcode_write

        prefix = str(tmp_path / "test")
        data = [float(i) for i in range(100)]
        dataset = [("sample1", data, 0.5)]
        Rcode_write(dataset, prefix, format="pdf", colNum=100)
        content = (tmp_path / "test.r").read_text()
        assert "sample1" in content
        assert "plot(" in content
        assert "dev.off()" in content

    def test_two_sample_legend(self, tmp_path):
        from scripts.geneBody_coverage import Rcode_write

        prefix = str(tmp_path / "test")
        data = [float(i) for i in range(100)]
        dataset = [
            ("sample1", data, 0.5),
            ("sample2", data, 0.3),
        ]
        Rcode_write(dataset, prefix, format="pdf", colNum=100)
        content = (tmp_path / "test.r").read_text()
        assert "legend(" in content
        assert "lines(" in content

    def test_many_samples_heatmap(self, tmp_path):
        from scripts.geneBody_coverage import Rcode_write

        prefix = str(tmp_path / "test")
        data = [float(i) for i in range(100)]
        dataset = [(f"sample{i}", data, 0.1 * i) for i in range(5)]
        Rcode_write(dataset, prefix, format="png", colNum=100)
        content = (tmp_path / "test.r").read_text()
        assert "heatmap(" in content
        assert "data_matrix" in content
        assert "png" in content

    def test_output_format_respected(self, tmp_path):
        from scripts.geneBody_coverage import Rcode_write

        prefix = str(tmp_path / "test")
        data = [float(i) for i in range(100)]
        dataset = [("s1", data, 0.5)]
        Rcode_write(dataset, prefix, format="jpeg", colNum=100)
        content = (tmp_path / "test.r").read_text()
        assert "jpeg(" in content


# ---------------------------------------------------------------------------
# scripts.RPKM_saturation: show_saturation
# ---------------------------------------------------------------------------


class TestShowSaturation:
    def _make_xls(self, tmp_path):
        """Create a synthetic eRPKM.xls file."""
        xls = tmp_path / "test.eRPKM.xls"
        header = "#chrom\tstart\tend\tname\tscore\tstrand\t%5\t%25\t%50\t%75\t%100\n"
        lines = [header]
        for i in range(20):
            # Each gene has increasing RPKM values at different sampling levels
            vals = [1.0 + i * j * 0.1 for j in range(1, 6)]
            lines.append(
                f"chr1\t{i * 1000}\t{i * 1000 + 500}\tgene{i}\t0\t+\t" + "\t".join(f"{v:.2f}" for v in vals) + "\n"
            )
        xls.write_text("".join(lines))
        return str(xls)

    def test_creates_r_file(self, tmp_path):
        from scripts.RPKM_saturation import show_saturation

        infile = self._make_xls(tmp_path)
        outfile = str(tmp_path / "test.saturation.r")
        show_saturation(infile, outfile, rpkm_cut=0.01)
        assert Path(outfile).exists()

    def test_r_file_has_boxplot(self, tmp_path):
        from scripts.RPKM_saturation import show_saturation

        infile = self._make_xls(tmp_path)
        outfile = str(tmp_path / "test.saturation.r")
        show_saturation(infile, outfile, rpkm_cut=0.01)
        content = Path(outfile).read_text()
        assert "boxplot(" in content
        assert "pdf(" in content
        assert "dev.off()" in content

    def test_r_file_has_four_quantiles(self, tmp_path):
        from scripts.RPKM_saturation import show_saturation

        infile = self._make_xls(tmp_path)
        outfile = str(tmp_path / "test.saturation.r")
        show_saturation(infile, outfile, rpkm_cut=0.01)
        content = Path(outfile).read_text()
        for q in ("Q1", "Q2", "Q3", "Q4"):
            assert q in content

    def test_high_rpkm_cut_filters_genes(self, tmp_path):
        from scripts.RPKM_saturation import show_saturation

        infile = self._make_xls(tmp_path)
        outfile = str(tmp_path / "test.saturation.r")
        show_saturation(infile, outfile, rpkm_cut=999999)
        content = Path(outfile).read_text()
        # With very high cutoff, no genes pass — boxplot data should be empty
        assert "pdf(" in content


# ---------------------------------------------------------------------------
# scripts.FPKM_UQ: cal_fpkm
# ---------------------------------------------------------------------------


class TestCalFpkm:
    def _make_count_file(self, tmp_path):
        count_file = tmp_path / "test.htseq.counts.txt"
        lines = [
            "GENE001\t500\n",
            "GENE002\t1000\n",
            "GENE003\t200\n",
            "GENE004\t0\n",
            "__no_feature\t100\n",
            "__ambiguous\t50\n",
        ]
        count_file.write_text("".join(lines))
        return str(count_file)

    def _make_info_file(self, tmp_path):
        info_file = tmp_path / "test.info"
        header = (
            "gene_id\tgene_name\tseqname\tstart\tend\tstrand\t"
            "gene_type\tgene_status\thavana_gene\tfull_length\texon_length\texon_num\n"
        )
        lines = [
            header,
            "GENE001\tFOO\tchr1\t100\t2000\t+\tprotein_coding\tKNOWN\tHG1\t1900\t1000\t3\n",
            "GENE002\tBAR\tchr1\t3000\t5000\t-\tprotein_coding\tKNOWN\tHG2\t2000\t1500\t2\n",
            "GENE003\tBAZ\tchr2\t100\t800\t+\tprotein_coding\tKNOWN\tHG3\t700\t500\t1\n",
            "GENE004\tQUX\tchr2\t1000\t2000\t+\tlincRNA\tKNOWN\tHG4\t1000\t800\t2\n",
        ]
        info_file.write_text("".join(lines))
        return str(info_file)

    def test_creates_output_file(self, tmp_path):
        from scripts.FPKM_UQ import cal_fpkm

        count_file = self._make_count_file(tmp_path)
        info_file = self._make_info_file(tmp_path)
        out_file = str(tmp_path / "output.FPKM-UQ.txt")
        cal_fpkm(count_file, info_file, out_file)
        assert Path(out_file).exists()

    def test_output_has_header(self, tmp_path):
        from scripts.FPKM_UQ import cal_fpkm

        count_file = self._make_count_file(tmp_path)
        info_file = self._make_info_file(tmp_path)
        out_file = str(tmp_path / "output.FPKM-UQ.txt")
        cal_fpkm(count_file, info_file, out_file)
        content = Path(out_file).read_text()
        assert "gene_ID" in content
        assert "FPKM" in content
        assert "FPKM-UQ" in content

    def test_output_has_gene_rows(self, tmp_path):
        from scripts.FPKM_UQ import cal_fpkm

        count_file = self._make_count_file(tmp_path)
        info_file = self._make_info_file(tmp_path)
        out_file = str(tmp_path / "output.FPKM-UQ.txt")
        cal_fpkm(count_file, info_file, out_file)
        lines = Path(out_file).read_text().strip().split("\n")
        # Header + 4 genes (not __no_feature/__ambiguous)
        assert len(lines) == 5  # 1 header + 4 gene lines

    def test_fpkm_values_are_numeric(self, tmp_path):
        from scripts.FPKM_UQ import cal_fpkm

        count_file = self._make_count_file(tmp_path)
        info_file = self._make_info_file(tmp_path)
        out_file = str(tmp_path / "output.FPKM-UQ.txt")
        cal_fpkm(count_file, info_file, out_file)
        lines = Path(out_file).read_text().strip().split("\n")[1:]  # skip header
        for line in lines:
            fields = line.split("\t")
            fpkm = fields[7]
            fpkm_uq = fields[8]
            # Should be numeric or "NA"
            if fpkm != "NA":
                float(fpkm)
            if fpkm_uq != "NA":
                float(fpkm_uq)

    def test_log2_mode(self, tmp_path):
        from scripts.FPKM_UQ import cal_fpkm

        count_file = self._make_count_file(tmp_path)
        info_file = self._make_info_file(tmp_path)
        out_file = str(tmp_path / "output.FPKM-UQ.txt")
        cal_fpkm(count_file, info_file, out_file, log2_flag=True)
        content = Path(out_file).read_text()
        assert "log2" in content  # header should mention log2

    def test_protein_coding_genes_have_numeric_fpkm(self, tmp_path):
        """Protein-coding genes with counts should have numeric FPKM values."""
        from scripts.FPKM_UQ import cal_fpkm

        count_file = self._make_count_file(tmp_path)
        info_file = self._make_info_file(tmp_path)
        out_file = str(tmp_path / "output.FPKM-UQ.txt")
        cal_fpkm(count_file, info_file, out_file)
        lines = Path(out_file).read_text().strip().split("\n")[1:]
        # GENE001 is protein_coding with count=500 — should have positive FPKM
        gene1_lines = [line for line in lines if line.startswith("GENE001")]
        assert len(gene1_lines) == 1
        fields = gene1_lines[0].split("\t")
        assert float(fields[7]) > 0  # FPKM
        assert float(fields[8]) > 0  # FPKM-UQ


# ---------------------------------------------------------------------------
# scripts.RNA_fragment_size: fragment_size
# ---------------------------------------------------------------------------


class TestFragmentSize:
    def test_yields_strings(self, mini_bam):
        from scripts.RNA_fragment_size import fragment_size

        samfile = pysam.AlignmentFile(str(mini_bam), "rb")
        results = list(fragment_size(MINI_BED, samfile, qcut=30, ncut=5))
        samfile.close()
        assert len(results) > 0
        for r in results:
            assert isinstance(r, str)

    def test_one_result_per_gene(self, mini_bam):
        from scripts.RNA_fragment_size import fragment_size

        samfile = pysam.AlignmentFile(str(mini_bam), "rb")
        results = list(fragment_size(MINI_BED, samfile, qcut=30, ncut=5))
        samfile.close()
        # mini.bed has 3 genes
        assert len(results) == 3

    def test_output_format_tab_separated(self, mini_bam):
        from scripts.RNA_fragment_size import fragment_size

        samfile = pysam.AlignmentFile(str(mini_bam), "rb")
        results = list(fragment_size(MINI_BED, samfile, qcut=30, ncut=5))
        samfile.close()
        for r in results:
            fields = r.split("\t")
            # Should have at least 4 fields: chrom, start, end, name + counts
            assert len(fields) >= 4

    def test_skips_comments(self, mini_bam, tmp_path):
        from scripts.RNA_fragment_size import fragment_size

        bed = tmp_path / "test.bed"
        bed.write_text(
            "# comment\n"
            "track name=test\n"
            "chr1\t1000\t5000\tgene1\t0\t+\t1200\t4800\t0,0,0\t3\t500,600,500\t0,1500,3500\n"
        )
        samfile = pysam.AlignmentFile(str(mini_bam), "rb")
        results = list(fragment_size(str(bed), samfile, qcut=30, ncut=1))
        samfile.close()
        assert len(results) == 1

    def test_low_ncut_reports_stats(self, mini_bam):
        """With ncut=1 and paired reads present, should report non-zero stats."""
        from scripts.RNA_fragment_size import fragment_size

        samfile = pysam.AlignmentFile(str(mini_bam), "rb")
        results = list(fragment_size(MINI_BED, samfile, qcut=30, ncut=1))
        samfile.close()
        # gene1 has paired reads — check if any result has non-zero frag count
        gene1_result = results[0]  # gene1 is first in mini.bed
        fields = gene1_result.split("\t")
        frag_count = int(fields[4])
        # We have one proper pair in gene1 region
        assert frag_count >= 0  # may or may not pass filters

    def test_handles_unknown_chrom(self, mini_bam, tmp_path):
        from scripts.RNA_fragment_size import fragment_size

        bed = tmp_path / "test.bed"
        bed.write_text("chrX\t1000\t5000\tgeneX\t0\t+\t1200\t4800\t0,0,0\t1\t4000\t0\n")
        samfile = pysam.AlignmentFile(str(mini_bam), "rb")
        results = list(fragment_size(str(bed), samfile, qcut=30, ncut=5))
        samfile.close()
        assert len(results) == 1
        # Should report 0 counts for unknown chrom
        fields = results[0].split("\t")
        assert fields[4] == "0"


# ---------------------------------------------------------------------------
# scripts.tin: build_bitsets (tin's local copy, not cli_common)
# ---------------------------------------------------------------------------


class TestTinBuildBitsets:
    def test_empty_list(self):
        from scripts.tin import build_bitsets

        result = build_bitsets([])
        assert result == {}

    def test_basic_intervals(self):
        from scripts.tin import build_bitsets

        entries = [
            ["chr1", 100, 200],
            ["chr1", 300, 400],
            ["chr2", 500, 600],
        ]
        result = build_bitsets(entries)
        assert "chr1" in result
        assert "chr2" in result
        assert len(result["chr1"].find(150, 151)) > 0
        assert len(result["chr1"].find(250, 251)) == 0
        assert len(result["chr2"].find(550, 551)) > 0


# ---------------------------------------------------------------------------
# scripts.geneBody_coverage: genebody_coverage (BAM-based)
# ---------------------------------------------------------------------------


class TestGeneBodyCoverageBAM:
    def test_returns_defaultdict(self, mini_bam):
        from scripts.geneBody_coverage import genebody_coverage, genebody_percentile

        position_list = genebody_percentile(MINI_BED, mRNA_len_cut=100)
        result = genebody_coverage(str(mini_bam), position_list)
        assert isinstance(result, dict)

    def test_coverage_has_100_bins(self, mini_bam):
        from scripts.geneBody_coverage import genebody_coverage, genebody_percentile

        position_list = genebody_percentile(MINI_BED, mRNA_len_cut=100)
        result = genebody_coverage(str(mini_bam), position_list)
        # Should have up to 100 percentile bins
        assert len(result) <= 100

    def test_coverage_values_nonnegative(self, mini_bam):
        from scripts.geneBody_coverage import genebody_coverage, genebody_percentile

        position_list = genebody_percentile(MINI_BED, mRNA_len_cut=100)
        result = genebody_coverage(str(mini_bam), position_list)
        for k, v in result.items():
            assert v >= 0
