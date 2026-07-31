# De novo Iso-Seq Transcript Pipeline — Full Guide

This guide walks through the entire pipeline from raw PacBio data to a final,
usable set of gene-family transcripts, written for someone who has never run
any of this before. It covers **two different clustering approaches** we
tried (Cogent, and a CD-HIT-EST/MMseqs2 fallback), why we needed both, and
exactly how a successful Cogent run looks in practice (demonstrated on a
chemoreceptor-focused subset).

---

## 1. The big picture

**Starting point:** PacBio Kinnex long-read RNA sequencing data for an
insect species with no reference genome.

**Goal:** turn millions of individual sequencing reads into a manageable set
of "this is probably one real gene" transcript groups, so they can be
BLASTed and manually annotated in Geneious — replacing what Trinity used to
do for short-read data.

**The pipeline has two phases:**

1. **A shared first phase** (always the same, regardless of which clustering
   approach you use afterward): raw reads → primer removal → error
   correction → one FASTA file of clean, consensus transcripts.
2. **A clustering phase**, where you group those transcripts into gene
   families. We tried two different tools for this:
   - **Cogent** — reconstructs an actual gene structure per family. More
     biologically informative, but computationally expensive on very large,
     highly-repetitive datasets.
   - **MMseqs2** (with CD-HIT-EST-style settings) — a fast, pure
     sequence-similarity clustering approach. Less structurally aware, but
     dramatically faster and scales to millions of sequences without
     getting stuck.

We started with Cogent on the *full* dataset, hit a real scaling wall
(explained in Section 4), and switched to MMseqs2 for full-scale runs while
keeping Cogent for smaller, targeted subsets (like a curated chemoreceptor
gene list) where it works cleanly and gives richer output.

---

## 2. Naming conventions used throughout

Understanding these up front will make every later section easier to
follow.

| Pattern | Meaning |
|---|---|
| `<sample>.fq` | Raw exported FASTQ for one sample (e.g. `Elaabr3MA.fq`). **Do not use this as pipeline input** — see Section 3. |
| `<sample>.bam` inside `demultiplexed_files/` | The real starting point. Already past `skera` and `lima` (done by the sequencing core before you received the data). |
| `refine_<sample>.bam` | Output of `isoseq refine` — Full-Length Non-Chimeric (FLNC) reads with poly-A tails confirmed. |
| `clustered_<sample>.bam` / `.fasta` | Output of `isoseq cluster2` — consensus isoforms (sequencing errors corrected, redundant reads collapsed). This is the file every downstream clustering approach actually starts from. |
| `transcript/N` | The naming format for every consensus isoform (e.g. `transcript/5419`). This ID is preserved through every downstream step, so you can always trace a sequence back to its origin. |
| `is:i:N` (a BAM tag, not a filename) | Attached to every consensus isoform — tells you how many raw reads were collapsed into it. Higher = more confidently a real transcript, not noise. |
| `fullDB*` | MMseqs2's internal database files (binary, not human-readable). Always paired with a `.tsv` for the readable version. |
| `fullDB_clu.tsv` | Readable clustering result: two columns, `representative-transcript` and `member-transcript`. |
| `fullDB_clu_rep.fasta` | One representative sequence per cluster — the actual usable output of the MMseqs2 approach. |
| `isoseq_flnc.fasta` | Cogent's required input filename (we rename whatever FASTA we're feeding it to match this, since Cogent's scripts expect it). |
| `cogent2.renamed.fasta` (inside a family folder) | Cogent's reconstructed gene locus for that one family. This is the real output of the Cogent approach. |
| `pathN` (inside a Cogent header) | When a family has more than one biologically plausible reconstructed structure, each gets its own `pathN` label. Not a copy — a genuinely different possible gene structure. |
| `COGENT.DONE` | An empty marker file Cogent writes when a family finishes successfully. Used to check progress/completion. |

---

## 3. Environment setup (one-time)

Two separate conda environments are used:

- **`isoseq-denovo`** — PacBio's official tools (`lima`, `isoseq`, `skera`,
  `samtools`), plus `cd-hit`, `mmseqs2`, and `git`.
- **`cogent`** — Cogent itself, plus its older, pinned dependencies
  (Python 3.9, specific versions of `networkx`, `numpy`, `scikit-image`,
  and `parasail` — Cogent is an older tool and breaks with current default
  versions of these libraries).

The full, tested setup script lives in this repo as
`setup_isoseq_environment.sh`. Run it once on a fresh machine; it builds
both environments and applies all the version pins we discovered were
necessary through trial and error.

**Switching between environments:**
```bash
mamba activate isoseq-denovo   # for refine, cluster2, cd-hit-est, mmseqs2
mamba activate cogent          # for run_mash.py, process_kmer_to_graph.py, reconstruct_contig.py
```

---

## 4. Shared first phase: raw data → consensus transcripts

### 4.1 Confirm your starting point

**Important lesson learned:** always check whether `skera` and `lima` have
already been run on your data before assuming you need to run them
yourself. Check the BAM's processing history:

```bash
samtools view -H demultiplexed_files/<sample>.bam | grep "^@PG"
```

If you see `skera-split` and `lima`/`lima.1` entries in that output, this
work is already done — start at Step 4.2 directly. (In our case, the
sequencing core had already done this for every sample.)

### 4.2 `isoseq refine` — extract genuine, complete transcripts

**What it does:** looks at each read and checks three things — does it have
the correct primers in the right places, does it have a real poly-A tail,
and is it free of chimeric (accidentally-fused) artifacts? Only reads that
pass all three become "FLNC" (Full-Length Non-Chimeric) reads.

```bash
isoseq refine demultiplexed_files/<sample>.bam \
  IsoSeq_v2_primers_12.fasta \
  refine_<sample>.bam --require-polya
```

**Output:** `refine_<sample>.bam`, plus a `filter_summary.report.json`
showing exactly how many reads passed each filter. In our runs, pass rates
were consistently 99.9%+.

### 4.3 `isoseq cluster2` — build consensus isoforms

**What it does:** groups together FLNC reads that are near-identical copies
of the same molecule (accounting for sequencing errors), producing one
clean, error-corrected "consensus" sequence per real transcript.

```bash
isoseq cluster2 refine_<sample>.bam clustered_<sample>.bam -j 32 --singletons
```

**Critical flag:** `--singletons`. Without it, any transcript supported by
only 1 read is silently dropped from the output — which matters a lot if
you care about comparing high-confidence "master" transcripts against
low-confidence partial ones (see the `is:i:N` tag in Section 2).

**Output:** `clustered_<sample>.bam` — this is **the real starting point**
for every clustering approach described below. Convert it to FASTA:

```bash
samtools fasta clustered_<sample>.bam > clustered_<sample>.fasta
```

And build a weights file (mapping each transcript to its read-support
count, pulled from the `is:i:N` BAM tag):

```bash
samtools view clustered_<sample>.bam | awk '{
  match($0, /is:i:[0-9]+/);
  split(substr($0, RSTART, RLENGTH), a, ":");
  print $1"\t"a[3]
}' > clustered_<sample>.weights
```

---

## 5. Approach A: Cogent (structure-aware clustering)

Cogent doesn't just ask "are these sequences similar" — it reconstructs an
actual plausible gene locus from a family of transcripts, letting it
distinguish real alternative-splice isoforms from unrelated genes that just
happen to share some sequence. This matters a lot for large, closely-related
gene families (like chemosensory receptors).

### 5.1 The three Cogent steps

```bash
# Step 1: compare every transcript against every other one
run_mash.py isoseq_flnc.fasta --cpus <N>

# Step 2: sort transcripts into gene families based on that comparison
process_kmer_to_graph.py -c isoseq_flnc.weights \
  isoseq_flnc.fasta isoseq_flnc.fasta.s1000k30.dist \
  isoseq_flnc isoseq_flnc

# Step 3: reconstruct a gene locus for every family found
generate_batch_cmd_for_Cogent_reconstruction.py isoseq_flnc > reconstruction_cmds.txt
cat reconstruction_cmds.txt | xargs -P <N> -I {} bash -c '{}'
```

(Cogent's scripts expect the input file to be literally named
`isoseq_flnc.fasta` and the weights file `isoseq_flnc.weights` — rename or
symlink your real file to match before running.)

### 5.2 The scaling problem we hit (and why it matters)

`run_mash.py` compares **every sequence against every other sequence** —
the amount of work grows with the *square* of your sequence count. On our
full dataset (1.87 million transcripts), this step alone was measured to
take over 100 days unparallelized, and with careful use of `--cpus` and a
CD-HIT-EST redundancy pre-filter, we still only got it down to ~2.3 days.

Worse: `process_kmer_to_graph.py` (the family-sorting step) got
structurally stuck — **92% of all transcripts turned out to be one single,
massive, interconnected group**, confirmed independently two different ways
(a direct graph analysis of the raw distance file, and a repeated attempt
with `run_preCluster.py`'s pre-binning approach, which hit the same
86%-in-one-bin pattern). This is a real property of highly duplicated gene
families, not a bug — but it made the algorithm computationally impractical
at full genome-wide scale.

**The fix that worked:** `run_preCluster.py --dun_use_partial` (a
large-dataset pre-binning step, documented for >20,000-sequence inputs)
excludes "partial hit" alignments that were letting unrelated transcripts
chain together indirectly. This dropped the largest bin from 86% of the
data down to ~0.2%, and the full pipeline completed successfully at full
scale afterward — see `RUNBOOK.md` and `COMMAND_LOG.md` for the exact
troubleshooting history (including patching several Python 2-only code
paths in the legacy `cDNA_Cupcake` toolkit that `run_preCluster.py` depends
on).

### 5.3 Where Cogent works cleanly: smaller, targeted datasets

For a smaller, focused input (under a few thousand sequences — well under
the 20,000-sequence threshold where pre-binning becomes necessary), Cogent
runs the standard 3-step process directly, quickly, and without any of the
scaling issues above.

**Worked example — chemoreceptor gene subset (`Elaabr3MA_chemoreceptors.fasta`, 1,197 sequences):**

```bash
mkdir -p cogent_chemoreceptor_run
cd cogent_chemoreceptor_run

# Bring in the input, renamed to Cogent's expected filename
cp ~/Downloads/Elaabr3MA_chemoreceptors.fasta isoseq_flnc.fasta

# Build the matching weights file (only for the transcripts in this subset)
grep "^>" isoseq_flnc.fasta | sed 's/>//' > chemo_names.txt
samtools view <path_to>/clustered_Elaabr3MA.bam | awk 'NR==FNR{names[$1];next} $1 in names {
  match($0, /is:i:[0-9]+/);
  split(substr($0, RSTART, RLENGTH), a, ":");
  print $1"\t"a[3]
}' chemo_names.txt - > isoseq_flnc.weights

# Run all three Cogent steps
mamba activate cogent
run_mash.py isoseq_flnc.fasta --cpus 12
process_kmer_to_graph.py -c isoseq_flnc.weights isoseq_flnc.fasta \
  isoseq_flnc.fasta.s1000k30.dist isoseq_flnc isoseq_flnc
generate_batch_cmd_for_Cogent_reconstruction.py isoseq_flnc > reconstruction_cmds.txt
cat reconstruction_cmds.txt | xargs -P 12 -I {} bash -c '{}'
```

**Result:** 1,197 input transcripts → 178 gene families → 465 reconstructed
sequences (some families produced more than one plausible `pathN`
structure — genuinely different possible gene structures, not duplicates).
The entire run completed in well under a minute, with zero errors and no
repeat of the giant-cluster problem — because the dataset was small and
targeted rather than genome-wide.

**Note on hidden formatting issues:** if a FASTA file comes from a
Windows-originated source (e.g. a professor's exported file), check for
hidden `^M` carriage-return characters before running anything —
`head -3 file.fasta | cat -A` will reveal them. They silently break
exact-text matching (like our weights lookup) without producing any error.
Fix with `sed -i 's/\r$//' file.fasta`.

### 5.4 Combining the final Cogent output

```bash
find isoseq_flnc -name "cogent2.renamed.fasta" -exec cat {} \; > FINAL_reconstructed.fasta
```

This walks every family folder and combines all reconstructed sequences
into one file — the genuine final deliverable of a Cogent run.

---

## 6. Approach B: MMseqs2 (fast, sequence-similarity clustering)

For full-scale, genome-wide datasets where Cogent's structural approach
becomes impractical, MMseqs2 provides a fast, sequence-identity-based
alternative. Settings and rationale are documented in the official
[MMseqs2 User Guide](https://www.mmseqs.com/latest/userguide.pdf).

### 6.1 The four steps

```bash
mamba activate isoseq-denovo

# Step 1: convert FASTA into MMseqs2's internal database format
mmseqs createdb clustered_<sample>.fasta fullDB

# Step 2: cluster — settings chosen specifically to let short/partial
# transcripts attach to full-length ones without becoming representatives
# themselves (guide, "Target coverage" section)
mkdir -p tmp
mmseqs cluster fullDB fullDB_clu tmp \
  --cov-mode 1 -c 0.95 --cluster-mode 2 --min-seq-id 0.97 --threads 40

# Step 3: convert to a readable table
mmseqs createtsv fullDB fullDB fullDB_clu fullDB_clu.tsv

# Step 4: extract one representative sequence per cluster
mmseqs createsubdb fullDB_clu fullDB fullDB_clu_rep
mmseqs convert2fasta fullDB_clu_rep fullDB_clu_rep.fasta
```

**What `--cov-mode 1 --cluster-mode 2` actually does:** `--cluster-mode 2`
sorts sequences by length and builds each cluster around the *longest*
sequence, absorbing everything connected to it (this mirrors CD-HIT's
classic algorithm). `--cov-mode 1` requires alignments to cover a strong
fraction of the *shorter/target* sequence specifically — letting partial
fragments attach to full-length matches without needing to cover the full
length of the long sequence themselves.

### 6.2 Real results from both full-scale samples

| Sample | FLNC reads | Consensus isoforms | MMseqs2 clusters | Largest cluster |
|---|---|---|---|---|
| Elaabr3MA | 19,683,798 | 1,869,503 | 463,029 | 3,756 (0.2%) |
| Elaabr5FA | 19,230,031 | 2,926,280 | 667,280 | 3,090 (0.1%) |

Both runs completed in minutes (not days), with a healthy, well-distributed
cluster size spread — no repeat of the giant-blob problem Cogent hit at
this scale.

### 6.3 Compressing / renaming final output

```bash
cp fullDB_clu_rep.fasta fullDB_<sample>.fasta
gzip fullDB_<sample>.fasta
```

Produces `fullDB_<sample>.fasta.gz` — a compressed, clearly-named final
deliverable.

---

## 7. Deciding which approach to use

| Situation | Recommended approach |
|---|---|
| Small, curated/targeted gene set (hundreds to a few thousand sequences) | **Cogent** — full structural detail, fast at this scale |
| Full transcriptome, genome-wide (hundreds of thousands to millions of sequences) | **MMseqs2** — scales cleanly; use Cogent's `--dun_use_partial` pre-binning approach only if structural detail is essential and you can accept a much longer runtime |
| Uncertain / exploratory | Start with a small test subset in Cogent to validate biology, then decide based on full dataset size |

---

## 8. Where everything lives (folder reference)

```
workFlow/
├── original_input_file/          Untouched original data (BAM, FASTA)
├── original_pipeline_output/      refine + cluster2 outputs
├── cogent_option_a_abandoned/     First (unfixed) Cogent attempt on full data — kept for reference
├── cogent_chemoreceptor_run/       Successful small-scale Cogent run
├── mmseqs_full_run/                Full-scale MMseqs2 run (Elaabr3MA)
│   ├── cluster_step2/              Raw clustering results
│   └── step_3/                     Extracted representative sequences
├── mmseqs_5fa_run/                 Full-scale MMseqs2 run (Elaabr5FA)
│   ├── cluster_step2/
│   └── step_3/
└── isoseq_cdhit_grouping/          Earlier CD-HIT-EST-based grouping (superseded by MMseqs2 approach)
```

---

## 9. Known issues and fixes (quick reference)

| Issue | Fix |
|---|---|
| bioconda package is `pbskera`, not `skera` | Install `pbskera`; the command itself is still `skera` |
| Cogent import errors (numpy/networkx/scikit-image) | Pin: `networkx<2.7 numpy<2 scikit-image<0.20 parasail` |
| `isoseq refine` rejects FASTQ input | Always use BAM — check for a pre-processed BAM in `demultiplexed_files/` before assuming you need raw FASTQ |
| Hidden `^M` characters break exact-text matching | `sed -i 's/\r$//' file` on any externally-sourced FASTA/text file |
| `run_preCluster.py` not in bioconda `cdna_cupcake` package | Build from source, specifically the `tofu2_v21` git tag (not `master` — the script was removed from later versions) |
| `run_preCluster.py` Python 2 syntax errors | Patch with `lib2to3`, plus manual fixes for `cPickle`, integer division (`/` → `//`), and old-style `raise` statements |
| Duplicate background processes writing to the same file | Always run one command per paste — never paste multiple `nohup ... &` commands together; verify with `ps aux` before stepping away |
