# De novo Iso-Seq transcriptome assembly pipeline

Replacing Trinity + manual Geneious annotation with a PacBio Iso-Seq
long-read pipeline for reference-genome-free transcript assembly, built
for the Mitchell Lab.

## Why this exists

Previously, transcript assembly for gene annotation (odorant, gustatory,
and ionotropic receptor families) relied on Trinity assembling short
Illumina reads, followed by manual annotation in Geneious. This pipeline
replaces that with PacBio Kinnex long reads, which capture each full
transcript in a single read - removing the need to computationally
stitch fragments together.

## Pipeline stages

1. `skera` - segment Kinnex arrays into individual reads (S-reads) -
   already done by the sequencing core for all samples used so far
2. `lima` - remove cDNA library primers - also already done by the
   sequencing core (confirmed via BAM `@PG` history)
3. `isoseq refine` - trim poly-A tails, filter chimeric reads -> FLNC reads
4. `isoseq cluster2` - error-correct and cluster into consensus isoforms
5. **Gene family clustering** - two approaches used depending on dataset
   scale (see below)
6. (later, per lab instruction) `oarfish` - expression quantification

## Gene family clustering: two approaches

See `PIPELINE_GUIDE.md` for the full walkthrough, explanations, and
worked examples of both.

- **Cogent** (`run_mash.py` -> `process_kmer_to_graph.py` ->
  `reconstruct_contig.py`) - reconstructs actual gene loci per family.
  Works cleanly and quickly on smaller, targeted datasets (validated on a
  1,197-transcript chemoreceptor subset -> 178 gene families, 465
  reconstructed sequences). At full genome-wide scale, this required a
  large-dataset pre-binning fix (`run_preCluster.py --dun_use_partial`)
  to avoid a structural bottleneck where the majority of transcripts
  formed one interconnected group - see `RUNBOOK.md` for the full
  investigation.
- **MMseqs2** (`createdb` -> `cluster` -> `createsubdb`/`convert2fasta`) -
  fast sequence-similarity clustering, used for full-scale, genome-wide
  runs. Settings (`--cov-mode 1 --cluster-mode 2`) chosen per the official
  MMseqs2 User Guide's recommendation for letting partial/fragment
  transcripts attach to full-length representatives.

## Environment

Two separate conda environments:

- `isoseq-denovo` - lima, isoseq, skera, samtools, cd-hit, mmseqs2
- `cogent` - Cogent, gmap, minimap2, py-spy

See `setup_isoseq_environment.sh` for the full build (includes the legacy
dependency version pins Cogent requires).

## Repo contents

- `README.md` - this file, high-level overview
- `PIPELINE_GUIDE.md` - complete, beginner-friendly step-by-step guide to
  every stage, including naming conventions and a full output/folder
  reference
- `RUNBOOK.md` - exact commands and real results for the first full-scale
  sample
- `COMMAND_LOG.md` - full chronological command history across the whole
  project, including dead ends and fixes
- `REFERENCES.md` - credits and citations for Cogent, MMseqs2, CD-HIT-EST,
  and the published methods that informed this pipeline's design
- `mmseqs2_userguide.pdf` - official MMseqs2 documentation, provided for
  reference alongside `REFERENCES.md`
- `setup_isoseq_environment.sh` - full environment build script

## Progress log

- [x] Environment setup (both conda environments built and verified,
      locally and on Roar Collab)
- [x] isoseq refine + cluster2 completed for two samples:
  - Elaabr3MA: 19,683,798 FLNC reads -> 1,869,503 consensus isoforms
  - Elaabr5FA: 19,230,031 FLNC reads -> 2,926,280 consensus isoforms
- [x] Cogent - full pipeline validated successfully on a targeted
      chemoreceptor subset (1,197 transcripts -> 178 gene families -> 465
      reconstructed sequences)
- [x] Cogent - full genome-wide scaling issue investigated and resolved
      (`--dun_use_partial` pre-binning fix); documented in `RUNBOOK.md`
- [x] MMseqs2 clustering completed at full scale for both samples:
  - Elaabr3MA: 463,029 gene family clusters
  - Elaabr5FA: 667,280 gene family clusters
- [ ] Additional sample(s) - in progress
- [ ] Add read-support (`count=X`) labels to final representative FASTA
      headers
- [ ] BLAST final gene family representatives against reference database
- [ ] Manual annotation in Geneious

## Known issues / decisions log

- GPU (CUDA) acceleration investigated for mash/CD-HIT-EST/mmseqs2
  clustering - not applicable for our workload; these tools either lack
  real GPU code paths for this use case, or (in mmseqs2's case) only
  support GPU acceleration for protein search, not nucleotide clustering.
- Cogent's family-finding step hit a real structural wall at full genome
  scale (92% of transcripts in one interconnected group, confirmed two
  independent ways). Resolved using `run_preCluster.py --dun_use_partial`
  (large-dataset pre-binning), which required building `cDNA_Cupcake`
  from source (`tofu2_v21` tag) and patching several Python 2-only code
  paths - full details in `RUNBOOK.md` and `COMMAND_LOG.md`.
- CD-HIT-EST and MMseqs2 both tested as fast, sequence-similarity-based
  alternatives to Cogent for full-scale clustering. Settled on MMseqs2
  with `--cov-mode 1 --cluster-mode 2` (documented in the official
  MMseqs2 User Guide for exactly this fragment-absorption use case) as
  the production approach for full genome-wide datasets.
- Roar Collab (PSU HPC) set up in parallel as a scaling option - blocked on account billing minutes (default "open" account has ~2.86 credits,
  insufficient for a multi-day job). Would need a PI-sponsored account
  from Professor Mitchell to use for future large-scale runs.

  ___________________________________________________________________________________________________
  Gemini and Claude.ai were used for debugging and understanding command prompt structures. 
