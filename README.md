# De novo Iso-Seq transcriptome assembly pipeline

Replacing Trinity + manual Geneious annotation with a PacBio Iso-Seq
long-read pipeline for reference-genome-free transcript assembly,
built for the Mitchell Lab.

## Why this exists

Previously, transcript assembly for gene annotation (odorant, gustatory,
and ionotropic receptor families) relied on Trinity assembling short
Illumina reads, followed by manual annotation in Geneious. This pipeline
replaces that with PacBio Kinnex long reads, which capture each full
transcript in a single read - removing the need to computationally
stitch fragments together.

## Pipeline stages

1. skera - segment Kinnex arrays into individual reads (S-reads)
2. lima - remove cDNA library primers
3. isoseq refine - trim poly-A tails, filter chimeric reads -> FLNC reads
4. isoseq cluster2 - error-correct and cluster into consensus isoforms
5. Cogent - group isoforms into gene families via reconstructed loci
6. (later) oarfish - expression quantification

## Environment

Two separate conda environments:
- isoseq-denovo - lima, isoseq, skera, samtools
- cogent - Cogent, gmap, minimap2

See setup_isoseq_environment.sh for the full build.

## Progress log

## Progress log

- [x] Environment setup (both conda environments built and verified, locally and on Roar Collab)
- [x] First lima/refine run on real data (Elaabr3MA: 19,683,798 FLNC reads from 19,696,797 full-length reads, 99.96% pass rate)
- [x] isoseq cluster2 run (1,869,503 consensus isoforms, --singletons included)
- [x] CD-HIT-EST redundancy pre-filter (1,869,503 -> 635,736 clusters, 66% redundancy removed)
- [x] Cogent run_mash.py full run (635,736 transcripts, ~24 hours, 11.7GB distance file)
- [ ] Cogent process_kmer_to_graph.py (family-finding) - in progress, one large subgraph taking longer than expected
- [ ] Cogent reconstruct_contig.py across all families
- [ ] Final alignment: map all transcripts back to reconstructed loci
- [ ] Manual annotation in Geneious

## Known issues / decisions log

- GPU (CUDA) acceleration investigated for mash/CD-HIT-EST/mmseqs2 clustering -
  not applicable; these tools lack real GPU code paths for this workload.
- mmseqs2 tested as a faster alternative to CD-HIT-EST - rejected due to
  unresolved coverage-mode discrepancies producing inconsistent cluster counts.
  Sticking with CD-HIT-EST (published precedent, Li et al. 2017, Cell Discovery).
- Roar Collab (PSU HPC) set up in parallel as a scaling option - blocked on
  account billing minutes (open/default account has ~2.86 credits, insufficient
  for a multi-day job). Needs PI account from Professor
