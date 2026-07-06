# isoseq-denovo_pipeline
cat > README.md << 'EOF'
# De novo Iso-Seq transcriptome assembly pipeline

Replacing Trinity + manual Geneious annotation with a PacBio Iso-Seq
long-read pipeline for reference-genome-free transcript assembly,
built for the McKenna Lab.

## Why this exists

Previously, transcript assembly for gene annotation (odorant, gustatory,
and ionotropic receptor families) relied on Trinity assembling short
Illumina reads, followed by manual annotation in Geneious. This pipeline
replaces that with PacBio Kinnex long reads, which capture each full
transcript in a single read - removing the need to computationally
stitch fragments together.

## Pipeline stages

1. `skera` - segment Kinnex arrays into individual reads (S-reads)
2. `lima` - remove cDNA library primers
3. `isoseq refine` - trim poly-A tails, filter chimeric reads -> FLNC reads
4. `isoseq cluster2` - error-correct and cluster into consensus isoforms
5. `Cogent` - group isoforms into gene families via reconstructed loci
6. (later) `oarfish` - expression quantification

## Environment

Two separate conda environments:
- `isoseq-denovo` - lima, isoseq, skera, samtools
- `cogent` - Cogent, gmap, minimap2

See `setup_isoseq_environment.sh` for the full build.

## Progress log

- [x] Environment setup (both conda environments built and verified)
- [ ] First lima/refine run on real data
- [ ] isoseq cluster2 run
- [ ] Cogent family grouping
- [ ] Manual annotation in Geneious
EOF
