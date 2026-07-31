# References and Credits

This pipeline relies heavily on the following open-source tools and their
documentation. Full credit to the original authors — this repo only
documents how we applied their work to our specific dataset.

## PacBio Iso-Seq toolkit

- `lima`, `isoseq` (refine, cluster2), `skera` — PacBio Biosciences.
  https://github.com/PacificBiosciences

## Cogent

Reference-free reconstruction of gene loci from Iso-Seq transcripts, used
for structurally-aware gene family clustering (Section 5 of
`PIPELINE_GUIDE.md`).

- **Repository:** https://github.com/Magdoll/Cogent
- **Author:** Elizabeth Tseng (Magdoll), PacBio Biosciences
- **Companion toolkit:** https://github.com/Magdoll/cDNA_Cupcake
  (specifically the `run_preCluster.py` large-dataset pre-binning script,
  built from the `tofu2_v21` tag — see `RUNBOOK.md` for why this specific
  tag was required)

## MMseqs2

Fast sequence clustering, used for full-scale gene family grouping when
Cogent's approach became computationally impractical (Section 6 of
`PIPELINE_GUIDE.md`).

- **Repository:** https://github.com/soedinglab/MMseqs2
- **User Guide (included in this repo as `mmseqs2_userguide.pdf`):**
  https://www.mmseqs.com/latest/userguide.pdf
- **Authors:** Martin Steinegger, Milot Mirdita, Eli Levy Karin, Lars von
  den Driesch, Clovis Galiez, Johannes Söding
- **Citation:** Steinegger M and Söding J. MMseqs2 enables sensitive
  protein sequence searching for the analysis of massive data sets.
  *Nature Biotechnology*, 2017.

The specific `--cov-mode 1` / `--cluster-mode 2` clustering settings used
throughout this pipeline (letting partial/fragment transcripts attach to
full-length representatives) are documented in the guide's "Target
coverage" and "Clustering modes" sections — see `PIPELINE_GUIDE.md`
Section 6.1 for the exact reasoning.

## CD-HIT-EST

Used as an earlier redundancy pre-filter approach before switching fully to
MMseqs2 for full-scale clustering.

- **Repository:** https://github.com/weizhongli/cdhit
- Parameter combination (`-c`, `-G 0`, `-aL`, `-aS`) informed by: Li A. et
  al., "Long read reference genome-free reconstruction of a full-length
  transcriptome from Astragalus membranaceus," *Cell Discovery*, 2017 —
  the same redundancy-pre-filter-before-family-reconstruction strategy
  applied here.
