#!/usr/bin/env bash
# De novo Iso-Seq transcript workflow - environment setup
# Two conda environments: isoseq-denovo (PacBio toolkit) and cogent (Cogent + its aligners)

set -euo pipefail

# 1. PacBio Iso-Seq toolkit: lima, isoseq (refine/cluster2), skera, samtools
mamba create -y -n isoseq-denovo -c bioconda -c conda-forge lima isoseq pbskera samtools

conda run -n isoseq-denovo lima --version
conda run -n isoseq-denovo isoseq --version
conda run -n isoseq-denovo skera --version

# 2. Cogent - separate environment (older pinned dependencies: Python 3.9, gmap for alignment)
mamba create -y -n cogent -c bioconda -c conda-forge python=3.9 pip gmap minimap2 samtools

cd ~
git clone https://github.com/Magdoll/Cogent.git
cd Cogent
conda run -n cogent pip install .
cd ~
conda run -n cogent pip install "networkx<2.7" "numpy<2" "scikit-image<0.20" parasail --force-reinstall

conda run -n cogent python -c "import Cogent; print('Cogent import OK')"

echo "Setup complete."
echo "  mamba activate isoseq-denovo  -> lima, isoseq, skera, samtools"
echo "  mamba activate cogent         -> Cogent"

