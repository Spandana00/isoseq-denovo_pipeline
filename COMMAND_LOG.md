# Command log

Chronological record of every command actually used to set up and run this
pipeline, across both the local lab workstation and Roar Collab. Kept as a
reference for reproducing this setup, and as a record of what was tried
(including dead ends) for the writeup to Professor Mitchell.

---

## 1. Environment setup (local workstation)

```bash
# Install Miniforge (conda + mamba)
wget -O ~/miniforge.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash ~/miniforge.sh -b -p "$HOME/miniforge3"
mamba shell init
source ~/.bashrc

# Core PacBio Iso-Seq toolkit
mamba create -y -n isoseq-denovo -c bioconda -c conda-forge lima isoseq pbskera samtools
# NOTE: bioconda package name is pbskera, not skera (binary is still `skera`)

conda run -n isoseq-denovo lima --version
conda run -n isoseq-denovo isoseq --version
conda run -n isoseq-denovo skera --version

# cd-hit and mmseqs2 (redundancy pre-filter tools)
mamba install -y -n isoseq-denovo -c bioconda -c conda-forge cd-hit mmseqs2

# git (needed for cloning Cogent - not preinstalled on this workstation)
mamba install -y -n base -c conda-forge git

# Cogent - separate environment (legacy pinned dependencies)
mamba create -y -n cogent -c bioconda -c conda-forge python=3.9 pip gmap minimap2 samtools

cd ~
git clone https://github.com/Magdoll/Cogent.git
cd Cogent
pip install .

# Fix legacy dependency versions (discovered via trial and error - see Section 5)
pip install "networkx<2.7" "numpy<2" "scikit-image<0.20" parasail --force-reinstall

python -c "import Cogent; print('Cogent import OK')"
```

---

## 2. GitHub repository setup

```bash
mamba install -y -n base -c conda-forge git   # git only in cogent env at first, fixed by installing to base

cd ~
git clone https://github.com/Spandana00/isoseq-denovo_pipeline.git
cd isoseq-denovo_pipeline

git config --global user.name "Spandana00"
git config --global user.email "spandana.b2007@gmail.com"
git config --global credential.helper store
git config --global pull.rebase false

# README and setup script written via nano (not heredoc - heredoc caused
# repeated issues with GitHub's web editor pasting literal `cat <<EOF` text)
nano README.md
nano setup_isoseq_environment.sh

git add .
git commit -m "Add setup script and update README with pipeline overview"
git push origin main
```

---

## 3. Locating and diagnosing the raw data

```bash
# Finding the file across mounted drives (permission issues on prof's drive)
sudo find / -iname "Elaabr3MA*" 2>/dev/null

# Copied to own account to avoid repeated sudo/auth prompts
# (done via GUI file browser, into ~/Desktop/workFlow/)

head -4 ~/Desktop/workFlow/Elaabr3MA.fq   # confirmed real CCS read, poly-A intact

# Confirming these were NOT raw - lima had already been run by the
# sequencing core (found via the demultiplexed_files/ folder + README from
# the sequencing office, and confirmed via @PG history in the BAM header)
samtools view -H demultiplexed_files/Elaabr3MA.bam | grep "^@PG"
```

---

## 4. Core pipeline run - sample Elaabr3MA (first pass, small-scale validation)

```bash
mamba activate isoseq-denovo
cd ~/Desktop/workFlow

# isoseq refine (on the already-demultiplexed BAM from the sequencing core)
isoseq refine demultiplexed_files/Elaabr3MA.bam IsoSeq_v2_primers_12.fasta \
  refine_Elaabr3MA.bam --require-polya
# Result: 19,683,798 FLNC reads (99.96% pass rate)

# isoseq cluster2 - NOTE: --singletons required or count=1 isoforms are
# silently dropped
isoseq cluster2 refine_Elaabr3MA.bam clustered_Elaabr3MA.bam -j 32 --singletons
# Result: 1,869,503 consensus isoforms

samtools view -c clustered_Elaabr3MA.bam
```

---

## 5. Cogent - small-scale test (1,000 transcripts, 176 families)

```bash
# Convert BAM to FASTA (Cogent requires FASTA, not BAM)
samtools fasta clustered_Elaabr3MA.bam > clustered_Elaabr3MA.fasta

# Build weights file from the is:i:N BAM tag
samtools view clustered_Elaabr3MA.bam | awk '{
  match($0, /is:i:[0-9]+/);
  split(substr($0, RSTART, RLENGTH), a, ":");
  print $1"\t"a[3]
}' > clustered_Elaabr3MA.weights

# Test subset
head -2000 clustered_Elaabr3MA.fasta > cogent_test.fasta

mamba activate cogent
run_mash.py cogent_test.fasta
# Hit missing dependency errors here - fixed with the pip pins in Section 1

process_kmer_to_graph.py -c cogent_test.weights cogent_test.fasta \
  cogent_test.fasta.s1000k30.dist cogent_test cogent_test
# Result: 176 gene family partitions

reconstruct_contig.py cogent_test/cogent_test_1
# Result: success, clean reconstructed locus

# Batch-run reconstruction across all 176 test partitions
generate_batch_cmd_for_Cogent_reconstruction.py cogent_test > cogent_test_cmds.txt
cat cogent_test_cmds.txt | xargs -P 32 -I {} bash -c "{}"
find cogent_test -name "COGENT.DONE" | wc -l
# Result: 176/176 successful
```

---

## 6. Scale investigation (50k-transcript benchmarks)

```bash
head -100000 clustered_Elaabr3MA.fasta > cogent_test50k.fasta   # 50,000 transcripts

# Single-threaded baseline
time run_mash.py cogent_test50k.fasta
# ~145 min projected -> 100+ days at full scale

# Parallelized
time run_mash.py cogent_test50k.fasta --cpus 32
# 6m24s -> ~6.3 days at full scale

# CD-HIT-EST pre-filter, then parallelized mash
cd-hit-est -i cogent_test50k.fasta -o cogent_test50k.cdhit.fasta \
  -c 0.99 -G 0 -aL 0.00 -aS 0.99 -AS 30 -M 0 -d 0 -p 1 -T 32
# 50,000 -> 29,498 clusters

time run_mash.py cogent_test50k.cdhit.fasta --cpus 32
# 2m19s -> ~2.3 days at full scale (chosen approach)
```

### mmseqs2 investigated as CD-HIT-EST alternative (not adopted)

```bash
mmseqs createdb cogent_test50k.fasta cogent_test50k_db --createdb-mode 2
mmseqs easy-cluster cogent_test50k.fasta mmseqs_gpu_out mmseqs_tmp \
  --min-seq-id 0.99 -c 0.99 --cov-mode 1 --gpu 1
# 12.87s, 33,998 clusters - but GPU confirmed NOT actually engaged
# (verified via nvidia-smi monitoring during the run)

mmseqs easy-cluster cogent_test50k.fasta mmseqs_fixed_out mmseqs_tmp \
  --min-seq-id 0.99 -c 0.99 --cov-mode 2
# 15.56s, 38,157 clusters
# Neither cov-mode reproduced CD-HIT-EST's -aS/-aL short-fragment-absorption
# behavior reliably. Decision: stick with CD-HIT-EST (published precedent).
```

### GPU investigation (not applicable)

```bash
nvidia-smi   # confirmed RTX 4500 Ada present, 24GB, mostly idle
# Researched: mash, CD-HIT-EST, minimap2, CBC solver = no CUDA support.
# mmseqs2 --gpu exists but only implemented for `search`, not `cluster`,
# and built for protein search, not nucleotide clustering. Not usable here.
```

---

## 7. Roar Collab (PSU HPC) setup - blocked on billing

```bash
ssh sjb7561@submit.hpc.psu.edu

# Storage check
ls -la /storage/work/$USER
sinfo   # partitions: basic, standard, himem, interactive, mgc-* (restricted)

# Miniforge + environments (same as Section 1, installed to /storage/work)
cd ~/work
wget -O miniforge.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash miniforge.sh -b -p ~/work/miniforge3
~/work/miniforge3/bin/mamba shell init
source ~/.bashrc

git clone https://github.com/Spandana00/isoseq-denovo_pipeline.git
cd isoseq-denovo_pipeline
bash setup_isoseq_environment.sh   # ran clean - all fixes already baked in

# Data transfer (from local workstation, not Roar)
# run from ~/Desktop/workFlow on the LAB machine:
rsync -avP clustered_Elaabr3MA.bam sjb7561@submit.hpc.psu.edu:/storage/work/sjb7561/isoseq-denovo_pipeline/

# Regenerated FASTA + weights on Roar
samtools fasta clustered_Elaabr3MA.bam > clustered_Elaabr3MA.fasta
grep "^>" clustered_Elaabr3MA.fasta | sed 's/>//' > all_transcript_names.txt
samtools view clustered_Elaabr3MA.bam | awk '{
  match($0, /is:i:[0-9]+/);
  split(substr($0, RSTART, RLENGTH), a, ":");
  print $1"\t"a[3]
}' > clustered_Elaabr3MA.weights

# Slurm batch job (cogent_full_job.slurm) - submitted but never ran
sbatch cogent_full_job.slurm
squeue -u sjb7561
# Status: PD (pending), reason AssocGrpBillingMinutes
# Account "open" has ~2.86 credits - insufficient for a multi-day job.
# BLOCKED - needs a PI-sponsored account from Professor Mitchell.
```

---

## 8. Full-scale local run (final approach)

```bash
mamba install -y -n base -c conda-forge tmux   # tmux behaved unreliably on
                                                 # this system - switched to nohup

cd ~/Desktop/workFlow
mamba activate isoseq-denovo

# Full-scale CD-HIT-EST (all 1,869,503 transcripts)
nohup cd-hit-est -i clustered_Elaabr3MA.fasta -o clustered_Elaabr3MA.cdhit.fasta \
  -c 0.99 -G 0 -aL 0.00 -aS 0.99 -AS 30 -M 0 -d 0 -p 1 -T 40 > cdhit_full.log 2>&1 &
# Result: 1,869,503 -> 635,736 clusters (66% redundancy removed), ~90 min

grep -c "^>" clustered_Elaabr3MA.cdhit.fasta

# Weights file for reduced set
grep "^>" clustered_Elaabr3MA.cdhit.fasta | sed 's/>//' > cdhit_names.txt
samtools view clustered_Elaabr3MA.bam | awk 'NR==FNR{names[$1];next} $1 in names {
  match($0, /is:i:[0-9]+/);
  split(substr($0, RSTART, RLENGTH), a, ":");
  print $1"\t"a[3]
}' cdhit_names.txt - > clustered_Elaabr3MA.cdhit.weights

# Disable sleep before leaving the machine unattended
gsettings set org.gnome.settings-daemon.plugins.power sleep-inactive-ac-type 'nothing'
gsettings set org.gnome.desktop.session idle-delay 0

# Full-scale run_mash.py
mamba activate cogent
nohup run_mash.py clustered_Elaabr3MA.cdhit.fasta --cpus 40 > mash_full.log 2>&1 &
# Result: ~24 hours, 636 chunks, 11.7GB distance file, completed successfully
```

---

## 9. Family-finding troubleshooting (in progress)

```bash
# First attempt - default --sim_threshold (0.05) - got stuck
nohup process_kmer_to_graph.py \
  -c clustered_Elaabr3MA.cdhit.weights \
  clustered_Elaabr3MA.cdhit.fasta \
  clustered_Elaabr3MA.cdhit.fasta.s1000k30.dist \
  cogent_full cogent_full > graph_full.log 2>&1 &
# Graph: 601,527 nodes, 22,963,944 edges
# Stuck on 1 subgraph for 20+ hours, memory grew to 92.8GB - killed

kill 2197855

# Second attempt - stricter --sim_threshold (0.15)
nohup process_kmer_to_graph.py \
  --sim_threshold 0.15 \
  -c clustered_Elaabr3MA.cdhit.weights \
  clustered_Elaabr3MA.cdhit.fasta \
  clustered_Elaabr3MA.cdhit.fasta.s1000k30.dist \
  cogent_full_v2 cogent_full_v2 > graph_full_v2.log 2>&1 &
# Graph: 574,394 nodes, 14,661,980 edges (fewer edges, as expected)
# Memory stable this time (~16GB) but still stuck on 1 subgraph after 1.5+ hrs

# Diagnostics run in parallel (did not affect the running job):

# 1. py-spy - inspect the running process's call stack directly
mamba install -y -n cogent -c conda-forge py-spy
sudo env "PATH=$PATH" py-spy dump --pid 2566304
# Confirmed: genuinely inside recursive ncut subdivision + scipy eigsh
# (eigenvalue solver) - real computation, not frozen

# 2. Independent connected-components check on the raw distance file
python3 -c "
import networkx as nx
G = nx.Graph()
with open('clustered_Elaabr3MA.cdhit.fasta.s1000k30.dist') as f:
    for i, line in enumerate(f):
        parts = line.strip().split()
        if len(parts) >= 3 and float(parts[2]) < 0.15:
            G.add_edge(parts[0], parts[1])
        if i % 500000 == 0:
            print(f'{i} lines processed', flush=True)
print('Total nodes:', G.number_of_nodes())
components = sorted(nx.connected_components(G), key=len, reverse=True)
print('Number of connected components:', len(components))
print('Largest component size:', len(components[0]))
print('Second largest:', len(components[1]) if len(components) > 1 else 0)
"
# RESULT: 635,736 total nodes, 24,542 connected components,
#         largest component = 586,581 (92% of all transcripts),
#         second largest = only 60
#
# CONCLUSION: 92% of transcripts fall into one massive interconnected
# component even at the stricter 0.15 threshold. This is a real, likely
# biological finding, not a bug - needs Professor Mitchell's input on how
# to proceed (further threshold tightening vs. accepting this structure
# vs. other approaches).
```
