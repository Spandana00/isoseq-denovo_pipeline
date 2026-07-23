# Pipeline runbook - Elaabr3MA sample

Exact commands used, in order, with real results from this run.

## 1. isoseq refine

    isoseq refine demultiplexed_files/Elaabr3MA.bam IsoSeq_v2_primers_12.fasta \
      refine_Elaabr3MA.bam --require-polya

Result: 19,683,798 FLNC reads (99.96% pass rate from 19,696,797 full-length reads)

## 2. isoseq cluster2

    isoseq cluster2 refine_Elaabr3MA.bam clustered_Elaabr3MA.bam -j 32 --singletons

Note: --singletons is required, or isoforms with only 1 supporting read are
silently dropped from the output.

Result: 1,869,503 consensus isoforms

## 3. Convert to FASTA (Cogent requires FASTA, not BAM)

    samtools fasta clustered_Elaabr3MA.bam > clustered_Elaabr3MA.fasta

## 4. Build weights file (maps transcript ID -> supporting read count)

    samtools view clustered_Elaabr3MA.bam | awk '{
      match($0, /is:i:[0-9]+/);
      split(substr($0, RSTART, RLENGTH), a, ":");
      print $1"\t"a[3]
    }' > clustered_Elaabr3MA.weights

Note: the "is:i:N" BAM tag holds the read-support count (equivalent to the
"count=N" field mentioned in FASTA-based workflows).

## 5. CD-HIT-EST redundancy pre-filter

    cd-hit-est -i clustered_Elaabr3MA.fasta -o clustered_Elaabr3MA.cdhit.fasta \
      -c 0.99 -G 0 -aL 0.00 -aS 0.99 -AS 30 -M 0 -d 0 -p 1 -T 40

Parameters from Li et al. 2017 (Cell Discovery), same combination used for
a reference-genome-free Iso-Seq transcriptome.

Result: 1,869,503 -> 635,736 clusters (66% redundancy removed), ~90 min at
full scale (40 threads).

## 6. Rebuild weights file for the reduced set

    grep "^>" clustered_Elaabr3MA.cdhit.fasta | sed 's/>//' > cdhit_names.txt
    samtools view clustered_Elaabr3MA.bam | awk 'NR==FNR{names[$1];next} $1 in names {
      match($0, /is:i:[0-9]+/);
      split(substr($0, RSTART, RLENGTH), a, ":");
      print $1"\t"a[3]
    }' cdhit_names.txt - > clustered_Elaabr3MA.cdhit.weights

## 7. Cogent run_mash.py

    run_mash.py clustered_Elaabr3MA.cdhit.fasta --cpus 40

Result: ~24 hours on 635,736 sequences (636 chunks, ~202,000 pairwise
comparisons), 11.7GB distance file.

## 8. Cogent process_kmer_to_graph.py

    process_kmer_to_graph.py -c clustered_Elaabr3MA.cdhit.weights \
      clustered_Elaabr3MA.cdhit.fasta clustered_Elaabr3MA.cdhit.fasta.s1000k30.dist \
      cogent_full cogent_full

Status: in progress. Graph has 601,527 nodes / 22,963,944 edges. No --cpus
flag on this step - single-threaded, runtime unpredictable. One large
subgraph has taken 4+ hours without resolving as of this writing. If this
does not converge within ~24h, fallback plan: rerun with a stricter
--sim_threshold (tested default 0.05; considering 0.15) to force smaller,
faster-converging subgraphs.

## Environment fixes required (see setup_isoseq_environment.sh)

- bioconda package name is pbskera, not skera (binary is still called skera)
- Cogent requires: networkx<2.7, numpy<2, scikit-image<0.20, parasail
  (installed separately - not declared correctly in Cogent's own setup.py)
- git required manual install into base conda env on this workstation
- lima/isoseq refine require BAM/FASTA/FASTQ (isoseq refine specifically
  rejects plain FASTQ with a ".fastq" extension issue - use BAM)

