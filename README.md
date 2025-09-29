# Genomic diversity of SARS-CoV-2 in Pakistan during the fourth wave of pandemic

Reproducible code that mirrors the Journal of Medical Virology article (2022). DOI: 10.1002/jmv.27957

## Program summary
This repository provides an end to end pipeline that matches the study design. Steps run under one Snakemake workflow and can be executed on a workstation or small server.

1) Inputs
   - Paired end FASTQ from amplicon libraries sequenced on Illumina iSeq 2x150.

2) Quality control and trimming
   - FastQC for initial QC.
   - Trimmomatic with Q30 sliding window and minimum length 50.

3) Reference based mapping and consensus
   - Fetch Wuhan Hu 1 reference (NC_045512.2) or provide a local FASTA.
   - Map with BWA MEM, sort and index with SAMtools. Duplicate marking with Picard is optional.
   - Call variants with bcftools. Mask sites below a depth threshold and write a consensus per sample.

4) Lineage assignment
   - Pangolin on combined consensus to assign PANGO lineages.

5) Phylogeny and context
   - Fetch context sequences by accession using Entrez. If you have GISAID, place exported FASTA into `data-private/` and set paths in config.
   - MAFFT alignment and IQ-TREE with 1000 ultrafast bootstraps. Optional ModelFinder can be toggled.
   - Optional Nextstrain build with Augur and Auspice.

6) Outputs
   - Per sample consensus: `results/consensus/<sample>.fa`
   - Combined consensus: `results/consensus/all_consensus.fasta`
   - Pangolin lineages: `results/pangolin/lineage_report.csv`
   - Alignment: `results/aln/wg_alignment.fasta`
   - Phylogeny: `results/iqtree/wg.treefile`

7) Repro and compliance
   - MIT LICENSE and CITATION.cff included.
   - Do not commit clinical or restricted data. Set NCBI_EMAIL for Entrez calls.

## Requirements
- Python 3.11 or newer
- Option A: pip and virtualenv
- Option B: conda or mamba
- Snakemake for the full workflow
- Nextstrain tools for optional time resolved builds

### NCBI usage note
Set a contact email once per shell for E-utilities. Optional API key improves rate limits.
```bash
export NCBI_EMAIL="you@example.com"
export NCBI_API_KEY="xxxxxxxxxxxxxxxxxxxxxxxxxxxx"   # optional
```

## Quick verification
```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r env/requirements.txt
python analysis/scripts/example_qc_plot.py --in data-example/example_counts.tsv --out results-example/example_plot.png
```

## One command end to end run
```bash
export NCBI_EMAIL="you@example.com"
conda env create -f env/environment.yml
conda activate sarscov2-fourthwave-env
snakemake -s workflow/Snakefile -c 4 --printshellcmds
```

## Configuration
Edit `config/config.yaml`. Example:
```yaml
pairs:
  - sample: ISL_001
    r1: data-private/ISL_001_R1.fastq.gz
    r2: data-private/ISL_001_R2.fastq.gz

reference:
  acc: NC_045512.2

context_accessions:
  - NC_045512.2
  - MW599079.1

pangolin:
  enable: true

augur:
  enable: true

params:
  threads: 4
  min_depth_consensus: 10
  min_qual: 20
  iqtree_model_wg: GTR+G+I
  bootstrap: 1000
  use_model_finder: false
  filter_min_length: 29500
  subsampling: 250
  clock_rate: 0.0008
  clock_std_dev: 0.0004
  root: Wuhan/Hu-1/2019
```

## Nextstrain build
You can produce an Auspice JSON for interactive viewing.

1. Prepare `data-private/metadata.tsv` with columns: `strain`, `date`, `region`, `division`, `location`. Strain must match FASTA headers.
2. Run:
   ```bash
   make nextstrain
   ```
3. Start Auspice and open `results/nextstrain/auspice/sarscov2_fourth_wave.json`.

## How to cite
- Paper: Umair M, Ikram A, Rehman Z, Haider SA, Ammar M, Badar N, Ali Q, Rana MS, Salman M. Genomic diversity of SARS-CoV-2 in Pakistan during the fourth wave of pandemic. Journal of Medical Virology. 2022. https://doi.org/10.1002/jmv.27957
- Software: Haider SA. SARS-CoV-2 fourth wave Pakistan analysis. Version 1.0.0. GitHub repository. Include commit hash when available.

## References
- Andrews S. 2010. FastQC. Babraham Bioinformatics.
- Bolger AM, Lohse M, Usadel B. 2014. Trimmomatic. Bioinformatics 30:2114-2120.
- Li H. 2013. BWA-MEM. arXiv:1303.3997.
- Li H, et al. 2009. SAMtools. Bioinformatics 25:2078-2079.
- Danecek P, et al. 2021. BCFtools. GigaScience 10:giab008.
- Katoh K, Standley DM. 2013. MAFFT. Mol Biol Evol 30:772-780.
- Minh BQ, et al. 2020. IQ-TREE 2. Mol Biol Evol 37:1530-1534.
- O'Toole Á, et al. 2021. Pangolin. Virus Evol 7:veab064.
- Huddleston J, et al. 2021. Augur. JOSS 6:2906.
- Sagulenko P, Puller V, Neher RA. 2018. TreeTime. Virus Evol 4:vex042.
- Shu Y, McCauley J. 2017. GISAID. Euro Surveill 22:30494.
- Camacho C, et al. 2009. BLAST+. BMC Bioinformatics 10:421.
- Kans J. Entrez E-utilities Help. NCBI.
- Köster J, Rahmann S. 2012. Snakemake. Bioinformatics 28:2520-2522.
