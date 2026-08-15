# DogMAG

DogMAG is a reusable canine faecal metagenome resource comprising long-read
assemblies, quality-controlled and dereplicated metagenome-assembled genome
(MAG) catalogues, taxonomic metadata, viral/proviral candidate metadata, and
reproducible workflow scripts.

This repository accompanies the DogMAG Scientific Data manuscript. It contains
the manuscript, figures, supplementary tables, public accession mappings, and
the scripts required to reproduce the reported processing and summaries.

## Key catalogue outputs

- 30,556 polished or reassembled candidate versions evaluated by BASALT.
- 11,276 BASALT-selected bin/version records.
- 3,418 reselected medium-quality-or-better MAG candidates.
- 792 representatives at 99% ANI and 135 representatives at 95% ANI.
- GTDB-Tk taxonomy for all 792 strain-like representatives.
- 22,068 viral/proviral candidate rows before final quality filtering.

The sequence files themselves are not stored in Git. Raw reads, assemblies,
MAG FASTAs, and deposited genome records are available through the ENA projects
and article data record described in the manuscript and `accessions/`.

## Repository map

- `manuscript/`: reconciled Markdown and Word manuscript.
- `figures/`: publication figures and the source tables used for Figures 5 and 6.
- `supplementary_tables/`: the nine manuscript supplementary tables.
- `scripts/`: public workflow and reporting scripts grouped by analysis stage.
- `accessions/`: public ENA sample, read, and MAG accession mappings.
- `docs/DATA_SHARING_SCOPE.md`: what is and is not distributed here.
- `checksums/SHA256SUMS.tsv`: SHA-256 inventory of repository files.

## Reuse

The 95% ANI catalogue is intended for species/SGB-like catalogue analyses and
the 99% ANI catalogue for strain-like analyses. DogMAG can also be used with
miniTax to obtain rapid, detailed taxonomic profiles from canine metagenomes.

## Citation

Please cite the DogMAG data descriptor once published. Interim citation
metadata are provided in `CITATION.cff`.
