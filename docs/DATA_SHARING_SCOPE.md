# Data sharing scope

## Included in this Git repository

- Final manuscript source and rendered DOCX.
- Publication figures and figure source tables.
- Public supplementary metadata tables.
- Public ENA accession mappings.
- Workflow, catalogue, taxonomy, viral, mapping, and deposition scripts.

## Available through ENA or the article data record

- Raw sequencing reads.
- The 81 long-read assemblies and linked assembly products.
- The 11,276-member final BASALT binset.
- The 3,418 reselected medium-quality-or-better candidate MAG FASTAs.
- The 792 dRep99 and 135 dRep95 representative MAG FASTAs.
- Viral/proviral candidate sequence FASTAs and depletion-panel FASTAs.

## Deliberately excluded from Git

- Raw FASTQ, BAM, SAM, and intermediate mapping files.
- The full BASALT working directory and all 30,556 candidate FASTA versions.
- Bowtie2, minimap2, DIAMOND, GTDB, CheckM2, and other databases or indexes.
- Temporary work directories, logs, checkpoints, and restart snapshots.
- Private dog-name crosswalks, internal barcode mappings, credentials, and test receipts.

Supplementary Table 2 preserves candidate-level provenance but uses portable,
repository-relative sequence locations. Supplementary Table 5 likewise uses
portable assembly paths rather than local server paths.
