# Data sharing scope

## Included in this Git repository

- Stable data-access links and deposition-status tables.
- Data-release figures and figure source tables, pending a versioned Figshare release.
- Public supplementary metadata tables, pending a versioned Figshare release.
- Public ENA accession mappings.
- Workflow, catalogue, taxonomy, viral, mapping, and deposition scripts.

## Available or planned through ENA

- Source reads are linked through PRJEB75753, PRJEB82125, and PRJEB85420.
- DogMAG reads are registered under PRJEB115259 but remain held/private.
- All 41 final dog-wise primary metagenome assemblies have accepted analysis
  accessions ERZ29880033-ERZ29880073 and linked coassembly BioSamples
  ERS31153171-ERS31153211.
- Of the 135 dRep95 representative MAGs, 122 multi-contig assemblies have
  accepted ERZ accessions and 13 single-contig assemblies remain pending.

See `../data_access/` for the authoritative link and status index.

## Planned versioned data release

Figures, supplementary metadata tables, viral/proviral candidate sequences,
depletion panels, and any sequence products not deposited in ENA may be
released through Figshare. Until a DOI is assigned and verified, repository
copies of the figures and tables are retained and no placeholder DOI is used.

## Deliberately excluded from Git

- Raw FASTQ, BAM, SAM, and intermediate mapping files.
- The full BASALT working directory and all 30,556 candidate FASTA versions.
- Bowtie2, minimap2, DIAMOND, GTDB, CheckM2, and other databases or indexes.
- Temporary work directories, logs, checkpoints, and restart snapshots.
- Private dog-name crosswalks, internal barcode mappings, credentials, and test receipts.
- Manuscript drafts, editorial files, journal submission forms, and operational
  ENA action notes.

Supplementary Table 2 preserves candidate-level provenance but uses portable,
repository-relative sequence locations. Supplementary Table 5 likewise uses
portable assembly paths rather than local server paths.
