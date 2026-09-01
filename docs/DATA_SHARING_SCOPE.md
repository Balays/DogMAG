# Data sharing scope

## Included in this Git repository

- Stable data-access links and deposition-status tables.
- Data-release figures and figure source tables, pending a versioned Figshare release.
- Public supplementary metadata tables, pending a versioned Figshare release.
- Public ENA accession mappings.
- Workflow, catalogue, taxonomy, viral, mapping, and deposition scripts.

## Available through ENA

- Source reads are linked through PRJEB75753, PRJEB82125, and PRJEB85420.
- PRJEB115259 and all 84 DogMAG raw-read experiments and runs are public.
- All 41 final dog-wise primary metagenome assemblies are public under
  ERZ29880033-ERZ29880073 with linked coassembly BioSamples.
- All 135 dRep95 representative MAG BioSamples are public under
  SAMEA123300624-SAMEA123300758.
- All 122 submitted multi-contig MAG assembly analyses are public with ERZ
  accessions.
- The 13 single-contig MAG BioSamples are public; their FASTAs remain in the
  associated article data package because ENA's CONTIG route did not accept
  single-sequence assemblies.

These layers were publicly verified on 2026-09-01. Some MAG ERS aliases may
still return 404 even though the corresponding SAMEA BioSample records are
public. See `../data_access/` for the authoritative link and status index.

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
