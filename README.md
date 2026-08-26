# DogMAG

DogMAG is a reusable canine faecal metagenome data resource comprising
dog-wise assemblies, quality-controlled and dereplicated metagenome-assembled
genome (MAG) catalogues, taxonomic metadata, viral/proviral candidate metadata,
and reproducible workflow scripts.

This is the public data and workflow repository from which the DogMAG data
descriptor was prepared. It is not the journal submission repository. Article
drafts, editorial material, submission checklists, private crosswalks, ENA
receipts, credentials, and operational submission notes are intentionally not
distributed here.

## Key catalogue outputs

- 30,556 polished or reassembled candidate versions evaluated by BASALT.
- 11,276 BASALT-selected bin/version records.
- 3,418 reselected medium-quality-or-better MAG candidates.
- 792 representatives at 99% ANI and 135 representatives at 95% ANI.
- GTDB-Tk taxonomy for all 792 strain-like representatives.
- 22,068 viral/proviral candidate rows before final quality filtering.

The sequence files themselves are not stored in Git. Stable access points and
current deposition status are listed in `data_access/`:

- source and DogMAG reads are linked through their ENA studies;
- all 41 final dog-wise primary metagenome assemblies have accepted `ERZ`
  accessions (`ERZ29880033`-`ERZ29880073`) and linked coassembly BioSamples;
- all 135 representative MAG BioSamples are registered in ENA;
- 122 multi-contig MAG assemblies have accepted `ERZ` accessions; and
- 13 single-contig MAG submissions remain pending an ENA-approved route.

The DogMAG ENA study, source-sample records, and 41 primary metagenome
assemblies are public. Newly submitted DogMAG read experiments and MAG-derived
sample/assembly records remain held at the time of this repository release;
their accession mappings are retained here so the links become resolvable when
ENA releases those record layers.

## Repository map

- `data_access/`: stable links and deposition status for reads, assemblies, and MAGs.
- `figures/`: current data-release figures and source tables, retained until a
  versioned Figshare release is available.
- `supplementary_tables/`: current data-release metadata tables, retained until
  a versioned Figshare release is available.
- `scripts/`: public workflow and reporting scripts grouped by analysis stage.
- `accessions/`: public ENA sample, read, primary-assembly, and MAG accession mappings.
- `docs/DATA_SHARING_SCOPE.md`: what is and is not distributed here.
- `checksums/SHA256SUMS.tsv`: SHA-256 inventory of repository files.

## Reuse

The 95% ANI catalogue is intended for species/SGB-like catalogue analyses and
the 99% ANI catalogue for strain-like analyses. DogMAG can also be used with
miniTax to obtain rapid, detailed taxonomic profiles from canine metagenomes.

## Citation

Please cite the DogMAG data descriptor and the versioned data release once
published. Interim citation metadata are provided in `CITATION.cff`.
