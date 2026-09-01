# Public accession mappings

This directory contains public-safe mappings exported from ENA submission work.

- `long_read_manifest.tsv`: public coded sample and FASTQ manifest for PRJEB115259.
- `source_sample_accessions.txt`: production ENA source-sample accessions.
- `primary_metagenome_assembly_accessions.tsv`: 41 public dog-wise primary
  metagenome assemblies, their coassembly BioSamples, BASALT-derived coverage,
  depth-table provenance, and ERZ analysis accessions.
- `mag_sample_accessions.tsv`: production MAG BioSample accession mapping. Its
  `PRIVATE` and `holdUntilDate` columns preserve the original registration
  receipt state; all 135 corresponding SAMEA BioSamples were public when
  rechecked on 2026-09-01.
- `mag_assembly_accessions_pass122.tsv`: 122 public multi-contig MAG analyses.
- Supplementary Table 9 records all 135 representatives, including the 13
  single-contig MAGs whose BioSamples are public but whose FASTAs are retained
  in the associated article data package.

For user-facing access links and current public status, see
`../data_access/`. Submission receipts, test-service records, credentials,
private crosswalks, and operational ENA action notes are intentionally excluded.
