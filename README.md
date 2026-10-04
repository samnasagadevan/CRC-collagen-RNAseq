# Input data

Raw data are not tracked in git. Place these files here before running the script:

| File | Description |
|---|---|
| `Combined_4PDO_rawdata.csv` | PDO raw gene counts. First column = Ensembl gene ID, one column per sample. |
| `sample_metadata.csv` | One row per PDO sample (16 = 4 PDO lines x 2 matrices x 2 treatments). Required columns: `Sample` (must match the count-matrix column names), `Group` and `Batch` (= PDO line). `PDO_line`, `Matrix` and `Treatment` are descriptive. See `sample_metadata_template.csv`. |
| `TCGA_COAD_read_counts.csv` | TCGA-COAD RNA-seq read counts. Column `Entrez_Gene_Id` (and optionally `Hugo_Symbol`), then one column per patient. |
| `Int_epi_genelist.csv` | Intestinal marker genes, one HGNC symbol per row in a column named `Genes`. |

Groups used by the analysis:

| Group | Matrix | Treatment |
|---|---|---|
| `MGcontrol` | Matrigel | none |
| `MGCOLcontrol` | Matrigel + collagen 1 | none |
| `MGtreated` | Matrigel | 5FU |
| `MGCOLtreated` | Matrigel + collagen 1 | 5FU |

Data availability:
- PDO RNA-seq: **[GEO / ArrayExpress accession]**
- TCGA-COAD: **[source and release, e.g. cBioPortal TCGA-COAD]**
