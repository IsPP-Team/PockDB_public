# Pock-DB public binding-site preparation workflow

This repository provides a public companion workflow for preparing and
finalizing Pock-DB-style protein-ligand binding-site records with matched
ligand-based (LB) and geometry-based (GB) pocket representations.

The workflow is designed for external users who want to reproduce, inspect, or
validate individual PDB biological assemblies using the same file-naming and
annotation conventions as Pock-DB. It does not require or distribute a local
PockDrug installation. Instead, it prepares cleaned biological assembly files
that users submit manually to the public PockDrug web server, then parses the
downloaded PockDrug server ZIP archive to build the final structural files and
annotation table.

This public workflow is intended for per-entry or small-subset processing. It is
not intended to rerun the complete large-scale production workflow used to
generate the full Pock-DB 2026 release.

## Relation to Pock-DB

Pock-DB is a standardized dataset of protein-ligand binding-site entries derived
from experimentally determined PDB structures. Each retained binding-site entry
is associated with two complementary structural representations:

- an LB (ligand-based) pocket, defined from the protein environment surrounding an
  experimentally observed ligand;
- a GB (geometry-based) pocket, detected independently from cavity geometry.

This repository implements a public server-compatible workflow for generating
Pock-DB-style outputs for a selected PDB biological assembly. The retained GB
pocket is selected as the geometry-based cavity with the highest strictly
positive SO overlap with the corresponding LB pocket.

## Repository contents

```text
pockdb_public/
├── prepare_pdb_bank.py
├── requirements.txt
├── README.md
├── .gitignore
├── data/
│   └── Table1_additive_like_or_non_primary_ligands_CCD_HET_codes.csv
└── examples/
    └── output_pockdrug/
        └── 1rm8_assembly1/
            └── result_fpocket_prox5_5.zip

```

`examples/output_pockdrug/1rm8_assembly1/result_fpocket_prox5_5.zip` is an example PockDrug web-server output for PDB ID `1rm8`, biological assembly `1`.

## What the workflow does

The script has two stages:

1. `prepare-server-inputs`
   - retrieves and processes one PDB entry;
   - prepares cleaned biological assembly files;
   - writes the cleaned assembly PDB/CIF files to the final output directory;
   - prints instructions for submitting the cleaned assembly PDB file(s) to the PockDrug web server.

2. `finalize`
   - reads a PockDrug web-server ZIP result;
   - remaps PockDrug server ligand identifiers to the cleaned mmCIF identifiers;
   - selects the geometry-based pocket with the highest positive SO value;
   - exports the retained ligand, ligand environment, LB pocket, GB pocket, and final CSV table.

The final retained-entry filters are:

```text
C_RESIDUES_LB > 3
SO_GB > 0
```

`C_RESIDUES_GB` is exported as a descriptor of the matched geometry-based pocket, but it is **not** used as a retained-entry filter.

## Requirements

The workflow requires Python 3.11 or later and the Python packages listed in `requirements.txt`.

PyMOL is required for structure preparation and file conversion. 
The recommended installation uses the provided Conda environment file.

```bash
conda env create -f environment.yml
conda activate pockdb
```

The script also queries public RCSB/PDBe web services, so an internet connection is required.

## Quick start with the included 1RM8 example

Run the following commands from the root of the repository:

```bash
cd pockdb_public
```

### Step 1: prepare cleaned assemblies for the PockDrug server

```bash
python prepare_pdb_bank.py prepare-server-inputs 1rm8 \
  --clean-output-dir examples/clean_assemblies
```

This creates the server-input files under:

```text
examples/clean_assemblies/1rm8/pockdrug_server_input/
```

For the included example, the file to upload to the PockDrug server is expected to be:

```text
examples/clean_assemblies/1rm8/pockdrug_server_input/1rm8_1_cleaned.pdb
```

The corresponding cleaned CIF file is also kept for traceability:

```text
examples/clean_assemblies/1rm8/pockdrug_server_input/1rm8_1_cleaned.cif
```

### Step 2: run the PockDrug web server

Open the PockDrug web server:

```text
https://pockdrug.rpbs.univ-paris-diderot.fr/cgi-bin/index.py?page=Druggability
```

Use the following settings:

1. Go to **Druggability Prediction using protein(s)**.
2. In **Protein(s) information**, choose **upload your PDB file** and upload the cleaned assembly PDB file, for example `1rm8_1_cleaned.pdb`.
3. Under **Pocket estimation method(s)**, select both:
   - `fpocket`
   - `prox`
4. Set **Ligand proximity threshold between 4 Å and 12 Å** to:

```text
5.5 Å
```

5. Submit the job.
6. Download the resulting ZIP archive.

For the worked example, the downloaded server archive is already included here:

```text
examples/output_pockdrug/1rm8_assembly1/result_fpocket_prox5_5.zip
```

### Step 3: finalize the entry from the PockDrug server ZIP

For PDB ID `1rm8`, assembly `1`, run:

```bash
python prepare_pdb_bank.py finalize 1rm8 \
  --num-assembly 1 \
  --clean-output-dir examples/results \
  --pockdrug-server-results examples/output_pockdrug/1rm8_assembly1/result_fpocket_prox5_5.zip
```

Expected terminal message for the included example:

```text
For PDB ID 1rm8, assembly 1: 1 validated binding-site entries were retained after applying filters.
```

## Output structure

After the two steps, the output directory should contain:

```text
pockdb_public/
└── examples/
    └── output_pockdrug/
    │   └── 1rm8_assembly1/
    │       └── result_fpocket_prox5_5.zip
    │
    └── clean_assemblies/
    │   └── 1rm8/
    │       └── pockdrug_server_input/
    │           └── 1rm8_1_cleaned.cif
    │           └── 1rm8_1_cleaned.pdb
    └── results/
        └── 1rm8_1/
            └── 1rm8_1_BAT_F_800__GB_pocket_atm.pdb
            └── 1rm8_1_BAT_F_800__GB_pocket_res.pdb
            └── 1rm8_1_BAT_F_800__LB_pocket_atm.pdb
            └── 1rm8_1_BAT_F_800__LB_pocket_res.pdb
            └── 1rm8_1_BAT_F_800__ligand.pdb
            └── 1rm8_1_BAT_F_800__ligand_environment.cif
            └── 1rm8_1_cleaned.cif
            └── 1rm8__final_website.csv
```

The final folder name includes the assembly number:

```text
<pdb_id>_<num_assembly>/
```

This makes it easier to distinguish outputs when several assemblies are tested independently on the PockDrug server.

## Several assemblies for one PDB ID

If the preparation stage generates several cleaned assemblies, for example:

```text
examples/clean_assemblies/<pdb_id>/pockdrug_server_input/<pdb_id>_1_cleaned.pdb
examples/clean_assemblies/<pdb_id>/pockdrug_server_input/<pdb_id>_2_cleaned.pdb
```

submit each assembly PDB file separately to the PockDrug web server.

Then run `finalize` once per assembly, each time with the matching `--num-assembly` value and the corresponding server ZIP file:

```bash
python prepare_pdb_bank.py finalize <pdb_id> \
  --num-assembly 1 \
  --clean-output-dir results \
  --pockdrug-server-results path/to/assembly1_result.zip
```

```bash
python prepare_pdb_bank.py finalize <pdb_id> \
  --num-assembly 2 \
  --clean-output-dir results \
  --pockdrug-server-results path/to/assembly2_result.zip
```

The outputs will be written to separate folders:

```text
results/<pdb_id>_1/
results/<pdb_id>_2/
```

## Temporary files

Users only need to provide the final output directory with:

```text
--clean-output-dir
```

Temporary files are created automatically under:

```text
<clean-output-dir>/_temporary_work/
```

This temporary directory is removed automatically at the end of each command. It is not part of the final output.

## Final CSV

The final CSV is written as:

```text
results/<pdb_id>_<num_assembly>/<pdb_id>_<num_assembly>__final_website.csv
```

It contains structure-level, ligand-level, ligand-environment, and pocket-level information. The final table includes descriptors for both retained representations:

- `_LB`: ligand-based pocket descriptors;
- `_GB`: geometry-based pocket descriptors.

The retained geometry-based pocket is the geometry-based cavity with the highest strictly positive SO value for the corresponding ligand-based pocket.

## Additive like or non primary ligand annotation table

The repository includes the Additive like or non primary ligand CCD/HET-code table:

```text
data/Table1_additive_like_or_non_primary_ligands_CCD_HET_codes.csv
```

An alternative additive like or non primary ligand table can be lied with:

```bash
--additives-table path/to/custom_table.csv
```

The custom table must contain a `CCD/HET code` column.

## Useful commands

Show help:

```bash
python prepare_pdb_bank.py --help
```

Show help for the preparation stage:

```bash
python prepare_pdb_bank.py prepare-server-inputs --help
```

Show help for the finalization stage:

```bash
python prepare_pdb_bank.py finalize --help
```

## Troubleshooting

### `--num-assembly` is missing

The `finalize` command requires an assembly number. Use for example:

```bash
--num-assembly 1
```

### Zero retained binding-site entries

A result with zero retained entries means that no candidate passed the final filters:

```text
C_RESIDUES_LB > 3
SO_GB > 0
```

It may also happen if the PockDrug ZIP file does not correspond to the selected assembly number.

### The PockDrug server ligand name differs from the final ligand name

This can happen because PockDrug server names and cleaned mmCIF identifiers may differ. For example, in the 1RM8 example, the server ligand identifier is remapped to the cleaned structure identifier used in the final output:

```text
BAT_A_1 -> BAT_F_800
```

The final exported files use the cleaned structure identifier.

### Paths with spaces

Use quotes if your paths contain spaces:

```bash
python prepare_pdb_bank.py finalize 1rm8 \
  --num-assembly 1 \
  --clean-output-dir "my results" \
  --pockdrug-server-results "examples/1rm8_assembly1/result_fpocket_prox5_5.zip"
```
