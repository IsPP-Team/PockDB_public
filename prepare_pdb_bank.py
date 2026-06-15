# SPDX-License-Identifier: BSD-3-Clause
#!/usr/bin/env python3

"""Prepare Pock-DB-style protein-ligand binding-site entries from PDB data.

This public workflow is a companion implementation for Pock-DB 2026. It
prepares cleaned biological assemblies for submission to the public PockDrug
web server, then finalizes matched ligand-based (LB) and geometry-based (GB)
pocket files from the downloaded server ZIP output.

The script is intended for per-entry processing, inspection, and validation of
individual PDB biological assemblies. It does not reproduce the complete
large-scale production workflow used to build the full Pock-DB 2026 release.

Only the SO overlap score is used to match ligand-based and geometry-based
pocket representations. Retained entries are exported using Pock-DB-compatible
file naming and annotation conventions, including the composite binding-site
identifier:

    <PDB_ID>_<num_assembly>_<Ligand_ID>_<Chain_Asym_ID_-_ligand>_<num_ligand>

The implementation keeps the scientific logic of the original workflow while
making the code suitable for public reuse: 
the workflow is executed through a command-line interface.
"""

from __future__ import annotations

import argparse
import ast
import logging
import math
import os
import re
import shutil
import time
import zipfile
import io
from collections import Counter, defaultdict
from glob import glob
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import requests

try:
    from pymol import cmd
except ModuleNotFoundError:
    cmd = None

from scipy.spatial import cKDTree

LOGGER = logging.getLogger(__name__)
PUBLIC_RELEASE_VERSION = "v1.0.0-server-so-single-workdir"


DEFAULT_DISTANCE_THRESHOLD = "5.5"
DEFAULT_SERVER_TMP_SUBDIR = "pockdrug_server_tmp"
DEFAULT_ADDITIVES_TABLE = (
    Path(__file__).resolve().parent
    / "data"
    / "Table1_additive_like_or_non_primary_ligands_CCD_HET_codes.csv"
)

# The public Git version ships with a curated additive like or non primary ligand table. Users can
# override it with --additives-table if they want to use a different list.
ADDITIVES_TABLE_PATH: str | None = str(DEFAULT_ADDITIVES_TABLE)

# Runtime globals are configured from command-line arguments before processing.
tmp_dir = ""
tmp_root = ""
dir_files = ""

# PyMOL is used as a structural preprocessing engine. Disable verbose feedback
# to keep batch logs readable.
if cmd is not None:
    cmd.feedback("disable", "all", "details")

DISTANCE_THRESHOLD = "5.5"
RCSB_SESSION = requests.Session()
RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"

# ---------------------------------------------------------------------------
# RCSB and PDBe metadata retrieval
# ---------------------------------------------------------------------------


def query_rcsb_assembly(pdb_id):
    """
    Query the RCSB GraphQL API for a given PDB identifier.

    Parameters
    ----------
    pdb_id : str
        The PDB identifier.

    Returns
    -------
    dict or None
        Parsed entry dictionary if successful, otherwise None.
    """
    try:
        query = f'\n        {{\n          entry(entry_id: "{pdb_id}") {{\n            rcsb_id\n            assemblies {{rcsb_assembly_container_identifiers {{assembly_id}}\n                         rcsb_struct_symmetry {{oligomeric_state}} }}\n            em_3d_reconstruction {{resolution}}\n            exptl {{method}}\n            exptl_crystal_grow {{method \n                                pdbx_details}}\n            pdbx_nmr_ensemble {{conformers_submitted_total_number}}\n            pdbx_nmr_sample_details {{contents}}\n            rcsb_accession_info {{deposit_date}}\n            rcsb_entry_info {{deposited_polymer_monomer_count\n                             resolution_combined\n                             structure_determination_methodology}}\n            rcsb_primary_citation {{pdbx_database_id_DOI\n                                   title}}\n            struct {{title}}\n            struct_keywords {{pdbx_keywords}}\n            polymer_entities {{polymer_entity_instances {{rcsb_polymer_entity_instance_container_identifiers {{asym_id \n                                                                                                               auth_asym_id}}}}\n                              rcsb_cluster_membership {{cluster_id identity}}\n                              rcsb_entity_source_organism {{ncbi_scientific_name}}\n                              rcsb_polymer_entity {{formula_weight \n                                                   pdbx_description}}\n                              rcsb_polymer_entity_container_identifiers {{entity_id reference_sequence_identifiers {{database_accession \n                                                                                                                     database_name}}}}}}\n            nonpolymer_entities {{nonpolymer_comp {{chem_comp {{formula \n                                                               formula_weight \n                                                               id name}}\n                                                    rcsb_chem_comp_descriptor {{InChI\n                                                                                SMILES}}\n                                                    rcsb_chem_comp_related {{resource_accession_code resource_name}}}}\n                                 rcsb_nonpolymer_entity_annotation {{type}}\n                                 rcsb_nonpolymer_entity_container_identifiers {{entity_id}}\n                                 nonpolymer_entity_instances {{rcsb_nonpolymer_entity_instance_container_identifiers {{asym_id \n                                                                                                                       auth_asym_id}}}}}}\n            branched_entities {{rcsb_branched_entity_container_identifiers {{chem_comp_monomers}}\n                                branched_entity_instances {{rcsb_branched_struct_conn {{connect_target {{label_asym_id\n                                                                                                         label_comp_id}}}}\n            }}}}\n        }}\n        }}\n        '
        headers = {"Content-Type": "application/json"}
        for attempt in range(10):
            try:
                response = RCSB_SESSION.post(
                    RCSB_GRAPHQL_URL, json={"query": query}, headers=headers, timeout=20
                )
                data = response.json()
                entry = data["data"]["entry"]
                return entry
            except Exception:
                time.sleep(1)
        print(f"Error: {pdb_id}: max retries exceeded")
        return None
    except Exception as e:
        print(f"Error: {pdb_id}: {e}")
        return None


def assembly_oligomeric(entry):
    """
    Extract oligomeric state information for each assembly
    in a given RCSB entry dictionary.

    Parameters
    ----------
    entry: dict
        Entry dictionary returned by the RCSB API.

    Returns
    -------
    dict
        Dictionary mapping assembly_id to oligomeric state(s).
    """
    dict_ass_oligomeric = {}
    for el in entry["assemblies"]:
        oligomeric = []
        assembly = el["rcsb_assembly_container_identifiers"]["assembly_id"]
        if el["rcsb_struct_symmetry"]:
            for sym in el["rcsb_struct_symmetry"]:
                oligomeric.append(sym["oligomeric_state"])
        else:
            oligomeric.append(None)
        if len(oligomeric) == 1:
            dict_ass_oligomeric[assembly] = oligomeric[0]
        else:
            dict_ass_oligomeric[assembly] = oligomeric
    return dict_ass_oligomeric


def process_pdb(entry):
    """
    Process a PDB entry retrieved from RCSB and construct
    structured pandas DataFrames combining general information,
    ligand data, and polymer chain data.

    Parameters
    ----------
    entry : dict
        Entry dictionary returned by the RCSB API.

    Returns
    -------
    pandas.DataFrame
        Final merged DataFrame containing structural metadata.
    """
    df_general = pd.DataFrame(
        [
            {
                "pdb_id": entry.get("rcsb_id"),
                "Experimental Method": (
                    entry.get("exptl", {})[0].get("method")
                    if entry.get("exptl")
                    else None
                ),
                "Refinement Resolution (Å)": (
                    (entry.get("rcsb_entry_info", {}).get("resolution_combined") or [])[
                        0
                    ]
                    if len(
                        entry.get("rcsb_entry_info").get("resolution_combined") or []
                    )
                    == 1
                    else (
                        entry.get("rcsb_entry_info").get("resolution_combined") or []
                        if entry.get("rcsb_entry_info").get("resolution_combined") or []
                        else None
                    )
                ),
                "Crystal Growth Procedure": (
                    entry.get("exptl_crystal_grow", [{}])[0].get("pdbx_details")
                    if entry.get("exptl_crystal_grow")
                    else None
                ),
                "Deposition Date": entry.get("rcsb_accession_info", {})
                .get("deposit_date")
                .split("T")[0],
                "Structure Determination Methodology": entry.get(
                    "rcsb_entry_info", {}
                ).get("structure_determination_methodology"),
                "DOI": (
                    entry.get("rcsb_primary_citation", {}).get("pdbx_database_id_DOI")
                    if entry.get("rcsb_primary_citation")
                    else None
                ),
                "Title": (
                    entry.get("rcsb_primary_citation", {}).get("title")
                    if entry.get("rcsb_primary_citation")
                    else None
                ),
                "Structure Title": entry.get("struct", {}).get("title"),
                "Stucture Keywords": (
                    entry.get("struct_keywords", {}).get("pdbx_keywords")
                    if entry.get("struct_keywords")
                    else None
                ),
                "List of Unique Monosaccharides": entry.get("branched_entities", {}),
            }
        ]
    )
    rows = []
    if entry.get("nonpolymer_entities", []):
        for item in entry.get("nonpolymer_entities"):
            base = {
                "entity_id": item.get(
                    "rcsb_nonpolymer_entity_container_identifiers", {}
                ).get("entity_id"),
                "Ligand Formula": (
                    item.get("nonpolymer_comp", {}).get("chem_comp", {}).get("formula")
                    if item.get("nonpolymer_comp")
                    else None
                ),
                "Ligand MW": (
                    item.get("nonpolymer_comp", {})
                    .get("chem_comp", {})
                    .get("formula_weight")
                    if item.get("nonpolymer_comp")
                    else None
                ),
                "Ligand ID": (
                    item.get("nonpolymer_comp", {}).get("chem_comp", {}).get("id")
                    if item.get("nonpolymer_comp")
                    else None
                ),
                "Ligand Name": (
                    item.get("nonpolymer_comp", {}).get("chem_comp", {}).get("name")
                    if item.get("nonpolymer_comp")
                    else None
                ),
                "InChI": (
                    item.get("nonpolymer_comp", {})
                    .get("rcsb_chem_comp_descriptor", {})
                    .get("InChI")
                    .split("InChI=")[1]
                    if item.get("nonpolymer_comp", {})
                    and item.get("nonpolymer_comp", {}).get(
                        "rcsb_chem_comp_descriptor", {}
                    )
                    and item.get("nonpolymer_comp", {})
                    .get("rcsb_chem_comp_descriptor", {})
                    .get("InChI")
                    else None
                ),
                "Ligand SMILES": (
                    item.get("nonpolymer_comp", {})
                    .get("rcsb_chem_comp_descriptor", {})
                    .get("SMILES")
                    if item.get("nonpolymer_comp", {})
                    and item.get("nonpolymer_comp", {}).get(
                        "rcsb_chem_comp_descriptor", {}
                    )
                    and item.get("nonpolymer_comp", {})
                    .get("rcsb_chem_comp_descriptor", {})
                    .get("SMILES")
                    else None
                ),
            }
            grouped = defaultdict(list)
            if (
                item.get("nonpolymer_comp", {}).get("rcsb_chem_comp_related", {})
                if item.get("nonpolymer_comp")
                else None
            ):
                for item2 in item.get("nonpolymer_comp", {}).get(
                    "rcsb_chem_comp_related", {}
                ):
                    if item2["resource_name"] == "ChEBI":
                        grouped[item2["resource_name"]].append(
                            item2["resource_accession_code"].split("CHEBI:")[1]
                        )
                    else:
                        grouped[item2["resource_name"]].append(
                            item2["resource_accession_code"]
                        )
            grouped = dict(grouped)
            grouped = {f"Ligand-{k}-Accession Code(s)": v for k, v in grouped.items()}
            base.update(grouped)
            if item.get("nonpolymer_entity_instances", {}):
                for instance in item.get("nonpolymer_entity_instances", {}):
                    ids = instance.get(
                        "rcsb_nonpolymer_entity_instance_container_identifiers", {}
                    )
                    row = base.copy()
                    row["Asym ID - nonpolymer"] = ids.get("asym_id")
                    row["Auth Asym ID - nonpolymer"] = ids.get("auth_asym_id")
                    rows.append(row)
    df_lgds = pd.DataFrame(rows)
    rows_chains = []
    for item in entry.get("polymer_entities", []):
        rcsb_entity_source_organism = item.get("rcsb_entity_source_organism") or [{}]
        rcsb_polymer_entity = item.get("rcsb_polymer_entity") or [{}]
        rcsb_polymer_entity_container_identifiers = item.get(
            "rcsb_polymer_entity_container_identifiers", {}
        )
        reference_sequence_identifiers = rcsb_polymer_entity_container_identifiers.get(
            "reference_sequence_identifiers"
        ) or [{}]
        for polymer_entity_instance in item.get("polymer_entity_instances", []):
            rcsb_polymer_entity_instance_container_identifiers = (
                polymer_entity_instance.get(
                    "rcsb_polymer_entity_instance_container_identifiers", {}
                )
            )
            base = {
                "Asym ID": rcsb_polymer_entity_instance_container_identifiers.get(
                    "asym_id"
                ),
                "Auth Asym ID": rcsb_polymer_entity_instance_container_identifiers.get(
                    "auth_asym_id"
                ),
                "Source Organism": [
                    ref.get("ncbi_scientific_name")
                    for ref in rcsb_entity_source_organism
                    if ref.get("ncbi_scientific_name")
                ],
                "Molecular Weight (Entity)(KDa)": rcsb_polymer_entity.get(
                    "formula_weight"
                ),
                "Macromolecule Name": rcsb_polymer_entity.get("pdbx_description"),
                "Entity ID_polymer": rcsb_polymer_entity_container_identifiers.get(
                    "entity_id"
                ),
                "Accession Code(s)": [
                    ref.get("database_accession")
                    for ref in reference_sequence_identifiers
                    if ref.get("database_accession")
                ],
                "Database Name": [
                    ref.get("database_name")
                    for ref in reference_sequence_identifiers
                    if ref.get("database_name")
                ],
            }
            database = base["Database Name"]
            if len(set(database)) == 1:
                base[f"{database[0]}-Accession Code(s)"] = base["Accession Code(s)"]
                base.pop("Accession Code(s)", None)
                base.pop("Database Name", None)
            elif base["Accession Code(s)"] == [] and base["Database Name"] == []:
                base.pop("Accession Code(s)", None)
                base.pop("Database Name", None)
            rows_chains.append(base)
    df_chains = pd.DataFrame(rows_chains)
    df_prot_repeated = pd.concat([df_general] * len(df_chains), ignore_index=True)
    df_chain_merged = pd.concat(
        [df_prot_repeated, df_chains.reset_index(drop=True)], axis=1
    )
    if len(df_lgds) != 0:
        df_final = df_chain_merged.merge(
            df_lgds,
            left_on="Auth Asym ID",
            right_on="Auth Asym ID - nonpolymer",
            how="left",
        )
    else:
        df_final = df_chain_merged.copy()
    return df_final


def list_saccharide_ligands(infos_fromrcsb):
    """
    Extract ligand IDs and unique monosaccharides
    from the processed RCSB dataframe.

    Parameters
    ----------
    infos_fromrcsb : pandas.DataFrame
        DataFrame returned by process_pdb.

    Returns
    -------
    tuple
        (list_ligands, list_unique_monosaccharides)
    """
    list_sacc = []
    for row in infos_fromrcsb["List of Unique Monosaccharides"]:
        if row is not None:
            for row_dict in row:
                for info_sacc in row_dict["rcsb_branched_entity_container_identifiers"][
                    "chem_comp_monomers"
                ]:
                    list_sacc.append(info_sacc)
    list_Monosaccharides = list(set(list_sacc))
    list_ligands = list(infos_fromrcsb["Ligand ID"])
    return (list_ligands, list_Monosaccharides)


def monosaccharides_chain(entry):
    """
    Extract monosaccharide chain identifiers from an RCSB entry.

    For each branched entity present in the PDB entry,
    retrieves the connected monosaccharide components
    and formats them as:

        <label_comp_id>_<label_asym_id>

    Parameters
    ----------
    entry : dict
        RCSB assembly entry returned by GraphQL query.

    Returns
    -------
    list
        Unique list of monosaccharide chain identifiers.
    """
    sacc_chain = []
    if entry["branched_entities"] is not None:
        for branch in entry["branched_entities"]:
            for sacc in branch["branched_entity_instances"]:
                for connect_target in sacc["rcsb_branched_struct_conn"]:
                    sacc_chain.append(
                        f"{connect_target['connect_target']['label_comp_id']}_{connect_target['connect_target']['label_asym_id']}"
                    )
    return list(set(sacc_chain))


def ensure_list(x):
    if x is None:
        return []
    if isinstance(x, float) and pd.isna(x):
        return []
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, list):
        return x
    if isinstance(x, str):
        try:
            val = ast.literal_eval(x)
            return val if isinstance(val, list) else [val]
        except Exception:
            return [x]
    return [x]


LIGAND_CACHE = {}


def generate_pdbe_ligand_output(ligand_id):
    """
    Retrieve ligand summary information from the PDBe Graph API
    with caching and retry.

    Parameters
    ----------
    ligand_id : str
        Ligand three-letter code or PDBe compound identifier.

    Returns
    -------
    dict or None
        Ligand summary dictionary if request succeeds,
        otherwise None.
    """
    if ligand_id in LIGAND_CACHE:
        return LIGAND_CACHE[ligand_id]
    url_pdbe = f"https://www.ebi.ac.uk/pdbe/graph-api/compound/summary/{ligand_id}"
    for attempt in range(10):
        try:
            response = requests.get(url_pdbe, timeout=10)
            if response.status_code == 200:
                ligand_data = response.json()
                if ligand_id in ligand_data:
                    result = ligand_data[ligand_id][0]
                else:
                    result = None
                LIGAND_CACHE[ligand_id] = result
                return result
            else:
                LIGAND_CACHE[ligand_id] = None
                return None
        except Exception:
            time.sleep(5)
    return None


def generate_pdbe_ligand_dict_crosslinks(output_pdbe):
    """
    Extract selected cross-link resources from PDBe ligand output.

    Filters cross-link entries based on `list_keep_resource`
    and formats them into a standardized dictionary:

        'Ligand-<Resource>-Accession Code(s)' : [ids]

    Special case:
        'CCDC' is renamed to 'CCDC/CSD'.

    Parameters
    ----------
    output_pdbe : dict
        Dictionary returned by PDBe ligand summary API.

    Returns
    -------
    dict
        Dictionary of filtered cross-link accession codes.
        Empty dictionary if no cross-links are present.
    """
    dict_pdbe = defaultdict(list)
    if output_pdbe["cross_links"]:
        for item in output_pdbe["cross_links"]:
            if item["resource"] == "CCDC":
                dict_pdbe["CCDC/CSD"].append(item["resource_id"])
            elif item["resource"] == "ChEBI":
                dict_pdbe[item["resource"]].append(item["resource_id"])
            else:
                dict_pdbe[item["resource"]].append(item["resource_id"])
        dict_pdbe = dict(dict_pdbe)
        dict_pdbe = {f"Ligand-{k}-Accession Code(s)": v for k, v in dict_pdbe.items()}
        return dict_pdbe
    else:
        return dict_pdbe


def add_ligands_pdbe_infos(infos_fromrcsb):
    """
    Enrich RCSB ligand dataframe with extended PDBe information (cross-links // physicochemical properties // cofactors)

    For each ligand:
    - Retrieve PDBe summary data
    - Extract selected cross-links
    - Add selected cross-links as new columns if necessary
    - Add physicochemical properties
    - Determine cofactor-like status (PDBe API)
    - Evaluate Lipinski rule compliance

    Parameters
    ----------
    infos_fromrcsb : pandas.DataFrame
        DataFrame containing at least a "Ligand ID" column.

    Returns
    -------
    pandas.DataFrame
        Updated dataframe including:
        - Cross-link accession codes
        - Physicochemical properties
        - Cofactor annotation
        - Lipinski compliance flag
    """
    list_of_desc = []
    for idx in infos_fromrcsb.index:
        ligand_id = infos_fromrcsb.at[idx, "Ligand ID"]
        if ligand_id:
            output_pdbe = generate_pdbe_ligand_output(ligand_id)
            if output_pdbe is not None:
                ligand_dict = generate_pdbe_ligand_dict_crosslinks(output_pdbe)
                for col, new_values in ligand_dict.items():
                    if col not in infos_fromrcsb.columns:
                        infos_fromrcsb[col] = [[] for _ in range(len(infos_fromrcsb))]
                    infos_fromrcsb.at[idx, col] = list(
                        set(ensure_list(infos_fromrcsb.at[idx, col]) + new_values)
                    )
                if output_pdbe["phys_chem_properties"] is not None:
                    desc_dict = output_pdbe["phys_chem_properties"]
                else:
                    desc_dict = {}
            else:
                desc_dict = {}
            url_cofactor = f"https://www.ebi.ac.uk/pdbe/api/v2/pdb/compound/cofactors/het/{ligand_id}"
            response = requests.get(url_cofactor)
            if response.status_code == 200:
                desc_dict["PDBe-Cofactor_like"] = True
            elif response.status_code == 404:
                desc_dict["PDBe-Cofactor_like"] = False
            else:
                desc_dict["PDBe-Cofactor_like"] = None
            required_keys = [
                "lipinski_hbd",
                "lipinski_hba",
                "Ligand MW",
                "crippen_clog_p",
            ]
            if all(
                (k in desc_dict and desc_dict[k] is not None for k in required_keys)
            ):
                criteria = [
                    desc_dict["lipinski_hbd"] <= 5,
                    desc_dict["lipinski_hba"] <= 10,
                    desc_dict["Ligand MW"] <= 500,
                    desc_dict["crippen_clog_p"] <= 5,
                ]
                desc_dict["Lipinski"] = sum(criteria) >= 3
            else:
                desc_dict["Lipinski"] = None
            list_additives = load_additive_codes()
            if ligand_id in list_additives:
                desc_dict["Additive_like_or_non_primary_ligand"] = True
            else:
                desc_dict["Additive_like_or_non_primary_ligand"] = False
            list_of_desc.append(desc_dict)
    desc_df = pd.DataFrame(list_of_desc)
    infos_fromrcsb = pd.concat([infos_fromrcsb.reset_index(drop=True), desc_df], axis=1)
    return infos_fromrcsb


# ---------------------------------------------------------------------------
# Biological assembly selection and structure preparation
# ---------------------------------------------------------------------------


def assembly_ids_single_pdb(pdb_id, query_data):
    """
    Retrieve all assembly IDs for a single PDB entry.

    Parameters
    ----------
    pdb_id : str
        PDB identifier.
    query_data : dict
        GraphQL response data for the PDB entry.

    Returns
    -------
    list or tuple
        List of assembly IDs if successful.
        In case of failure, returns ({}, pdb_id).
    """
    try:
        assembly_ids = [
            a["rcsb_assembly_container_identifiers"]["assembly_id"]
            for a in query_data["assemblies"]
        ]
        return assembly_ids
    except Exception:
        return ({}, pdb_id)


def find_all_minimal_combinations(pdb_id, assembly_ids):
    """
    Find minimal combinations of assemblies whose total atom count
    matches the ASU atom count.

    Parameters
    ----------
    pdb_id : str
        PDB identifier.
    assembly_ids : list
        List of assembly identifiers.

    Returns
    -------
    list
        List of valid assembly combinations matching ASU atom count.
    """
    cmd.cd(tmp_dir)
    cmd.reinitialize()
    asu_obj = "asu_ref"
    cmd.fetch(pdb_id, asu_obj, type="cif", async_=0)
    if asu_obj not in cmd.get_object_list():
        return []
    asu_atoms = cmd.count_atoms(asu_obj)
    assembly_objects = {}
    atom_counts = {}
    for assembly_number in assembly_ids:
        obj_name = f"asm_{assembly_number}"
        cmd.set("assembly", assembly_number)
        try:
            cmd.fetch(pdb_id, obj_name, async_=0)
        except Exception:
            continue
        if obj_name in cmd.get_object_list():
            assembly_objects[assembly_number] = obj_name
            atom_counts[assembly_number] = cmd.count_atoms(obj_name)
    if not assembly_objects:
        cmd.delete("all")
        return []
    assembly_keys = sorted(
        assembly_objects.keys(), key=lambda x: atom_counts[x], reverse=True
    )
    for r in range(1, len(assembly_keys) + 1):
        valid_combos_count = []
        for combo in combinations(assembly_objects.keys(), r):
            test_atoms = sum((atom_counts[a] for a in combo))
            if test_atoms == asu_atoms:
                valid_combos_count.append(list(combo))
        if valid_combos_count:
            break
    cmd.reinitialize()
    return valid_combos_count


def process_single_pdb(pdb_id, entry, assembly_ids, result):
    """
    Process a single PDB entry for assembly analysis.

    Parameters
    ----------
    pdb_id : str
        PDB identifier.
    entry : dict
        Entry data from RCSB.
    assembly_ids : list
        List of assembly identifiers.
    result : any
        Previously computed result to return.

    Returns
    -------
    dict or tuple
        Processed result dictionary if successful.
        In case of failure, returns ({}, pdb_id).
    """
    try:
        if not assembly_ids:
            return ({pdb_id: []}, None)
        filename = f"{pdb_id}.cif".lower()
        if os.path.exists(filename):
            os.remove(filename)
        return result
    except Exception:
        return ({}, pdb_id)


def list_hetatm(obj_name, list_ligands, list_Monosaccharides):
    """
    Retrieve HETATM residues excluding known ligands and monosaccharides.

    Parameters
    ----------
    obj_name : str
        Name of the PyMOL object.

    Returns
    -------
    set or None
        Set of (chain, resn, resi) tuples if remaining HETATM are found,
        otherwise None.
    """
    cmd.remove(f"{obj_name} and resn HOH")
    het_selection = f"{obj_name} and hetatm"
    het_residues = set()
    cmd.iterate(
        het_selection,
        "het_residues.add((chain, resn, resi))",
        space={"het_residues": het_residues},
    )
    if het_residues:
        for chain, resn, resi in sorted(het_residues):
            if resn in list_ligands or resn in list_Monosaccharides:
                het_residues.remove((chain, resn, resi))
    if het_residues == set():
        return None
    else:
        return het_residues


def list_hetatm_ligands_rcsb(obj_name):
    """
    Retrieve all HETATM residues (ligands) from a PyMOL object.

    Parameters
    ----------
    obj_name : str
        Name of the PyMOL object.

    Returns
    -------
    set
        Set of (chain, resn, resi) tuples corresponding to HETATM.
    """
    cmd.remove(f"{obj_name} and resn HOH")
    het_selection = f"{obj_name} and hetatm"
    het_residues = set()
    cmd.iterate(
        het_selection,
        "het_residues.add((chain, resn, resi))",
        space={"het_residues": het_residues},
    )
    return het_residues


def get_used_residues():
    """
    Retrieve all numeric residue identifiers currently used in the structure.

    Returns
    -------
    set
        Set of integer residue numbers.
    """
    used = set()
    for atom in cmd.get_model("all").atom:
        resi = atom.resi.strip()
        numeric = ""
        for c in resi:
            if c.isdigit():
                numeric += c
            else:
                break
        if numeric != "":
            used.add(int(numeric))
    return used


def renumber_ligand(ligand, chain, old_res):
    """
    Renumber a ligand residue to avoid conflicts with existing residues.

    Parameters
    ----------
    ligand : str
        Ligand residue name (e.g., HEM, NAD).
    chain : str
        Chain identifier.
    old_res : str
        Original residue number.
    """
    used_residues = get_used_residues()
    new_resi = max(used_residues) + 1
    cmd.alter(
        f"resn {ligand} and chain {chain} and resi {old_res}", f"resi='{new_resi}'"
    )
    return {f"{ligand}_{chain}": [old_res, f"{new_resi}"]}


def has_insertion_code(resi):
    """
    Detect PyMOL residues carrying an insertion code.
    In PyMOL, pdbx_PDB_ins_code is generally represented in resi,
    e.g. 100A, 100B.
    """
    resi = str(resi).strip()
    return re.match("^-?\\d+[A-Za-z]+$", resi) is not None


def numeric_part_resi(resi):
    """
    Extract numeric part of a residue identifier.
    Example: 100A -> 100
    """
    resi = str(resi).strip()
    m = re.match("^-?\\d+", resi)
    return int(m.group(0)) if m else 10**9


def renumber_inserted_residues(obj_name):
    """
    Renumber ATOM residues carrying an insertion code.

    This removes technical ambiguity between residues such as:
        100 and 100A

    Strategy:
        - keep the same chain
        - assign a new residue number above the current maximum
        - do not modify HETATM ligands, because they are already handled
          by renumber_ligand()
    """
    used_residues = get_used_residues()
    if not used_residues:
        return {}
    inserted_residues = set()
    cmd.iterate(
        f"{obj_name} and not hetatm and not solvent",
        "inserted_residues.add((chain, resn, resi)) if has_insertion_code(resi) else None",
        space={
            "inserted_residues": inserted_residues,
            "has_insertion_code": has_insertion_code,
        },
    )
    if not inserted_residues:
        return {}
    mapping = {}
    new_resi = max(used_residues) + 1
    for chain, resn, old_resi in sorted(
        inserted_residues, key=lambda x: (x[0], numeric_part_resi(x[2]), x[2], x[1])
    ):
        selection = f"{obj_name} and not hetatm and chain {chain} and resi {old_resi}"
        cmd.alter(selection, f"resi='{new_resi}'")
        mapping[f"INS_{resn}_{chain}_{old_resi}"] = [old_resi, str(new_resi)]
        new_resi += 1
    cmd.sort()
    return mapping


def enforce_single_altloc(obj_name):
    """
    Enforce a single alternative location (altloc) per residue.

    Ligand altloc is chosen based on occupancy.
    Protein residues follow ligand choice if available,
    otherwise 'A', otherwise alphabetical order.

    Parameters
    ----------
    obj_name : str
        Name of the PyMOL object.

    Returns
    -------
    tuple
        (residue_alts, ligand_choice, chosen_alt)
    """
    chosen_alt = []
    residue_alts = defaultdict(set)
    space = {"residue_alts": residue_alts}
    cmd.iterate(
        obj_name,
        "residue_alts.setdefault(chain+'_'+resi+'_'+resn, set()).add(alt) if alt != '' else None",
        space=space,
    )
    ligand_occ = defaultdict(float)
    space2 = {"ligand_occ": ligand_occ}
    cmd.iterate(
        f"{obj_name} and hetatm and not alt ''",
        "ligand_occ[alt] += q if alt != '' else None",
        space=space2,
    )
    ligand_choice = None
    if ligand_occ:
        ligand_choice = sorted(ligand_occ.items(), key=lambda x: (-x[1], x[0]))[0][0]
    for key, alts in residue_alts.items():
        alts = sorted([a for a in alts if a != ""])
        if not alts:
            continue
        if len(alts) == 1:
            chosen = alts[0]
        elif ligand_choice and ligand_choice in alts:
            chosen = ligand_choice
        elif "A" in alts:
            chosen = "A"
        else:
            chosen = alts[0]
        chain, resi, resn = key.split("_")
        chosen_alt.append(chosen)
        cmd.remove(
            f"{obj_name} and chain {chain} and resi {resi} and not alt ''+{chosen}"
        )
        cmd.alter(
            f"{obj_name} and chain {chain} and resi {resi} and alt {chosen}", "alt=''"
        )
    return (residue_alts, ligand_choice, chosen_alt)


def write_info_header(file_path):
    """
    Write the header section of the information summary file.

    Parameters
    ----------
    file_path : str
        Path to the output information file.
    """
    with open(file_path, "w") as f:
        f.write("# ============================================\n")
        f.write("# FILE INFORMATION SUMMARY\n")
        f.write("# ============================================\n")
        f.write("# Format: TAB separated table\n")
        f.write("# One row = one assembly\n")
        f.write("#\n")
        f.write("# Index description:\n")
        f.write("# PDB              : PDB identifier\n")
        f.write("# ASM_START_N      : Number of assemblies before filtering\n")
        f.write(
            "# ASM_COMB         : Minimal valid assembly combinations (atom count pre-filter) (list of list)\n"
        )
        f.write("# ASM_END_N        : Number of assemblies kept after filtering\n")
        f.write(
            "# ASM_LIG_N        : Number of assemblies containing at least one ligand\n"
        )
        f.write("# ASM_ID           : Assembly identifier currently processed\n")
        f.write(
            "# HET_TO_ATOM      : Modified residues converted from HETATM to ATOM\n"
        )
        f.write("# LIG_RCSB         : Ligands confirmed by RCSB (dictionary)\n")
        f.write("# tMONOSACCHARIDES : Monosaccharides confirmed by RCSB (dictionary)\n")
        f.write(
            "# LIG_5CH_MAP      : Mapping of 5-character ligand names to 3-character PDB-compatible names (dictionary)\n"
        )
        f.write("# ALT_RES          : Altloc detected per residue (dictionary)\n")
        f.write("# ALT_LIG_CHOICE   : Chosen ligand alternative\n")
        f.write("# ALT_STRUCT       : Final altloc kept in structure (list)\n")
        f.write(
            "# RENUMBER_LIG     : Renumbered ligands and ATOM residues with insertion codes (dictionary) {identifier:[old_number,new_number]}\n"
        )
        f.write("# NUMBER_H         : Number of hydrogen atoms removed (int)\n")
        f.write("# ============================================\n\n")
        f.write(
            "PDB\tASM_ID\tHET_TO_ATOM\tLIG_RCSB\tMONOSACCHARIDES\tLIG_5CH_MAP\tALT_RES\tALT_LIG_CHOICE\tALT_STRUCT\tRENUMBER_LIG\tNUMBER_H\n"
        )


def write_assembly_row(
    file_path,
    pdb_id,
    asm_id,
    het_to_atom,
    lig_rcsb,
    monosaccharides_chain,
    lig_5ch_map,
    alt_res,
    alt_choice,
    alt_struct,
    renumber_lig,
    number_h,
):
    """
    Append one assembly row to the information file.

    Parameters
    ----------
    file_path : str
        Path to the output file.
    pdb_id : str
        PDB identifier.
    asm_id : str
        Assembly identifier.
    het_to_atom : any
        Residues converted from HETATM to ATOM.
    lig_rcsb : any
        Ligands confirmed by RCSB.
    monosaccharides_chain : any
        Monosaccharides confirmed by RCSB.
    lig_5ch_map : dict or None
        Mapping of long ligand names to 3-character names.
    alt_res : dict or None
        Altloc per residue.
    alt_choice : str or None
        Selected ligand altloc.
    alt_struct : list or None
        Final altloc(s) kept in structure.
    renumber_lig : dict or None
        Renumber a ligand residue
    number_h : int
        Number of hydrogen atoms removed

    """
    with open(file_path, "a") as f:
        f.write(
            f"{pdb_id}\t{asm_id}\t{het_to_atom}\t{lig_rcsb}\t{monosaccharides_chain}\t{lig_5ch_map}\t{alt_res}\t{alt_choice}\t{alt_struct}\t{renumber_lig}\t{number_h}\n"
        )


def write_global_info(file_path, start_n, comb, end_n, lig_n):
    """
    Append global assembly statistics at the end of the file.

    Parameters
    ----------
    file_path : str
        Path to the output file.
    start_n : int
        Number of assemblies before filtering.
    comb : any
        Minimal assembly combinations.
    end_n : int
        Number of assemblies kept after filtering.
    lig_n : int
        Number of assemblies containing at least one ligand.
    """
    with open(file_path, "a") as f:
        f.write("\n# GLOBAL INFO\n")
        f.write(f"# ASM_START_N\t{start_n}\n")
        f.write(f"# ASM_COMB\t{comb}\n")
        f.write(f"# ASM_END_N\t{end_n}\n")
        f.write(f"# ASM_LIG_N\t{lig_n}\n")


def prepare_files_id(
    pdb_id, entry, dir_files, list_ligands, list_monosaccharides_chain
):
    """
    Prepare directories and processed structure files for a given PDB ID.

    Parameters
    ----------
    pdb_id : str
        PDB identifier.
    entry : dict
        RCSB entry data.
    dir_files : str
        Base directory where results are written.

    Returns
    -------
    int
        Number of assemblies containing at least one ligand.
    """
    start_assembly_ids = assembly_ids_single_pdb(pdb_id, entry)
    result = find_all_minimal_combinations(pdb_id, start_assembly_ids)
    if process_single_pdb(pdb_id, entry, start_assembly_ids, result) != []:
        assembly_ids = process_single_pdb(pdb_id, entry, start_assembly_ids, result)[0]
    else:
        assembly_ids = start_assembly_ids
    tot_assembly_ids = 0
    pdb_dir = f"{dir_files}{pdb_id}"
    if os.path.exists(pdb_dir):
        shutil.rmtree(pdb_dir)
    os.makedirs(pdb_dir)
    info_file = f"{pdb_dir}/infos_0.txt"
    write_info_header(info_file)
    list_Monosaccharides = []
    for sacc in list_monosaccharides_chain:
        list_Monosaccharides.append(sacc.split("_")[0])
    dict_assembly_infos_tmp = {}
    for assembly_number in assembly_ids:
        cmd.reinitialize()
        all_dict_renumber_ligand_all = {}
        assembly_dir = f"{pdb_dir}/{assembly_number}"
        os.makedirs(assembly_dir)
        obj_name = f"{pdb_id}_{assembly_number}"
        cmd.set("assembly", assembly_number)
        cmd.fetch(pdb_id, obj_name, async_=0)
        cmd.remove("resn HOH")
        cmd.remove("resn DOD")
        n_H = cmd.count_atoms(f"{obj_name} and elem H")
        cmd.remove(f"{obj_name} and elem H")
        cmd.remove(f"{obj_name} and elem D")
        list_hetatm_toatom = list_hetatm(obj_name, list_ligands, list_Monosaccharides)
        if list_hetatm_toatom:
            for chain, resn, resi in sorted(list_hetatm_toatom):
                selection = (
                    f"{obj_name} and resn {resn} and resi {resi} and chain {chain}"
                )
                cmd.set("pdb_hetatm_guess", 0)
                cmd.alter(selection, 'type="ATOM"')
        dict_renumber_inserted = renumber_inserted_residues(obj_name)
        all_dict_renumber_ligand_all.update(dict_renumber_inserted)
        list_hetatm_ligands = list_hetatm_ligands_rcsb(obj_name)
        lig_5ch_map = {}
        for chain, resn, resi in list_hetatm_ligands:
            if len(resn) > 4:
                truncated = resn[:3]
                lig_5ch_map[resn] = truncated
        if not lig_5ch_map:
            lig_5ch_map = None
        list_resi_num = [resi for chain, resn, resi in list_hetatm_ligands]
        if len(set(list_resi_num)) != len(list_hetatm_ligands):
            for chain, ligand, old_resi in sorted(list_hetatm_ligands):
                dict_renumber_ligand_all = renumber_ligand(ligand, chain, old_resi)
                all_dict_renumber_ligand_all.update(dict_renumber_ligand_all)
            list_hetatm_ligands = list_hetatm_ligands_rcsb(obj_name)
        ligand_selection = f"{obj_name} and hetatm"
        ligand_atoms = cmd.count_atoms(ligand_selection)
        alt_res = None
        alt_choice = None
        alt_struct = None
        if ligand_atoms > 0:
            tot_assembly_ids += 1
            alt_data = enforce_single_altloc(obj_name)
            alt_res = dict(alt_data[0])
            alt_choice = alt_data[1]
            alt_struct = list(set(alt_data[2]))
            cmd.sort()
            if lig_5ch_map is not None:
                cmd.save(f"{assembly_dir}/{obj_name}.cif", state=-1)
                for lig5 in lig_5ch_map.keys():
                    lig3 = lig_5ch_map[lig5]
                    cmd.alter(f"resn {lig5} and hetatm", f"resn='{lig3}'")
                cmd.save(f"{assembly_dir}/{obj_name}.pdb", state=-1)
            else:
                cmd.save(f"{assembly_dir}/{obj_name}.cif", state=-1)
                cmd.save(f"{assembly_dir}/{obj_name}.pdb", state=-1)
        if dict_assembly_infos_tmp != {}:
            difference = 0
            for key_dict_assembly_infos_tmp in dict_assembly_infos_tmp.keys():
                if dict_assembly_infos_tmp[key_dict_assembly_infos_tmp] != (
                    pdb_id,
                    list_hetatm_toatom,
                    list_hetatm_ligands,
                    list_monosaccharides_chain,
                    lig_5ch_map,
                    alt_res,
                    alt_choice,
                    alt_struct,
                    all_dict_renumber_ligand_all,
                    n_H,
                ):
                    difference += 1
                else:
                    pass
            if difference == len(dict_assembly_infos_tmp.keys()):
                write_assembly_row(
                    info_file,
                    pdb_id,
                    assembly_number,
                    list_hetatm_toatom,
                    list_hetatm_ligands,
                    list_monosaccharides_chain,
                    lig_5ch_map,
                    alt_res,
                    alt_choice,
                    alt_struct,
                    all_dict_renumber_ligand_all,
                    n_H,
                )
            else:
                tot_assembly_ids -= 1
                assembly_ids.remove(f"{assembly_number}")
                if os.path.exists(f"{assembly_dir}"):
                    shutil.rmtree(f"{assembly_dir}")
        else:
            write_assembly_row(
                info_file,
                pdb_id,
                assembly_number,
                list_hetatm_toatom,
                list_hetatm_ligands,
                list_monosaccharides_chain,
                lig_5ch_map,
                alt_res,
                alt_choice,
                alt_struct,
                all_dict_renumber_ligand_all,
                n_H,
            )
        dict_assembly_infos_tmp[assembly_number] = (
            pdb_id,
            list_hetatm_toatom,
            list_hetatm_ligands,
            list_monosaccharides_chain,
            lig_5ch_map,
            alt_res,
            alt_choice,
            alt_struct,
            all_dict_renumber_ligand_all,
            n_H,
        )
        cmd.delete(obj_name)
        filename = f"{pdb_id}.cif".lower()
        if os.path.exists(filename):
            os.remove(filename)
        cmd.reinitialize()
    write_global_info(
        info_file, len(start_assembly_ids), result, len(assembly_ids), tot_assembly_ids
    )
    os.rename(info_file, f"{pdb_dir}/{pdb_id}__infos_assemblies.txt")
    return tot_assembly_ids


PROTEIN_RESIDUES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}
RNA_RESIDUES = {"A", "U", "C", "G", "I"}
DNA_RESIDUES = {"DA", "DT", "DG", "DC"}
NUCLEIC_RESIDUES = RNA_RESIDUES | DNA_RESIDUES
WATER_NAMES = {"HOH", "H2O", "DOD"}
IONS = {"NA", "K", "MG", "CA", "ZN", "FE", "CU", "MN", "NI", "CO", "CL", "BR", "I"}
ATOMIC_WEIGHTS = {
    "H": 1.008,
    "HE": 4.0026,
    "LI": 7.0,
    "BE": 9.012183,
    "B": 10.81,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.99840316,
    "NE": 20.18,
    "NA": 22.99840316,
    "MG": 24.305,
    "AL": 26.981538,
    "SI": 28.085,
    "P": 30.973762,
    "S": 32.07,
    "CL": 35.45,
    "AR": 39.9,
    "K": 39.0983,
    "CA": 40.08,
    "SC": 44.95591,
    "TI": 47.867,
    "V": 50.9415,
    "CR": 51.996,
    "MN": 54.93804,
    "FE": 55.84,
    "CO": 58.93319,
    "NI": 58.693,
    "CU": 63.55,
    "ZN": 65.4,
    "GA": 69.723,
    "GE": 72.63,
    "AS": 74.92159,
    "SE": 78.97,
    "BR": 79.9,
    "KR": 83.8,
    "RB": 85.468,
    "SR": 87.62,
    "Y": 88.90584,
    "ZR": 91.22,
    "NB": 92.90637,
    "MO": 95.95,
    "TC": 96.90636,
    "RU": 101.1,
    "RH": 102.9055,
    "PD": 106.42,
    "AG": 107.868,
    "CD": 112.41,
    "IN": 114.818,
    "SN": 118.71,
    "SB": 121.76,
    "TE": 127.6,
    "I": 126.9045,
    "XE": 131.29,
    "CS": 132.905452,
    "BA": 137.33,
    "LA": 138.9055,
    "CE": 140.116,
    "PR": 140.90766,
    "ND": 144.24,
    "PM": 132.905452,
    "SM": 150.4,
    "EU": 151.964,
    "GD": 157.25,
    "TB": 158.92535,
    "DY": 162.5,
    "HO": 164.93033,
    "ER": 167.26,
    "TM": 168.93422,
    "YB": 173.05,
    "LU": 174.9667,
    "HF": 178.49,
    "TA": 180.9479,
    "W": 183.84,
    "RE": 186.207,
    "OS": 190.2,
    "IR": 192.22,
    "PT": 195.08,
    "AU": 196.96657,
    "HG": 200.59,
    "TL": 204.383,
    "PB": 207,
    "BI": 208.9804,
    "PO": 208.98243,
    "AT": 209.98715,
    "RN": 222.01758,
    "D": 1.008,
}
ION_CHARGES = {
    "NA": +1,
    "K": +1,
    "MG": +2,
    "CA": +2,
    "ZN": +2,
    "FE": +2,
    "CU": +2,
    "MN": +2,
    "NI": +2,
    "CO": +2,
    "CL": -1,
    "BR": -1,
    "I": -1,
}

# ---------------------------------------------------------------------------
# mmCIF parsing and ligand-neighbor analysis
# ---------------------------------------------------------------------------


def normalize_element(elem: str) -> str:
    """Normalize element symbol; treat deuterium as hydrogen."""
    if elem.upper() == "D":
        return "H"
    return elem.upper()


def compute_formula(elements: list[str]) -> str:
    """Compute a chemical formula string from a list of elements."""
    normalized = [normalize_element(e) for e in elements]
    counts = Counter(normalized)
    return "".join(
        (
            f"{elem}{(counts[elem] if counts[elem] > 1 else '')}"
            for elem in sorted(counts)
        )
    )


def compute_molecular_weight(elements: list[str]) -> float:
    """Compute molecular weight from a list of elements."""
    return sum((ATOMIC_WEIGHTS.get(normalize_element(e), 0) for e in elements))


def estimate_charge(resname: str, elements: list[str]) -> int:
    """Estimate the ionic charge of a residue if known."""
    return ION_CHARGES.get(resname, 0)


def distance(coord1: tuple, coord2: tuple) -> float:
    """Compute Euclidean distance between two 3D coordinates."""
    return math.sqrt(sum(((a - b) ** 2 for a, b in zip(coord1, coord2))))


def group_atoms_by_molecule(atoms: list) -> tuple[dict, dict]:
    """Group atoms by molecule (resname, chain, resid, ins_code)."""
    molecules = defaultdict(list)
    mol_tags = {}
    for resname, chain, chainauth, resid, ins_code, element, coord, tag in atoms:
        key = (resname, chain, chainauth, resid, ins_code)
        molecules[key].append((element, coord, tag))
        mol_tags.setdefault(key, set()).add(tag)
    return (molecules, mol_tags)


def parse_mmcif_atoms(filepath: str) -> list:
    """Parse atoms from a mmCIF file and return a list of tuples."""
    with open(filepath, "r") as f:
        lines = f.readlines()
    start = None
    for i, line in enumerate(lines):
        if line.strip().startswith("_atom_site."):
            start = i
            break
    if start is None:
        return []
    columns = []
    i = start
    while i < len(lines) and lines[i].startswith("_atom_site."):
        columns.append(lines[i].strip())
        i += 1
    col_idx = {col: idx for idx, col in enumerate(columns)}
    data = []
    for line in lines[i:]:
        if line.startswith("#") or line.strip() == "":
            continue
        parts = re.findall("'[^']*'|\\\"[^\\\"]*\\\"|\\S+", line.strip())
        if len(parts) < len(columns):
            continue
        data.append(parts)
    atoms = []
    seen_atoms = set()
    group_idx = col_idx.get("_atom_site.group_PDB")
    resname_idx = col_idx.get("_atom_site.label_comp_id")
    chain_idx = col_idx.get("_atom_site.label_asym_id")
    chainauth_idx = col_idx.get("_atom_site.auth_asym_id")
    auth_seq_idx = col_idx.get("_atom_site.label_seq_id")
    atom_name_idx = col_idx.get("_atom_site.label_atom_id")
    type_symbol_idx = col_idx.get("_atom_site.type_symbol")
    ins_code_idx = col_idx.get("_atom_site.pdbx_PDB_ins_code")
    x_idx = col_idx.get("_atom_site.Cartn_x")
    y_idx = col_idx.get("_atom_site.Cartn_y")
    z_idx = col_idx.get("_atom_site.Cartn_z")
    for d in data:
        try:
            group_val = d[group_idx] if group_idx is not None else d[0]
            group_val = group_val.strip().upper()
            if group_val not in ("ATOM", "HETATM"):
                continue
            if (
                resname_idx is None
                or chain_idx is None
                or auth_seq_idx is None
                or (atom_name_idx is None)
                or (x_idx is None)
                or (y_idx is None)
                or (z_idx is None)
            ):
                continue
            resname = d[resname_idx].upper()
            chain = d[chain_idx]
            chainauth = d[chainauth_idx]
            auth_resid = d[auth_seq_idx]
            atom_name = d[atom_name_idx]
            element = None
            if type_symbol_idx is not None and type_symbol_idx < len(d):
                element = d[type_symbol_idx].upper()
            else:
                element = atom_name[0].upper() if atom_name else ""
            ins_code = (
                d[ins_code_idx]
                if ins_code_idx is not None and ins_code_idx < len(d)
                else "."
            )
            x = float(d[x_idx])
            y = float(d[y_idx])
            z = float(d[z_idx])
            atom_key = (resname, chain, chainauth, auth_resid, ins_code, atom_name)
            if atom_key not in seen_atoms:
                tag = "atom" if group_val == "ATOM" else "hetatm"
                atoms.append(
                    (
                        resname,
                        chain,
                        chainauth,
                        auth_resid,
                        ins_code,
                        element,
                        (x, y, z),
                        tag,
                    )
                )
                seen_atoms.add(atom_key)
        except Exception:
            continue
    return atoms


def count_nearby_residues(
    molecule_atoms: list,
    all_atoms: list,
    threshold: float = float(DISTANCE_THRESHOLD),
    nucleic: bool = False,
) -> tuple[int, list[str]]:
    """Count nearby residues within a distance threshold for protein or nucleic acids."""
    reference_set = NUCLEIC_RESIDUES if nucleic else PROTEIN_RESIDUES
    neighbors = set()
    for e1, c1, r1, id1, ch1 in molecule_atoms:
        for e2, c2, r2, id2, ch2, chauth2, tag2 in all_atoms:
            if (r1, id1, ch1) == (r2, id2, ch2):
                continue
            if (
                distance(c1, c2) <= threshold
                and r2 in reference_set
                and (tag2 == "atom")
            ):
                neighbors.add((r2, id2, ch2))
    neighbor_list = sorted(
        [f"{res}:{resid}:{chain}" for res, resid, chain in neighbors]
    )
    return (len(neighbors), neighbor_list)


def analyze_pdb(
    pdb_id: str, list_assemblies: list, pdb_root_dir: str, tot_assembly_ids: str
) -> list[dict]:
    """
    Analyze mmCIF files for a given PDB ID and return ligand/hetatm results.

    Args:
        pdb_id (str): PDB identifier.
        list_assemblies (int): List of assembly IDs to process.
        pdb_root_dir (str): Root directory where mmCIF files are stored.

    Returns:
        list[dict]: A list of dictionaries with ligand information, neighbors,
                    molecular weights, formulas, and nearby residues.
    """
    results = []
    for assembly in list_assemblies:
        assembly_dir = os.path.join(pdb_root_dir, pdb_id, str(assembly))
        cif_files = glob(os.path.join(assembly_dir, "*.cif"))
        for filepath in cif_files:
            atoms = parse_mmcif_atoms(filepath)
            if not atoms:
                continue
            molecules, mol_tags = group_atoms_by_molecule(atoms)
            all_atoms_full = [
                (elem, coord, resname, resid, chain, chainauth, tag)
                for resname, chain, chainauth, resid, ins_code, elem, coord, tag in atoms
            ]
            protein_chains = {
                chain
                for (
                    resname,
                    chain,
                    chainauth,
                    resid,
                    ins_code,
                ), tags in mol_tags.items()
                if "atom" in tags
            }
            resname_counts = Counter(
                (
                    resname
                    for resname, chain, chainauth, resid, ins_code in molecules
                    if "hetatm"
                    in mol_tags.get((resname, chain, chainauth, resid, ins_code), set())
                )
            )
            mol_info = {}
            for mol_id, atom_data in molecules.items():
                resname, chain, chainauth, resid, ins_code = mol_id
                elements = [e for e, _, _ in atom_data]
                coords = [c for _, c, _ in atom_data]
                tag_final = (
                    "hetatm" if "hetatm" in mol_tags.get(mol_id, set()) else "atom"
                )
                mol_info[mol_id] = {
                    "elements": elements,
                    "coords": coords,
                    "atoms": atom_data,
                    "resname": resname,
                    "chain": chain,
                    "resid": resid,
                    "occurrence": resname_counts[resname],
                    "tag": tag_final,
                }
            protein_name = os.path.splitext(os.path.basename(filepath))[0]
            chain_info_str = (
                f"{{{len(protein_chains)}:{','.join(sorted(protein_chains))}}}"
            )
            coords_all = [
                xyz for atom in mol_info.values() for _, xyz, _ in atom["atoms"]
            ]
            meta_all = [
                (res, ch, chauth, resid, elem, tag)
                for resid_info, atom in mol_info.items()
                for elem, xyz, tag in atom["atoms"]
                for res, ch, chauth, resid, _ in [resid_info]
            ]
            tree = cKDTree(coords_all)
            for mol_id, info in mol_info.items():
                resname, chain, chainauth, resid, ins_code = mol_id
                if (
                    info["tag"] == "atom"
                    or resname in WATER_NAMES
                    or resname in NUCLEIC_RESIDUES
                ):
                    continue
                elements = info["elements"]
                atoms_list = info["atoms"]
                charge = estimate_charge(resname, elements)
                lig_coords = [xyz for _, xyz, _ in atoms_list]
                idx_in_sphere = set(
                    (
                        idx
                        for pt in lig_coords
                        for idx in tree.query_ball_point(
                            pt, r=float(DISTANCE_THRESHOLD)
                        )
                    )
                )
                neighbors_dict = defaultdict(dict)
                neighbors_dict_auth = defaultdict(dict)
                peptide_neighbor_chains = set()
                for idx in idx_in_sphere:
                    (
                        neigh_res,
                        neigh_chain,
                        neigh_chain_auth,
                        neigh_resid,
                        neigh_elem,
                        neigh_tag,
                    ) = meta_all[idx]
                    neigh_coord = coords_all[idx]
                    if (neigh_res, neigh_chain, neigh_resid) == (resname, chain, resid):
                        continue
                    for lig_elem, lig_xyz, lig_tag in atoms_list:
                        dist = distance(lig_xyz, neigh_coord)
                        if dist <= float(DISTANCE_THRESHOLD):
                            if neigh_tag == "hetatm" and neigh_res not in WATER_NAMES:
                                prev = neighbors_dict[neigh_res].get(neigh_chain, dist)
                                neighbors_dict[neigh_res][
                                    f"{neigh_chain}_{neigh_resid}"
                                ] = min(prev, dist)
                                prev = neighbors_dict_auth[neigh_res].get(
                                    neigh_chain_auth, dist
                                )
                                neighbors_dict_auth[neigh_res][
                                    f"{neigh_chain_auth}_{neigh_resid}"
                                ] = min(prev, dist)
                            elif neigh_tag == "atom":
                                if (
                                    neigh_res in PROTEIN_RESIDUES
                                    or neigh_res in NUCLEIC_RESIDUES
                                ):
                                    peptide_neighbor_chains.add(neigh_chain)
                                elif neigh_res not in WATER_NAMES:
                                    prev = neighbors_dict[neigh_res].get(
                                        neigh_chain, dist
                                    )
                                    neighbors_dict[neigh_res][
                                        f"{neigh_chain}_{neigh_resid}"
                                    ] = min(prev, dist)
                                    prev = neighbors_dict_auth[neigh_res].get(
                                        neigh_chain_auth, dist
                                    )
                                    neighbors_dict_auth[neigh_res][
                                        f"{neigh_chain_auth}_{neigh_resid}"
                                    ] = min(prev, dist)
                molecule_atoms_full = [
                    (elem, coord, mol_id[0], mol_id[2], mol_id[1])
                    for elem, coord, _ in atoms_list
                ]
                nb_protein, list_protein = count_nearby_residues(
                    molecule_atoms_full, all_atoms_full
                )
                nb_nucleic, list_nucleic = count_nearby_residues(
                    molecule_atoms_full, all_atoms_full, nucleic=True
                )
                neighbors_str = []
                for res, chains in neighbors_dict.items():
                    for cid in chains:
                        neighbors_str.append(res)
                neighbors_str_chain = []
                for res, chains in neighbors_dict.items():
                    for cid in chains:
                        neighbors_str_chain.append(f"{res}_{cid}")
                neighbors_dist = []
                for res, chains in neighbors_dict.items():
                    for cid, dist in chains.items():
                        neighbors_dist.append(f"{res}_{cid}__{dist:.2f}")
                neighbors_str_auth = []
                for res, chains in neighbors_dict_auth.items():
                    for cid in chains:
                        neighbors_str_auth.append(f"{res}_{cid}")
                neighbors_dist_auth = []
                for res, chains in neighbors_dict_auth.items():
                    for cid, dist in chains.items():
                        neighbors_dist_auth.append(f"{res}_{cid}__{dist:.2f}")
                results.append(
                    {
                        "pdb_id": pdb_id,
                        "pdb_assembly_assemblies": f"{protein_name}_{tot_assembly_ids}",
                        "assembly": protein_name.split("_")[1],
                        "Ligand ID": resname,
                        "occurrence": info["occurrence"],
                        "Asym ID - nonpolymer": chain,
                        "Auth Asym ID - nonpolymer": chainauth,
                        "resid": resid,
                        "charge": charge,
                        "neighbors": neighbors_str,
                        "neighbors_chain_resid": neighbors_str_chain,
                        "neighbors_dist_A": neighbors_dist,
                        "neighbors_authchain": neighbors_str_auth,
                        "neighbors_authchain_dist_A": neighbors_dist_auth,
                        "protein_chains": chain_info_str,
                        "residues_near_A": nb_protein,
                        "residues_near_list_A": "; ".join(list_protein),
                        "nucleic_near_A": nb_nucleic,
                        "nucleic_near_list_A": "; ".join(list_nucleic),
                        "ligand_chain": f"{resname}_{chain}",
                        "ligand_chain_resid": f"{resname}_{chain}_{resid}",
                    }
                )
    return results


def create_files_neighbors(df_from_cif: "pd.DataFrame", pdb_id: str, pdb_root_dir: str):
    """
    Create mmCIF files filtered to keep only ligand neighbors.

    Args:
        df_from_cif (pd.DataFrame): DataFrame with columns:
            - pdb_assembly
            - resid, chain, ligand
            - neighbors (neighbor list string)
        pdb_id (str): PDB identifier.
        pdb_root_dir (str): Root directory where mmCIF files are stored.
    """
    list_files_assemblies = list(set(df_from_cif["assembly"]))
    for assembly in list_files_assemblies:
        df_assembly = df_from_cif[df_from_cif["assembly"] == assembly]
        cif_path = os.path.join(
            pdb_root_dir, pdb_id, assembly, f"{pdb_id}_{assembly}.cif"
        )
        if not os.path.exists(cif_path):
            LOGGER.warning("File not found: %s", cif_path)
            continue
        for _, row in df_assembly.iterrows():
            cmd.reinitialize()
            cmd.load(cif_path, "obj")
            selections = []
            resi_lig = str(row["resid"])
            resn_lig = str(row["Ligand ID"])
            chain_lig = str(row["Asym ID - nonpolymer"])
            voisin = row["neighbors_authchain"]
            if voisin and voisin != "":
                neighbors_list = [v for v in voisin]
                for v in neighbors_list:
                    resn, chain, resi = v.split("_")
                    sel = f"(obj and hetatm and chain {chain} and resi {resi} and resn {resn})"
                    selections.append(sel)
                keep_selection = " or ".join(selections)
                cmd.remove(f"obj and not ({keep_selection})")
                out_dir = os.path.join(pdb_root_dir, pdb_id, assembly)
                out_file = os.path.join(
                    out_dir,
                    f"{resn_lig}_{chain_lig}_{resi_lig}",
                    f"{pdb_id}_{assembly}_{resn_lig}_{chain_lig}_{resi_lig}__neighbors.cif",
                )
                cmd.save(out_file, "obj")
            if not voisin:
                continue
            cmd.reinitialize()


# ---------------------------------------------------------------------------
# PockDrug server result extraction
# ---------------------------------------------------------------------------


def get_cif_centers(cif_path):
    """
    Load a CIF file and compute geometric centers of non-polymeric ligands.

    This function parses the `_atom_site` section of a CIF file and
    calculates mean coordinates for HETATM groups.

    Parameters
    ----------
    cif_path : str
        Path to CIF file.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing ligand centers.
    """
    cif_cache = {}
    atoms, headers = ([], [])
    with open(cif_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("_atom_site."):
                headers.append(line)
            elif headers and (not line.startswith("_")):
                parts = line.split()
                if parts:
                    atoms.append(parts)
    if not headers or not atoms:
        return pd.DataFrame()
    n_cols = len(headers)
    clean_atoms = [parts[:n_cols] + [None] * (n_cols - len(parts)) for parts in atoms]
    df = pd.DataFrame(clean_atoms, columns=headers)
    nonpol = df[df["_atom_site.group_PDB"] == "HETATM"].copy()
    for col in ["_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z"]:
        nonpol[col] = pd.to_numeric(nonpol[col], errors="coerce")
    centers = (
        nonpol.groupby(
            [
                "_atom_site.label_comp_id",
                "_atom_site.label_asym_id",
                "_atom_site.auth_asym_id",
                "_atom_site.label_seq_id",
            ]
        )[["_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z"]]
        .mean()
        .reset_index()
    )
    cif_cache[cif_path] = centers
    return centers


def ligand_center_from_pdb(pdb_path):
    """
    Compute the geometric center of ligand atoms from a PDB file.

    Parameters
    ----------
    pdb_path : str
        Path to ligand PDB file.

    Returns
    -------
    numpy.ndarray or None
        Mean XYZ coordinates of ligand atoms, or None if unavailable.
    """
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    coords.append(
                        [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                    )
                except Exception:
                    continue
    return np.mean(coords, axis=0) if coords else None


def find_best_match(lig_center, centers):
    """
    Find the closest CIF ligand to a given ligand center.

    Parameters
    ----------
    lig_center : numpy.ndarray
        XYZ center of ligand from PDB.
    centers : pandas.DataFrame
        CIF ligand centers.

    Returns
    -------
    tuple
        label_asym_id, auth_asym_id,
        ligand_chain_resid, ligand_authchain_resid
    """
    centers["distance"] = np.linalg.norm(
        centers[
            ["_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z"]
        ].values
        - lig_center,
        axis=1,
    )
    best = centers.loc[centers["distance"].idxmin()]
    return (
        best["_atom_site.label_asym_id"],
        best["_atom_site.auth_asym_id"],
        best["ligand_chain_resid"],
        best["ligand_authchain_resid"],
    )


def load_matrix(path):
    """
    Load a tabulated SO matrix file if valid.

    The file must:
    - Exist
    - Contain more than one line
    - Not contain 'NA' in second line

    Parameters
    ----------
    path : str
        Path to matrix file.

    Returns
    -------
    pandas.DataFrame or None
        Loaded matrix or None if invalid.
    """
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        lines = f.readlines()
    if len(lines) <= 1 or "NA" in lines[1] or lines[0] == "\n":
        return None
    return pd.read_csv(path, sep="\t")


def safe_max_idx(df):
    """
    Compute idxmax and max while handling all-zero columns.

    If a column contains only zeros, its max and idxmax are set to None.

    Parameters
    ----------
    df : pandas.DataFrame

    Returns
    -------
    tuple
        idxmax_vals : pandas.Series
        max_vals : pandas.Series
    """
    max_vals = df.max()
    idxmax_vals = df.idxmax()
    idxmax_vals[max_vals == 0] = None
    max_vals[max_vals == 0] = None
    return (idxmax_vals, max_vals)



def pockdrug_result_columns() -> list[str]:
    """Return the SO-only PockDrug result columns used internally."""
    return [
        "ID",
        "Ligands",
        "Pockets_SOmax",
        "ScoreSO_PocketSOmax",
        "Asym_ID",
        "Auth_Asym_ID",
        "ligand_chain_resid",
        "ligand_authchain_resid",
    ]


def empty_pockdrug_result() -> pd.DataFrame:
    """Return an empty SO-only PockDrug result table."""
    return pd.DataFrame(columns=pockdrug_result_columns())


def add_ligand_identifier_columns(centers: pd.DataFrame) -> pd.DataFrame:
    """Add label-chain and auth-chain ligand identifiers to CIF centers."""
    if centers.empty:
        centers["ligand_chain_resid"] = []
        centers["ligand_authchain_resid"] = []
        return centers

    centers = centers.copy()
    centers["ligand_chain_resid"] = [
        f"{row['_atom_site.label_comp_id']}_{row['_atom_site.label_asym_id']}_{row['_atom_site.label_seq_id']}"
        for _, row in centers.iterrows()
    ]
    centers["ligand_authchain_resid"] = [
        f"{row['_atom_site.label_comp_id']}_{row['_atom_site.auth_asym_id']}_{row['_atom_site.label_seq_id']}"
        for _, row in centers.iterrows()
    ]
    return centers


def ligand_name_without_proximity_suffix(ligand_key: str) -> str:
    """Return the ligand identifier without the PockDrug proximity suffix."""
    suffix = f"_prox{DISTANCE_THRESHOLD.replace('.', '_')}"
    return str(ligand_key).split(suffix)[0]


def analyze_pockdrug_results(tmp_path_results, path_file_cif, pdb_id):
    """Analyze SO-only PockDrug output and map ligands to detected pockets.

    Only the SO overlap matrix is used. A ligand is considered unmatched when
    its maximum SO score is missing or equal to zero.
    """
    so_path = os.path.join(
        tmp_path_results,
        "0000",
        f"prox{DISTANCE_THRESHOLD.replace('.', '_')}",
        "SO.txt",
    )
    df_so = load_matrix(so_path)
    if not isinstance(df_so, pd.DataFrame):
        return empty_pockdrug_result(), {}, [], []

    so_idxmax, so_max = safe_max_idx(df_so)
    result = pd.DataFrame(
        {
            "ID": pdb_id,
            "Ligands": df_so.columns,
            "Pockets_SOmax": so_idxmax.values,
            "ScoreSO_PocketSOmax": so_max.values,
        }
    )

    centers = add_ligand_identifier_columns(get_cif_centers(path_file_cif))
    asym_list: list[str | None] = []
    auth_list: list[str | None] = []
    ligand_asym_list: list[str | None] = []
    ligand_auth_list: list[str | None] = []
    clean_ligand_names: list[str] = []

    for ligand_key in result["Ligands"]:
        ligand_name = ligand_name_without_proximity_suffix(str(ligand_key))
        clean_ligand_names.append(ligand_name)

        # The PockDrug server may rename ligand-proximity pockets using its own
        # internal numbering (for example BAT_A_1), while the cleaned mmCIF keeps
        # canonical residue identifiers (for example BAT_F_800). Therefore the
        # exact ligand file extracted from the cleaned assembly may be unavailable.
        # In that case, the ligand-based proximity pocket is used as a spatial
        # proxy to recover the closest ligand in the cleaned mmCIF.
        ligand_pdb = os.path.join(tmp_path_results, "0000", f"{ligand_name}.pdb")
        prox_pdb = os.path.join(
            tmp_path_results,
            "0000",
            f"prox{DISTANCE_THRESHOLD.replace('.', '_')}",
            f"{ligand_name}_prox{DISTANCE_THRESHOLD.replace('.', '_')}.pdb",
        )

        if centers.empty:
            asym_list.append(None)
            auth_list.append(None)
            ligand_asym_list.append(None)
            ligand_auth_list.append(None)
            continue

        ligand_center = None
        if os.path.isfile(ligand_pdb):
            ligand_center = ligand_center_from_pdb(ligand_pdb)

        if ligand_center is None and os.path.isfile(prox_pdb):
            ligand_center = ligand_center_from_pdb(prox_pdb)

        if ligand_center is None:
            asym_list.append(None)
            auth_list.append(None)
            ligand_asym_list.append(None)
            ligand_auth_list.append(None)
            continue

        asym, auth, ligand_asym, ligand_auth = find_best_match(ligand_center, centers)
        asym_list.append(asym)
        auth_list.append(auth)
        ligand_asym_list.append(ligand_asym)
        ligand_auth_list.append(ligand_auth)

    result["Ligands"] = clean_ligand_names
    result["Asym_ID"] = asym_list
    result["Auth_Asym_ID"] = auth_list
    result["ligand_chain_resid"] = ligand_asym_list
    result["ligand_authchain_resid"] = ligand_auth_list

    ligands_without_so: list[str] = []
    ligand_to_so_pocket: dict[str, list[str]] = {}
    for _, row in result.iterrows():
        ligand_id = row["ligand_chain_resid"]
        pocket_id = row["Pockets_SOmax"]
        so_score = row["ScoreSO_PocketSOmax"]

        if ligand_id is None or pd.isna(ligand_id) or pocket_id is None or pd.isna(pocket_id):
            if ligand_id is not None and not pd.isna(ligand_id):
                ligands_without_so.append(ligand_id)
            continue

        so_score = float(so_score)
        if so_score <= 0:
            ligands_without_so.append(ligand_id)
            continue

        ligand_to_so_pocket[ligand_id] = [f"{pocket_id}___SOmax__SO_{so_score:.2f}"]

    return result, ligand_to_so_pocket, ligands_without_so, []


def descriptor_name_map() -> dict[str, str]:
    """Map server descriptor names to the compact PockDrug column names."""
    return {
        "Ne2 atom": "p_NE2_atom",
        "Diameter hull": "DIAMETER_HULL",
        "Polar residues": "p_polar_residues",
        "Smallest size": "SMALLEST_SIZE",
        "Nlys atom": "p_Nlys_atom",
        "Ntrp atom": "p_Ntrp_atom",
        "Aromatic residues": "p_aromatic_residues",
        "Volume hull": "VOLUME_HULL",
        "Otyr atom": "p_Otyr_atom",
        "Nb RES": "C_RESIDUES",
        "Surface hull": "SURFACE_HULL",
        "Ooh atom": "p_Ooh_atom",
        "Hydrophobic kyte": "hydrophobic_kyte",
        "Radius cylinder": "RADIUS_CYLINDER",
        "Aliphatic residues": "p_aliphatic_residues",
        "Nd1 atom": "p_ND1_atom",
        "Hydrophobic residues": "p_hydrophobic_residues",
        "Score Drugg": "Score_Drugg",
        "Confidence": "Confidence_Drugg",
    }


def read_descriptor_file(path: str) -> pd.DataFrame:
    """Read a two-column PockDrug descriptor file into one-row dataframe."""
    if not os.path.isfile(path):
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", header=None).set_index(0).T
    df.columns = df.columns.str.replace("pocket_", "", regex=False)
    return df


def copy_if_exists(source: str, destination: str) -> bool:
    """Copy a file when it exists and return whether the copy was done."""
    if not os.path.isfile(source):
        LOGGER.warning("Missing expected file: %s", source)
        return False
    os.makedirs(os.path.dirname(destination), exist_ok=True)
    shutil.copyfile(source, destination)
    return True


def remove_intermediate_pdb_directory(intermediate_root: str, pdb_id: str) -> None:
    """Delete the intermediate working directory for one processed PDB entry."""
    pdb_dir = os.path.join(ensure_trailing_separator(intermediate_root), pdb_id.lower())
    if os.path.isdir(pdb_dir):
        shutil.rmtree(pdb_dir)


def remove_empty_directory(path: str | os.PathLike[str]) -> None:
    """Remove a directory when it exists and is empty."""
    path = str(path)
    if os.path.isdir(path) and not os.listdir(path):
        os.rmdir(path)


def extract_ligand_pdb_from_mmcif(
    cif_path: str | os.PathLike[str], ligand_chain_resid: str, output_pdb: str | os.PathLike[str]
) -> bool:
    """Extract one ligand from a cleaned mmCIF file and write a minimal PDB file.

    The ligand identifier must use the canonical label-based form
    ``<resname>_<label_asym_id>_<label_seq_id>``, for example ``BAT_F_800``.
    This fallback is used when the PockDrug server has renamed the ligand
    internally and no exact server ligand PDB file is available.
    """
    cif_path = Path(cif_path)
    output_pdb = Path(output_pdb)

    if not cif_path.is_file():
        return False

    try:
        ligand_id, label_chain, label_resid = str(ligand_chain_resid).split("_", 2)
    except ValueError:
        return False

    with open(cif_path, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    start = None
    for index, line in enumerate(lines):
        if line.strip().startswith("_atom_site."):
            start = index
            break
    if start is None:
        return False

    columns = []
    index = start
    while index < len(lines) and lines[index].strip().startswith("_atom_site."):
        columns.append(lines[index].strip())
        index += 1

    col_idx = {column: i for i, column in enumerate(columns)}
    required = [
        "_atom_site.group_PDB",
        "_atom_site.label_atom_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_seq_id",
        "_atom_site.type_symbol",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
    ]
    if any(column not in col_idx for column in required):
        return False

    output_lines = []
    serial = 1
    for line in lines[index:]:
        if line.startswith("#") or not line.strip():
            continue
        parts = re.findall(r"'[^']*'|\"[^\"]*\"|\S+", line.strip())
        if len(parts) < len(columns):
            continue

        group = parts[col_idx["_atom_site.group_PDB"]].strip().upper()
        resname = parts[col_idx["_atom_site.label_comp_id"]].strip().upper()
        chain = parts[col_idx["_atom_site.label_asym_id"]].strip()
        resid = parts[col_idx["_atom_site.label_seq_id"]].strip()

        if group != "HETATM" or resname != ligand_id or chain != label_chain or resid != label_resid:
            continue

        atom_name = parts[col_idx["_atom_site.label_atom_id"]].strip().strip("'\"")
        element = parts[col_idx["_atom_site.type_symbol"]].strip().strip("'\"").upper()[:2]
        try:
            x = float(parts[col_idx["_atom_site.Cartn_x"]])
            y = float(parts[col_idx["_atom_site.Cartn_y"]])
            z = float(parts[col_idx["_atom_site.Cartn_z"]])
        except ValueError:
            continue

        pdb_chain = chain[:1] if chain else "A"
        try:
            pdb_resid = int(float(resid))
        except ValueError:
            pdb_resid = 1

        output_lines.append(
            f"HETATM{serial:5d} {atom_name[:4]:<4s} {resname[:3]:>3s} {pdb_chain:1s}"
            f"{pdb_resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"  1.00  0.00          {element:>2s}\n"
        )
        serial += 1

    if not output_lines:
        return False

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with open(output_pdb, "w", encoding="utf-8") as handle:
        handle.writelines(output_lines)
        handle.write("END\n")
    return True


def copy_essential_pockdrug_files(
    final_df, ligand_to_so_pocket, assembly, tmp_path_results, pdb_id, reference_cif=None
):
    """Copy SO-selected PockDrug files and descriptor tables.

    The exported intermediate files keep the legacy names expected by the
    downstream cleaner, but contain only the ligand-based pocket, the SO-best
    geometry-based pocket, the ligand file, and the two descriptor tables.
    """
    prox_tag = f"prox{DISTANCE_THRESHOLD.replace('.', '_')}"

    for _, row in final_df.iterrows():
        ligand_key = row["Ligands"]
        ligand_chain_resid = row["ligand_chain_resid"]
        if ligand_chain_resid not in ligand_to_so_pocket:
            continue

        ligand_name, ligand_chain, ligand_resid = str(ligand_chain_resid).split("_", 2)
        ligand_dir = os.path.join(dir_files, pdb_id, str(assembly), str(ligand_chain_resid))
        os.makedirs(ligand_dir, exist_ok=True)

        for pocket_info in ligand_to_so_pocket[ligand_chain_resid]:
            geom_pocket_id = pocket_info.split("___", 1)[0]
            so_match = re.search(r"SO_([0-9.]+)", pocket_info)
            so_score = float(so_match.group(1)) if so_match else None

            prox_atm = os.path.join(
                tmp_path_results,
                "0000",
                prox_tag,
                f"{ligand_key}_{prox_tag}.pdb",
            )
            prox_res = os.path.join(
                tmp_path_results,
                "0000",
                prox_tag,
                f"{ligand_key}_{prox_tag}_res.pdb",
            )
            prox_des = os.path.join(
                tmp_path_results,
                "0000",
                prox_tag,
                f"{ligand_key}_{prox_tag}.des",
            )
            geom_atm = os.path.join(
                tmp_path_results,
                "0000",
                "protein_out",
                "pockets",
                f"{geom_pocket_id}.pdb",
            )
            geom_res = os.path.join(
                tmp_path_results,
                "0000",
                "protein_out",
                "pockets",
                f"{geom_pocket_id}_res.pdb",
            )
            geom_des = os.path.join(
                tmp_path_results,
                "0000",
                "protein_out",
                "pockets",
                f"{geom_pocket_id}.des",
            )

            copy_if_exists(
                prox_atm,
                os.path.join(
                    ligand_dir,
                    f"{pdb_id}_{assembly}_{ligand_chain_resid}__pocket_prox_atm.pdb",
                ),
            )
            copy_if_exists(
                prox_res,
                os.path.join(
                    ligand_dir,
                    f"{pdb_id}_{assembly}_{ligand_chain_resid}__pocket_prox_res.pdb",
                ),
            )
            overlap = f"SOmax___SO_{so_score:.2f}" if so_score is not None else "SOmax___SO_NA"
            copy_if_exists(
                geom_atm,
                os.path.join(
                    ligand_dir,
                    f"{pdb_id}_{assembly}_{ligand_chain_resid}__{overlap}__pocket_geom_atm.pdb",
                ),
            )
            copy_if_exists(
                geom_res,
                os.path.join(
                    ligand_dir,
                    f"{pdb_id}_{assembly}_{ligand_chain_resid}__{overlap}__pocket_geom_res.pdb",
                ),
            )

            df_prox_des = read_descriptor_file(prox_des)
            if not df_prox_des.empty:
                df_prox_des["estimation"] = "prox"
                df_prox_des["pdb_id"] = str(pdb_id)
                df_prox_des["assembly"] = assembly
                df_prox_des["Ligand ID"] = str(ligand_name)
                df_prox_des["Asym ID - nonpolymer"] = ligand_chain
                df_prox_des["resid"] = ligand_resid
                df_prox_des.to_csv(
                    os.path.join(
                        ligand_dir,
                        f"{pdb_id}_{assembly}_{ligand_chain_resid}__pocket_prox_des.csv",
                    ),
                    sep="\t",
                    index=False,
                )

            df_geom_des = read_descriptor_file(geom_des)
            if not df_geom_des.empty:
                df_geom_des["estimation"] = "geom"
                df_geom_des["pdb_id"] = str(pdb_id)
                df_geom_des["assembly"] = assembly
                df_geom_des["Ligand ID"] = str(ligand_name)
                df_geom_des["Asym ID - nonpolymer"] = ligand_chain
                df_geom_des["resid"] = ligand_resid
                df_geom_des["SO"] = so_score
                df_geom_des.to_csv(
                    os.path.join(
                        ligand_dir,
                        f"{pdb_id}_{assembly}_{ligand_chain_resid}__{overlap}__pocket_geom_des.csv",
                    ),
                    sep="\t",
                    index=False,
                )

            source_ligand = os.path.join(tmp_path_results, "0000", f"{ligand_key}.pdb")
            destination_ligand = os.path.join(ligand_dir, f"{ligand_chain_resid}.pdb")

            # In server-output mode, the standalone ligand PDB is not always
            # provided under the SO column name (for example BAT_A_1.pdb).
            # This is expected: the prox pocket file is sufficient for spatial
            # remapping, and the final ligand file can be reconstructed from
            # the cleaned assembly CIF using the remapped identifier
            # (for example BAT_F_800). Avoid logging a false warning before
            # trying this intended fallback.
            if os.path.isfile(source_ligand):
                os.makedirs(os.path.dirname(destination_ligand), exist_ok=True)
                shutil.copyfile(source_ligand, destination_ligand)
            elif reference_cif and extract_ligand_pdb_from_mmcif(
                reference_cif, ligand_chain_resid, destination_ligand
            ):
                continue
            else:
                source_ligand_cif = os.path.join(tmp_path_results, "0000", f"{ligand_key}.cif")
                destination_ligand_cif = os.path.join(ligand_dir, f"{ligand_chain_resid}.cif")
                if not copy_if_exists(source_ligand_cif, destination_ligand_cif):
                    LOGGER.warning(
                        "Could not create ligand file for %s from server output or cleaned CIF.",
                        ligand_chain_resid,
                    )


def delete_tmp_pockdrug(tmpdir):
    """Delete a temporary directory used to normalize one server result."""
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)



def read_server_descriptor_table(path: str) -> pd.DataFrame:
    """Read a PockDrug server descriptor table with an unlabeled first column."""
    with open(path, encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if not lines:
        return pd.DataFrame()

    header = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:]]
    if rows and len(rows[0]) == len(header) + 1:
        header = ["identifier", *header]
    return pd.DataFrame(rows, columns=header)


def write_descriptor_table_as_des(table: pd.DataFrame, output_dir: str) -> None:
    """Convert a server descriptor table into one ``.des`` file per row."""
    if table.empty:
        return

    mapping = descriptor_name_map()
    id_col = table.columns[0]
    for _, row in table.iterrows():
        identifier = str(row[id_col])
        des_path = os.path.join(output_dir, f"{identifier}.des")
        with open(des_path, "w", encoding="utf-8") as handle:
            for source_col, target_col in mapping.items():
                if source_col not in table.columns:
                    continue
                handle.write(f"pocket_{target_col}\t{row[source_col]}\n")


def build_residue_file_from_atom_file(
    reference_structure: str,
    atom_file: str,
    output_file: str,
    distance_cutoff: float = 0.05,
) -> None:
    """Build a residue-level pocket file from an atom-level pocket file.

    If PyMOL cannot recover the corresponding residues, the atom-level file is
    copied as a conservative fallback so that downstream file names remain
    available.
    """
    if cmd is None or not os.path.isfile(reference_structure) or not os.path.isfile(atom_file):
        if os.path.isfile(atom_file):
            shutil.copyfile(atom_file, output_file)
        return

    try:
        cmd.reinitialize()
        cmd.load(reference_structure, "reference_structure")
        cmd.load(atom_file, "atom_pocket")
        selection = f"byres (reference_structure within {distance_cutoff} of atom_pocket)"
        if cmd.count_atoms(selection) == 0:
            shutil.copyfile(atom_file, output_file)
        else:
            cmd.save(output_file, selection, state=-1)
    except Exception as exc:
        LOGGER.warning("Could not build residue-level file for %s: %s", atom_file, exc)
        if os.path.isfile(atom_file):
            shutil.copyfile(atom_file, output_file)
    finally:
        if cmd is not None:
            cmd.reinitialize()


def extract_ligand_file_from_structure(
    reference_structure: str, ligand_key: str, output_file: str
) -> None:
    """Extract a ligand PDB file from the cleaned assembly structure."""
    if cmd is None or not os.path.isfile(reference_structure):
        return

    ligand_id, chain_id, residue_id = str(ligand_key).split("_", 2)
    try:
        cmd.reinitialize()
        cmd.load(reference_structure, "reference_structure")
        selection = (
            f"reference_structure and hetatm and resn {ligand_id} "
            f"and chain {chain_id} and resi {residue_id}"
        )
        if cmd.count_atoms(selection) > 0:
            cmd.save(output_file, selection, state=-1)
    except Exception as exc:
        LOGGER.warning("Could not extract ligand %s: %s", ligand_key, exc)
    finally:
        if cmd is not None:
            cmd.reinitialize()


def normalize_server_pockdrug_results(
    server_zip: str | os.PathLike[str],
    normalized_root: str | os.PathLike[str],
    reference_structure: str,
) -> str:
    """Normalize a PockDrug server ZIP to the internal working layout.

    The server provides compact archives. This function recreates only the files
    required downstream: ``SO.txt``, ligand-based pockets, geometry-based
    pockets, descriptor ``.des`` files, residue-level pocket files, and ligand
    PDB files extracted from the cleaned assembly.
    """
    server_zip = str(server_zip)
    normalized_root = str(normalized_root)
    prox_tag = f"prox{DISTANCE_THRESHOLD.replace('.', '_')}"
    root_0000 = os.path.join(normalized_root, "0000")
    prox_dir = os.path.join(root_0000, prox_tag)
    geom_dir = os.path.join(root_0000, "protein_out", "pockets")
    os.makedirs(prox_dir, exist_ok=True)
    os.makedirs(geom_dir, exist_ok=True)

    with zipfile.ZipFile(server_zip) as archive:
        archive.extract("SO.txt", prox_dir)

        with zipfile.ZipFile(io.BytesIO(archive.read("pocket_PPE.zip"))) as geom_zip:
            geom_zip.extractall(geom_dir)

        with zipfile.ZipFile(io.BytesIO(archive.read(f"pocket_{prox_tag}.zip"))) as prox_zip:
            prox_zip.extractall(prox_dir)

    geom_table_path = os.path.join(geom_dir, "pocket_PPE.txt")
    prox_table_path = os.path.join(prox_dir, f"pocket_{prox_tag}.txt")

    if os.path.isfile(geom_table_path):
        geom_table = read_server_descriptor_table(geom_table_path)
        write_descriptor_table_as_des(geom_table, geom_dir)

    ligand_keys: list[str] = []
    if os.path.isfile(prox_table_path):
        prox_table = read_server_descriptor_table(prox_table_path)
        write_descriptor_table_as_des(prox_table, prox_dir)
        first_col = prox_table.columns[0]
        for full_key in prox_table[first_col].astype(str):
            ligand_key = ligand_name_without_proximity_suffix(full_key)
            ligand_keys.append(ligand_key)
            extract_ligand_file_from_structure(
                reference_structure,
                ligand_key,
                os.path.join(root_0000, f"{ligand_key}.pdb"),
            )

    for atom_file in glob(os.path.join(prox_dir, "*.pdb")):
        if atom_file.endswith("_res.pdb"):
            continue
        output_file = atom_file.replace(".pdb", "_res.pdb")
        build_residue_file_from_atom_file(reference_structure, atom_file, output_file)

    for atom_file in glob(os.path.join(geom_dir, "*_atm.pdb")):
        output_file = atom_file.replace(".pdb", "_res.pdb")
        build_residue_file_from_atom_file(reference_structure, atom_file, output_file)

    return normalized_root


def resolve_server_results_path(
    server_results: str | os.PathLike[str] | None, pdb_id: str, assembly: str
) -> str | None:
    """Resolve a server ZIP path for one assembly."""
    if server_results is None:
        return None

    path = Path(server_results)
    if path.is_file():
        return str(path)

    if not path.is_dir():
        raise FileNotFoundError(f"PockDrug server result path does not exist: {path}")

    patterns = [
        f"{pdb_id}_{assembly}*.zip",
        f"{pdb_id.upper()}_{assembly}*.zip",
        f"*{assembly}*.zip",
        "*.zip",
    ]
    for pattern in patterns:
        matches = sorted(path.glob(pattern))
        if len(matches) == 1:
            return str(matches[0])

    raise FileNotFoundError(
        f"Could not identify a unique server ZIP for {pdb_id} assembly {assembly} in {path}."
    )

def ensure_trailing_separator(path: str) -> str:
    """Return a path string ending with the OS-specific path separator."""
    path = str(path)
    return path if path.endswith(os.sep) else path + os.sep


def initialize_done_file(output_dir: str) -> None:
    """Create the pipeline completion summary file with its header if needed."""
    done_file = os.path.join(output_dir, "pdb_done.txt")
    if os.path.exists(done_file):
        return

    with open(done_file, "w", encoding="utf-8") as handle:
        handle.write(
            "PDB_ID\tlist_validated_ligands\tnumber_validated_ligands\t"
            "ligand_nucleic_near_A\tligands_SOnull\t"
            "PockDrug_failures\tall_ligands_fromrcsb\n"
        )


def load_additive_codes(
    additive_table_path: str | os.PathLike[str] | None = None,
) -> list[str]:
    """Load additive like or non primary ligand CCD/HET codes.

    By default, the function reads the bundled data table stored in
    ``data/Table1_additive_like_or_non_primary_ligands_CCD_HET_codes.csv``. The table may
    be semicolon-separated and contain either a ``CCD/HET code`` column or a
    ``HET code`` column. Missing files are tolerated so that the rest of the
    pipeline remains usable.
    """
    path = Path(additive_table_path or ADDITIVES_TABLE_PATH or DEFAULT_ADDITIVES_TABLE)

    if not path.exists():
        LOGGER.warning("Additive like or non primary ligand table not found: %s", path)
        return []

    try:
        table = pd.read_csv(
            path,
            sep=";",
            encoding="utf-8-sig",
            keep_default_na=False,
            dtype=str,
        )
    except Exception as exc:
        LOGGER.warning("Could not read semicolon-separated additive table %s: %s", path, exc)
        return []

    code_column = None
    for candidate in ("CCD/HET code", "HET code", "het_code", "code"):
        if candidate in table.columns:
            code_column = candidate
            break

    if code_column is None:
        code_column = table.columns[0]
        LOGGER.warning(
            "No explicit additive-code column found in %s. Using first column: %s",
            path,
            code_column,
        )

    additive_codes = (
        table[code_column]
        .dropna()
        .astype(str)
        .str.strip()
        .str.upper()
    )
    return sorted({code for code in additive_codes if code})



# ---------------------------------------------------------------------------
# Final CSV and clean-file export
# ---------------------------------------------------------------------------

COLS_WEBSITE = [
    "PDB_ID",
    "num_assembly",
    "Deposition_Date",
    "DOI",
    "Title",
    "Structure_Title",
    "Experimental_Method",
    "Crystal_Growth_Procedure",
    "Refinement_Resolution_(Å)",
    "Source_Organism",
    "Macromolecule_Name",
    "Structure_Keywords",
    "Molecular_Weight_(Entity)(KDa)",
    "oligomeric",
    "Entity_ID_polymer",
    "Auth_Asym_ID",
    "Asym_ID",
    "List_of_Unique_Monosaccharides",
    "UniProt-Accession_Code(s)",
    "Ligand_ID",
    "Ligand_Name",
    "Ligand_Formula",
    "Ligand_MW",
    "Ligand_SMILES",
    "Ligand_InChI",
    "entity_id",
    "Chain_Asym_ID_-_ligand",
    "Chain_Auth_Asym_ID_-_ligand",
    "num_ligand",
    "occurrence",
    "Ligand-DrugBank-Accession_Code(s)",
    "Ligand-PubChem-Accession_Code(s)",
    "Ligand-ChEBI-Accession_Code(s)",
    "Ligand-ChEMBL-Accession_Code(s)",
    "Ligand-BindingDb-Accession_Code(s)",
    "Ligand-SureChEMBL-Accession_Code(s)",
    "Ligand-CCDC/CSD-Accession_Code(s)",
    "Ligand-ZINC-Accession_Code(s)",
    "crippen_clog_p",
    "num_rotatable_bonds",
    "num_aromatic_rings",
    "num_hbd",
    "num_hba",
    "tpsa",
    "lipinski_hba",
    "lipinski_hbd",
    "PDBe-Cofactor_like",
    "Lipinski",
    "Additive_like_or_non_primary_ligand",
    "protein_chains_environment",
    "neighbors",
    "hydrophobic_kyte_LB",
    "hydrophobic_kyte_GB",
    "p_hydrophobic_residues_LB",
    "p_hydrophobic_residues_GB",
    "p_polar_residues_LB",
    "p_polar_residues_GB",
    "p_aromatic_residues_LB",
    "p_aromatic_residues_GB",
    "p_aliphatic_residues_LB",
    "p_aliphatic_residues_GB",
    "p_Otyr_atom_LB",
    "p_Otyr_atom_GB",
    "p_NE2_atom_LB",
    "p_NE2_atom_GB",
    "p_Nlys_atom_LB",
    "p_Nlys_atom_GB",
    "p_Ntrp_atom_LB",
    "p_Ntrp_atom_GB",
    "p_Ooh_atom_LB",
    "p_Ooh_atom_GB",
    "p_ND1_atom_LB",
    "p_ND1_atom_GB",
    "SURFACE_HULL_LB",
    "SURFACE_HULL_GB",
    "DIAMETER_HULL_LB",
    "DIAMETER_HULL_GB",
    "VOLUME_HULL_LB",
    "VOLUME_HULL_GB",
    "SMALLEST_SIZE_LB",
    "SMALLEST_SIZE_GB",
    "RADIUS_CYLINDER_LB",
    "RADIUS_CYLINDER_GB",
    "C_RESIDUES_LB",
    "C_RESIDUES_GB",
    "SO_GB",
]

CSV_RENAME_MAP = {
    "pdb_id": "PDB_ID",
    "PDB ID": "PDB_ID",
    "assembly": "num_assembly",
    "num assembly": "num_assembly",
    "Deposition Date": "Deposition_Date",
    "Experimental Method": "Experimental_Method",
    "Crystal Growth Procedure": "Crystal_Growth_Procedure",
    "Refinement Resolution (Å)": "Refinement_Resolution_(Å)",
    "Source Organism": "Source_Organism",
    "Macromolecule Name": "Macromolecule_Name",
    "Structure Title": "Structure_Title",
    "Stucture Keywords": "Structure_Keywords",
    "Structure Keywords": "Structure_Keywords",
    "Molecular Weight (Entity)(KDa)": "Molecular_Weight_(Entity)(KDa)",
    "Entity ID_polymer": "Entity_ID_polymer",
    "List of Unique Monosaccharides": "List_of_Unique_Monosaccharides",
    "UniProt-Accession Code(s)": "UniProt-Accession_Code(s)",
    "Ligand ID": "Ligand_ID",
    "Ligand Name": "Ligand_Name",
    "Ligand Formula": "Ligand_Formula",
    "Ligand MW": "Ligand_MW",
    "Ligand SMILES": "Ligand_SMILES",
    "InChI": "Ligand_InChI",
    "Ligand InChI": "Ligand_InChI",
    "Asym ID - nonpolymer": "Chain_Asym_ID_-_ligand",
    "Auth Asym ID - nonpolymer": "Chain_Auth_Asym_ID_-_ligand",
    "Chain Asym ID - ligand": "Chain_Asym_ID_-_ligand",
    "Chain Auth Asym ID - ligand": "Chain_Auth_Asym_ID_-_ligand",
    "resid": "num_ligand",
    "num ligand": "num_ligand",
    "protein_chains": "protein_chains_environment",
    "Ligand-DrugBank-Accession Code(s)": "Ligand-DrugBank-Accession_Code(s)",
    "Ligand-PubChem-Accession Code(s)": "Ligand-PubChem-Accession_Code(s)",
    "Ligand-ChEBI-Accession Code(s)": "Ligand-ChEBI-Accession_Code(s)",
    "Ligand-ChEMBL-Accession Code(s)": "Ligand-ChEMBL-Accession_Code(s)",
    "Ligand-BindingDb-Accession Code(s)": "Ligand-BindingDb-Accession_Code(s)",
    "Ligand-SureChEMBL-Accession Code(s)": "Ligand-SureChEMBL-Accession_Code(s)",
    "Ligand-CCDC/CSD-Accession Code(s)": "Ligand-CCDC/CSD-Accession_Code(s)",
    "Ligand-ZINC-Accession Code(s)": "Ligand-ZINC-Accession_Code(s)",
}

KEY_COLUMNS = [
    "PDB_ID",
    "num_assembly",
    "Ligand_ID",
    "Chain_Asym_ID_-_ligand",
    "num_ligand",
]


def standardize_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """Rename columns to the final website naming convention."""
    df = df.rename(columns=CSV_RENAME_MAP).copy()
    df.columns = [CSV_RENAME_MAP.get(col, col.replace(" ", "_")) for col in df.columns]
    return df


def normalize_merge_keys(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize key columns used to merge final-info and descriptors."""
    df = df.copy()
    for col in ["PDB_ID", "Ligand_ID", "Chain_Asym_ID_-_ligand"]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip().str.upper()
    for col in ["num_assembly", "num_ligand"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def read_final_infos_for_single_pdb(path: str) -> pd.DataFrame:
    """Read and standardize the metadata table produced for one PDB."""
    df = pd.read_csv(path, sep="/", dtype={"pdb_id": str, "Ligand ID": str})
    df = standardize_column_names(df)
    return normalize_merge_keys(df)


def read_descriptor_tables_for_single_pdb(
    pdb_dir: str,
    suffix: str,
    column_suffix: str,
    must_contain: str | None = None,
) -> pd.DataFrame:
    """Read one-pDB descriptor tables and suffix non-key descriptor columns."""
    paths = sorted(glob(os.path.join(pdb_dir, "*", "*", f"*{suffix}")))
    if must_contain is not None:
        paths = [path for path in paths if must_contain in os.path.basename(path)]

    frames = []
    for path in paths:
        try:
            df = pd.read_csv(path, sep="\t", dtype={"pdb_id": str, "Ligand ID": str})
        except Exception as exc:
            LOGGER.warning("Could not read descriptor table %s: %s", path, exc)
            continue

        df = standardize_column_names(df)
        df = normalize_merge_keys(df)
        rename = {
            col: col if col in KEY_COLUMNS else f"{col}_{column_suffix}"
            for col in df.columns
        }
        frames.append(df.rename(columns=rename))

    if not frames:
        return pd.DataFrame(columns=KEY_COLUMNS)
    return pd.concat(frames, ignore_index=True)


def clean_scalar_values(df: pd.DataFrame) -> pd.DataFrame:
    """Clean scalar values in the final website table."""
    df = df.copy()
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].replace({np.nan: None, "[]": None, "NAN": "NA"})

    for col in ["Chain_Asym_ID_-_ligand", "Asym_ID"]:
        if col in df.columns:
            df[col] = df[col].replace(["NAN"], pd.NA).fillna("NA").astype(str)

    for col in ["Refinement_Resolution_(Å)", "Molecular_Weight_(Entity)(KDa)"]:
        if col in df.columns:
            if col == "Refinement_Resolution_(Å)":
                df[col] = df[col].apply(mean_if_list_string)
            df[col] = pd.to_numeric(df[col], errors="coerce").round(2)
            df[col] = df[col].astype(object).where(df[col].notna(), None)

    float_columns = [
        "Ligand_MW",
        "crippen_clog_p",
        "tpsa",
        "hydrophobic_kyte_LB",
        "hydrophobic_kyte_GB",
        "p_hydrophobic_residues_LB",
        "p_hydrophobic_residues_GB",
        "p_polar_residues_LB",
        "p_polar_residues_GB",
        "p_aromatic_residues_LB",
        "p_aromatic_residues_GB",
        "p_aliphatic_residues_LB",
        "p_aliphatic_residues_GB",
        "p_Otyr_atom_LB",
        "p_Otyr_atom_GB",
        "p_NE2_atom_LB",
        "p_NE2_atom_GB",
        "p_Nlys_atom_LB",
        "p_Nlys_atom_GB",
        "p_Ntrp_atom_LB",
        "p_Ntrp_atom_GB",
        "p_Ooh_atom_LB",
        "p_Ooh_atom_GB",
        "p_ND1_atom_LB",
        "p_ND1_atom_GB",
        "SURFACE_HULL_LB",
        "SURFACE_HULL_GB",
        "DIAMETER_HULL_LB",
        "DIAMETER_HULL_GB",
        "VOLUME_HULL_LB",
        "VOLUME_HULL_GB",
        "SMALLEST_SIZE_LB",
        "SMALLEST_SIZE_GB",
        "RADIUS_CYLINDER_LB",
        "RADIUS_CYLINDER_GB",
        "SO_GB",
    ]
    for col in float_columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").round(2)
            df[col] = df[col].astype(object).where(df[col].notna(), None)

    integer_columns = [
        "num_rotatable_bonds",
        "num_aromatic_rings",
        "num_hba",
        "num_hbd",
        "lipinski_hba",
        "lipinski_hbd",
        "C_RESIDUES_LB",
        "C_RESIDUES_GB",
        "Entity_ID_polymer",
        "entity_id",
        "num_assembly",
        "num_ligand",
        "occurrence",
    ]
    for col in integer_columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
            df[col] = df[col].apply(lambda value: None if pd.isna(value) else int(value))

    return df


def mean_if_list_string(value):
    """Return the mean when a value is a stringified numeric list."""
    if isinstance(value, str) and value.startswith("["):
        try:
            values = ast.literal_eval(value)
            return float(sum(values) / len(values))
        except Exception:
            return value
    return value


def build_final_csv_for_single_pdb(pdb_id: str, pdb_dir: str) -> pd.DataFrame:
    """Build the final one-PDB CSV from metadata and pocket descriptors."""
    final_infos_path = os.path.join(pdb_dir, f"{pdb_id}__final_infos.csv")
    if not os.path.isfile(final_infos_path):
        LOGGER.warning("Final metadata table not found: %s", final_infos_path)
        return pd.DataFrame(columns=COLS_WEBSITE)

    df_info = read_final_infos_for_single_pdb(final_infos_path)
    df_lb = read_descriptor_tables_for_single_pdb(
        pdb_dir, suffix="pocket_prox_des.csv", column_suffix="LB"
    )
    df_gb = read_descriptor_tables_for_single_pdb(
        pdb_dir,
        suffix="pocket_geom_des.csv",
        column_suffix="GB",
        must_contain="SOmax",
    )

    merged = df_info.merge(df_lb, on=KEY_COLUMNS, how="left")
    merged = merged.loc[:, ~merged.columns.duplicated()]
    merged = merged.merge(df_gb, on=KEY_COLUMNS, how="left")
    merged = merged.loc[:, ~merged.columns.duplicated()]


    filter_columns = ["C_RESIDUES_LB", "SO_GB"]
    diagnostic_filter_columns = ["C_RESIDUES_LB", "C_RESIDUES_GB", "SO_GB"]
    for col in diagnostic_filter_columns:
        if col not in merged.columns:
            merged[col] = np.nan
        merged[col] = pd.to_numeric(merged[col], errors="coerce")

    valid_mask = (merged["C_RESIDUES_LB"] > 3) & (merged["SO_GB"] > 0)
    rejected = merged.loc[~valid_mask].copy()
    valid = merged.loc[valid_mask].copy()
    if valid.empty and not merged.empty:
        diagnostic_cols = [
            col for col in KEY_COLUMNS + diagnostic_filter_columns if col in merged.columns
        ]
        LOGGER.warning(
            "No row passed the final filters for %s. Candidate rows before filtering:\n%s",
            pdb_id,
            merged[diagnostic_cols].to_string(index=False),
        )
    remove_rejected_intermediate_files(pdb_dir, rejected, valid)

    valid = standardize_column_names(valid)
    for col in COLS_WEBSITE:
        if col not in valid.columns:
            valid[col] = None
    valid = clean_scalar_values(valid[COLS_WEBSITE])
    return valid


def ligand_prefix_from_row(row: pd.Series) -> str:
    """Return the file prefix for one ligand row."""
    assembly = int(row["num_assembly"]) if pd.notna(row["num_assembly"]) else row["num_assembly"]
    ligand_num = int(row["num_ligand"]) if pd.notna(row["num_ligand"]) else row["num_ligand"]
    return (
        f"{str(row['PDB_ID']).lower()}_{assembly}_"
        f"{str(row['Ligand_ID']).upper()}_"
        f"{str(row['Chain_Asym_ID_-_ligand']).upper()}_{ligand_num}"
    )


def ligand_dir_from_row(pdb_dir: str, row: pd.Series) -> str:
    """Return the intermediate ligand directory for one row."""
    assembly = int(row["num_assembly"]) if pd.notna(row["num_assembly"]) else row["num_assembly"]
    ligand_num = int(row["num_ligand"]) if pd.notna(row["num_ligand"]) else row["num_ligand"]
    ligand_dir_name = (
        f"{str(row['Ligand_ID']).upper()}_"
        f"{str(row['Chain_Asym_ID_-_ligand']).upper()}_{ligand_num}"
    )
    return os.path.join(pdb_dir, str(assembly), ligand_dir_name)


def remove_rejected_intermediate_files(
    pdb_dir: str, rejected: pd.DataFrame, valid: pd.DataFrame
) -> None:
    """Remove intermediate ligand folders and empty assembly folders rejected by filters."""
    if rejected.empty:
        return

    for _, row in rejected.iterrows():
        ligand_dir = ligand_dir_from_row(pdb_dir, row)
        if os.path.isdir(ligand_dir):
            shutil.rmtree(ligand_dir)

    valid_assemblies = set(valid["num_assembly"].dropna().astype(int).astype(str))
    for assembly_dir in glob(os.path.join(pdb_dir, "*")):
        if not os.path.isdir(assembly_dir):
            continue
        assembly_id = os.path.basename(assembly_dir)
        if assembly_id not in valid_assemblies:
            shutil.rmtree(assembly_dir)


def clean_final_output_directory(destination_dir: str) -> None:
    """Clean final output files while preserving server-input assemblies.

    The ``pockdrug_server_input`` directory is intentionally kept because users
    may need to submit several cleaned assemblies to the PockDrug web server,
    rerun a server job, or keep the exact uploaded structures for traceability.
    """
    preserved_names = {"pockdrug_server_input"}
    os.makedirs(destination_dir, exist_ok=True)

    for name in os.listdir(destination_dir):
        if name in preserved_names:
            continue
        path = os.path.join(destination_dir, name)
        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.remove(path)


def final_assembly_output_name(pdb_id: str, num_assembly: str | int | None) -> str:
    """Return the final output folder name for one finalized assembly."""
    pdb_id = str(pdb_id).lower()
    if num_assembly is None:
        return pdb_id
    return f"{pdb_id}_{str(num_assembly)}"


def export_clean_pdb_outputs(
    pdb_id: str,
    intermediate_root: str,
    clean_output_dir: str,
    final_df: pd.DataFrame,
    num_assembly: str | int | None = None,
) -> str:
    """Export only final clean files and the one-PDB website CSV."""
    pdb_id = pdb_id.lower()
    source_pdb_dir = os.path.join(intermediate_root, pdb_id)
    destination_dir = os.path.join(
        clean_output_dir, final_assembly_output_name(pdb_id, num_assembly)
    )

    clean_final_output_directory(destination_dir)

    final_csv = os.path.join(destination_dir, f"{pdb_id}__final_website.csv")
    if final_df.empty:
        pd.DataFrame(columns=COLS_WEBSITE).to_csv(final_csv, index=False, sep="	")
        return destination_dir

    for assembly in sorted(set(final_df["num_assembly"].dropna().astype(int).astype(str))):
        source_cif = os.path.join(source_pdb_dir, assembly, f"{pdb_id}_{assembly}.cif")
        destination_cif = os.path.join(destination_dir, f"{pdb_id}_{assembly}_cleaned.cif")
        copy_if_exists(source_cif, destination_cif)

    for _, row in final_df.iterrows():
        assembly = str(int(row["num_assembly"]))
        ligand = f"{row['Ligand_ID']}_{row['Chain_Asym_ID_-_ligand']}_{int(row['num_ligand'])}"
        prefix = f"{pdb_id}_{assembly}_{ligand}"
        ligand_dir = os.path.join(source_pdb_dir, assembly, ligand)
        if not os.path.isdir(ligand_dir):
            continue

        for filename in os.listdir(ligand_dir):
            source_file = os.path.join(ligand_dir, filename)
            if not os.path.isfile(source_file):
                continue

            new_name = None
            if "pocket_prox_atm" in filename:
                new_name = f"{prefix}__LB_pocket_atm.pdb"
            elif "pocket_prox_res" in filename:
                new_name = f"{prefix}__LB_pocket_res.pdb"
            elif "pocket_geom_atm" in filename and "SOmax" in filename:
                new_name = f"{prefix}__GB_pocket_atm.pdb"
            elif "pocket_geom_res" in filename and "SOmax" in filename:
                new_name = f"{prefix}__GB_pocket_res.pdb"
            elif filename.startswith(ligand) and filename.endswith((".pdb", ".cif")):
                extension = filename.rsplit(".", 1)[-1]
                new_name = f"{prefix}__ligand.{extension}"
            elif "neighbors.cif" in filename:
                new_name = f"{prefix}__ligand_environment.cif"

            if new_name is not None:
                shutil.copyfile(source_file, os.path.join(destination_dir, new_name))

    final_csv = os.path.join(destination_dir, f"{pdb_id}__final_website.csv")
    final_df.to_csv(final_csv, index=False, sep="\t")
    return destination_dir



def export_cleaned_assemblies_for_server(
    pdb_id: str,
    intermediate_root: str,
    clean_output_dir: str,
    assemblies: list[str | int],
) -> str:
    """Export cleaned assembly files that must be uploaded to the PockDrug server.

    Both PDB and CIF files are exported. The PockDrug server expects the PDB
    file, while the CIF file is kept as a traceable cleaned structural record.
    """
    pdb_id = pdb_id.lower()
    source_pdb_dir = os.path.join(intermediate_root, pdb_id)
    destination_dir = os.path.join(clean_output_dir, pdb_id, "pockdrug_server_input")

    if os.path.isdir(destination_dir):
        shutil.rmtree(destination_dir)
    os.makedirs(destination_dir, exist_ok=True)

    for assembly in sorted({str(assembly) for assembly in assemblies}, key=str):
        source_assembly_dir = os.path.join(source_pdb_dir, assembly)
        for extension in ("pdb", "cif"):
            source_file = os.path.join(
                source_assembly_dir,
                f"{pdb_id}_{assembly}.{extension}",
            )
            destination_file = os.path.join(
                destination_dir,
                f"{pdb_id}_{assembly}_cleaned.{extension}",
            )
            copy_if_exists(source_file, destination_file)

    return destination_dir


def print_pockdrug_server_submission_instructions(
    pdb_id: str,
    server_input_dir: str,
    assemblies: list[str | int],
    candidate_counts: dict[str, int],
) -> None:
    """Print the manual PockDrug server submission instructions."""
    print("\nPockDrug server submission required")
    print("==================================")
    print(f"PDB ID: {pdb_id.lower()}")
    print(f"Cleaned assembly files were written to: {server_input_dir}")
    for assembly in sorted({str(assembly) for assembly in assemblies}, key=str):
        count = candidate_counts.get(str(assembly), 0)
        print(
            f"For PDB ID {pdb_id.lower()}, assembly {assembly}: "
            f"{count} candidate ligand-binding entity/entries were found in the "
            "cleaned assembly."
        )
    print("\nUpload each cleaned assembly PDB file to:")
    print("https://pockdrug.rpbs.univ-paris-diderot.fr/cgi-bin/index.py?page=Druggability")
    print("\nServer settings:")
    print("1. Go to: Druggability Prediction using protein(s).")
    print("2. In Protein(s) information, choose: upload your PDB file > Browse.")
    print("3. Select both pocket estimation methods: fpocket and prox.")
    print(f"4. Set the ligand proximity threshold to {DISTANCE_THRESHOLD} Å.")
    print("5. Submit the job and download the resulting ZIP file.")
    print("\nThen run the finalization command with:")
    print(
        "python prepare_pdb_bank.py finalize "
        f"{pdb_id.lower()} --num-assembly <assembly_number> "
        "--clean-output-dir /path/to/results/ "
        "--pockdrug-server-results /path/to/server_result.zip"
    )


def report_validated_binding_site_counts(final_df: pd.DataFrame, pdb_id: str) -> None:
    """Print the number of retained binding-site entries per assembly."""
    print("\nFinal validated binding-site entries")
    print("====================================")
    if final_df.empty:
        print(
            f"For PDB ID {pdb_id.lower()}: 0 validated binding-site entries were retained "
            "after applying filters."
        )
        return

    counts = final_df.groupby("num_assembly").size().to_dict()
    for assembly, count in sorted(counts.items(), key=lambda item: int(item[0])):
        print(
            f"For PDB ID {pdb_id.lower()}, assembly {int(assembly)}: "
            f"{int(count)} validated binding-site entries were retained "
            "after applying filters."
        )


def prepare_pockdrug_server_inputs(
    pdb_id: str,
    distance_threshold: str,
    dir_files: str,
    tmp_dir: str,
    clean_output_dir: str,
) -> str | None:
    """Prepare cleaned assembly files for manual PockDrug server submission.

    This first-stage workflow stops before any pocket matching step. It writes
    cleaned assembly files, prints the server submission instructions, and then
    removes the intermediate working directory.
    """
    dir_files = ensure_trailing_separator(dir_files)
    tmp_dir = ensure_trailing_separator(tmp_dir)

    globals()["DISTANCE_THRESHOLD"] = str(distance_threshold)
    globals()["dir_files"] = dir_files
    globals()["tmp_dir"] = tmp_dir

    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(dir_files, exist_ok=True)

    entry = query_rcsb_assembly(pdb_id)
    if entry is None:
        LOGGER.error("No RCSB entry could be retrieved for %s.", pdb_id)
        return None

    infos_fromrcsb = process_pdb(entry)
    if "Ligand ID" not in infos_fromrcsb.columns:
        LOGGER.warning("No ligand was reported by RCSB for %s.", pdb_id)
        return None

    list_ligands, _ = list_saccharide_ligands(infos_fromrcsb)
    list_monosaccharides_chain = monosaccharides_chain(entry)
    prepare_files_id(
        pdb_id,
        entry,
        dir_files,
        list_ligands,
        list_monosaccharides_chain,
    )

    infos_path = os.path.join(dir_files, pdb_id, f"{pdb_id}__infos_assemblies.txt")
    infos_files = pd.read_csv(infos_path, sep="\t", comment="#")
    assemblies_with_ligands = list(
        infos_files[infos_files["LIG_RCSB"].astype(str) != "set()"]["ASM_ID"]
    )

    assemblies_without_ligands = list(
        infos_files[infos_files["LIG_RCSB"].astype(str) == "set()"]["ASM_ID"]
    )
    for assembly_nolig in assemblies_without_ligands:
        assembly_dir = os.path.join(dir_files, pdb_id, str(assembly_nolig))
        if os.path.isdir(assembly_dir):
            shutil.rmtree(assembly_dir)

    results = analyze_pdb(pdb_id, assemblies_with_ligands, dir_files, len(assemblies_with_ligands))
    cif_ligands = pd.DataFrame(results)
    if not cif_ligands.empty and "assembly" in cif_ligands.columns:
        candidate_counts = (
            cif_ligands.groupby(cif_ligands["assembly"].astype(str))
            .size()
            .astype(int)
            .to_dict()
        )
    else:
        candidate_counts = {str(assembly): 0 for assembly in assemblies_with_ligands}

    server_input_dir = export_cleaned_assemblies_for_server(
        pdb_id=pdb_id,
        intermediate_root=dir_files,
        clean_output_dir=clean_output_dir,
        assemblies=assemblies_with_ligands,
    )

    print_pockdrug_server_submission_instructions(
        pdb_id=pdb_id,
        server_input_dir=server_input_dir,
        assemblies=assemblies_with_ligands,
        candidate_counts=candidate_counts,
    )

    remove_intermediate_pdb_directory(dir_files, pdb_id)

    return server_input_dir


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def normalize_metadata_cif_merge_keys(*dataframes: pd.DataFrame) -> None:
    """Normalize metadata/CIF merge keys in-place.

    RCSB metadata may contain uppercase PDB IDs (for example, ``1RM8``),
    whereas user input and local file paths may use lowercase identifiers
    (for example, ``1rm8``). The merge between metadata and CIF-derived ligand
    information must therefore use normalized string keys.
    """
    string_columns = [
        "pdb_id",
        "Ligand ID",
        "ligand_chain",
        "Asym ID - nonpolymer",
        "Auth Asym ID - nonpolymer",
    ]

    for df in dataframes:
        if df is None or df.empty:
            continue

        for column in string_columns:
            if column not in df.columns:
                continue

            df[column] = (
                df[column]
                .astype(str)
                .str.strip()
                .str.upper()
                .replace({"NONE": None, "NAN": None, "": None})
            )


def launch_pdb_bank(
    pdb_id,
    distance_threshold,
    dir_files,
    tmp_dir,
    tmp_root,
    pockdrug_server_results,
    clean_output_dir,
    num_assembly,
):
    """
    Finalize one PDB entry from PockDrug server results.

    This function prepares the cleaned assemblies again, normalizes the
    downloaded PockDrug server ZIP file(s), keeps only SO-based pocket matches,
    applies the final filters, exports clean files, writes the website CSV, and
    deletes intermediate files.

    Parameters
    ----------
    pdb_id : str
        PDB identifier.
    distance_threshold : str
        Distance threshold used for ligand-proximity pocket detection.
    dir_files : str
        Intermediate working directory.
    tmp_dir : str
        Temporary directory used by PyMOL fetch operations.
    tmp_root : str
        Temporary root used to normalize PockDrug server ZIP outputs.
    pockdrug_server_results : str
        Path to one server ZIP result, or to a directory containing one ZIP per assembly.
    clean_output_dir : str
        Final output directory.
    num_assembly : str or int
        Assembly number to finalize. The final folder is named
        ``<pdb_id>_<num_assembly>``.

    Returns
    -------
    None
    """
    dir_files = ensure_trailing_separator(dir_files)
    tmp_dir = ensure_trailing_separator(tmp_dir)

    globals()["DISTANCE_THRESHOLD"] = str(distance_threshold)
    globals()["dir_files"] = dir_files
    globals()["tmp_dir"] = tmp_dir
    globals()["tmp_root"] = tmp_root

    if pockdrug_server_results is None:
        raise ValueError(
            "Finalization requires --pockdrug-server-results. "
            "Run the prepare-server-inputs command first, upload each cleaned "
            "assembly PDB to the PockDrug server, download the ZIP result(s), "
            "and then run finalize."
        )

    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(dir_files, exist_ok=True)
    initialize_done_file(dir_files)

    entry = query_rcsb_assembly(pdb_id)
    if entry is None:
        LOGGER.error("No RCSB entry could be retrieved for %s.", pdb_id)
        return None

    infos_fromrcsb = process_pdb(entry)
    ligands_chain_all_assemblies = []
    all_ligands_SOnull = []
    all_ligands_Fail_Pockdrug = []
    unknown_ligands = []
    validated_ligand_records = {}
    list_monosaccharides_chain = monosaccharides_chain(entry)
    if "Ligand ID" in infos_fromrcsb.columns:
        dict_ass_oligomeric = assembly_oligomeric(entry)
        infos_fromrcsb["ligand_chain"] = (
            infos_fromrcsb["Ligand ID"] + "_" + infos_fromrcsb["Asym ID - nonpolymer"]
        )
        all_ligands_fromrcsb = list(infos_fromrcsb["ligand_chain"])
        list_ligands, list_Monosaccharides = list_saccharide_ligands(infos_fromrcsb)
        tot_assembly_ids = prepare_files_id(
            pdb_id, entry, dir_files, list_ligands, list_monosaccharides_chain
        )
        infos_files = pd.read_csv(
            f"{dir_files}{pdb_id}/{pdb_id}__infos_assemblies.txt", sep="\t", comment="#"
        )
        list_assemblies = list(
            infos_files[infos_files["LIG_RCSB"].astype(str) != "set()"]["ASM_ID"]
        )
        requested_assembly = str(num_assembly)
        available_assemblies = {str(assembly) for assembly in list_assemblies}
        if requested_assembly not in available_assemblies:
            raise ValueError(
                f"Assembly {requested_assembly} is not available for PDB ID {pdb_id}. "
                f"Available assemblies with ligands are: {sorted(available_assemblies)}"
            )
        list_assemblies = [
            assembly for assembly in list_assemblies if str(assembly) == requested_assembly
        ]
        list_assemblies_nolig = list(
            infos_files[infos_files["LIG_RCSB"].astype(str) == "set()"]["ASM_ID"]
        )
        for assembly_nolig in list_assemblies_nolig:
            if os.path.exists(f"{dir_files}{pdb_id}/{assembly_nolig}"):
                shutil.rmtree(f"{dir_files}{pdb_id}/{assembly_nolig}")
        results = analyze_pdb(pdb_id, list_assemblies, dir_files, tot_assembly_ids)
        cif_ligands = pd.DataFrame(results)
        if "assembly" in list(cif_ligands.columns):
            cif_ligands["oligomeric"] = (
                cif_ligands["assembly"].astype(str).map(dict_ass_oligomeric)
            )
        else:
            cif_ligands["oligomeric"] = None
        if "ligand_chain" in list(cif_ligands.columns):
            list_monosaccharides_chain_resid = list(
                cif_ligands[
                    cif_ligands["ligand_chain"].isin(list_monosaccharides_chain)
                ]["ligand_chain_resid"]
            )
        else:
            list_monosaccharides_chain_resid = []
        if "nucleic_near_A" in list(cif_ligands.columns):
            ligand_nucleic_near_A = list(
                cif_ligands[cif_ligands["nucleic_near_A"] != 0]["ligand_chain_resid"]
            )
        else:
            ligand_nucleic_near_A = []
        for assembly in list_assemblies:
            path_file_cif = f"{dir_files}{pdb_id}/{assembly}/{pdb_id}_{assembly}.cif"
            path_file_pdb = f"{dir_files}{pdb_id}/{assembly}/{pdb_id}_{assembly}.pdb"

            server_zip = resolve_server_results_path(
                pockdrug_server_results, pdb_id, str(assembly)
            )
            timestamp = str(int(time.time() * 1000))
            tmpdir = os.path.join(
                tmp_root, f"pockdrug_server_{pdb_id}_{assembly}_{timestamp}"
            )
            tmp_path_results = os.path.join(tmpdir, "normalized")
            os.makedirs(tmp_path_results, exist_ok=True)
            normalize_server_pockdrug_results(server_zip, tmp_path_results, path_file_pdb)
            (
                final_df,
                dict_lig_pocketSO,
                ligands_SOnull,
                ligands_Fail_Pockdrug,
            ) = analyze_pockdrug_results(tmp_path_results, path_file_cif, pdb_id)

            # Some server outputs can contain ligand-based pocket columns that
            # cannot be mapped back to a ligand in the cleaned assembly, for
            # example when the ligand PDB file is missing or when the geometric
            # center cannot be matched to the mmCIF ligand centers. Such rows
            # cannot be exported as binding-site entries because their canonical
            # ligand identifier (<ligand>_<chain>_<residue>) is unavailable.
            if "ligand_chain_resid" in final_df.columns:
                missing_ligand_id = (
                    final_df["ligand_chain_resid"].isna()
                    | final_df["ligand_chain_resid"]
                    .astype(str)
                    .str.strip()
                    .isin(["", "None", "nan", "NaN", "<NA>"])
                )
                if missing_ligand_id.any():
                    failed_ligands = (
                        final_df.loc[missing_ligand_id, "Ligands"]
                        .dropna()
                        .astype(str)
                        .tolist()
                    )
                    all_ligands_Fail_Pockdrug.extend(
                        f"{ligand}:unmatched_to_cleaned_structure"
                        for ligand in failed_ligands
                    )
                    LOGGER.warning(
                        "PDB %s assembly %s: %d ligand-based PockDrug "
                        "entry/entries could not be mapped back to the "
                        "cleaned assembly and will be skipped: %s",
                        pdb_id,
                        assembly,
                        int(missing_ligand_id.sum()),
                        failed_ligands,
                    )
                    final_df = final_df.loc[~missing_ligand_id].copy()

            for ligand_SOnull in ligands_SOnull:
                if ligand_SOnull not in list_monosaccharides_chain_resid:
                    all_ligands_SOnull.append(ligand_SOnull)
            for ligand_Fail_Pockdrug in ligands_Fail_Pockdrug:
                if ligand_Fail_Pockdrug not in list_monosaccharides_chain_resid:
                    all_ligands_Fail_Pockdrug.append(ligand_Fail_Pockdrug)
            final_df = final_df[
                ~final_df["ligand_chain_resid"].isin(ligands_SOnull)
            ]
            copy_essential_pockdrug_files(
                final_df, dict_lig_pocketSO, assembly, tmp_path_results, pdb_id, path_file_cif
            )
            delete_tmp_pockdrug(tmpdir)
            if ligands_SOnull != []:
                for ligand_sonull in ligands_SOnull:
                    if os.path.exists(
                        f"{dir_files}{pdb_id}/{assembly}/{ligand_sonull}"
                    ):
                        shutil.rmtree(f"{dir_files}{pdb_id}/{assembly}/{ligand_sonull}")
            if ligand_nucleic_near_A != []:
                for ligand_contact in ligand_nucleic_near_A:
                    if os.path.exists(
                        f"{dir_files}{pdb_id}/{assembly}/{ligand_contact}"
                    ):
                        shutil.rmtree(
                            f"{dir_files}{pdb_id}/{assembly}/{ligand_contact}"
                        )
            if list_monosaccharides_chain_resid != []:
                for monosacc in list_monosaccharides_chain_resid:
                    if os.path.exists(f"{dir_files}{pdb_id}/{assembly}/{monosacc}"):
                        shutil.rmtree(f"{dir_files}{pdb_id}/{assembly}/{monosacc}")
            valid_ligand_ids = (
                final_df["ligand_chain_resid"]
                .dropna()
                .astype(str)
                .str.strip()
            )
            valid_ligand_ids = [
                lig
                for lig in valid_ligand_ids
                if lig and lig not in {"None", "nan", "NaN", "<NA>"}
            ]

            for _, valid_row in final_df.iterrows():
                label_ligand_resid = valid_row.get("ligand_chain_resid")
                if pd.isna(label_ligand_resid):
                    continue
                label_ligand_resid = str(label_ligand_resid).strip()
                if not label_ligand_resid or label_ligand_resid in {"None", "nan", "NaN", "<NA>"}:
                    continue

                label_parts = label_ligand_resid.split("_")
                if len(label_parts) < 3:
                    continue

                auth_ligand_resid = valid_row.get("ligand_authchain_resid")
                auth_ligand_resid = None if pd.isna(auth_ligand_resid) else str(auth_ligand_resid).strip()
                auth_parts = auth_ligand_resid.split("_") if auth_ligand_resid else []

                label_chain_key = label_ligand_resid.rsplit("_", 1)[0]
                auth_chain_key = (
                    auth_ligand_resid.rsplit("_", 1)[0]
                    if auth_ligand_resid and auth_ligand_resid not in {"None", "nan", "NaN", "<NA>"}
                    else label_chain_key
                )

                validated_ligand_records[label_chain_key] = {
                    "label_chain_key": label_chain_key,
                    "auth_chain_key": auth_chain_key,
                    "ligand_id": label_parts[0],
                    "label_chain": label_parts[1],
                    "auth_chain": auth_parts[1] if len(auth_parts) >= 3 else label_parts[1],
                }

            for lig in valid_ligand_ids:
                if lig.split("_")[0] in ["UNL", "UNX"]:
                    unknown_ligands.append(lig)
                    if os.path.exists(f"{dir_files}{pdb_id}/{assembly}/{lig}"):
                        shutil.rmtree(f"{dir_files}{pdb_id}/{assembly}/{lig}")
            for lig in valid_ligand_ids:
                if (
                    lig not in ligand_nucleic_near_A
                    and lig not in ligands_SOnull
                    and (lig not in list_monosaccharides_chain_resid)
                    and (lig not in unknown_ligands)
                ):
                    ligands_chain_all_assemblies.append(lig.rsplit("_", 1)[0])
        if ligands_chain_all_assemblies != []:
            valid_label_chain_set = set(ligands_chain_all_assemblies)
            valid_records = {
                key: value
                for key, value in validated_ligand_records.items()
                if key in valid_label_chain_set
            }
            auth_to_label_chain = {
                record["auth_chain_key"]: record["label_chain_key"]
                for record in valid_records.values()
            }
            all_accepted_chain_keys = valid_label_chain_set | set(auth_to_label_chain)

            cif_ligands = (
                cif_ligands[cif_ligands["ligand_chain"].isin(valid_label_chain_set)]
                .reset_index()
                .drop("index", axis=1)
            )

            infos_fromrcsb = infos_fromrcsb.copy()
            infos_fromrcsb["_label_chain_candidate"] = (
                infos_fromrcsb["Ligand ID"].astype(str)
                + "_"
                + infos_fromrcsb["Asym ID - nonpolymer"].astype(str)
            )
            infos_fromrcsb["_auth_chain_candidate"] = (
                infos_fromrcsb["Ligand ID"].astype(str)
                + "_"
                + infos_fromrcsb["Auth Asym ID - nonpolymer"].astype(str)
            )
            infos_fromrcsb = infos_fromrcsb[
                infos_fromrcsb["_label_chain_candidate"].isin(all_accepted_chain_keys)
                | infos_fromrcsb["_auth_chain_candidate"].isin(all_accepted_chain_keys)
            ].copy()

            for idx, row in infos_fromrcsb.iterrows():
                candidates = [
                    row.get("_label_chain_candidate"),
                    row.get("_auth_chain_candidate"),
                    row.get("ligand_chain"),
                ]
                matched_label_key = None
                for candidate in candidates:
                    if candidate in valid_label_chain_set:
                        matched_label_key = candidate
                        break
                    if candidate in auth_to_label_chain:
                        matched_label_key = auth_to_label_chain[candidate]
                        break

                if matched_label_key is None or matched_label_key not in valid_records:
                    continue

                record = valid_records[matched_label_key]
                infos_fromrcsb.at[idx, "ligand_chain"] = record["label_chain_key"]
                infos_fromrcsb.at[idx, "Asym ID - nonpolymer"] = record["label_chain"]
                infos_fromrcsb.at[idx, "Auth Asym ID - nonpolymer"] = record["auth_chain"]

            infos_fromrcsb = (
                infos_fromrcsb.drop(
                    columns=["_label_chain_candidate", "_auth_chain_candidate"],
                    errors="ignore",
                )
                .reset_index()
                .drop("index", axis=1)
            )

            infos_fromrcsb["List of Unique Monosaccharides"] = [
                list_monosaccharides_chain
            ] * len(infos_fromrcsb)
            if len(infos_fromrcsb) != len(set(ligands_chain_all_assemblies)):
                LOGGER.warning(
                    "Some ligands were not recovered in the RCSB metadata section for %s.",
                    pdb_id,
                )
            if len(cif_ligands) != len(set(ligands_chain_all_assemblies)):
                LOGGER.warning(
                    "Some ligands were not recovered in the CIF parsing section for %s.",
                    pdb_id,
                )
            create_files_neighbors(cif_ligands, pdb_id, dir_files)
            infos_fromrcsb = infos_fromrcsb.replace({np.nan: None})
            infos_fromrcsb = add_ligands_pdbe_infos(infos_fromrcsb)

            # Ensure metadata and CIF-derived tables use identical merge keys.
            normalize_metadata_cif_merge_keys(infos_fromrcsb, cif_ligands)

            strict_merge_keys = [
                "ligand_chain",
                "pdb_id",
                "Asym ID - nonpolymer",
                "Auth Asym ID - nonpolymer",
                "Ligand ID",
            ]
            infos_fromrcsb_cif = pd.merge(
                infos_fromrcsb,
                cif_ligands,
                on=strict_merge_keys,
                how="inner",
            )

            if infos_fromrcsb_cif.empty and not infos_fromrcsb.empty and not cif_ligands.empty:
                fallback_keys = ["ligand_chain", "pdb_id", "Ligand ID"]
                fallback = pd.merge(
                    infos_fromrcsb,
                    cif_ligands,
                    on=fallback_keys,
                    how="inner",
                    suffixes=("", "_cif"),
                )
                if not fallback.empty:
                    for col in ["Asym ID - nonpolymer", "Auth Asym ID - nonpolymer"]:
                        cif_col = f"{col}_cif"
                        if cif_col in fallback.columns:
                            fallback[col] = fallback[cif_col]
                            fallback = fallback.drop(columns=[cif_col])
                    infos_fromrcsb_cif = fallback

            infos_fromrcsb_cif = infos_fromrcsb_cif.drop(
                [
                    "residues_near_A",
                    "residues_near_list_A",
                    "nucleic_near_A",
                    "nucleic_near_list_A",
                    "ligand_chain",
                    "neighbors_chain_resid",
                    "neighbors_authchain",
                    "ligand_chain_resid",
                ],
                axis=1,
            )
            infos_fromrcsb_cif.to_csv(
                f"{dir_files}{pdb_id}/{pdb_id}__final_infos.csv", sep="/"
            )

            final_website_df = build_final_csv_for_single_pdb(
                pdb_id, os.path.join(dir_files, pdb_id)
            )
            export_clean_pdb_outputs(
                pdb_id=pdb_id,
                intermediate_root=dir_files,
                clean_output_dir=clean_output_dir,
                final_df=final_website_df,
                num_assembly=num_assembly,
            )
            report_validated_binding_site_counts(final_website_df, pdb_id)
            remove_intermediate_pdb_directory(dir_files, pdb_id)
            remove_empty_directory(tmp_root)

            with open(f"{dir_files}/pdb_done.txt", "a+") as f:
                f.write(
                    f"{pdb_id}\t{ligands_chain_all_assemblies}\t{len(ligands_chain_all_assemblies)}\t{ligand_nucleic_near_A}\t{all_ligands_SOnull}\t{all_ligands_Fail_Pockdrug}\t{all_ligands_fromrcsb}\n"
                )
        else:
            empty_final_df = pd.DataFrame(columns=COLS_WEBSITE)
            export_clean_pdb_outputs(
                pdb_id=pdb_id,
                intermediate_root=dir_files,
                clean_output_dir=clean_output_dir,
                final_df=empty_final_df,
                num_assembly=num_assembly,
            )
            report_validated_binding_site_counts(empty_final_df, pdb_id)
            remove_intermediate_pdb_directory(dir_files, pdb_id)
            remove_empty_directory(tmp_root)
            with open(f"{dir_files}/pdb_done.txt", "a+") as f:
                f.write(
                    f"{pdb_id}\t{ligands_chain_all_assemblies}\t{len(ligands_chain_all_assemblies)}\t{ligand_nucleic_near_A}\t{all_ligands_SOnull}\t{all_ligands_Fail_Pockdrug}\t{all_ligands_fromrcsb}\n"
                )
    else:
        empty_final_df = pd.DataFrame(columns=COLS_WEBSITE)
        export_clean_pdb_outputs(
            pdb_id=pdb_id,
            intermediate_root=dir_files,
            clean_output_dir=clean_output_dir,
            final_df=empty_final_df,
            num_assembly=num_assembly,
        )
        report_validated_binding_site_counts(empty_final_df, pdb_id)
        remove_intermediate_pdb_directory(dir_files, pdb_id)
        remove_empty_directory(tmp_root)
        with open(f"{dir_files}/pdb_done.txt", "a+") as f:
            f.write(
                f"{pdb_id}\t{ligands_chain_all_assemblies}\t{len(ligands_chain_all_assemblies)}\t{[]}\t{[]}\t{[]}\t{[]}\n"
            )


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments shared by the two workflow stages."""
    parser.add_argument("pdb_id", help="PDB identifier to process, for example 1ABC.")
    parser.add_argument(
        "--distance-threshold",
        default=DEFAULT_DISTANCE_THRESHOLD,
        help="Distance threshold used for ligand proximity and PockDrug server settings.",
    )
    parser.add_argument(
        "--clean-output-dir",
        required=True,
        help="Final directory where server inputs, clean files, and final CSV are written.",
    )
    parser.add_argument(
        "--additives-table",
        default=None,
        help=(
            "Optional additive like or non primary ligand table. Defaults to the bundled "
            "data/Table1_additive_like_or_non_primary_ligands_CCD_HET_codes.csv file."
        ),
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
        help="Logging verbosity.",
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for the two-stage server workflow."""
    parser = argparse.ArgumentParser(
        description=(
            "Prepare one PDB entry for the PockDrug server, then finalize the "
            "binding-site export from downloaded server ZIP results."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser(
        "prepare-server-inputs",
        help="Create cleaned assembly PDB files to upload to the PockDrug server.",
    )
    add_common_arguments(prepare_parser)

    finalize_parser = subparsers.add_parser(
        "finalize",
        help="Build final files and CSV from PockDrug server ZIP output(s).",
    )
    add_common_arguments(finalize_parser)
    finalize_parser.add_argument(
        "--num-assembly",
        default=None,
        help=(
            "Assembly number corresponding to the PockDrug server ZIP to finalize. "
            "The final output folder is named <pdb_id>_<num_assembly>."
        ),
    )
    finalize_parser.add_argument(
        "--pockdrug-server-results",
        required=True,
        help=(
            "Path to one PockDrug server ZIP result, or to a directory containing "
            "one ZIP result per assembly."
        ),
    )
    return parser.parse_args(argv)


def configure_runtime(args: argparse.Namespace) -> None:
    """Configure logging and the additive like or non primary ligand table path."""
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )

    global ADDITIVES_TABLE_PATH
    ADDITIVES_TABLE_PATH = args.additives_table or str(DEFAULT_ADDITIVES_TABLE)


def resolve_output_and_work_paths(clean_output_dir: str) -> tuple[str, dict[str, str]]:
    """Resolve the final output directory and internal temporary paths.

    The user-facing command only requires ``--clean-output-dir``. Relative
    paths are resolved against the current working directory. Temporary files
    are created below ``<clean-output-dir>/_temporary_work`` and the full
    temporary directory is removed at the end of each command.
    """
    output_dir = Path(clean_output_dir).expanduser().resolve()
    root = output_dir / "_temporary_work"
    work_paths = {
        "root": str(root),
        "intermediate_dir": str(root / "intermediate_work"),
        "tmp_fetch_dir": str(root / "tmp_fetch"),
        "server_tmp_dir": str(root / DEFAULT_SERVER_TMP_SUBDIR),
    }
    return str(output_dir), work_paths


def remove_internal_work_directories(work_paths: dict[str, str]) -> None:
    """Remove all internal temporary directories created by this workflow."""
    root = work_paths.get("root")
    if root and os.path.isdir(root):
        shutil.rmtree(root, ignore_errors=True)


def main(argv: list[str] | None = None) -> int:
    """Run either the preparation stage or the finalization stage."""
    args = parse_args(argv)
    configure_runtime(args)
    clean_output_dir, work_paths = resolve_output_and_work_paths(args.clean_output_dir)

    if args.command == "finalize" and args.num_assembly is None:
        print(
            "Which assembly number should be finalized? "
            "Please provide it with --num-assembly, for example: --num-assembly 1."
        )
        return 2

    if cmd is None:
        LOGGER.critical(
            "PyMOL is required to run this pipeline. Install PyMOL in the "
            "execution environment before processing structures."
        )
        return 1

    pdb_id = args.pdb_id.lower()

    try:
        if args.command == "prepare-server-inputs":
            prepare_pockdrug_server_inputs(
                pdb_id=pdb_id,
                distance_threshold=str(args.distance_threshold),
                dir_files=work_paths["intermediate_dir"],
                tmp_dir=work_paths["tmp_fetch_dir"],
                clean_output_dir=clean_output_dir,
            )
            return 0

        if args.command == "finalize":
            args.pockdrug_server_results = str(
                Path(args.pockdrug_server_results).expanduser().resolve()
            )
            launch_pdb_bank(
                pdb_id=pdb_id,
                distance_threshold=str(args.distance_threshold),
                dir_files=work_paths["intermediate_dir"],
                tmp_dir=work_paths["tmp_fetch_dir"],
                tmp_root=work_paths["server_tmp_dir"],
                pockdrug_server_results=args.pockdrug_server_results,
                clean_output_dir=clean_output_dir,
                num_assembly=str(args.num_assembly),
            )
            return 0

        raise ValueError(f"Unsupported command: {args.command}")
    finally:
        remove_internal_work_directories(work_paths)


if __name__ == "__main__":
    raise SystemExit(main())
