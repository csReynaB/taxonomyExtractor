"""
Author: Carlos S Reyna-Blanco
Email: carlos.reynablanco@meduniwien.ac.at
Date: 2025-07-02
Description: Script to retrieve organism taxonomy from uniprot and peptide metadata.
Dependencies: requests, xml.etree.ElementTree, pandas, tqdm
"""
import argparse
import requests
import xml.etree.ElementTree as ET
import pandas as pd
import time
import re
from tqdm import tqdm
import numpy as np
from urllib.parse import quote
import csv



def get_organism_from_uniprot(rep_id):
    url = f"https://www.uniprot.org/uniprot/{rep_id}.xml"
    try:
        response = requests.get(url)
        if response.status_code != 200:
            #print(f"Failed to retrieve {rep_id}: HTTP {response.status_code}")
            return None

        root = ET.fromstring(response.content)
        ns = {'up': 'http://uniprot.org/uniprot'}

        result = {'common': None, 'scientific': None, 'tax_id': None}

        for name in root.findall('.//up:organism/up:name', namespaces=ns):
            name_type = name.attrib.get('type', '').lower()
            if name_type == 'common':
                result['common'] = name.text
            elif name_type == 'scientific':
                result['scientific'] = name.text

        taxon = root.find('.//up:organism/up:dbReference[@type="NCBI Taxonomy"]', namespaces=ns)
        if taxon is not None:
            result['tax_id'] = taxon.attrib['id']

        return result if result['scientific'] else None
    except Exception as e:
        #print(f"Error processing {rep_id}: {e}")
        return None

def lookup_organism(rep):
    """
    Returns a dictionary with 'common', 'scientific', and 'tax_id' for the given rep ID.
    First tries UniProt API.
    If the rep ID contains underscores, tries splitting to find alternative candidates.
    """
    if not rep  or not isinstance(rep, str):
        return {'common': None, 'scientific': None, 'tax_id': None}

    # Try full rep ID
    result = get_organism_from_uniprot(rep)
    if result:
        return result

    # Decide on splitting logic
    if "_" in rep:
        parts = rep.split("_")
    else:
        parts = rep.split(", ")

    # Try parts of the rep ID (e.g., split by '_' or ', ')
    for candidate in map(str.strip, parts):#rep.split("_"):
        if not candidate:
            continue
        result = get_organism_from_uniprot(candidate)
        if result:
            return result
    #     elif candidate in manual_mapping:
    #         return {'common': None, 'scientific': manual_mapping[candidate], 'tax_id': None}
    #
    # # Fallback to manual mapping using full rep
    # if rep in manual_mapping:
    #     return {'common': None, 'scientific': manual_mapping[rep], 'tax_id': None}

    return {'common': None, 'scientific': None, 'tax_id': None}

def search_taxid_by_name(scientific_name):
    """
    Query NCBI for scientific_name first.
    If no hit, query UniProt taxonomy API and:
      - Try the primary scientificName field
      - Then try any synonyms or otherNames found in the UniProt entry
    Returns the first NCBI taxid found, or None.
    """
    def ncbi_lookup(name):
        term = quote(str(name))
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
            f"esearch.fcgi?db=taxonomy&retmode=json&term={term}[SCIN]"
        )
        try:
            r = requests.get(url, timeout=5)
            r.raise_for_status()
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            return ids[0] if ids else None
        except Exception:
            return None

    # Try NCBI on the given name
    taxid = ncbi_lookup(scientific_name)
    if taxid:
        return taxid

    # Fallback to UniProt taxonomy search
    up_url = (
        "https://rest.uniprot.org/taxonomy/search"
        f"?query={quote(str(scientific_name))}&format=json"
    )

    try:
        r = requests.get(up_url, timeout=5)
        r.raise_for_status()
        results = r.json().get("results", [])
        if not results:
            return None

        entry = results[0]
        # Try the primary scientificName from UniProt
        main_name = entry.get("scientificName")
        taxid = ncbi_lookup(main_name)
        if taxid:
            return taxid

        # Then try any synonyms or other names
        synonyms = entry.get("synonyms", []) + entry.get("otherNames", [])
        for syn in synonyms:
            taxid = ncbi_lookup(syn)
            if taxid:
                return taxid

    except Exception:
        pass

    return None


def get_lineage_from_taxid(tax_id):
    url = f"https://rest.uniprot.org/taxonomy/{tax_id}"
    try:
        response = requests.get(url)
        if response.status_code != 200:
            #print(f"Failed to retrieve taxonomy for {tax_id}")
            return {}

        data = response.json()
        lineage = {rank['rank'].lower(): rank['scientificName'] for rank in data.get('lineage', [])}
        lineage[data.get('rank').lower()] = data.get('scientificName')  # Add most specific name
        lineage['common'] = data.get('commonName', '')         # Add common name if available

        return lineage
    except Exception as e:
        #print(f"Error retrieving lineage: {e}")
        return {}

def get_taxid_from_protein_id(protein_id):
    """
    Given a protein accession (e.g. YP_009321702.1), fetch its TaxID and organism name from NCBI.
    """
    if not isinstance(protein_id, str) or " " in protein_id:
        return None

    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "protein",
        "id": protein_id,
        "retmode": "json"
    }

    try:
        response = requests.get(base_url, params=params, timeout=10)
        response.raise_for_status()
        data = response.json()
        uid = next(iter(data["result"]["uids"]))
        result = data["result"][uid]
        taxid = result.get("taxid")
        #organism = result.get("organism")
        return str(taxid) #{"tax_id": str(taxid), "organism": organism}
    except Exception as e:
        #print(f"[Error] Failed to retrieve TaxID for {protein_id}: {e}")
        return None

def extract_rep_id(description):
    """
    Extracts the repID from a description string.
    The function looks for a pattern like 'RepID=<rep_value>'.
    If no repID is found, returns None.
    """
    if not description:
        return None

    # Regular expression pattern: looks for "RepID=" followed by one or more alphanumeric or underscore characters.
    # Adjust the pattern if your repIDs can include other characters.
    match = re.search(r'RepID=([\w\-:]+)', description)
    if match:
        return match.group(1)
    else:
        return None


def clean_iedb_organism_name(name):
    if pd.isna(name):
        return None

    # Step 1: Split by '&' and take the first group
    first_block = name.split('&')[0].strip()

    # Step 2: Extract content inside the first set of square brackets
    match = re.search(r"\[([^\]]+)\]", first_block)
    if not match:
        return None

    content = match.group(1)

    # Step 3: Split by ';' and take the first entry
    content = content.split(';')[0].strip()

    # Step 4: Remove anything in parentheses
    content = re.sub(r"\s*\(.*?\)", "", content)

    # Step 5: Final strip
    return content.strip("'\" ")

# Apply to each row with fallback to phage_name
def get_organism_complete_name(row):
    primary = clean_iedb_organism_name(row['IEDB_organism_name'])
    if primary:
        return primary

    # Fallback to phage_name if valid
    phage_val = row.get('phage_name', None)
    if pd.notna(phage_val) and str(phage_val).strip().lower() not in ['false', '', 'nan']:
        return str(phage_val).strip()

    return None

def collect_valid_ids(row):
    rep_cols = ['allergome_uniprot', 'allergen_uniprot', 'iedb_uniprot', 'fummy_uniprot', 'gened_uniprot']
    valid_ids = []
    for col in rep_cols:
        val = row[col]
        if isinstance(val, str):
            cleaned = val.strip()
            if cleaned.lower() not in ['false', 'nan', '']:
                valid_ids.append(cleaned)
        elif pd.notna(val) and val is not False:
            valid_ids.append(str(val).strip())
    return ", ".join(valid_ids) if valid_ids else None


def get_first_peptide_id_valid(row):
    for val in row:
        if isinstance(val, str) and val.strip().lower() not in ['', 'false', 'nan']:
            return val.strip()
        elif pd.notna(val) and val is not False:
            return val
    return None

def extract_taxid(desc):
    match = re.search(r'TaxID=(\d+)', str(desc))
    return match.group(1) if match else None


def annotate_taxonomy(
    df,
    organism_col=None,
    rep_id_col=None,
    prot_id_col=None,
    taxid_col=None,
    method_priority=("organism", "rep_id", "prot_id", "tax_id"),  # explicit order
    outfile_path="taxonomy_output.csv",
    max_rows=None
):
    lineage_fields = ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'common']
    id_cache = {}

    # Map method name to column name
    method_cols = {
        "organism": organism_col,
        "rep_id": rep_id_col,
        "prot_id": prot_id_col,
        "tax_id": taxid_col
    }

    with open(outfile_path, "w", newline="", encoding="utf-8") as out_f:
        writer = csv.DictWriter(
            out_f,
            fieldnames=["peptide_name"] + lineage_fields
        )
        writer.writeheader()

        if max_rows:
            df = df.head(max_rows)

        for idx, row in tqdm(df.iterrows(), total=len(df), desc="Fetching taxonomy"):
            row_data = dict.fromkeys(lineage_fields, "")
            lineage = {}

            # ----- Prioritization logic -----
            for method in method_priority:
                col = method_cols.get(method)
                if not col or col not in row:
                    continue
                val = row[col]
                if pd.isna(val):
                    continue

                if method == "organism" and pd.notna(val) and isinstance(val, str) and val.strip(): #if method == "organism" and str(val).strip() and val is not pd.NA:
                    if val in id_cache:
                        tax_id = id_cache[val]
                        if tax_id in id_cache:
                            lineage = id_cache[tax_id]
                        else:
                            lineage = get_lineage_from_taxid(tax_id) or {}
                            id_cache[tax_id] = lineage
                    else:
                        tax_id = search_taxid_by_name(val)
                        if tax_id:
                            id_cache[val] = tax_id
                            lineage = get_lineage_from_taxid(tax_id) or {}
                            id_cache[tax_id] = lineage

                elif method == "rep_id" and pd.notna(val):
                    rep_data = lookup_organism(val)
                    row_data["common"] = rep_data.get("common", "")
                    tax_id = rep_data.get("tax_id")
                    if tax_id:
                        if tax_id in id_cache:
                            lineage = id_cache[tax_id]
                        else:
                            lineage = get_lineage_from_taxid(tax_id) or {}
                            id_cache[tax_id] = lineage

                elif method == "prot_id" and pd.notna(val):
                    if val in id_cache:
                        tax_id = id_cache[val]
                        if tax_id in id_cache:
                            lineage = id_cache[tax_id]
                        else:
                            lineage = get_lineage_from_taxid(tax_id) or {}
                            id_cache[tax_id] = lineage
                    else:
                        tax_id = get_taxid_from_protein_id(val)
                        if tax_id:
                            id_cache[val] = tax_id
                            lineage = get_lineage_from_taxid(tax_id) or {}
                            id_cache[tax_id] = lineage

                elif method == "tax_id" and pd.notna(val):
                    tax_id = val
                    if tax_id in id_cache:
                        lineage = id_cache[tax_id]
                    else:
                        lineage = get_lineage_from_taxid(tax_id) or {}
                        id_cache[tax_id] = lineage

                if lineage:
                    break
            # Write output
            for key in lineage_fields:
                if key in lineage:
                    row_data[key] = lineage[key]

            if not any(row_data.values()):
                continue

            writer.writerow({"peptide_name": idx, **row_data})
            out_f.flush()

def parse_priority_list(val):
    return [v.strip() for v in val.split(",") if v.strip()]

def main():
    parser = argparse.ArgumentParser(description="Annotate a peptide metadata CSV with taxonomy lineages.")

    parser.add_argument("--input", required=True, help="Path to the input CSV file")
    parser.add_argument("--output", default="taxonomy_output.csv", help="Path to save the annotated CSV file")
    parser.add_argument("--organism_col", default=None, help="Column name for scientific names")
    parser.add_argument("--rep_id_col", default=None, help="Column name for representative IDs")
    parser.add_argument("--prot_id_col", default=None, help="Column name for protein IDs")
    parser.add_argument("--taxid_col", default=None, help="Column name for TaxIDs")

    #parser.add_argument("--method_priority", default="organism,rep_id,tax_id", type=list_of_strings,
    #                    help="Comma-separated list of method priority (possible valies: organism,rep_id,tax_id)")
    parser.add_argument("--method_priority", type=parse_priority_list, default=["organism", "rep_id", "prot_id", "tax_id"],
                        help="Comma-separated list of method priority (e.g. 'organism,tax_id')")

    parser.add_argument("--max_rows", type=int, default=None, help="Limit number of rows processed (for testing)")

    args = parser.parse_args()

    dtype_map = {}
    if args.taxid_col:
        dtype_map[args.taxid_col] = str

    df = pd.read_csv(args.input, index_col=0, low_memory=False, dtype=dtype_map)
#    if args.taxid_col and args.taxid_col in df.columns:
#        df[args.taxid_col] = df[args.taxid_col].apply(lambda x: str(int(x)) if pd.notna(x) else None)

    annotate_taxonomy(
        df=df,
        organism_col=args.organism_col,
        rep_id_col=args.rep_id_col,
        prot_id_col=args.prot_id_col,
        taxid_col=args.taxid_col,
        method_priority=args.method_priority,
        outfile_path=args.output,
        max_rows=args.max_rows
    )

if __name__ == "__main__":
    main()
