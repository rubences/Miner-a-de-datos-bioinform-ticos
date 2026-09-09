from __future__ import annotations

from pathlib import Path
from typing import Optional
from urllib.parse import quote

import pandas as pd
import requests
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from rdkit import Chem
from rdkit.Chem import Descriptors


PDB_IDS = ["1tup", "2xyz", "3def", "4ogq", "5jkl", "6mno", "7pqr", "8stu", "9vwx", "10yza"]


def _safe_get_json(url: str, timeout: int = 30) -> Optional[dict]:
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        return response.json()
    except requests.RequestException as error:
        print(f"[ERROR] No se pudo consultar {url}: {error}")
        return None
    except ValueError as error:
        print(f"[ERROR] La respuesta no es JSON válido para {url}: {error}")
        return None


def download_pdb_cif_files(pdb_ids: list[str], output_dir: Path) -> tuple[list[str], list[str]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded, failed = [], []

    for pdb_id in pdb_ids:
        pdb_id = pdb_id.lower()
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        target = output_dir / f"{pdb_id}.cif"
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            target.write_bytes(response.content)
            downloaded.append(pdb_id)
            print(f"[OK] Descargado {target}")
        except requests.RequestException as error:
            failed.append(pdb_id)
            print(f"[ERROR] No se pudo descargar {pdb_id}.cif: {error}")

    return downloaded, failed


def get_uniprot_id_from_pdb(pdb_id: str) -> Optional[str]:
    entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.lower()}"
    entry = _safe_get_json(entry_url)
    if not entry:
        return None

    entity_ids = (
        entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
    )

    for entity_id in entity_ids:
        poly_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id.lower()}/{entity_id}"
        polymer_data = _safe_get_json(poly_url)
        if not polymer_data:
            continue
        references = (
            polymer_data.get("rcsb_polymer_entity_container_identifiers", {}).get(
                "reference_sequence_identifiers", []
            )
        )
        for ref in references:
            if ref.get("database_name") == "UniProt":
                return ref.get("database_accession")
    return None


def get_uniprot_dataframe(uniprot_id: str) -> pd.DataFrame:
    data = _safe_get_json(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json")
    if not data:
        return pd.DataFrame(
            columns=[
                "Uniprot_id",
                "Fecha_publicacion",
                "Fecha_modificacion",
                "Revisado",
                "Nombre_del_gen",
                "Sinónimos",
                "Organismo",
                "PDB_ids",
            ]
        )

    audit = data.get("entryAudit", {})
    genes = data.get("genes", [])
    gene_name = ""
    gene_synonyms: list[str] = []
    if genes:
        gene_name = genes[0].get("geneName", {}).get("value", "")
        gene_synonyms = [syn.get("value", "") for syn in genes[0].get("synonyms", []) if syn.get("value")]

    pdb_ids = [
        xref.get("id")
        for xref in data.get("uniProtKBCrossReferences", [])
        if xref.get("database") == "PDB" and xref.get("id")
    ]

    row = {
        "Uniprot_id": data.get("primaryAccession", uniprot_id),
        "Fecha_publicacion": audit.get("firstPublicDate", ""),
        "Fecha_modificacion": audit.get("lastAnnotationUpdateDate", ""),
        "Revisado": "Swiss-Prot" if "reviewed" in data.get("entryType", "").lower() else "TrEMBL",
        "Nombre_del_gen": gene_name,
        "Sinónimos": ", ".join(gene_synonyms),
        "Organismo": data.get("organism", {}).get("scientificName", ""),
        "PDB_ids": ", ".join(sorted(set(pdb_ids))),
    }
    return pd.DataFrame([row])


def get_cofactor_pubchem_dataframe(uniprot_id: str) -> pd.DataFrame:
    data = _safe_get_json(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json")
    columns = ["Compuesto", "Pubchem_id", "Peso_molecular", "Inchi", "Inchikey", "Iupac_name"]
    if not data:
        return pd.DataFrame(columns=columns)

    cofactors: list[str] = []
    for comment in data.get("comments", []):
        if comment.get("commentType") != "COFACTOR":
            continue
        for cofactor in comment.get("cofactors", []):
            if cofactor.get("name"):
                cofactors.append(cofactor["name"])

    rows = []
    for cofactor_name in sorted(set(cofactors)):
        cid_data = _safe_get_json(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quote(cofactor_name)}/cids/JSON"
        )
        cids = (cid_data or {}).get("IdentifierList", {}).get("CID", [])
        if not cids:
            print(f"[WARN] No se encontró CID en PubChem para: {cofactor_name}")
            continue

        cid = cids[0]
        props = _safe_get_json(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            f"{cid}/property/ExactMass,InChI,InChIKey,IUPACName/JSON"
        )
        prop_row = ((props or {}).get("PropertyTable", {}).get("Properties") or [{}])[0]
        rows.append(
            {
                "Compuesto": cofactor_name,
                "Pubchem_id": cid,
                "Peso_molecular": prop_row.get("ExactMass"),
                "Inchi": prop_row.get("InChI"),
                "Inchikey": prop_row.get("InChIKey"),
                "Iupac_name": prop_row.get("IUPACName"),
            }
        )

    return pd.DataFrame(rows, columns=columns)


def get_heteromolecules_from_structure(cif_file: Path) -> list[str]:
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("target", str(cif_file))
    hetero = set()
    for residue in structure.get_residues():
        hetfield = residue.id[0]
        if hetfield.strip() and hetfield != "W":
            hetero.add(residue.resname)
    return sorted(hetero)


def heteromolecule_dataframe_from_mmcif_dict(cif_file: Path) -> pd.DataFrame:
    mmcif_data = MMCIF2Dict(str(cif_file))
    comp_ids = mmcif_data.get("_pdbx_entity_nonpoly.comp_id", [])
    names = mmcif_data.get("_pdbx_entity_nonpoly.name", [])

    if isinstance(comp_ids, str):
        comp_ids = [comp_ids]
    if isinstance(names, str):
        names = [names]

    max_len = min(len(comp_ids), len(names))
    rows = [
        {"ID_tres_letras": comp_ids[i], "Nombre_heteromolecula": names[i]}
        for i in range(max_len)
    ]
    return pd.DataFrame(rows)


def add_pubchem_smiles(df: pd.DataFrame) -> pd.DataFrame:
    smiles_values = []
    for _, row in df.iterrows():
        name = str(row.get("Nombre_heteromolecula", "")).strip()
        comp_id = str(row.get("ID_tres_letras", "")).strip()
        smiles = None

        for query in [name, comp_id]:
            if not query:
                continue
            data = _safe_get_json(
                "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
                f"{quote(query)}/property/CanonicalSMILES/JSON"
            )
            props = (data or {}).get("PropertyTable", {}).get("Properties", [])
            if props and props[0].get("CanonicalSMILES"):
                smiles = props[0]["CanonicalSMILES"]
                break
        smiles_values.append(smiles)

    result = df.copy()
    result["SMILES"] = smiles_values
    return result


def write_sdf_with_molecular_weight(df: pd.DataFrame, output_sdf: Path) -> int:
    writer = Chem.SDWriter(str(output_sdf))
    written = 0
    for _, row in df.iterrows():
        smiles = row.get("SMILES")
        if not smiles:
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        mol.SetProp("_Name", str(row.get("Nombre_heteromolecula", "")))
        mol.SetProp("Molecular_weight", str(Descriptors.MolWt(mol)))
        writer.write(mol)
        written += 1
    writer.close()

    supplier = Chem.SDMolSupplier(str(output_sdf), removeHs=False)
    valid = sum(1 for mol in supplier if mol is not None)
    if valid != written:
        raise ValueError("El archivo SDF contiene moléculas no válidas para RDKit.")
    return written


def main() -> None:
    base_dir = Path(__file__).resolve().parent
    data_dir = base_dir / "data"

    print("\n=== 1A) Descarga de CIF desde PDB ===")
    downloaded, failed = download_pdb_cif_files(PDB_IDS, data_dir)
    print(f"Descargados: {downloaded}")
    print(f"No descargados: {failed}")

    print("\n=== 1B) Obtención de datos UniProt para 1tup ===")
    uniprot_id = get_uniprot_id_from_pdb("1tup")
    if not uniprot_id:
        print("[ERROR] No se pudo determinar un ID de UniProt para 1tup.")
        return

    print(f"UniProt ID para 1tup: {uniprot_id}")
    uniprot_df = get_uniprot_dataframe(uniprot_id)
    print(uniprot_df.to_string(index=False))
    sequence = ""
    uniprot_json = _safe_get_json(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json")
    if uniprot_json:
        sequence = uniprot_json.get("sequence", {}).get("value", "")
        protein_name = (
            uniprot_json.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value", "")
        )
        print(f"\nNombre completo de proteína: {protein_name}")
        print(f"Secuencia (1 letra): {sequence}")

    print("\n=== 1C) Cofactores en UniProt y datos de PubChem ===")
    print("En UniProt, el cofactor se identifica en comments[].commentType == 'COFACTOR'.")
    cofactor_df = get_cofactor_pubchem_dataframe(uniprot_id)
    print(cofactor_df.to_string(index=False))

    print("\n=== 2A/2B/2C/2D) Análisis de 4ogq.cif con Biopython y RDKit ===")
    cif_file = data_dir / "4ogq.cif"
    if not cif_file.exists():
        print("[ERROR] No existe 4ogq.cif para continuar con las tareas 2A-2D.")
        return

    hetero_list = get_heteromolecules_from_structure(cif_file)
    print(f"Heteromoléculas (sin agua): {hetero_list}")

    hetero_df = heteromolecule_dataframe_from_mmcif_dict(cif_file)
    hetero_df = add_pubchem_smiles(hetero_df)
    print(hetero_df.to_string(index=False))

    sdf_path = base_dir / "heteromoleculas_4ogq.sdf"
    total = write_sdf_with_molecular_weight(hetero_df, sdf_path)
    print(f"Moléculas guardadas en SDF: {total} ({sdf_path})")


if __name__ == "__main__":
    main()
