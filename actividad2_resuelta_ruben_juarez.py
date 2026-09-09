# -*- coding: utf-8 -*-
"""
MUBIO07 - Programación en Python
Actividad 2: Minería de datos bioinformáticos

Autor: Rubén Juárez
Descripción:
    Resolución reproducible de la actividad mediante las API REST de RCSB PDB,
    UniProt y PubChem, junto con Biopython, Pandas y RDKit.

La ejecución completa:
    1) Descarga los mmCIF indicados por el enunciado, controlando errores.
    2) Mapea PDB 1TUP -> UniProt mediante la API de ID Mapping.
    3) Extrae los metadatos requeridos de UniProt y los guarda en DataFrame.
    4) Localiza el cofactor de p53 y consulta sus propiedades en PubChem.
    5) Parsea 4OGQ.cif con MMCIFParser y MMCIF2Dict.
    6) Recupera SMILES de las heteromoléculas mediante PubChem.
    7) Genera un SDF multipropiedad y verifica que RDKit puede releerlo.

Uso:
    python actividad2_resuelta_ruben_juarez.py

Opcional:
    python actividad2_resuelta_ruben_juarez.py --data-dir data/pdb --output-dir resultados
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import time
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import quote

import pandas as pd
import requests
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from rdkit import Chem
from rdkit.Chem import Descriptors
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


# ---------------------------------------------------------------------------
# Configuración general
# ---------------------------------------------------------------------------

RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
UNIPROT_API = "https://rest.uniprot.org"
PUBCHEM_API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

PDB_IDS = [
    "1tup", "2xyz", "3def", "4ogq", "5jkl",
    "6mno", "7pqr", "8stu", "9vwx", "10yza",
]

POLL_INTERVAL_SECONDS = 1.0
PUBCHEM_PAUSE_SECONDS = 0.22
HTTP_TIMEOUT_SECONDS = 30

LOGGER = logging.getLogger("mubio07_act2")


# ---------------------------------------------------------------------------
# Utilidades HTTP
# ---------------------------------------------------------------------------

def crear_sesion_http() -> requests.Session:
    """
    Crea una sesión HTTP con reintentos exponenciales para errores transitorios.

    No se reintentan los 404, ya que en esta actividad un identificador inexistente
    debe registrarse como fallo y no bloquear el resto de la ejecución.
    """
    retry = Retry(
        total=4,
        connect=4,
        read=4,
        backoff_factor=0.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET", "POST"}),
        raise_on_status=False,
    )
    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": (
                "MUBIO07-Actividad2/1.0 "
                "(academic bioinformatics exercise; Python requests)"
            )
        }
    )
    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


def comprobar_respuesta(response: requests.Response, contexto: str) -> None:
    """Lanza una excepción HTTP con un mensaje contextual más informativo."""
    try:
        response.raise_for_status()
    except requests.HTTPError as exc:
        raise requests.HTTPError(
            f"{contexto}: HTTP {response.status_code} - {response.url}"
        ) from exc


def a_lista(value: Any) -> list[Any]:
    """Normaliza valores escalares/listas de MMCIF2Dict a una lista."""
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


# ===========================================================================
# ACTIVIDAD 1A - Descarga de estructuras PDB en formato mmCIF
# ===========================================================================

def es_pdb_id_clasico_valido(pdb_id: str) -> bool:
    """
    Comprueba el formato clásico de un identificador PDB: 4 caracteres y
    comienzo numérico. Es suficiente para la lista proporcionada en el enunciado.
    """
    return bool(re.fullmatch(r"[0-9][A-Za-z0-9]{3}", pdb_id.strip()))


def descargar_cif_pdb(
    pdb_ids: Iterable[str],
    destino: Path,
    session: requests.Session,
) -> pd.DataFrame:
    """
    Descarga una lista de archivos mmCIF desde RCSB PDB.

    Devuelve un DataFrame con el estado de cada descarga. Un error individual
    nunca detiene las restantes descargas.
    """
    destino.mkdir(parents=True, exist_ok=True)
    estados: list[dict[str, Any]] = []

    for raw_id in pdb_ids:
        pdb_id = raw_id.strip().upper()
        registro: dict[str, Any] = {
            "PDB_id": pdb_id,
            "Descargado": False,
            "Archivo": "",
            "Estado": "",
        }

        if not es_pdb_id_clasico_valido(pdb_id):
            registro["Estado"] = (
                "ID no válido en formato PDB clásico (se esperaban 4 caracteres)"
            )
            estados.append(registro)
            LOGGER.warning("%s -> %s", pdb_id, registro["Estado"])
            continue

        url = RCSB_CIF_URL.format(pdb_id=pdb_id)
        archivo = destino / f"{pdb_id}.cif"
        temporal = archivo.with_suffix(".cif.part")

        try:
            response = session.get(
                url,
                stream=True,
                timeout=HTTP_TIMEOUT_SECONDS,
            )
            comprobar_respuesta(response, f"Descarga de {pdb_id}")

            with temporal.open("wb") as fh:
                for chunk in response.iter_content(chunk_size=64 * 1024):
                    if chunk:
                        fh.write(chunk)

            # Comprobación ligera de integridad del mmCIF.
            if temporal.stat().st_size < 20:
                raise ValueError("El archivo descargado está vacío o es demasiado pequeño")

            with temporal.open("rb") as fh:
                cabecera = fh.read(256).lower()

            if b"data_" not in cabecera:
                raise ValueError("El contenido recibido no parece un archivo mmCIF")

            temporal.replace(archivo)
            registro.update(
                {
                    "Descargado": True,
                    "Archivo": str(archivo),
                    "Estado": f"OK ({archivo.stat().st_size:,} bytes)",
                }
            )
            LOGGER.info("%s -> descargado correctamente", pdb_id)

        except (requests.RequestException, OSError, ValueError) as exc:
            temporal.unlink(missing_ok=True)
            registro["Estado"] = f"ERROR: {exc}"
            LOGGER.warning("%s -> %s", pdb_id, registro["Estado"])

        estados.append(registro)

    return pd.DataFrame(estados)


# ===========================================================================
# ACTIVIDAD 1B - PDB 1TUP -> UniProt y extracción de metadatos
# ===========================================================================

def mapear_pdb_a_uniprot(
    pdb_id: str,
    session: requests.Session,
    max_intentos: int = 30,
) -> str:
    """
    Mapea un PDB ID a UniProtKB usando el servicio oficial UniProt ID Mapping.

    Se utiliza:
        POST /idmapping/run
        GET  /idmapping/status/{jobId}
        GET  /idmapping/details/{jobId}
        GET  redirectURL de resultados
    """
    pdb_id = pdb_id.upper()

    submit = session.post(
        f"{UNIPROT_API}/idmapping/run",
        data={"from": "PDB", "to": "UniProtKB", "ids": pdb_id},
        timeout=HTTP_TIMEOUT_SECONDS,
    )
    comprobar_respuesta(submit, "Envío del trabajo de ID Mapping de UniProt")

    job_id = submit.json().get("jobId")
    if not job_id:
        raise RuntimeError("UniProt no devolvió un jobId para el ID Mapping")

    terminado = False

    for _ in range(max_intentos):
        status_response = session.get(
            f"{UNIPROT_API}/idmapping/status/{job_id}",
            timeout=HTTP_TIMEOUT_SECONDS,
            allow_redirects=False,
        )

        # UniProt puede responder 303 cuando los resultados ya están disponibles.
        if status_response.status_code == 303:
            terminado = True
            break

        comprobar_respuesta(status_response, "Consulta de estado del ID Mapping")
        status_data = status_response.json()
        estado = status_data.get("jobStatus")

        if estado in {"RUNNING", "NEW"}:
            time.sleep(POLL_INTERVAL_SECONDS)
            continue

        if estado == "FINISHED" or "results" in status_data or "failedIds" in status_data:
            terminado = True
            break

        if estado in {"ERROR", "FAILED"}:
            raise RuntimeError(f"El ID Mapping de UniProt terminó con estado {estado}")

        time.sleep(POLL_INTERVAL_SECONDS)

    if not terminado:
        raise TimeoutError("El ID Mapping de UniProt excedió el tiempo de espera")

    details = session.get(
        f"{UNIPROT_API}/idmapping/details/{job_id}",
        timeout=HTTP_TIMEOUT_SECONDS,
    )
    comprobar_respuesta(details, "Obtención de detalles del ID Mapping")

    redirect_url = details.json().get("redirectURL")
    if not redirect_url:
        # Fallback al endpoint genérico documentado.
        redirect_url = f"{UNIPROT_API}/idmapping/results/{job_id}"

    results_response = session.get(
        redirect_url,
        params={"format": "json", "size": 500},
        timeout=HTTP_TIMEOUT_SECONDS,
    )
    comprobar_respuesta(results_response, "Obtención de resultados del ID Mapping")

    results_data = results_response.json()
    results = results_data.get("results", [])
    if not results:
        failed = results_data.get("failedIds", [])
        raise LookupError(
            f"No se encontró correspondencia UniProt para {pdb_id}. "
            f"failedIds={failed}"
        )

    # Si hubiese varias correspondencias, se prioriza una entrada revisada.
    accessions: list[tuple[str, bool]] = []

    for item in results:
        target = item.get("to")
        if isinstance(target, dict):
            accession = (
                target.get("primaryAccession")
                or target.get("id")
                or target.get("uniProtkbId")
            )
            entry_type = str(target.get("entryType", "")).lower()
            reviewed = "reviewed" in entry_type and "unreviewed" not in entry_type
        else:
            accession = str(target) if target else None
            reviewed = False

        if accession:
            accessions.append((accession, reviewed))

    if not accessions:
        raise LookupError(f"UniProt devolvió resultados sin accession para {pdb_id}")

    accessions.sort(key=lambda x: x[1], reverse=True)
    return accessions[0][0]


def obtener_entrada_uniprot(
    accession: str,
    session: requests.Session,
) -> dict[str, Any]:
    """Recupera una entrada UniProtKB completa en JSON."""
    response = session.get(
        f"{UNIPROT_API}/uniprotkb/{accession}.json",
        timeout=HTTP_TIMEOUT_SECONDS,
    )
    comprobar_respuesta(response, f"Consulta de UniProt {accession}")
    return response.json()


def nombre_proteina_uniprot(data: dict[str, Any]) -> str:
    """Obtiene el nombre recomendado de proteína con fallbacks."""
    description = data.get("proteinDescription", {})

    recommended = description.get("recommendedName", {})
    full_name = recommended.get("fullName", {})
    if isinstance(full_name, dict) and full_name.get("value"):
        return str(full_name["value"])

    for key in ("submissionNames", "alternativeNames"):
        values = description.get(key, [])
        if values:
            candidate = values[0].get("fullName", {})
            if isinstance(candidate, dict) and candidate.get("value"):
                return str(candidate["value"])

    return ""


def extraer_info_uniprot(data: dict[str, Any]) -> pd.DataFrame:
    """
    Construye el DataFrame solicitado por el enunciado.

    Nota:
    El texto del ejercicio también solicita el nombre completo de la proteína y
    su secuencia, aunque esas dos variables no aparecen en la lista de columnas
    que se imprime en el enunciado. Para no perder información se incluyen como
    columnas adicionales, manteniendo todas las columnas obligatorias.
    """
    audit = data.get("entryAudit", {})
    entry_type = str(data.get("entryType", ""))

    reviewed = (
        "Sí (Swiss-Prot)"
        if ("reviewed" in entry_type.lower() and "unreviewed" not in entry_type.lower())
        else "No (TrEMBL)"
    )

    genes = data.get("genes", [])
    gene_name = ""
    synonyms: list[str] = []

    if genes:
        gene_name = genes[0].get("geneName", {}).get("value", "")
        synonyms = [
            synonym.get("value", "")
            for synonym in genes[0].get("synonyms", [])
            if synonym.get("value")
        ]

    organism_data = data.get("organism", {})
    organism = organism_data.get("scientificName", "")
    common_name = organism_data.get("commonName")
    if common_name:
        organism = f"{organism} ({common_name})"

    pdb_ids = sorted(
        {
            ref.get("id", "")
            for ref in data.get("uniProtKBCrossReferences", [])
            if ref.get("database") == "PDB" and ref.get("id")
        }
    )

    sequence = data.get("sequence", {}).get("value", "")

    row = {
        "Uniprot_id": data.get("primaryAccession", ""),
        "Fecha_publicacion": audit.get("firstPublicDate", ""),
        "Fecha_modificacion": audit.get("lastAnnotationUpdateDate", ""),
        "Revisado": reviewed,
        "Nombre_del_gen": gene_name,
        "Sinónimos": ", ".join(synonyms),
        "Organismo": organism,
        "PDB_ids": ", ".join(pdb_ids),
        # Dos columnas adicionales requeridas por el texto descriptivo.
        "Nombre_proteina": nombre_proteina_uniprot(data),
        "Secuencia": sequence,
    }

    df = pd.DataFrame([row])

    columnas_rubrica = [
        "Uniprot_id",
        "Fecha_publicacion",
        "Fecha_modificacion",
        "Revisado",
        "Nombre_del_gen",
        "Sinónimos",
        "Organismo",
        "PDB_ids",
    ]
    assert set(columnas_rubrica).issubset(df.columns)

    return df


# ===========================================================================
# ACTIVIDAD 1C - Cofactor UniProt y consulta PubChem
# ===========================================================================

def localizar_cofactores_uniprot(
    data: dict[str, Any],
) -> list[dict[str, str]]:
    """
    Busca los comentarios estructurados de tipo COFACTOR de UniProt.

    También conserva la referencia ChEBI y la nota explicativa para poder
    justificar explícitamente dónde se identifica el cofactor.
    """
    encontrados: list[dict[str, str]] = []

    for comment in data.get("comments", []):
        if comment.get("commentType") != "COFACTOR":
            continue

        note = "; ".join(
            text.get("value", "")
            for text in comment.get("note", {}).get("texts", [])
            if text.get("value")
        )

        for cofactor in comment.get("cofactors", []):
            xref = cofactor.get("cofactorCrossReference", {})
            encontrados.append(
                {
                    "Nombre": cofactor.get("name", ""),
                    "Base_referencia": xref.get("database", ""),
                    "ID_referencia": xref.get("id", ""),
                    "Nota": note,
                    "Ubicacion_JSON": "comments -> commentType='COFACTOR'",
                }
            )

    return encontrados


def buscar_cids_pubchem(
    nombre: str,
    session: requests.Session,
) -> list[int]:
    """Obtiene CIDs de PubChem por nombre/sinónimo."""
    encoded = quote(nombre, safe="")
    response = session.get(
        f"{PUBCHEM_API}/compound/name/{encoded}/cids/JSON",
        timeout=HTTP_TIMEOUT_SECONDS,
    )

    if response.status_code == 404:
        return []

    comprobar_respuesta(response, f"Búsqueda PubChem para '{nombre}'")
    return [
        int(cid)
        for cid in response.json().get("IdentifierList", {}).get("CID", [])
    ]


def propiedades_pubchem(
    cid: int,
    properties: str,
    session: requests.Session,
) -> dict[str, Any]:
    """Recupera propiedades químicas de un CID mediante PUG REST."""
    response = session.get(
        f"{PUBCHEM_API}/compound/cid/{cid}/property/{properties}/JSON",
        timeout=HTTP_TIMEOUT_SECONDS,
    )
    comprobar_respuesta(response, f"Propiedades PubChem CID {cid}")

    rows = response.json().get("PropertyTable", {}).get("Properties", [])
    if not rows:
        raise LookupError(f"PubChem no devolvió propiedades para CID {cid}")
    return rows[0]


def consultar_cofactor_pubchem(
    cofactor_name: str,
    session: requests.Session,
) -> pd.DataFrame:
    """
    Consulta en PubChem el cofactor obtenido previamente desde UniProt.

    'Peso_molecular' contiene ExactMass, porque el enunciado pide de forma
    explícita el peso/masa molecular exacta.
    """
    cids = buscar_cids_pubchem(cofactor_name, session)

    # Para Zn(2+) PubChem reconoce también "zinc ion"; se añade como fallback
    # semántico sin fijar el CID manualmente.
    if not cids and "zn" in cofactor_name.lower():
        cids = buscar_cids_pubchem("zinc ion", session)

    if not cids:
        raise LookupError(
            f"No se encontró el cofactor '{cofactor_name}' en PubChem"
        )

    cid = cids[0]
    props = propiedades_pubchem(
        cid,
        "ExactMass,InChI,InChIKey,IUPACName",
        session,
    )

    return pd.DataFrame(
        [
            {
                "Compuesto": cofactor_name,
                "Pubchem_id": cid,
                "Peso_molecular": props.get("ExactMass"),
                "Inchi": props.get("InChI", ""),
                "Inchikey": props.get("InChIKey", ""),
                "Iupac_name": props.get("IUPACName", ""),
            }
        ],
        columns=[
            "Compuesto",
            "Pubchem_id",
            "Peso_molecular",
            "Inchi",
            "Inchikey",
            "Iupac_name",
        ],
    )


# ===========================================================================
# ACTIVIDAD 2A - MMCIFParser: lista de heteromoléculas de 4OGQ
# ===========================================================================

def heteromoleculas_con_mmcifparser(cif_path: Path) -> list[Any]:
    """
    Parsea 4OGQ con MMCIFParser y devuelve todos los residuos hetero no acuosos.

    En Biopython:
        residue.id[0] == " " -> residuo estándar
        residue.id[0] == "W" -> agua
        otros valores         -> heteromolécula
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("4OGQ", str(cif_path))

    heteromoleculas = [
        residue
        for residue in structure.get_residues()
        if residue.id[0] not in {" ", "W"}
    ]

    return heteromoleculas


# ===========================================================================
# ACTIVIDAD 2B - MMCIF2Dict: DataFrame de heteromoléculas
# ===========================================================================

def dataframe_heteromoleculas_con_mmcif2dict(cif_path: Path) -> pd.DataFrame:
    """
    Extrae `_pdbx_entity_nonpoly.name` y `_pdbx_entity_nonpoly.comp_id`.

    Se excluye HOH porque el apartado A especifica que no deben incluirse aguas,
    y se elimina cualquier posible duplicado por comp_id.
    """
    mmcif = MMCIF2Dict(str(cif_path))

    names = a_lista(mmcif.get("_pdbx_entity_nonpoly.name"))
    comp_ids = a_lista(mmcif.get("_pdbx_entity_nonpoly.comp_id"))

    if not names or not comp_ids:
        raise KeyError(
            "No se encontraron las claves _pdbx_entity_nonpoly.name/comp_id"
        )

    if len(names) != len(comp_ids):
        raise ValueError(
            "Las columnas _pdbx_entity_nonpoly de nombre e ID tienen "
            "longitudes distintas"
        )

    df = pd.DataFrame(
        {
            "Compuesto": names,
            "PDB_comp_id": comp_ids,
        }
    )

    df = (
        df[df["PDB_comp_id"].astype(str).str.upper() != "HOH"]
        .drop_duplicates(subset=["PDB_comp_id"])
        .reset_index(drop=True)
    )

    return df


# ===========================================================================
# ACTIVIDAD 2C - PubChem: SMILES de cada heteromolécula
# ===========================================================================

def resolver_compuesto_pubchem(
    nombre: str,
    comp_id: str,
    session: requests.Session,
) -> tuple[int | None, str | None, str]:
    """
    Resuelve una heteromolécula contra PubChem.

    Estrategia:
        1. Buscar por el nombre químico extraído del propio mmCIF, tal como pide
           el enunciado.
        2. Si no hay resultado, usar el PDB comp_id como sinónimo de respaldo.
        3. Recuperar IsomericSMILES/CanonicalSMILES y conservar el primero
           disponible, dando prioridad a la estereoquímica.
    """
    intentos = [nombre]
    if comp_id and comp_id.upper() != nombre.upper():
        intentos.append(comp_id)

    for termino in intentos:
        try:
            cids = buscar_cids_pubchem(termino, session)
        except requests.RequestException:
            cids = []

        if not cids:
            continue

        cid = cids[0]
        try:
            props = propiedades_pubchem(
                cid,
                "IsomericSMILES,CanonicalSMILES",
                session,
            )
        except requests.RequestException:
            continue

        smiles = (
            props.get("IsomericSMILES")
            or props.get("CanonicalSMILES")
            or props.get("SMILES")
            or props.get("ConnectivitySMILES")
        )

        if smiles:
            return cid, str(smiles), termino

    return None, None, ""


def agregar_smiles_pubchem(
    df: pd.DataFrame,
    session: requests.Session,
) -> pd.DataFrame:
    """Añade CID, SMILES y término de resolución PubChem al DataFrame."""
    rows: list[dict[str, Any]] = []

    for row in df.to_dict(orient="records"):
        cid, smiles, termino = resolver_compuesto_pubchem(
            str(row["Compuesto"]),
            str(row["PDB_comp_id"]),
            session,
        )

        enriched = dict(row)
        enriched["PubChem_CID"] = cid
        enriched["SMILES"] = smiles
        enriched["Consulta_PubChem"] = termino
        rows.append(enriched)

        if cid is None:
            LOGGER.warning(
                "PubChem: no se pudo resolver %s (%s)",
                row["PDB_comp_id"],
                row["Compuesto"],
            )
        else:
            LOGGER.info(
                "PubChem: %s -> CID %s",
                row["PDB_comp_id"],
                cid,
            )

        # Pausa prudente para no enviar ráfagas innecesarias al servicio público.
        time.sleep(PUBCHEM_PAUSE_SECONDS)

    return pd.DataFrame(rows)


# ===========================================================================
# ACTIVIDAD 2D - SDF + RDKit + Molecular_weight
# ===========================================================================

def obtener_sdf_pubchem(
    cid: int,
    session: requests.Session,
) -> str:
    """Descarga un registro SDF 2D de PubChem para un CID."""
    response = session.get(
        f"{PUBCHEM_API}/compound/cid/{cid}/SDF",
        params={"record_type": "2d"},
        timeout=HTTP_TIMEOUT_SECONDS,
    )
    comprobar_respuesta(response, f"SDF de PubChem CID {cid}")
    return response.text


def mol_desde_sdf_o_smiles(
    sdf_text: str,
    smiles: str | None,
) -> Chem.Mol | None:
    """
    Intenta crear un Mol de RDKit desde SDF.

    Si una estructura poco convencional no pudiera sanitizarse directamente,
    se realiza un segundo intento y, como último respaldo, se usa el SMILES
    recuperado en el apartado C.
    """
    mol = Chem.MolFromMolBlock(sdf_text, sanitize=True, removeHs=False)
    if mol is not None:
        return mol

    mol = Chem.MolFromMolBlock(sdf_text, sanitize=False, removeHs=False)
    if mol is not None:
        try:
            Chem.SanitizeMol(mol)
            return mol
        except Exception:
            pass

    if smiles:
        return Chem.MolFromSmiles(smiles)

    return None


def generar_sdf_heteromoleculas(
    df: pd.DataFrame,
    output_path: Path,
    session: requests.Session,
) -> pd.DataFrame:
    """
    Genera un único SDF con una molécula por heteromolécula única.

    Añade las propiedades:
        - PDB_comp_id
        - PubChem_CID
        - Molecular_weight

    Tras escribir el fichero, lo vuelve a abrir con SDMolSupplier para verificar
    que cada registro puede convertirse en un objeto Mol de RDKit.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    writer = Chem.SDWriter(str(output_path))
    if writer is None:
        raise OSError(f"No se pudo abrir {output_path} para escritura")

    registros: list[dict[str, Any]] = []
    escritos = 0

    try:
        for row in df.to_dict(orient="records"):
            comp_id = str(row["PDB_comp_id"])
            nombre = str(row["Compuesto"])
            cid = row.get("PubChem_CID")
            smiles = row.get("SMILES")

            status = {
                "PDB_comp_id": comp_id,
                "Compuesto": nombre,
                "PubChem_CID": cid,
                "Escrito_SDF": False,
                "Molecular_weight": None,
                "Error": "",
            }

            if pd.isna(cid):
                status["Error"] = "Sin CID PubChem; no es posible solicitar SDF"
                registros.append(status)
                continue

            try:
                cid_int = int(cid)
                sdf_text = obtener_sdf_pubchem(cid_int, session)
                mol = mol_desde_sdf_o_smiles(
                    sdf_text,
                    None if pd.isna(smiles) else str(smiles),
                )

                if mol is None:
                    raise ValueError("RDKit no pudo construir un objeto Mol")

                molecular_weight = float(Descriptors.MolWt(mol))

                mol.SetProp("_Name", f"{comp_id} | {nombre}")
                mol.SetProp("PDB_comp_id", comp_id)
                mol.SetProp("PubChem_CID", str(cid_int))
                mol.SetProp("Molecular_weight", f"{molecular_weight:.4f}")

                writer.write(mol)
                escritos += 1

                status["Escrito_SDF"] = True
                status["Molecular_weight"] = molecular_weight

            except (requests.RequestException, ValueError, OSError) as exc:
                status["Error"] = str(exc)

            registros.append(status)
            time.sleep(PUBCHEM_PAUSE_SECONDS)

    finally:
        writer.close()

    # Validación independiente: todos los registros realmente escritos deben
    # poder releerse como objetos Mol.
    supplier = Chem.SDMolSupplier(str(output_path), removeHs=False)
    reread = list(supplier)

    if len(reread) != escritos:
        raise RuntimeError(
            f"Validación SDF: se escribieron {escritos} registros, "
            f"pero RDKit releyó {len(reread)}"
        )

    invalidos = [i for i, mol in enumerate(reread) if mol is None]
    if invalidos:
        raise RuntimeError(
            "Validación SDF: RDKit devolvió Mol=None en registros "
            f"{invalidos}"
        )

    LOGGER.info(
        "SDF validado: %d/%d registros escritos y releídos correctamente",
        escritos,
        len(df),
    )

    return pd.DataFrame(registros)


# ===========================================================================
# Presentación y ejecución completa
# ===========================================================================

def imprimir_titulo(texto: str) -> None:
    print("\n" + "=" * 79)
    print(texto)
    print("=" * 79)


def guardar_dataframe(
    df: pd.DataFrame,
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8-sig")


def ejecutar(data_dir: Path, output_dir: Path) -> None:
    """Ejecuta de forma secuencial todos los apartados de la actividad."""
    session = crear_sesion_http()
    data_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ 1A
    imprimir_titulo("1A. Descarga de archivos .cif desde RCSB PDB")
    df_downloads = descargar_cif_pdb(PDB_IDS, data_dir, session)
    print(df_downloads.to_string(index=False))
    guardar_dataframe(df_downloads, output_dir / "01_descargas_pdb.csv")

    # ------------------------------------------------------------------ 1B
    imprimir_titulo("1B. Identificación de 1TUP en UniProt y extracción de datos")
    uniprot_id = mapear_pdb_a_uniprot("1TUP", session)
    print(f"1TUP -> UniProtKB: {uniprot_id}")

    uniprot_json = obtener_entrada_uniprot(uniprot_id, session)
    df_uniprot = extraer_info_uniprot(uniprot_json)
    print(df_uniprot.drop(columns=["Secuencia"]).to_string(index=False))
    print(f"\nLongitud de la secuencia: {len(df_uniprot.loc[0, 'Secuencia'])} aa")
    guardar_dataframe(df_uniprot, output_dir / "02_uniprot_1TUP.csv")

    # ------------------------------------------------------------------ 1C
    imprimir_titulo("1C. Cofactor de 1TUP/p53 y propiedades PubChem")
    cofactores = localizar_cofactores_uniprot(uniprot_json)
    if not cofactores:
        raise LookupError("No se encontraron comentarios COFACTOR en la entrada UniProt")

    for c in cofactores:
        print(
            f"Cofactor: {c['Nombre']} | "
            f"{c['Base_referencia']} {c['ID_referencia']} | "
            f"{c['Nota']} | Ubicación: {c['Ubicacion_JSON']}"
        )

    # La entrada de p53 contiene un único tipo de cofactor: Zn(2+).
    df_cofactor = consultar_cofactor_pubchem(cofactores[0]["Nombre"], session)
    print(df_cofactor.to_string(index=False))
    guardar_dataframe(df_cofactor, output_dir / "03_cofactor_pubchem.csv")

    # ------------------------------------------------------------------ 2A
    cif_4ogq = data_dir / "4OGQ.cif"
    if not cif_4ogq.exists():
        raise FileNotFoundError(
            "No existe 4OGQ.cif. Revise la descarga del apartado 1A "
            "o coloque el archivo en el directorio de datos."
        )

    imprimir_titulo("2A. Heteromoléculas de 4OGQ con MMCIFParser")
    hetero_residues = heteromoleculas_con_mmcifparser(cif_4ogq)
    print(f"Número de instancias heteromoleculares sin agua: {len(hetero_residues)}")
    print(
        "IDs únicos detectados:",
        sorted({res.get_resname() for res in hetero_residues}),
    )

    # ------------------------------------------------------------------ 2B
    imprimir_titulo("2B. _pdbx_entity_nonpoly con MMCIF2Dict")
    df_hetero = dataframe_heteromoleculas_con_mmcif2dict(cif_4ogq)
    print(df_hetero.to_string(index=False))
    guardar_dataframe(df_hetero, output_dir / "04_heteromoleculas_4OGQ.csv")

    # ------------------------------------------------------------------ 2C
    imprimir_titulo("2C. SMILES de heteromoléculas mediante PubChem")
    df_smiles = agregar_smiles_pubchem(df_hetero, session)
    print(df_smiles.to_string(index=False))
    guardar_dataframe(df_smiles, output_dir / "05_heteromoleculas_smiles.csv")

    # ------------------------------------------------------------------ 2D
    imprimir_titulo("2D. Generación y validación del archivo SDF")
    sdf_path = output_dir / "heteromoleculas_4OGQ.sdf"
    df_sdf = generar_sdf_heteromoleculas(df_smiles, sdf_path, session)
    print(df_sdf.to_string(index=False))
    guardar_dataframe(df_sdf, output_dir / "06_validacion_sdf.csv")

    total = len(df_sdf)
    ok = int(df_sdf["Escrito_SDF"].sum()) if total else 0
    print(f"\nSDF final: {sdf_path}")
    print(f"Registros válidos RDKit: {ok}/{total}")

    imprimir_titulo("ACTIVIDAD COMPLETADA")


def construir_parser_argumentos() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="MUBIO07 - Actividad 2 de minería de datos bioinformáticos"
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("data/pdb"),
        help="Directorio de archivos PDB/mmCIF (por defecto: data/pdb)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("resultados"),
        help="Directorio de resultados (por defecto: resultados)",
    )
    return parser


def main() -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s | %(message)s",
    )

    args = construir_parser_argumentos().parse_args()

    try:
        ejecutar(args.data_dir, args.output_dir)
    except KeyboardInterrupt:
        LOGGER.error("Ejecución cancelada por el usuario")
        return 130
    except Exception as exc:
        LOGGER.exception("La actividad terminó con error: %s", exc)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
