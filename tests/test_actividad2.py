from pathlib import Path
import importlib.util

import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[1] / "actividad2_resuelta_ruben_juarez.py"
SPEC = importlib.util.spec_from_file_location("actividad2", MODULE_PATH)
actividad2 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(actividad2)


def test_validacion_pdb_ids():
    assert actividad2.es_pdb_id_clasico_valido("1TUP")
    assert actividad2.es_pdb_id_clasico_valido("9vwx")
    assert not actividad2.es_pdb_id_clasico_valido("10yza")
    assert not actividad2.es_pdb_id_clasico_valido("ABC")


def test_extraer_info_uniprot_cubre_rubrica_y_campos_omitidos():
    data = {
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "primaryAccession": "P04637",
        "entryAudit": {
            "firstPublicDate": "1987-08-13",
            "lastAnnotationUpdateDate": "2026-09-02",
        },
        "organism": {
            "scientificName": "Homo sapiens",
            "commonName": "Human",
        },
        "proteinDescription": {
            "recommendedName": {
                "fullName": {"value": "Cellular tumor antigen p53"}
            }
        },
        "genes": [
            {
                "geneName": {"value": "TP53"},
                "synonyms": [{"value": "P53"}],
            }
        ],
        "sequence": {"value": "MEEPQ"},
        "uniProtKBCrossReferences": [
            {"database": "PDB", "id": "1TUP"},
            {"database": "PDB", "id": "2OCJ"},
            {"database": "GO", "id": "GO:0003677"},
        ],
    }

    df = actividad2.extraer_info_uniprot(data)

    required = {
        "Uniprot_id",
        "Fecha_publicacion",
        "Fecha_modificacion",
        "Revisado",
        "Nombre_del_gen",
        "Sinónimos",
        "Organismo",
        "PDB_ids",
        "Nombre_proteina",
        "Secuencia",
    }
    assert required.issubset(df.columns)
    assert df.loc[0, "Uniprot_id"] == "P04637"
    assert df.loc[0, "Nombre_del_gen"] == "TP53"
    assert df.loc[0, "Sinónimos"] == "P53"
    assert df.loc[0, "Secuencia"] == "MEEPQ"
    assert "1TUP" in df.loc[0, "PDB_ids"]


def test_localizar_cofactor():
    data = {
        "comments": [
            {
                "commentType": "COFACTOR",
                "cofactors": [
                    {
                        "name": "Zn(2+)",
                        "cofactorCrossReference": {
                            "database": "ChEBI",
                            "id": "CHEBI:29105",
                        },
                    }
                ],
                "note": {
                    "texts": [
                        {"value": "Binds 1 zinc ion per subunit."}
                    ]
                },
            }
        ]
    }

    cofactors = actividad2.localizar_cofactores_uniprot(data)
    assert len(cofactors) == 1
    assert cofactors[0]["Nombre"] == "Zn(2+)"
    assert cofactors[0]["ID_referencia"] == "CHEBI:29105"
    assert "COFACTOR" in cofactors[0]["Ubicacion_JSON"]


def test_a_lista():
    assert actividad2.a_lista("X") == ["X"]
    assert actividad2.a_lista(["X", "Y"]) == ["X", "Y"]
    assert actividad2.a_lista(None) == []
