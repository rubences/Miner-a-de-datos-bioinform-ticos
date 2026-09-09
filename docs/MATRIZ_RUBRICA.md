# Matriz de cobertura de la rúbrica

La actividad se evalúa sobre 10 puntos: 6,5 puntos para obtención de información mediante API y 3,5 puntos para manipulación bioinformática/química.

| Apartado | Puntos | Requisito | Implementación principal | Evidencia generada |
|---|---:|---|---|---|
| 1A | 2,5 | Descargar los `.cif` de la lista PDB y controlar fallos sin detener el programa | `descargar_cif_pdb()` | `resultados/01_descargas_pdb.csv` |
| 1B | 2,5 | Resolver `1TUP` en UniProt y extraer metadatos, gen, organismo, secuencia y PDB asociados | `mapear_pdb_a_uniprot()`, `obtener_entrada_uniprot()`, `extraer_info_uniprot()` | `resultados/02_uniprot_1TUP.csv` |
| 1C | 1,5 | Identificar el cofactor en UniProt y recuperar CID, masa exacta, InChI, InChIKey e IUPAC desde PubChem | `localizar_cofactores_uniprot()`, `consultar_cofactor_pubchem()` | `resultados/03_cofactor_pubchem.csv` |
| 2A | 0,5 | Parsear `4OGQ.cif` con `MMCIFParser()` y listar heteromoléculas sin agua | `heteromoleculas_con_mmcifparser()` | salida de consola + comprobación de IDs |
| 2B | 0,5 | Usar `MMCIF2Dict()` y `_pdbx_entity_nonpoly` para ID/nombre | `dataframe_heteromoleculas_con_mmcif2dict()` | `resultados/04_heteromoleculas_4OGQ.csv` |
| 2C | 1,0 | Consultar PubChem usando el nombre y añadir SMILES al DataFrame | `resolver_compuesto_pubchem()`, `agregar_smiles_pubchem()` | `resultados/05_heteromoleculas_smiles.csv` |
| 2D | 1,5 | Generar un SDF, añadir `Molecular_weight` con RDKit y verificar relectura como `Mol` | `generar_sdf_heteromoleculas()` | `heteromoleculas_4OGQ.sdf` + `06_validacion_sdf.csv` |
| **Total** | **10,0** |  |  |  |

## Controles adicionales que exceden el mínimo

- Sesión HTTP común con reintentos exponenciales para 429/5xx.
- `timeout` en todas las peticiones externas.
- Validación ligera del contenido descargado antes de aceptar un `.cif`.
- Priorización de una entrada UniProt revisada si un mapeo devuelve más de una candidata.
- Conservación de `Nombre_proteina` y `Secuencia` aunque el listado de columnas de la rúbrica las omita, porque el texto sí exige ambas.
- Estrategia PubChem por nombre químico con `PDB_comp_id` como fallback.
- Relectura independiente del SDF final mediante `Chem.SDMolSupplier`.
- Tests unitarios sin dependencia de red y CI automatizada.
