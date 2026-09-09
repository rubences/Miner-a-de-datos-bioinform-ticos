# MUBIO07 · Programación en Python · Actividad 2

## Minería de datos bioinformáticos

Resolución reproducible de la segunda actividad de **Programación en Python**, centrada en minería de datos bioinformáticos mediante servicios públicos y bibliotecas científicas.

### Tecnologías

- Python 3.11+
- `requests`
- `pandas`
- Biopython (`MMCIFParser`, `MMCIF2Dict`)
- RDKit
- RCSB Protein Data Bank
- UniProt REST API / ID Mapping API
- PubChem PUG REST

## Qué resuelve

El script `actividad2_resuelta_ruben_juarez.py` implementa los dos bloques de la rúbrica:

1. **Obtención de información sobre proteínas y compuestos químicos**
   - Descarga tolerante a fallos de los mmCIF indicados.
   - PDB `1TUP` → UniProt mediante el ID Mapping oficial.
   - Extracción de fechas, revisión Swiss-Prot/TrEMBL, gen, sinónimos, organismo, proteína, secuencia y PDB asociados.
   - Identificación estructurada del cofactor de p53.
   - Consulta de CID, masa exacta, InChI, InChIKey e IUPAC en PubChem.

2. **Manipulación de datos biológicos**
   - Parsing de `4OGQ.cif` mediante `MMCIFParser`.
   - Extracción de `_pdbx_entity_nonpoly` mediante `MMCIF2Dict`.
   - Obtención de SMILES mediante PubChem.
   - Generación de un SDF multipropiedad con `Molecular_weight`.
   - Relectura final del SDF con RDKit para comprobar que todos los registros escritos son objetos `Mol` válidos.

## Estructura

```text
.
├── actividad2_resuelta_ruben_juarez.py
├── environment.yml
├── requirements.txt
├── .gitignore
├── README.md
├── data/
│   └── pdb/
├── resultados/
├── docs/
│   └── RESULTADOS_ESPERADOS.md
└── tests/
    └── test_actividad2.py
```

Los archivos `.cif` y los resultados generados se crean durante la ejecución y no es necesario versionarlos.

## Instalación recomendada con Conda

```bash
conda env create -f environment.yml
conda activate mubio07-act2
```

## Ejecución

Desde la raíz del repositorio:

```bash
python actividad2_resuelta_ruben_juarez.py
```

También pueden personalizarse los directorios:

```bash
python actividad2_resuelta_ruben_juarez.py \
  --data-dir data/pdb \
  --output-dir resultados
```

## Resultados generados

La ejecución crea:

```text
resultados/
├── 01_descargas_pdb.csv
├── 02_uniprot_1TUP.csv
├── 03_cofactor_pubchem.csv
├── 04_heteromoleculas_4OGQ.csv
├── 05_heteromoleculas_smiles.csv
├── 06_validacion_sdf.csv
└── heteromoleculas_4OGQ.sdf
```

## Pruebas

Las pruebas incluidas validan la lógica pura de extracción de datos sin depender de los servicios remotos:

```bash
pytest -q
```

## Diseño y decisiones técnicas

- Las excepciones de una descarga PDB se aíslan por identificador.
- Se usan reintentos únicamente para errores HTTP transitorios.
- Se aplica `timeout` a todas las llamadas remotas.
- El ID Mapping de UniProt se ejecuta mediante el flujo asíncrono oficial.
- El DataFrame de UniProt conserva **todas las columnas expresamente indicadas por la rúbrica** y añade `Nombre_proteina` y `Secuencia`, porque el texto también exige ambas aunque no aparezcan en la lista de columnas del enunciado.
- La consulta a PubChem parte primero del **nombre químico extraído del mmCIF**, tal como exige la actividad, y usa el `PDB_comp_id` solo como respaldo.
- El peso del apartado del cofactor se interpreta como **masa exacta (`ExactMass`)**, no como peso molecular medio.
- El SDF final se vuelve a abrir con `Chem.SDMolSupplier`; escribir un fichero no se considera suficiente si RDKit no puede releerlo.

## Nota sobre GitHub

Para una actividad evaluable es recomendable mantener el repositorio **privado hasta que la calificación sea definitiva**, para evitar que terceros reutilicen el código antes del cierre de la entrega.
