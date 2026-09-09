# MUBIO07 · Programación en Python · Actividad 2

[![CI](https://github.com/rubences/Miner-a-de-datos-bioinform-ticos/actions/workflows/ci.yml/badge.svg)](https://github.com/rubences/Miner-a-de-datos-bioinform-ticos/actions/workflows/ci.yml)

## Minería de datos bioinformáticos

Resolución reproducible de la actividad centrada en **minería de datos bioinformáticos**, integración de API científicas y manipulación de estructuras biomoleculares.

El objetivo del repositorio es doble: cumplir literalmente la rúbrica académica y mantener una implementación técnicamente defendible, reproducible y comprobable.

## Tecnologías y fuentes

- Python 3.11+
- `requests`
- `pandas`
- Biopython: `MMCIFParser` y `MMCIF2Dict`
- RDKit
- [RCSB Protein Data Bank](https://www.rcsb.org/)
- [UniProt REST API](https://rest.uniprot.org/)
- [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest)

## Cobertura de la actividad

El script principal es **`actividad2_resuelta_ruben_juarez.py`**. Implementa los dos bloques de la rúbrica:

### 1. Obtención de información sobre proteínas y compuestos químicos — 6,5 puntos

- Descarga tolerante a fallos de los mmCIF indicados en el enunciado.
- Control individual de errores: un PDB inexistente o un ID inválido no interrumpe el resto de la ejecución.
- Mapeo `1TUP` → UniProtKB mediante el servicio oficial de **UniProt ID Mapping**.
- Extracción de fechas, estado Swiss-Prot/TrEMBL, gen, sinónimos, organismo, nombre de proteína, secuencia y PDB asociados.
- Identificación del cofactor de p53 a partir del comentario estructurado `COFACTOR` de UniProt.
- Consulta de CID, masa exacta, InChI, InChIKey e IUPAC mediante PubChem.

### 2. Manipulación de datos biológicos — 3,5 puntos

- Parsing de `4OGQ.cif` mediante `MMCIFParser`.
- Exclusión explícita de moléculas de agua.
- Extracción de `_pdbx_entity_nonpoly.name` y `_pdbx_entity_nonpoly.comp_id` mediante `MMCIF2Dict`.
- Obtención de SMILES de las heteromoléculas mediante PubChem.
- Generación de un único SDF con la propiedad `Molecular_weight`.
- Relectura final mediante `Chem.SDMolSupplier` para verificar que cada registro escrito es un objeto `Mol` válido de RDKit.

La correspondencia exacta entre rúbrica, funciones y resultados se documenta en [`docs/MATRIZ_RUBRICA.md`](docs/MATRIZ_RUBRICA.md).

## Estructura

```text
.
├── .github/
│   └── workflows/
│       └── ci.yml
├── actividad2_resuelta_ruben_juarez.py
├── environment.yml
├── requirements.txt
├── .gitignore
├── README.md
├── data/
│   └── pdb/
├── resultados/
├── docs/
│   ├── MATRIZ_RUBRICA.md
│   └── RESULTADOS_ESPERADOS.md
└── tests/
    └── test_actividad2.py
```

Los `.cif`, CSV y SDF generados son reproducibles y se excluyen del control de versiones.

## Instalación recomendada con Conda

```bash
conda env create -f environment.yml
conda activate mubio07-act2
```

Alternativamente, con un entorno Python ya creado:

```bash
python -m pip install -r requirements.txt
```

## Ejecución

Desde la raíz del repositorio:

```bash
python actividad2_resuelta_ruben_juarez.py
```

Directorios personalizados:

```bash
python actividad2_resuelta_ruben_juarez.py \
  --data-dir data/pdb \
  --output-dir resultados
```

## Resultados generados

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

## Pruebas y CI

Las pruebas unitarias se han diseñado para validar la lógica determinista sin depender de que RCSB, UniProt o PubChem estén disponibles en ese instante:

```bash
python -m py_compile actividad2_resuelta_ruben_juarez.py
pytest -q
```

GitHub Actions ejecuta automáticamente esas comprobaciones en cada `push` a `main` y en cada pull request. Las consultas a API públicas se mantienen fuera de CI para evitar falsos negativos por disponibilidad, rate limiting o cambios transitorios de red.

## Decisiones técnicas relevantes

- Todas las llamadas remotas aplican `timeout`.
- Los errores HTTP transitorios usan reintentos con backoff; los 404 no se reintentan de forma inútil.
- El DataFrame de UniProt conserva **todas las columnas indicadas literalmente en el enunciado** y añade `Nombre_proteina` y `Secuencia`, porque el texto también solicita ambos datos.
- `Peso_molecular` del apartado del cofactor contiene `ExactMass`, ya que el enunciado pide el valor molecular exacto.
- La información de las heteromoléculas se extrae del propio mmCIF; no se codifica manualmente una lista de compuestos.
- El archivo SDF no se considera válido solo por haberse escrito: RDKit debe poder releer todos los registros generados.

## Controles biológicos de referencia

Las fuentes oficiales permiten comprobar que:

- [RCSB 1TUP](https://www.rcsb.org/structure/1TUP) vincula la proteína p53 humana con UniProtKB `P04637`.
- [UniProt P04637](https://www.uniprot.org/uniprotkb/P04637/entry) corresponde a `TP53`, está revisada (Swiss-Prot), contiene 393 aa y describe unión a Zn²⁺.
- [PubChem Zinc(2+)](https://pubchem.ncbi.nlm.nih.gov/compound/32051) corresponde al CID `32051`.
- [RCSB 4OGQ](https://www.rcsb.org/structure/4OGQ) informa actualmente de 17 ligandos únicos en la unidad asimétrica.

Los valores dinámicos —especialmente fechas de modificación y referencias cruzadas— se recuperan siempre de las API durante la ejecución y no se sustituyen por constantes.

## Entrega académica

La plataforma de evaluación solicita un **script Python `.py`**. El fichero que debe utilizarse para la entrega es:

```text
actividad2_resuelta_ruben_juarez.py
```

El resto del repositorio aporta reproducibilidad, validación y trazabilidad, pero no es necesario adjuntarlo si la plataforma limita la entrega a un único script.
