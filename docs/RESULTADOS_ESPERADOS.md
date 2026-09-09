# Resultados de referencia y comprobaciones

Este documento contiene **controles de coherencia**, no datos que el programa deba tener codificados. Los valores dinámicos se consultan siempre en las API oficiales.

## 1TUP → p53 / UniProt

La ficha oficial de [RCSB 1TUP](https://www.rcsb.org/structure/1TUP) describe el complejo del supresor tumoral p53 con ADN y asocia la proteína humana con UniProtKB `P04637`.

Valores de referencia:

- **UniProt accession:** `P04637`
- **Entry name:** `P53_HUMAN`
- **Estado:** reviewed / Swiss-Prot
- **Gen:** `TP53`
- **Sinónimo:** `P53`
- **Organismo:** `Homo sapiens (Human)`
- **Nombre recomendado:** `Cellular tumor antigen p53`
- **Longitud de la secuencia canónica:** `393 aa`

Fuente: [UniProt P04637](https://www.uniprot.org/uniprotkb/P04637/entry).

Las fechas de publicación/modificación y el conjunto de referencias PDB pueden cambiar con nuevas versiones de UniProt; por ello no se fijan como constantes.

## Cofactor

UniProt identifica el cofactor de forma estructurada en:

```text
comments -> commentType = "COFACTOR"
```

Control esperado:

- Cofactor: `Zn(2+)`
- Referencia ChEBI: `CHEBI:29105`
- Nota: `Binds 1 zinc ion per subunit.`
- PubChem CID: `32051`
- Exact Mass: aproximadamente `63.929142 Da`
- InChI: `InChI=1S/Zn/q+2`
- InChIKey: `PTFCDOFLOPIGGS-UHFFFAOYSA-N`
- IUPAC: `zinc(2+)`
- SMILES de referencia: `[Zn+2]`

Fuente química: [PubChem CID 32051](https://pubchem.ncbi.nlm.nih.gov/compound/32051).

## 4OGQ

[RCSB 4OGQ](https://www.rcsb.org/structure/4OGQ) corresponde a *Internal Lipid Architecture of the Hetero-Oligomeric Cytochrome b6f Complex*. La ficha actual indica **17 ligandos únicos**.

En la **unidad asimétrica**, las instancias mostradas por RCSB suman **31 heteromoléculas no acuosas**. La asamblea biológica dimérica duplica esas posiciones; por eso no debe confundirse el número de instancias de la unidad asimétrica con el de la asamblea. `MMCIFParser` aplicado al archivo depositado trabaja sobre las coordenadas de la unidad asimétrica.

Los 17 `comp_id` únicos de referencia son:

```text
CLA, OPC, SQD, 1O2, 2WM, 3WM, HEC, 2WA, 7PH,
BCR, UMQ, 2WD, 8K6, MYS, FES, OCT, CD
```

El programa **no depende de esta lista**: obtiene tanto el identificador como el nombre desde `_pdbx_entity_nonpoly` del propio mmCIF. Esta enumeración existe únicamente como control independiente de coherencia.
