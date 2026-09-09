# Resultados de referencia y comprobaciones

Este documento no sustituye a la ejecución de las API; sirve para comprobar que
el programa está obteniendo información biológicamente coherente.

## 1TUP

La estructura **1TUP** corresponde al complejo del supresor tumoral p53 con ADN.
La proteína de p53 se asocia a la entrada UniProtKB:

- **UniProt accession:** `P04637`
- **Entry name:** `P53_HUMAN`
- **Estado:** reviewed / Swiss-Prot
- **Gen:** `TP53`
- **Sinónimo:** `P53`
- **Organismo:** `Homo sapiens (Human)`
- **Nombre recomendado:** `Cellular tumor antigen p53`
- **Longitud de la secuencia canónica:** 393 aa

La API de UniProt es la fuente autoritativa durante la ejecución, por lo que
fechas de modificación y listas de referencias pueden cambiar con nuevas
versiones de la base de datos.

## Cofactor

La entrada de UniProt identifica el cofactor en el comentario estructurado:

```text
comments -> commentType = "COFACTOR"
```

Información esperada:

- Cofactor: `Zn(2+)`
- Referencia UniProt/ChEBI: `CHEBI:29105`
- Nota: `Binds 1 zinc ion per subunit.`
- PubChem CID de zinc(2+): `32051`
- Exact Mass: aproximadamente `63.929142 Da`
- InChI: `InChI=1S/Zn/q+2`
- InChIKey: `PTFCDOFLOPIGGS-UHFFFAOYSA-N`
- IUPAC: `zinc(2+)`

## 4OGQ

4OGQ es el complejo fotosintético hetero-oligomérico cytochrome b6f.

El mmCIF contiene 62 instancias químicas/no poliméricas no acuosas,
correspondientes a **17 tipos únicos** de heteromolécula:

| PDB comp_id | Nombre |
|---|---|
| HEC | HEME C |
| UMQ | UNDECYL-MALTOSIDE |
| 7PH | (1R)-2-(dodecanoyloxy)-1-[(phosphonooxy)methyl]ethyl tetradecanoate |
| 8K6 | OCTADECANE |
| 2WM | (1S,8E)-1-{[(2S)-3-hydroxy-2-{[(1S)-1-hydroxyoctadecyl]oxy}propyl]oxy}octadec-8-en-1-ol |
| CLA | CHLOROPHYLL A |
| OPC | dioleoyl-phosphatidylcholine-like ligand (PDB OPC) |
| CD | CADMIUM ION |
| MYS | PENTADECANE |
| SQD | SULFOQUINOVOSYLDIACYLGLYCEROL |
| FES | FE2/S2 (INORGANIC) CLUSTER |
| 2WD | ether lipid ligand 2WD |
| 3WM | ether lipid ligand 3WM |
| 2WA | ether lipid ligand 2WA |
| OCT | N-OCTANE |
| 1O2 | galactolipid ligand 1O2 |
| BCR | BETA-CAROTENE |

El DataFrame del programa extrae los nombres completos del propio mmCIF, por lo
que no depende de esta tabla abreviada para resolver la actividad.
