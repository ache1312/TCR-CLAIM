# Reporte critico del benchmark de definicion de clonotipos TCR

Fecha de este reporte: 2026-06-01  
Carpeta revisada: `/mnt/data1/Andres/clonotype_definition_benchmark`  
Resultado consolidado principal usado: `results/cross_dataset_with_io_direct_source_20260504`

## Resumen ejecutivo

Este trabajo compara dos formas de representar clonotipos TCR en datos
scRNA-seq/scTCR-seq de cancer:

- `ct_strict`: clonotipo estricto basado en TCR alfa y beta dominantes, con genes
  V/J/C/D y CDR3 nucleotidico pareado.
- `ct_vgene`: representacion relajada basada solo en el par `TRAV-TRBV`.

La conclusion principal no es que `TRAV-TRBV` reemplace al clonotipo estricto.
La conclusion defendible es mas acotada:

`TRAV-TRBV` puede servir como una representacion de cribado para clones
expandidos, especialmente en CD8 tumorales, y para generar candidatos de sharing
entre tejidos. No es suficientemente especifica para medir diversidad total del
repertorio ni para llamar sharing clonal final sin confirmacion por CDR3 pareado.

En terminos practicos, esto sirve para reducir el espacio de busqueda en analisis
downstream biologico. Permite priorizar grupos TCR que probablemente contienen
clones relevantes, cruzarlos con fenotipo, tejido, expansion, respuesta
terapeutica o ubicacion espacial, y luego exigir confirmacion estricta antes de
hacer afirmaciones biologicas.

## Para que sirve biologicamente

La utilidad real de este benchmark esta en separar tres escenarios:

1. **Screening de clones tumorales expandidos**

   En CD8 tumorales, los top clones estrictos se conservan muy bien bajo la
   representacion `TRAV-TRBV`. Esto permite usar `TRAV-TRBV` como llave rapida
   para priorizar grupos clonales grandes cuando el objetivo inicial es encontrar
   candidatos, no cerrar una conclusion definitiva.

2. **Priorizacion de sharing entre tejidos**

   Un grupo `TRAV-TRBV` compartido entre tumor, ganglio, sangre o tejido normal
   puede marcar un candidato de migracion, expansion comun o convergencia
   biologica. Pero el benchmark muestra que muchos de esos grupos compartidos son
   aparentes: comparten genes V alfa/beta, pero no necesariamente el mismo CDR3
   pareado. Por eso `TRAV-TRBV` sirve para priorizar, no para confirmar.

3. **Anotacion funcional de clones candidatos**

   Una vez definido un clon estricto, o un grupo `TRAV-TRBV` candidato, se puede
   cruzar con:

   - estado funcional: agotamiento, citotoxicidad, proliferacion, fenotipo
     progenitor-like;
   - compartimento tisular: tumor primario, metastasis/ganglio, sangre, normal;
   - magnitud de expansion;
   - persistencia entre muestras o tratamientos;
   - antigenos candidatos si luego se integra con pMHC, prediccion o datos
     funcionales.

La regla de uso deberia ser:

> usar `TRAV-TRBV` para encontrar y ordenar candidatos; usar `ct_strict` para
> afirmar clonotipo, expansion clonal o sharing real.

## Estructura del proyecto revisado

El proyecto esta escrito principalmente en R y Seurat. Los componentes centrales
son:

- `R/01_tcr_normalization.R`: lectura, normalizacion y construccion de clonotipos.
- `R/02_seurat_integration.R`: carga de objetos Seurat o matrices 10x,
  integracion de metadatos y union con TCR.
- `R/03_metrics.R`: metricas principales de QC, concordancia, dominancia y
  sharing.
- `R/05_phenotype_clonality.R`: asociacion entre expansion clonal estricta y
  firmas transcripcionales.
- `R/06_mait_control.R`: control MAIT-like basado en `TRAV1-2` y genes TRAJ
  canonicos.
- `scripts/run_pipeline.R`: pipeline principal por manifest.
- `scripts/summarize_all_results.R`: sintesis cross-dataset.
- `tests/testthat/`: tests sinteticos de definiciones y metricas.

Los tests locales fueron ejecutados y pasan:

```bash
Rscript tests/testthat.R
```

Resultado: todos los tests pasaron.

## Como funciona el pipeline

### Entrada requerida

Cada muestra debe tener:

- matriz RNA o Seurat object;
- tabla de contigs TCR, preferentemente estilo 10x `filtered_contig_annotations.csv`;
- columnas equivalentes a barcode, chain, productive, V/J, CDR3 nucleotidico y
  soporte por UMI/read si existe.

El README define la politica del proyecto: datasets sin TRA/TRB productivo
pareado y CDR3 nucleotidico no deben degradarse silenciosamente a TRB-only o a
CDR3 aminoacidico.

### Normalizacion TCR

El archivo `R/01_tcr_normalization.R` hace lo siguiente:

1. estandariza nombres de columnas;
2. conserva solo cadenas `TRA` y `TRB`;
3. filtra contigs productivos, full-length y celulares si esas columnas estan
   disponibles;
4. para cada celula y cadena, elige una cadena dominante;
5. el criterio de dominancia es UMI descendente, luego reads descendente, luego
   desempate deterministico por genes y CDR3;
6. construye una tabla wide por celula.

Esto es importante: si una celula tiene multiples alfas o multiples betas
productivas, el pipeline no modela todas las combinaciones posibles. Escoge una
alfa y una beta dominantes. La tabla conserva flags `has_multi_tra` y
`has_multi_trb`, pero la definicion de clonotipo usa solo el par dominante.

### Definicion de clonotipos

El clonotipo estricto se construye solo si hay:

- `TRAV`;
- `TRBV`;
- CDR3 nucleotidico alfa;
- CDR3 nucleotidico beta.

La formula conceptual es:

```text
TRAV.TRAJ.TRAC;TRA_CDR3_NT_TRBV.TRBD.TRBJ.TRBC;TRB_CDR3_NT
```

La representacion relajada es:

```text
TRAV_TRBV
```

Por construccion, muchos clonotipos estrictos pueden caer dentro del mismo grupo
`TRAV-TRBV`. Esa relacion many-to-one es el punto biologico y tecnico que el
benchmark cuantifica.

### Celulas primarias del benchmark

El benchmark primario usa solo celulas:

- anotadas como `CD4` o `CD8`;
- con `ct_strict` no faltante;
- con `ct_vgene` no faltante.

Las metricas se calculan por:

```text
dataset_id, cancer_type, donor_id, sample_id, tissue_type, cell_class
```

Por lo tanto, las unidades de evidencia no son celulas individuales sino estratos
donor/sample/tejido/CD4-CD8.

## Datasets y elegibilidad

El resultado consolidado mas amplio incluye estos result sets:

| result_id | celulas RNA/QC | TCR pareado | CD4/CD8 pareado primario |
|---|---:|---:|---:|
| `gse121637` | 25,688 | 7,829 | 4,682 |
| `gse139555` | 200,626 | 60,847 | 55,933 |
| `gse148190` | 8,848 | 4,481 | 3,413 |
| `gse181061` | 47,390 | 29,637 | 23,611 |
| `gse185206` | 217,101 | 123,601 | 0 |
| `gse193371` | 68,934 | 38,502 | 38,502 |
| `gse200996` | 17,694,720 | 41,006 | 29,236 |
| `gse201425` | 84,190 | 24,250 | 0 |
| `gse242477` | 527,342 | 31,059 | 0 |
| `io_dataset` | 35,059 | 19,775 | 18,741 |

Tres datasets (`gse185206`, `gse201425`, `gse242477`) tienen TCR pareado pero
cero celulas CD4/CD8 primarias bajo la anotacion usada en este pipeline. No deben
contarse como evidencia primaria para las conclusiones CD4/CD8, aunque pueden
servir para QC o para trabajo futuro si se corrige la anotacion celular.

`gse200996` requiere cautela especial: la tabla QC reporta muchos millones de
celulas RNA/QC y solo decenas de miles con TCR pareado. Esto sugiere que el
universo RNA original no esta filtrado de la misma forma que el universo TCR.
Para downstream, cualquier conclusion debe basarse en las celulas TCR pareadas y
no en el total bruto.

## Resultados principales

### Concordancia de conteo de clones

La metrica principal compara:

```text
(numero de grupos TRAV-TRBV - numero de clonotipos estrictos) /
numero de clonotipos estrictos
```

Valores negativos indican compresion: `TRAV-TRBV` junta varios clonotipos
estrictos. Valores positivos indican que, bajo thresholds de expansion, el
agrupamiento relajado puede crear grupos que superan el umbral aunque los
clonotipos estrictos individuales no lo hagan.

#### Tumor CD8

| set analitico | unidades | diferencia relativa mediana | dentro de 10% | dentro de 20% | interpretacion |
|---|---:|---:|---:|---:|---|
| todos, `>=1` | 39 | -28.5% | 2.6% | 12.8% | diversidad comprimida |
| `>=2` | 39 | +5.3% | 61.5% | 87.2% | conteos similares |
| `>=5` | 39 | +9.7% | 48.7% | 69.2% | razonable para expandidos |
| `>=10` | 39 | +18.7% | 31.6% | 52.6% | util pero no exacto |
| top 10 estrictos | 39 | 0.0% | 100.0% | 100.0% | alta concordancia |

Interpretacion: en CD8 tumoral, `TRAV-TRBV` funciona bien para capturar los
clones dominantes y bastante bien para expansion moderada. No funciona para medir
riqueza total del repertorio.

#### Tumor CD4

| set analitico | unidades | diferencia relativa mediana | dentro de 10% | dentro de 20% | interpretacion |
|---|---:|---:|---:|---:|---|
| todos, `>=1` | 30 | -13.1% | 33.3% | 63.3% | diferencia significativa |
| `>=2` | 30 | +62.5% | 3.7% | 3.7% | no confiable |
| `>=5` | 30 | +135.4% | 27.3% | 31.8% | no confiable |
| `>=10` | 30 | +33.3% | 46.7% | 46.7% | inestable |
| top 10 estrictos | 30 | 0.0% | 96.7% | 96.7% | concordante en top clones |

Interpretacion: en CD4, `TRAV-TRBV` solo es defendible para top clones. No
deberia usarse como proxy general de expansion ni diversidad.

### Diversidad clonal

La comparacion de diversidad confirma que `TRAV-TRBV` comprime la riqueza:

| celula | tejido | ratio de riqueza relajada/estricta | diferencia relativa |
|---|---|---:|---:|
| CD4 | tumor | 0.676 | -32.4% |
| CD8 | tumor | 0.460 | -54.0% |
| CD4 | sangre | 0.687 | -31.3% |
| CD8 | sangre | 0.769 | -23.1% |
| CD4 | normal/adyacente | 0.740 | -26.0% |
| CD8 | normal/adyacente | 0.627 | -37.3% |

La compresion es especialmente fuerte en CD8 tumorales. Esto es esperable:
muchos CDR3 distintos pueden usar el mismo par de genes V alfa/beta.

Conclusion: `TRAV-TRBV` no debe usarse para estimar riqueza del repertorio,
diversidad Shannon efectiva ni diversidad Simpson como reemplazo del clonotipo
estricto.

### Top N clonotipos

Los analisis Top N son mas favorables:

| celula | top N | diferencia relativa mediana | dentro de 10% | dentro de 20% |
|---|---:|---:|---:|---:|
| CD8 | 10 | 0.0% | 100.0% | 100.0% |
| CD8 | 20 | -5.0% | 91.7% | 100.0% |
| CD8 | 50 | -7.0% | 70.8% | 95.8% |
| CD8 | 100 | -12.0% | 29.2% | 83.3% |
| CD4 | 10 | 0.0% | 100.0% | 100.0% |
| CD4 | 20 | -5.0% | 85.7% | 100.0% |
| CD4 | 50 | -6.0% | 80.9% | 85.7% |
| CD4 | 100 | -11.0% | 42.9% | 66.7% |

Interpretacion: mientras mas se expande el universo desde top 10 hacia top 100,
mas aparece la perdida de especificidad. El uso mas defendible es top clones, no
repertorio completo.

### Dominancia dentro de grupos TRAV-TRBV

La pregunta aqui es: cuando un grupo `TRAV-TRBV` contiene muchas celulas, esta
dominado por un solo clonotipo estricto o mezcla muchos clonotipos?

En tumor CD8:

| minimo de celulas por grupo TRAV-TRBV | contextos | top strict fraction mediana ponderada por celulas |
|---:|---:|---:|
| 1 | 13 | 83.9% |
| 5 | 13 | 83.2% |
| 10 | 13 | 86.0% |
| 50 | 12 | 90.2% |

En tumor CD4:

| minimo de celulas por grupo TRAV-TRBV | contextos | top strict fraction mediana ponderada por celulas |
|---:|---:|---:|
| 1 | 10 | 89.3% |
| 5 | 8 | 84.8% |
| 10 | 6 | 72.5% |
| 50 | 1 | 86.4% |

Interpretacion: grupos grandes `TRAV-TRBV`, especialmente CD8 tumorales,
frecuentemente estan dominados por un clonotipo estricto. Esto justifica usar
`TRAV-TRBV` para priorizar clones candidatos grandes. Pero no implica que cada
grupo sea monoclonal.

### Sharing entre tejidos

La metrica de sharing separa:

- clonotipos estrictos compartidos entre tejidos;
- grupos `TRAV-TRBV` compartidos;
- grupos `TRAV-TRBV` compartidos que tienen respaldo por al menos un clon
  estricto compartido;
- sharing aparente: grupo `TRAV-TRBV` compartido sin clon estricto compartido.

#### Tumor vs no tumor

| celula | set | strict shared | TRAV-TRBV shared | TRAV-TRBV respaldado | aparente only | precision-like |
|---|---|---:|---:|---:|---:|---:|
| CD8 | `>=1` | 1,580 | 2,526 | 1,398 | 1,128 | 55.3% |
| CD8 | `>=5` | 830 | 1,278 | 760 | 518 | 59.5% |
| CD8 | `>=10` | 471 | 682 | 440 | 242 | 64.5% |
| CD4 | `>=1` | 630 | 3,703 | 588 | 3,115 | 15.9% |
| CD4 | `>=5` | 156 | 1,393 | 151 | 1,242 | 10.8% |
| CD4 | `>=10` | 59 | 423 | 58 | 365 | 13.7% |

#### Tumor vs tumor

| celula | set | strict shared | TRAV-TRBV shared | TRAV-TRBV respaldado | aparente only | precision-like |
|---|---|---:|---:|---:|---:|---:|
| CD8 | `>=1` | 1,384 | 1,937 | 1,066 | 871 | 55.0% |
| CD8 | `>=5` | 792 | 1,248 | 662 | 586 | 53.0% |
| CD8 | `>=10` | 490 | 734 | 430 | 304 | 58.6% |
| CD4 | `>=1` | 8 | 22 | 8 | 14 | 36.4% |

Interpretacion: `TRAV-TRBV` tiene alta sensibilidad por construccion, pero
precision limitada. La aparicion de muchos grupos compartidos aparentes impide
usar `TRAV-TRBV` como llamada final de sharing clonal.

Una advertencia tecnica importante: la metrica de recall estricto por
`TRAV-TRBV` es 100% en muchos contextos porque `TRAV-TRBV` es una proyeccion del
clonotipo estricto. Si un clon estricto esta compartido, su par `TRAV-TRBV`
tambien lo estara. Por eso el recall alto no es evidencia independiente de
validez biologica; la precision-like es la metrica mas informativa.

## Analisis fenotipo-clonalidad

El analisis de fenotipo usa `ct_strict` como definicion primaria de expansion.
Las firmas se calculan como expresion media log-normalizada y luego se
z-normalizan dentro de contexto dataset/cancer/tejido/CD4-CD8.

Hallazgos tumorales CD8:

| estado | expansion `>=2`, delta z | expansion `>=5`, delta z | expansion `>=10`, delta z |
|---|---:|---:|---:|
| citotoxico | +0.58 | +0.65 | +0.68 |
| agotado | +0.34 | +0.41 | +0.44 |
| progenitor | -0.64 | -0.71 | -0.75 |
| progenitor-exhausted | -0.39 | -0.39 | -0.41 |
| proliferativo | +0.14 | +0.15 | +0.15 |
| Tc-lineage | +0.49 | +0.56 | +0.63 |

Interpretacion: los clones CD8 tumorales expandidos tienden a tener mayor
citotoxicidad y agotamiento, y menor firma progenitor/naive-like. Esto es
biologicamente plausible, pero debe considerarse asociacion, no causalidad ni
anotacion celular definitiva.

En CD4 tumorales, los resultados son menos limpios. Hay aumento de firma
citotoxica en clones expandidos, descenso de firma helper, y Treg variable. La
interpretacion CD4 requiere anotacion celular mas curada, especialmente para
Treg.

## Control MAIT-like

El control MAIT-like usa:

- `trav1_2`: dominante alfa `TRAV1-2`;
- `mait_like_alpha`: `TRAV1-2` con `TRAJ33`, `TRAJ12` o `TRAJ20`;
- `mait_like_alpha_cdr3`: lo anterior mas CDR3 aminoacidico alfa que empieza
  con `CAV`.

La utilidad de este control es mostrar un caso donde compartir biologia de V/J
alfa puede colapsar clonotipos estrictos distintos. En tumor CD8, el resumen
indica que los grupos MAIT-like alpha tienen menor dominancia por clon estricto
que los no `TRAV1-2`: aproximadamente 64.4% vs 84.1% de top strict fraction
ponderada por celulas, segun el resumen del reporte MAIT.

Esto respalda la cautela central del benchmark: compartir genes V no equivale a
compartir clonotipo.

## Limitaciones principales

### 1. `TRAV-TRBV` es una proyeccion many-to-one

Esta es la limitacion biologica fundamental. Muchos CDR3 alfa/beta distintos
pueden usar los mismos genes `TRAV` y `TRBV`. Por eso `TRAV-TRBV` comprime
riqueza y produce sharing aparente.

### 2. El clonotipo estricto depende de llamadas de genes y nomenclatura

`ct_strict` incluye genes V/J/C/D ademas de CDR3 nucleotidico. Esto puede ser
util para trazabilidad, pero hace que diferencias de nomenclatura, genes
faltantes o llamadas ambiguas separen clones que quizas biologicamente deberian
agruparse. Una libreria downstream deberia versionar esta politica.

### 3. La seleccion de cadena dominante elimina informacion

En celulas con multiples TRA o TRB productivos, el pipeline elige una cadena
dominante. Esto reduce ambiguedad operacional, pero puede ocultar:

- doble TCR real;
- multiplets;
- ruido de ensamblaje;
- clonotipos alternativos.

El flag multi-chain debe entrar en QC y, para analisis sensibles, puede ser
motivo de exclusion o analisis separado.

### 4. CD4/CD8 no siempre esta curado

Cuando no hay anotacion, el pipeline infiere CD4/CD8 desde expresion cruda de
`CD4`, `CD8A` y `CD8B`. Esto es fragil, especialmente en datos tumorales con
dropout, celulas activadas, dobles positivos o anotaciones incompletas.

### 5. Algunos datasets no aportan al analisis primario

Tres result sets del resumen amplio tienen cero celulas CD4/CD8 pareadas
primarias. No deben contribuir a conclusiones CD4/CD8 hasta corregir o agregar
metadatos de celula T.

### 6. Los resultados son asociaciones, no validacion antigenica

El benchmark no demuestra especificidad antigenica. Un clon expandido o
compartido puede ser relevante, pero se necesita validacion adicional: antigeno,
pMHC, TCR reconstruction, perturbacion, datos funcionales o evidencia clinica.

## Que se puede afirmar

Se puede afirmar:

- `TRAV-TRBV` conserva muy bien los top clones estrictos, especialmente en CD8
  tumorales.
- `TRAV-TRBV` comprime de forma importante la diversidad total del repertorio.
- `TRAV-TRBV` genera muchas llamadas aparentes de sharing, especialmente en CD4
  y en comparaciones tumor vs no tumor.
- Para clones grandes CD8 tumorales, muchos grupos `TRAV-TRBV` estan dominados
  por un clonotipo estricto.
- Los clones CD8 tumorales expandidos por `ct_strict` se asocian con firmas de
  citotoxicidad y agotamiento.

## Que no se debe afirmar

No se debe afirmar:

- que `TRAV-TRBV` define clonotipos equivalentes a CDR3 pareado;
- que `TRAV-TRBV` mide diversidad de repertorio de forma valida;
- que sharing por `TRAV-TRBV` es sharing clonal real;
- que los resultados demuestran especificidad antigenica;
- que las asociaciones fenotipicas son causales;
- que los datasets sin CD4/CD8 primario aportan evidencia al benchmark CD4/CD8.

## Propuesta de libreria Python downstream

La libreria no deberia venderse como una herramienta para reemplazar clonotipos
estrictos. Deberia ser una herramienta para:

1. normalizar TCR;
2. construir representaciones estrictas y relajadas;
3. cuantificar riesgo de colapso;
4. priorizar candidatos biologicos;
5. integrarse con objetos single-cell en Python.

### Nombre conceptual

`clonobench` o `tcrclonekit`.

### Objetivo

Permitir que un biologo computacional tome contigs TCR y metadatos celulares y
obtenga:

- clonotipos estrictos reproducibles;
- grupos relajados `TRAV-TRBV`;
- QC de elegibilidad;
- metricas de expansion;
- dominancia de CDR3 dentro de grupos V-gene;
- sharing estricto y sharing aparente;
- tablas listas para AnnData, scanpy, scirpy, modelos estadisticos o reportes.

### API minima propuesta

```python
import clonobench as cb

contigs = cb.read_10x_contigs("filtered_contig_annotations.csv")
tcr = cb.normalize_contigs(
    contigs,
    dataset_id="GSE139555",
    donor_id="patient_1",
    sample_id="tumor_1",
    tissue_type="tumor",
)

tcr = cb.define_clonotypes(tcr, strict_policy="paired_nt_vjc")
qc = cb.qc_summary(tcr)

agreement = cb.clone_count_agreement(
    tcr,
    groupby=["dataset_id", "donor_id", "sample_id", "tissue_type", "cell_class"],
    thresholds=[1, 2, 5, 10],
)

dominance = cb.vgene_dominance(tcr)
sharing = cb.tissue_sharing(tcr, tissue_col="tissue_type")
```

### Modulos sugeridos

```text
clonobench/
  io.py              # lectores 10x, AIRR, TSV/CSV wide
  normalize.py       # estandarizacion, cadenas dominantes, QC
  definitions.py     # ct_strict, ct_vgene, politicas alternativas
  metrics.py         # conteos, diversidad, dominancia, sharing
  phenotype.py       # joins con scores/fenotipos externos
  anndata.py         # integracion con AnnData/obs/obsm
  reporting.py       # tablas resumen y markdown/html
  tests/
```

### Principios de diseno

- Cada clonotipo debe tener un campo `definition_policy`.
- Nunca degradar silenciosamente de paired TRA/TRB a TRB-only.
- Guardar flags de ambiguedad: multi TRA, multi TRB, missing CDR3, missing V/J,
  barcode no mapeado.
- Separar claramente `candidate_group` de `confirmed_clone`.
- Hacer que todas las metricas puedan agruparse por donor/sample/tissue/cell
  state.
- Integrar con `AnnData.obs` sin requerir Seurat.

### Downstream biologico directo

Una libreria asi serviria para:

1. **Priorizar clones tumorales relevantes**

   Ordenar clones por expansion estricta, por dominancia dentro de `TRAV-TRBV`,
   y por enriquecimiento tumoral.

2. **Buscar candidatos de migracion o recirculacion**

   Detectar grupos relajados compartidos entre tejidos y luego confirmar por
   `ct_strict`.

3. **Relacionar clonotipo con fenotipo**

   Asociar expansion clonal con firmas de agotamiento, citotoxicidad,
   proliferacion, memoria/progenitor-like o Treg.

4. **Construir cohortes comparables**

   Estandarizar TCR entre datasets publicos y privados con QC explicito.

5. **Preparar candidatos para validacion experimental**

   Exportar clonotipos estrictos, CDR3 alfa/beta, genes V/J y tamaño clonal para
   reconstruccion TCR, pMHC screening o ensayos funcionales.

## Recomendaciones inmediatas

1. Mantener `ct_strict` como verdad operacional para clonotipo.
2. Usar `ct_vgene` solo como agrupacion relajada/candidata.
3. En reportes, cambiar lenguaje de "clonotipo `TRAV-TRBV`" a "grupo
   `TRAV-TRBV`" cuando no haya confirmacion por CDR3.
4. Auditar y corregir datasets con `n_primary_cd4_cd8_paired = 0`.
5. Revisar `GSE200996` para asegurar que el universo celular usado no mezcle raw
   droplets con celulas filtradas.
6. Agregar QC por porcentaje de barcodes TCR mapeados a RNA.
7. Agregar una metrica explicita por grupo `TRAV-TRBV`: numero de CDR3
   estrictos, clon dominante, fraccion dominante y entropia.
8. Para downstream, exportar una tabla canonica por celula y una tabla canonica
   por clon estricto.

## Veredicto

El trabajo es valioso, pero la conclusion debe formularse con precision.

La version fuerte, incorrecta, seria:

> `TRAV-TRBV` define clonotipos equivalentes a CDR3 pareado.

La version correcta es:

> `TRAV-TRBV` es una representacion relajada util para priorizar clones
> expandidos y candidatos de sharing, sobre todo en CD8 tumorales, pero comprime
> la diversidad y genera sharing aparente. Las conclusiones biologicas finales
> deben confirmarse con clonotipos estrictos basados en CDR3 nucleotidico
> pareado.

Esa distincion hace que el resultado sea util: no como sustituto de la biologia
del TCR, sino como capa de priorizacion y control de calidad para analisis
downstream.

## Archivos fuente principales

- `README.md`
- `R/01_tcr_normalization.R`
- `R/02_seurat_integration.R`
- `R/03_metrics.R`
- `R/05_phenotype_clonality.R`
- `R/06_mait_control.R`
- `results/cross_dataset_with_io_direct_source_20260504/summary.md`
- `results/cross_dataset_with_io_direct_source_20260504/cross_dataset_clone_agreement_summary.csv`
- `results/cross_dataset_with_io_direct_source_20260504/cross_dataset_tissue_sharing_summary.csv`
- `results/cross_dataset_with_io_direct_source_20260504/cross_dataset_vgene_dominance_by_size.csv`
- `results/share/clonal_diversity_comparison/clonal_diversity_summary.csv`
- `results/share/topn_clonotype_comparison/topn_strict_vs_relaxed_summary.csv`
- `results/phenotype_clonality_with_io_direct_source_20260504/summary.md`
- `results/mait_control_with_io_direct_source_20260504/summary.md`
