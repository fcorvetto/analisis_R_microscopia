# Análisis estadístico de datos de microscopía en R

Script de R para el análisis estadístico de datos de microscopía de fluorescencia, incluyendo modelos mixtos generalizados (GLMM), contrastes personalizados y visualizaciones con ggplot2.

## Descripción

Este script analiza el número y tamaño de condensados biomoleculares (BMCs/MLOs) por célula en función del genotipo y tratamiento, a partir de datos cuantificados en ImageJ.

## Pipeline de análisis

1. **Carga y preprocesamiento** de datos desde Excel
2. **Estadística descriptiva** por genotipo y tratamiento
3. **Modelado estadístico** con GLMM (distribución binomial negativa)
4. **Diagnóstico del modelo** con DHARMa (residuos simulados)
5. **Contrastes personalizados** con emmeans
6. **Visualizaciones** con ggplot2 (dispersión, violin plots)

## Librerías requeridas

```r
install.packages(c("readxl", "dplyr", "tidyr", "ggplot2",
                   "glmmTMB", "DHARMa", "car", "emmeans", "ggeffects"))
```

## Uso

1. Cloná el repositorio
2. Abrí el script en RStudio
3. Modificá la línea de `setwd()` con la ruta a tu carpeta de datos
4. Asegurate de tener el archivo `.xlsx` con los datos en esa carpeta
5. Ejecutá el script sección por sección

## Estructura del script

| Sección | Descripción |
|---|---|
| Importar datos | Carga el Excel y define factores |
| Descriptiva | Calcula media, SD, min y max por grupo |
| Modelo GLMM | Ajusta modelo con binomial negativa y efectos aleatorios |
| Diagnósticos | Verifica supuestos con residuos simulados |
| EMMEANS | Estima medias marginales del modelo |
| Contrastes | Comparaciones personalizadas con corrección de Holm |
| Gráficos | Plots de dispersión y violin plots por condición |

## Outputs generados

- `resumendatos_experimento_02_05_23.csv` — estadística descriptiva por grupo
- `tabla_contrastes_experimento_02_05_23.csv` — resultados de contrastes con p-valores ajustados
- Gráficos de dispersión (Número vs Tamaño de BMCs)
- Violin plots por tratamiento y genotipo

## Relacionado

- [Segmentación con Cellpose](https://github.com/fcorvetto/cell-pose) — genera las ROIs usadas para cuantificar en ImageJ
- [Macros de ImageJ](https://github.com/fcorvetto/imagej-macros) — cuantifica intensidad y área a partir de las ROIs
