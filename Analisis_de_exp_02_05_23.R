rm(list = ls())

# =========================
# DIRECTORIO
# =========================
setwd("ruta/a/tu/carpeta/de/datos")  # Modificar según tu sistema

# =========================
# LIBRERÍAS
# =========================
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(car)
library(emmeans)
library(ggeffects)

# =========================
# IMPORTAR DATOS
# =========================
Datos <- read_excel("02.05.23 mSmaug WT y SM Ba, Chx, Met y SAG para r.xlsx")

Datos <- Datos %>%
  mutate(
    Tratamiento = factor(Tratamiento, levels = c("Ba","Chx","Met","SAG")),
    Genotipo = factor(Genotipo, levels = c("Spmd","WT")),
    Cover = factor(Cover),
    Cover.Global = interaction(Cover, Tratamiento, Genotipo),
    Trat.Gen = interaction(Tratamiento, Genotipo),
    NumeroTamaño = Numero_focos * Tamaño_promedio_MLOs
  )

# =========================
# EXPLORACIÓN
# =========================
str(Datos)
summary(Datos)

# =========================
# DESCRIPTIVA
# =========================
variables <- c("Numero_focos","Area_celular","Int_celula_entera",
               "Intensidad_en_focos","Tamaño_promedio_MLOs","NumeroTamaño")

resumen <- Datos %>%
  group_by(Genotipo, Tratamiento) %>%
  summarise(
    n_celulas = n(),   # 👈 número de observaciones por grupo
    
    across(
      all_of(variables),
      list(
        media = mean,
        sd    = sd,
        min   = min,
        max   = max
      ),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(resumen)

write.csv(
  resumen,
  "resumendatos_experimento_02_05_23.csv",
  row.names = FALSE
)


# =========================
# MODELO GLM con offset)
# =========================
modelo2 <- glmmTMB(
  Numero_focos ~ Genotipo + Tratamiento + offset(log(Area_celular)),
  data = Datos,
  family = nbinom2()
)


modelo1 <- glmmTMB(Numero_focos ~ Genotipo * Tratamiento + offset(log(Area_celular)) + (1 | Cover.Global), 
                   data = Datos,
                   family = nbinom2())
#Aplique el modelo 1 porque es mejor ya que cuenta la interaccion
# =========================
# DIAGNÓSTICOS
# =========================
simulationOutput1 <- simulateResiduals(modelo1)
plot(simulationOutput1)
testDispersion(simulationOutput1)

# =========================
# RESUMEN MODELO
# =========================
summary(modelo1)
confint(modelo1)

drop1(modelo1, test = "Chi")
Anova(modelo1, type = 2)

# =========================
# EMMEANS
# =========================
emm1 <- emmeans(modelo1, ~ Tratamiento * Genotipo)

summary(emm1)
levels(emm1)

# Orden esperado:
# 1 Ba   Spmd
# 2 Chx  Spmd
# 3 Met  Spmd
# 4 SAG  Spmd
# 5 Ba   WT
# 6 Chx  WT
# 7 Met  WT
# 8 SAG  WT

# =========================
# CONTRASTES CUSTOM
# =========================
custom_contrastes <- list(
  "Spmd Ba vs WT Ba"    = c( 1,  0,  0,  0, -1,  0,  0,  0),
  "WT Ba vs WT Chx"     = c( 0,  0,  0,  0,  1, -1,  0,  0),
  "WT Ba vs WT Met"     = c( 0,  0,  0,  0,  1,  0, -1,  0),
  "WT Ba vs WT SAG"     = c( 0,  0,  0,  0,  1,  0,  0, -1),
  "Spmd Ba vs Spmd Chx" = c( 1, -1,  0,  0,  0,  0,  0,  0),
  "Spmd Ba vs Spmd Met" = c( 1,  0, -1,  0,  0,  0,  0,  0),
  "Spmd Ba vs Spmd SAG" = c( 1,  0,  0, -1,  0,  0,  0,  0)
)

# =========================
# APLICAR CONTRASTES
# =========================
contrastes_raw <- contrast(emm1, custom_contrastes)
contrastes_df <- as.data.frame(contrastes_raw)

# =========================
# TABLA FINAL
# =========================
tabla_contrastes <- contrastes_df %>%
  mutate(
    log_Ratio   = estimate,
    Ratio       = exp(estimate),
    p.value_adj = p.adjust(p.value, method = "holm")
  ) %>%
  select(
    contrast,
    Ratio,
    log_Ratio,
    SE,
    z.ratio,
    p.value,
    p.value_adj
  )

tabla_contrastes

# =========================
# PLOT EMMEANS
# =========================
plot(emm1, comparisons = TRUE)

# =========================
# PREDICCIONES
# =========================
pred1 <- ggpredict(modelo1, terms = c("Tratamiento", "Genotipo"))

plot(pred1) +
  labs(
    title = "Modelo exp 02.05.23",
    y = "Número estimado de focos por célula",
    x = "Tratamiento"
  ) +
  theme_minimal()

# =========================
# EXPORTAR TABLA
# =========================
write.csv(tabla_contrastes,
          "tabla_contrastes_experimento_02_05_23.csv",
          row.names = FALSE)



# --------- GRAFICOS De Numero vs Tamaño-------------------------------------------

# ----------------------------------------------------
# 1. ESTILOS GLOBALES
# ----------------------------------------------------

COLORES_TRATAMIENTO <- c(
  "Ba"  = "#f8776c",
  "Chx" = "#7bad04",
  "Met" = "#0fbec3",
  "SAG" = "#c77cff"
)

FORMAS_GENOTIPO <- c(
  "Spmd" = 17, # triángulo
  "WT"   = 16  # círculo
)

# ----------------------------------------------------
# 2. PROMEDIOS Y SE POR GRUPO
# ----------------------------------------------------

Datos_prom <- Datos %>%
  group_by(Genotipo, Tratamiento) %>%
  summarise(
    Count_prom   = mean(Numero_focos, na.rm = TRUE),
    AvgSize_prom = mean(Tamaño_promedio_MLOs, na.rm = TRUE),
    se_Count     = sd(Numero_focos, na.rm = TRUE) / sqrt(n()),
    se_Size      = sd(Tamaño_promedio_MLOs, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# ----------------------------------------------------
# 3. FUNCIÓN DE GRAFICACIÓN
# ----------------------------------------------------

graficar_dispersion <- function(
    datos_ind, datos_prom, titulo,
    x_lim = c(0, 100), y_lim = c(0, 2),
    leyenda = c("ambas", "solo_genotipo", "solo_tratamiento", "ninguna")
) {
  
  leyenda <- match.arg(leyenda)
  
  p <- ggplot(datos_ind, aes(x = Numero_focos, y = Tamaño_promedio_MLOs)) +
    
    # puntos individuales
    geom_point(
      aes(shape = Genotipo, color = Tratamiento),
      size = 0.8, alpha = 0.4
    ) +
    
    # barras de error horizontales
    geom_errorbarh(
      data = datos_prom,
      aes(
        x = Count_prom, y = AvgSize_prom,
        xmin = Count_prom - se_Count,
        xmax = Count_prom + se_Count,
        color = Tratamiento
      ),
      height = 0,
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    
    # barras de error verticales
    geom_errorbar(
      data = datos_prom,
      aes(
        x = Count_prom, y = AvgSize_prom,
        ymin = AvgSize_prom - se_Size,
        ymax = AvgSize_prom + se_Size,
        color = Tratamiento
      ),
      width = 0,
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    
    # puntos promedio
    geom_point(
      data = datos_prom,
      aes(
        x = Count_prom, y = AvgSize_prom,
        shape = Genotipo, color = Tratamiento
      ),
      size = 2,
      stroke = 1.2
    ) +
    
    scale_color_manual(values = COLORES_TRATAMIENTO, name = "Tratamiento") +
    scale_shape_manual(values = FORMAS_GENOTIPO, name = "Genotipo") +
    
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    
    labs(
      title = titulo,
      x = "Número de BMCs por célula",
      y = expression(paste("Tamaño promedio de BMCs (", mu, m^2, ")"))
    ) +
    
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      
      # quitar grillas
      #panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # reforzar ejes
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black")
    )
  
  if (leyenda == "solo_genotipo") {
    p <- p + guides(color = "none")
  } else if (leyenda == "solo_tratamiento") {
    p <- p + guides(shape = "none")
  } else if (leyenda == "ninguna") {
    p <- p + guides(shape = "none", color = "none")
  }
  
  return(p)
}

# ----------------------------------------------------
# 4. GENERACIÓN DE GRÁFICOS
# ----------------------------------------------------

# 4.1 Basales (comparación genotipo)
datos_Ba <- Datos %>% filter(Tratamiento == "Ba")
prom_Ba  <- Datos_prom %>% filter(Tratamiento == "Ba")

p_basales <- graficar_dispersion(
  datos_Ba, prom_Ba,
  "mSmaug1 WT vs mSmaug1 Spmd (Basal)",
  x_lim = c(0, 60),
  y_lim = c(0, 2),
  leyenda = "solo_genotipo"
)

print(p_basales)

# 4.2 WT (tratamientos)
datos_WT <- Datos %>% filter(Genotipo == "WT")
prom_WT  <- Datos_prom %>% filter(Genotipo == "WT")

p_WT <- graficar_dispersion(
  datos_WT, prom_WT,
  "mSmaug1 WT",
  x_lim = c(0, 60),
  y_lim = c(0, 2),
  leyenda = "solo_tratamiento"
)

print(p_WT)

# 4.3 Spmd (tratamientos)
datos_Spmd <- Datos %>% filter(Genotipo == "Spmd")
prom_Spmd  <- Datos_prom %>% filter(Genotipo == "Spmd")

p_Spmd <- graficar_dispersion(
  datos_Spmd, prom_Spmd,
  "mSmaug1 Spmd",
  x_lim = c(0, 60),
  y_lim = c(0, 2),
  leyenda = "solo_tratamiento"
)

print(p_Spmd)


#----Violin Plots-------

COLORES_TRATAMIENTO <- c(
  "Ba"  = "#f8776c",
  "Chx" = "#7bad04",
  "Met" = "#0fbec3",
  "SAG" = "#c77cff"
)

p_violin_media <- ggplot(Datos, aes(x = Tratamiento, y = Numero_focos, fill = Tratamiento)) +
  
  # violin
  geom_violin(
    alpha = 0.7,
    trim = FALSE,
    color = NA
  ) +
  
  # puntos individuales (fondo)
  geom_jitter(
    aes(color = Tratamiento),
    width = 0.15,
    size = 1.3,
    alpha = 0.25
  ) +
  
  # promedio por condición
  stat_summary(
    fun = mean,
    geom = "crossbar",
    size = 0.3,
    color = "black",
    width = 0.3
  ) +
  
  facet_wrap(~ Genotipo, nrow = 1, strip.position = "bottom") +
  
  scale_fill_manual(values = COLORES_TRATAMIENTO) +
  scale_color_manual(values = COLORES_TRATAMIENTO) +
  
  coord_cartesian(ylim = c(0, 60)) +
  
  labs(
    title = "Número de MLOs por condición",
    x = "Tratamiento",
    y = "Número de MLOs por célula"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    
    legend.position = "none",
    

    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

print(p_violin_media)
