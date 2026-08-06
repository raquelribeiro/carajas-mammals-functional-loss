########################################################################
############ Functional Diversity & Ecosystem Services ###########
# Rversion = 4.3.2 
# Gabriela Alves Ferreira (gabriela-alves77@hotmail.com)
# Adjusted for GitHub Repository Structure
#==============================================================================

# ------------------------------------------------------------------------------
# 0. Load Packages
# ------------------------------------------------------------------------------
library(FD)
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(janitor)
library(ggplot2)
library(stringr)
library(viridis)

# ------------------------------------------------------------------------------
# 1. Import Data
# ------------------------------------------------------------------------------
# Functional traits 
traits <- read_excel("data/traits.xlsx") %>%
  clean_names() %>%
  mutate(species = tolower(trimws(species)))

# Functional area data
dados <- read_excel("data/functional.xlsx") %>%
  clean_names() %>%
  mutate(species = tolower(trimws(species)))

# ------------------------------------------------------------------------------
# 2. Calculate base_value
# ------------------------------------------------------------------------------
traits_proc <- traits %>%
  mutate(
    grupo_score = case_when(
      functional_group %in% c("FRU", "HER") ~ 5,
      functional_group %in% c("ONI", "INS") ~ 3,
      functional_group == "CAR" ~ 2,
      TRUE ~ 1
    ),
    iucn_score = case_when(
      iucn_status == "CR" ~ 5,
      iucn_status == "EN" ~ 4,
      iucn_status == "VU" ~ 3,
      iucn_status == "NT" ~ 2,
      TRUE ~ 1
    ),
    estrato_score = case_when(
      grepl("^ar", foraging_stratum) ~ 5,
      grepl("^te", foraging_stratum) ~ 3,
      grepl("^su", foraging_stratum) ~ 1,
      TRUE ~ 2
    ),
    atividade_score = case_when(
      activity_cycle == 3 ~ 3,
      activity_cycle == 2 ~ 2,
      TRUE ~ 1
    )
  ) %>%
  mutate(
    base_value = rowMeans(select(., grupo_score, iucn_score, body_mass,
                                 dispersal_km, estrato_score, atividade_score), na.rm = TRUE)
  )

# ------------------------------------------------------------------------------
# 3. Calculate functional values & loss
# ------------------------------------------------------------------------------
base_value_df <- traits_proc %>%
  select(species, base_value) %>%
  group_by(species) %>%
  summarise(base_value = mean(base_value, na.rm = TRUE), .groups = "drop")

dados_funcional <- dados %>%
  select(-base_value) %>%
  inner_join(base_value_df, by = "species") %>%
  mutate(
    functional_value = base_value * current_area_km2,
    functional_area = base_value * future_area_km2,
    functional_loss = (functional_value - functional_area) / functional_value
  )

# ------------------------------------------------------------------------------
# 4. Functional diversity metrics (Community)
# ------------------------------------------------------------------------------
# Traits for FD analysis
traits_fd <- traits_proc %>%
  mutate(across(c(functional_group, iucn_status, foraging_stratum, activity_cycle), as.factor)) %>%
  column_to_rownames("species") %>%
  select(body_mass, dispersal_km, functional_group, iucn_status, foraging_stratum, activity_cycle)

# Presence-absence matrix by community
presenca_fd <- dados_funcional %>%
  mutate(presente = 1) %>%
  unite("comunidade", scenario_rcp, year, region, remove = FALSE) %>%
  select(comunidade, species, presente) %>%
  pivot_wider(names_from = species, values_from = presente, values_fill = 0) %>%
  column_to_rownames("comunidade")

# Calculate FD metrics
resultado_fd <- dbFD(x = traits_fd, a = presenca_fd, calc.FRic = TRUE)

# Final result formatting
metricas_fd <- as.data.frame(resultado_fd[c("FRic", "FDis")])
metricas_fd$redundancia <- resultado_fd$redundancy
metricas_fd$comunidade <- rownames(metricas_fd)

metricas_fd <- metricas_fd %>%
  separate(comunidade, into = c("scenario_rcp", "year", "region"), sep = "_") %>%
  mutate(year = as.integer(year))

# Export Table Results
write_csv(dados_funcional, "outputs/tables/valores_funcionais_multitraco.csv")
write_csv(metricas_fd, "outputs/tables/diversidade_funcional_multitraco.csv")


# ------------------------------------------------------------------------------
# 5. Functional Loss by Ecosystem Service (Visualization)
# ------------------------------------------------------------------------------
# Note: Ensure 'functional_values_services.csv' is in the data folder.
dados_service <- read_csv("data/functional_values_services.csv", show_col_types = FALSE)

# Calculate losses by service
perda_servico <- dados_service %>%
  group_by(Ecosystem_service, scenario_rcp, year, region) %>%
  summarise(
    functional_value = sum(functional_value, na.rm = TRUE),
    functional_area = sum(functional_area, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    perda_absoluta = functional_area - functional_value,
    perda_relativa = round(perda_absoluta / functional_value, 4)
  )

# Rank services by mean loss
ordem_servico <- perda_servico %>%
  group_by(Ecosystem_service) %>%
  summarise(media_perda = mean(perda_relativa, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(media_perda)) %>%
  pull(Ecosystem_service)

# Convert to ordered factor
perda_servico$Ecosystem_service <- factor(perda_servico$Ecosystem_service, levels = ordem_servico)

# Create Plot
viridis_colors <- viridis(2, option = "D")
colors <- c("PA" = viridis_colors[1], "UC" = viridis_colors[2])

grafico <- ggplot(perda_servico, aes(x = perda_relativa, y = Ecosystem_service, fill = region)) +
  geom_col(position = "dodge") +
  facet_wrap(~ scenario_rcp) +
  scale_fill_manual(values = colors) +
  labs(
    x = "Relative functional loss",
    y = "Ecosystem service"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    strip.text = element_text(size = 14)
  )

# Export Figures
formatos <- c("pdf", "png", "tiff")
for (formato in formatos) {
  ggsave(
    filename = paste0("outputs/figures/figure_functional_loss_by_service.", formato),
    plot = grafico,
    width = 12,
    height = 8,
    units = "in",
    dpi = 600
  )
}