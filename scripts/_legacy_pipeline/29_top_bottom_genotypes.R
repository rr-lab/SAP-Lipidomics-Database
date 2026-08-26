# Top and bottom genotypes by aggregate fitted lipid signal within each condition.

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

dir.create("fig/supp", recursive = TRUE, showWarnings = FALSE)
dir.create("table/supp", recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  CTL = "data/SPATS_fitted/non_normalized_intensities/Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv",
  LIN = "data/SPATS_fitted/non_normalized_intensities/Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv"
)

core_classes <- c(
  "TG", "PC", "PG", "PE", "PA", "PS", "LPC", "LPE",
  "SQDG", "DG", "MG", "DGDG", "MGDG"
)

lipid_class <- function(x) {
  case_when(
    str_detect(x, "^DGDG\\(") ~ "DGDG",
    str_detect(x, "^MGDG\\(") ~ "MGDG",
    str_detect(x, "^LPC\\(")  ~ "LPC",
    str_detect(x, "^LPE\\(")  ~ "LPE",
    str_detect(x, "^SQDG\\(") ~ "SQDG",
    str_detect(x, "^TG\\(")   ~ "TG",
    str_detect(x, "^DG\\(")   ~ "DG",
    str_detect(x, "^MG\\(")   ~ "MG",
    str_detect(x, "^PC\\(")   ~ "PC",
    str_detect(x, "^PG\\(")   ~ "PG",
    str_detect(x, "^PE\\(")   ~ "PE",
    str_detect(x, "^PA\\(")   ~ "PA",
    str_detect(x, "^PS\\(")   ~ "PS",
    TRUE ~ NA_character_
  )
}

summarize_condition <- function(path, condition) {
  dat <- vroom(path, show_col_types = FALSE)
  id_col <- names(dat)[1]
  names(dat)[1] <- "Genotype"

  candidate_cols <- setdiff(names(dat), c("Genotype", "PlotID", "row", "col"))
  class_map <- tibble(
    Lipid = candidate_cols,
    Class = lipid_class(candidate_cols)
  ) %>%
    filter(Class %in% core_classes)

  if (nrow(class_map) == 0) {
    stop("No core lipid columns identified in ", path)
  }

  class_totals <- lapply(core_classes, function(cls) {
    cols <- class_map$Lipid[class_map$Class == cls]
    if (length(cols) == 0) return(rep(0, nrow(dat)))
    rowSums(as.data.frame(dat[, cols, drop = FALSE]), na.rm = TRUE)
  })
  names(class_totals) <- paste0("Sum_", core_classes)

  bind_cols(
    tibble(
      Condition = condition,
      Genotype = as.character(dat$Genotype),
      SpeciesIncluded = nrow(class_map)
    ),
    as_tibble(class_totals)
  ) %>%
    mutate(
      AggregateLipidSignal = rowSums(across(starts_with("Sum_")), na.rm = TRUE),
      SourceFile = path,
      OriginalIDColumn = id_col
    )
}

all_scores <- bind_rows(
  summarize_condition(input_files[["CTL"]], "CTL"),
  summarize_condition(input_files[["LIN"]], "LIN")
) %>%
  filter(is.finite(AggregateLipidSignal), AggregateLipidSignal > 0)

extremes <- all_scores %>%
  group_by(Condition) %>%
  arrange(AggregateLipidSignal, .by_group = TRUE) %>%
  mutate(
    RankAscending = row_number(),
    RankDescending = n() - row_number() + 1,
    Group = case_when(
      RankAscending <= 10 ~ "Bottom 10",
      RankDescending <= 10 ~ "Top 10",
      TRUE ~ NA_character_
    ),
    RankWithinGroup = case_when(
      Group == "Bottom 10" ~ RankAscending,
      Group == "Top 10" ~ RankDescending,
      TRUE ~ NA_integer_
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(Group)) %>%
  arrange(Condition, factor(Group, levels = c("Top 10", "Bottom 10")), RankWithinGroup)

write.csv(
  extremes,
  "table/supp/SuppTable_Top_Bottom_Lipid_Genotypes.csv",
  row.names = FALSE
)

ranking_colors <- c(
  "CTL Top 10" = "#006D2C",
  "CTL Bottom 10" = "#74C476",
  "LIN Top 10" = "#D94801",
  "LIN Bottom 10" = "#FDAE6B"
)

class_colors <- c(
  PC = "#00441B", PA = "#1B7837", PE = "#41AB5D",
  LPC = "#66C2A4", LPE = "#2CA25F", PG = "#78C679",
  PS = "#C2E699", DG = "#54278F", DGDG = "#F768A1",
  MG = "#8941ED", MGDG = "#FBB4D9", SQDG = "#9D4D6C",
  TG = "#ED804A"
)

# Show the contribution of every core lipid-class sum for the selected extremes.
class_sum_long <- extremes %>%
  select(Condition, Genotype, Group, RankWithinGroup, starts_with("Sum_")) %>%
  pivot_longer(
    cols = starts_with("Sum_"),
    names_to = "Class",
    values_to = "ClassSum"
  ) %>%
  mutate(
    Class = str_remove(Class, "^Sum_"),
    Class = factor(Class, levels = core_classes),
    PanelGenotype = paste(Condition, Group, Genotype, sep = "__"),
    Panel = factor(
      paste(Condition, Group),
      levels = c("CTL Top 10", "LIN Top 10", "CTL Bottom 10", "LIN Bottom 10")
    )
  )

genotype_levels <- class_sum_long %>%
  distinct(Condition, Group, RankWithinGroup, PanelGenotype) %>%
  arrange(
    factor(Condition, levels = c("CTL", "LIN")),
    factor(Group, levels = c("Top 10", "Bottom 10")),
    desc(RankWithinGroup)
  ) %>%
  pull(PanelGenotype)

class_sum_long <- class_sum_long %>%
  mutate(
    PanelGenotype = factor(PanelGenotype, levels = genotype_levels),
    Group = factor(Group, levels = c("Top 10", "Bottom 10")),
    Condition = factor(Condition, levels = c("CTL", "LIN"))
  )

class_sum_figure <- ggplot(
  class_sum_long,
  aes(x = ClassSum, y = PanelGenotype, fill = Class)
) +
  geom_col(width = 0.76, color = "white", linewidth = 0.12) +
  facet_wrap(~Panel, ncol = 2, scales = "free") +
  scale_fill_manual(values = class_colors, drop = FALSE) +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_x_continuous(
    labels = label_number(scale_cut = cut_short_scale()),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(x = "Fitted lipid-class sum", y = NULL, fill = "Lipid class") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.45),
    axis.text = element_text(color = "black"),
    strip.text = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "grey94", color = "grey65"),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.margin = margin(8, 12, 8, 8)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ggsave(
  "fig/supp/SuppFig_Lipid_Class_Sums_Top_Bottom_Genotypes.png",
  class_sum_figure,
  width = 15,
  height = 12,
  dpi = 300,
  bg = "white"
)

# Rank genotypes separately for each class using only species from that class.
class_specific_extremes <- all_scores %>%
  select(Condition, Genotype, starts_with("Sum_")) %>%
  pivot_longer(
    cols = starts_with("Sum_"),
    names_to = "Class",
    values_to = "ClassSignal"
  ) %>%
  mutate(Class = str_remove(Class, "^Sum_")) %>%
  group_by(Condition, Class) %>%
  arrange(ClassSignal, .by_group = TRUE) %>%
  mutate(
    RankAscending = row_number(),
    RankDescending = n() - row_number() + 1,
    Group = case_when(
      RankAscending <= 10 ~ "Bottom 10",
      RankDescending <= 10 ~ "Top 10",
      TRUE ~ NA_character_
    ),
    RankWithinGroup = case_when(
      Group == "Bottom 10" ~ RankAscending,
      Group == "Top 10" ~ RankDescending,
      TRUE ~ NA_integer_
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(Group)) %>%
  mutate(
    RankingGroup = paste(Condition, Group),
    Panel = factor(
      RankingGroup,
      levels = c("CTL Top 10", "LIN Top 10", "CTL Bottom 10", "LIN Bottom 10")
    )
  )

write.csv(
  class_specific_extremes,
  "table/supp/SuppTable_ClassSpecific_Top_Bottom_Genotypes.csv",
  row.names = FALSE
)

dir.create(
  "fig/supp/class_specific_genotype_rankings",
  recursive = TRUE,
  showWarnings = FALSE
)

make_class_specific_figure <- function(class_name) {
  plot_df <- class_specific_extremes %>%
    filter(Class == class_name) %>%
    mutate(
      OrderValue = if_else(Group == "Top 10", ClassSignal, -ClassSignal),
      PanelGenotype = paste(Panel, Genotype, sep = "__")
    )

  genotype_order <- plot_df %>%
    distinct(Panel, PanelGenotype, OrderValue) %>%
    arrange(Panel, OrderValue) %>%
    pull(PanelGenotype)

  plot_df <- plot_df %>%
    mutate(PanelGenotype = factor(PanelGenotype, levels = genotype_order))

  ggplot(plot_df, aes(x = ClassSignal, y = PanelGenotype, fill = RankingGroup)) +
    geom_col(width = 0.72, color = "grey25", linewidth = 0.25) +
    geom_text(
      aes(label = label_number(scale_cut = cut_short_scale(), accuracy = 0.1)(ClassSignal)),
      hjust = 1.08,
      color = "white",
      fontface = "bold",
      size = 3.2
    ) +
    facet_wrap(~Panel, ncol = 2, scales = "free") +
    scale_fill_manual(values = ranking_colors, name = NULL) +
    scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
    scale_x_continuous(
      labels = label_number(scale_cut = cut_short_scale()),
      expand = expansion(mult = c(0, 0.04))
    ) +
    labs(x = paste0("Fitted ", class_name, " lipid sum"), y = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.45),
      axis.text = element_text(color = "black"),
      strip.text = element_text(face = "bold", size = 14),
      strip.background = element_rect(fill = "grey94", color = "grey65"),
      legend.position = "bottom",
      plot.margin = margin(8, 12, 8, 8)
    )
}

invisible(lapply(core_classes, function(class_name) {
  class_figure <- make_class_specific_figure(class_name)
  output_path <- file.path(
    "fig/supp/class_specific_genotype_rankings",
    paste0("SuppFig_ClassSpecific_", class_name, "_Genotypes.png")
  )
  ggsave(
    output_path,
    class_figure,
    width = 14,
    height = 11,
    dpi = 300,
    bg = "white"
  )
}))

# Consistently upper-ranked genotypes: top decile for a class in both conditions.
class_percentiles <- all_scores %>%
  select(Condition, Genotype, starts_with("Sum_")) %>%
  pivot_longer(
    cols = starts_with("Sum_"),
    names_to = "Class",
    values_to = "ClassSignal"
  ) %>%
  mutate(Class = str_remove(Class, "^Sum_")) %>%
  group_by(Condition, Class) %>%
  mutate(Percentile = percent_rank(ClassSignal)) %>%
  ungroup()

rank_stability <- class_percentiles %>%
  select(Condition, Genotype, Class, ClassSignal, Percentile) %>%
  pivot_wider(
    names_from = Condition,
    values_from = c(ClassSignal, Percentile),
    names_glue = "{Condition}_{.value}"
  ) %>%
  filter(!is.na(CTL_Percentile), !is.na(LIN_Percentile)) %>%
  mutate(
    MeanPercentile = (CTL_Percentile + LIN_Percentile) / 2,
    PercentileDifference = abs(CTL_Percentile - LIN_Percentile),
    ConsistentlyUpperRanked = CTL_Percentile >= 0.90 & LIN_Percentile >= 0.90,
    Class = factor(Class, levels = core_classes)
  )

consistent_high <- rank_stability %>%
  filter(ConsistentlyUpperRanked) %>%
  group_by(Genotype) %>%
  mutate(StableHighClassCount = n()) %>%
  ungroup() %>%
  arrange(desc(StableHighClassCount), desc(MeanPercentile), Genotype, Class)

repeated_consistent_high <- consistent_high %>%
  filter(StableHighClassCount >= 2)

write.csv(
  rank_stability,
  "table/supp/SuppTable_ClassSpecific_CTL_LIN_Percentile_Ranks.csv",
  row.names = FALSE
)
write.csv(
  consistent_high,
  "table/supp/SuppTable_Consistently_Upper_Ranked_Genotypes.csv",
  row.names = FALSE
)

candidate_counts <- consistent_high %>%
  distinct(Genotype, Class) %>%
  count(Genotype, name = "StableHighClassCount") %>%
  filter(StableHighClassCount >= 2) %>%
  arrange(desc(StableHighClassCount), Genotype)

top_candidate_summary <- repeated_consistent_high %>%
  arrange(Genotype, desc(MeanPercentile)) %>%
  group_by(Genotype) %>%
  summarise(
    StableHighClassCount = first(StableHighClassCount),
    Classes = paste(as.character(Class), collapse = "; "),
    MeanPercentileAcrossCalls = mean(MeanPercentile),
    MinimumCTLPercentile = min(CTL_Percentile),
    MinimumLINPercentile = min(LIN_Percentile),
    .groups = "drop"
  ) %>%
  arrange(desc(StableHighClassCount), desc(MeanPercentileAcrossCalls), Genotype)

write.csv(
  top_candidate_summary,
  "table/supp/SuppTable_Top_Consistently_Upper_Ranked_Genotypes.csv",
  row.names = FALSE
)

candidate_order <- rev(candidate_counts$Genotype)

consistency_matrix <- expand_grid(
  Genotype = candidate_counts$Genotype,
  Class = core_classes
) %>%
  left_join(
    consistent_high %>%
      select(Genotype, Class, MeanPercentile),
    by = c("Genotype", "Class")
  ) %>%
  mutate(
    Genotype = factor(Genotype, levels = candidate_order),
    Class = factor(Class, levels = core_classes)
  )

p_consistency_heatmap <- ggplot(
  consistency_matrix,
  aes(x = Class, y = Genotype, fill = MeanPercentile)
) +
  geom_tile(color = "grey82", linewidth = 0.35) +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0.90, 1),
    breaks = c(0.90, 0.95, 1),
    labels = label_percent(),
    na.value = "white",
    name = "Mean percentile"
  ) +
  labs(x = "Lipid class", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black", linewidth = 0.4),
    legend.position = "bottom",
    plot.margin = margin(8, 3, 8, 8)
  )

p_candidate_counts <- ggplot(
  candidate_counts %>% mutate(Genotype = factor(Genotype, levels = candidate_order)),
  aes(x = StableHighClassCount, y = Genotype)
) +
  geom_col(fill = "#277F8E", width = 0.72, color = "grey25", linewidth = 0.25) +
  geom_text(aes(label = StableHighClassCount), hjust = -0.25, size = 3.2) +
  scale_x_continuous(
    breaks = 0:max(candidate_counts$StableHighClassCount),
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(x = "Number of classes", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    plot.margin = margin(8, 8, 8, 0)
  )

p_consistency_matrix <- p_consistency_heatmap | p_candidate_counts +
  plot_layout(widths = c(5.2, 1.2))

ggsave(
  "fig/supp/SuppFig_Top_Consistently_Upper_Ranked_Genotypes.png",
  p_consistency_matrix,
  width = 15,
  height = 8,
  dpi = 300,
  bg = "white"
)

message("Saved: fig/supp/SuppFig_Lipid_Class_Sums_Top_Bottom_Genotypes.png")
message("Saved: fig/supp/SuppFig_Top_Consistently_Upper_Ranked_Genotypes.png")
message("Saved: table/supp/SuppTable_Top_Bottom_Lipid_Genotypes.csv")
message("Saved: table/supp/SuppTable_ClassSpecific_Top_Bottom_Genotypes.csv")
message("Saved: table/supp/SuppTable_ClassSpecific_CTL_LIN_Percentile_Ranks.csv")
message("Saved: table/supp/SuppTable_Consistently_Upper_Ranked_Genotypes.csv")
message("Saved: table/supp/SuppTable_Top_Consistently_Upper_Ranked_Genotypes.csv")
message("Saved class-specific figures under: fig/supp/class_specific_genotype_rankings/")
