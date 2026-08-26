# ==============================================================================
# Figure 2 -- Genotypic lipid variation is independent of race and structure
#   A  Top 10 high-variance lipid species, CTL
#   B  Top 10 high-variance lipid species, LIN
#   C  Variance explained by botanical race (6 groups)
#   D  Variance explained by marker-based genetic cluster (K = 6)
#
# Panels A/B reproduce the computation in 20_variance_species_overview.R
# (TIC-normalised, variance across genotypes, named species only).
# Panels C/D read the Kruskal-Wallis output of 33_race_structure_tests.py.
# ==============================================================================
suppressPackageStartupMessages({
  library(vroom); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(patchwork)
})

data_root <- Sys.getenv("SOLD_DATA", "/Users/nirwantandukar/Documents/Github/SoLD_paper/data")
race_csv  <- Sys.getenv("RACE_CSV",
  "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/table/race/race_structure_lipid_tests.csv")
out_file  <- Sys.getenv("FIG2_OUT",
  "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/fig/main/Figure2_Genotypic_Variation.png")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

# project palettes (identical to 20_variance_species_overview.R)
class_colors <- c(
  PC = "#00441B", PA = "#1B7837", PE = "#41AB5D", LPC = "#66C2A4", LPE = "#2CA25F",
  PG = "#78C679", PS = "#C2E699", DG = "#54278F", DGDG = "#F768A1", MG = "#8941ED",
  MGDG = "#FBB4D9", SQDG = "#9D4D6C", TG = "#ED804A"
)
condition_colors <- c(CTL = "#440154FF", LIN = "#FDE725FF")

plot_theme <- theme_minimal(base_size = 13) +
  theme(axis.text = element_text(colour = "black", size = 9),
        axis.title = element_text(face = "bold", size = 11),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        plot.tag = element_text(face = "bold", size = 16))

get_lipid_class <- function(x) {
  cls <- str_extract(x, "^[A-Za-z0-9]+(?=\\()")
  ifelse(cls %in% names(class_colors), cls, NA_character_)
}

# ---- A, B: high-variance species --------------------------------------------
spats <- file.path(data_root, "SPATS_fitted/non_normalized_intensities")
read_cond <- function(f) vroom(file.path(spats, f), show_col_types = FALSE) %>% dplyr::select(-c(2, 3, 4))
control_raw <- read_cond("Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv")
low_raw     <- read_cond("Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv")
names(control_raw)[1] <- names(low_raw)[1] <- "Compound_Name"

tic_normalize <- function(df) {
  X  <- df %>% dplyr::select(-Compound_Name) %>% dplyr::select(where(is.numeric)) %>% as.matrix()
  rs <- rowSums(X, na.rm = TRUE); rs[rs == 0 | is.na(rs)] <- NA_real_
  Xp <- sweep(X, 1, rs, FUN = "/") * 100; Xp[is.na(Xp)] <- 0
  bind_cols(df %>% dplyr::select(Compound_Name), as.data.frame(Xp))
}

top_var <- function(df_tic, cond, top_n = 10) {
  lc <- setdiff(names(df_tic), "Compound_Name"); lc <- lc[grepl("\\(", lc)]
  df_tic %>% dplyr::select(all_of(lc)) %>%
    mutate(across(everything(), ~suppressWarnings(as.numeric(.x)))) %>%
    summarise(across(everything(), ~var(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "Lipid", values_to = "Variance") %>%
    arrange(desc(Variance)) %>% mutate(Rank = row_number(), Condition = cond,
                                       Class = get_lipid_class(Lipid)) %>%
    filter(Rank <= top_n)
}

tv_c <- top_var(tic_normalize(control_raw), "CTL")
tv_l <- top_var(tic_normalize(low_raw),     "LIN")

var_panel <- function(d) {
  ggplot(d, aes(reorder(Lipid, Variance), Variance, fill = Class)) +
    geom_col(width = .7, colour = "black", linewidth = .3) + coord_flip() +
    scale_fill_manual(values = class_colors, na.value = "grey70") +
    scale_y_continuous(expand = expansion(mult = c(0, .06))) +
    labs(x = NULL, y = "Variance (% TIC)") +
    plot_theme + theme(panel.grid.major.y = element_blank(), legend.position = "right")
}
pA <- var_panel(tv_c); pB <- var_panel(tv_l)

# ---- C, D: variance explained by race and by genetic cluster ----------------
race <- vroom(race_csv, show_col_types = FALSE) %>%
  mutate(Condition = recode(Condition, Control = "CTL", LowInput = "LIN", CTL = "CTL", LIN = "LIN"),
         Condition = factor(Condition, c("CTL", "LIN")),
         sig = !is.na(q_BH) & q_BH < 0.05)

ord <- race %>% group_by(Trait) %>% summarise(m = max(epsilon2, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(Trait)
race <- race %>% mutate(Trait = factor(Trait, levels = rev(ord)))

grp_panel <- function(g) {
  d <- race %>% filter(Grouping == g)
  nsig <- sum(d$sig, na.rm = TRUE)
  ggplot(d, aes(Trait, epsilon2, fill = Condition)) +
    geom_col(position = position_dodge(.78), width = .7, colour = "black", linewidth = .25) +
    geom_text(data = subset(d, sig), aes(label = "***"), position = position_dodge(.78),
              hjust = -.2, size = 3.2) +
    coord_flip() +
    scale_fill_manual(values = condition_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, .12))) +
    annotate("text", x = 1, y = max(d$epsilon2, na.rm = TRUE),
             label = sprintf("%d trait%s significant after BH correction",
                             nsig, ifelse(nsig == 1, "", "s")),
             hjust = 1, vjust = -.2, size = 2.9, fontface = "italic", colour = "grey25") +
    labs(x = NULL, y = expression(bold(paste("Variance explained (", epsilon^2, ")")))) +
    plot_theme + theme(panel.grid.major.y = element_blank(), legend.position = "right")
}
pC <- grp_panel(unique(race$Grouping)[1])
pD <- grp_panel(unique(race$Grouping)[2])

row1 <- (pA | pB) + plot_layout(guides = "collect") & theme(legend.position = "right")
row2 <- (pC | pD) + plot_layout(guides = "collect") & theme(legend.position = "right")
fig  <- row1 / row2 + plot_annotation(tag_levels = "A")
ggsave(out_file, fig, width = 14, height = 10, dpi = 300, bg = "white")

cat("groupings:", paste(unique(race$Grouping), collapse = " | "), "\n")
cat("CTL top3:", paste(head(tv_c$Lipid, 3), collapse = ", "), "\n")
cat("LIN top3:", paste(head(tv_l$Lipid, 3), collapse = ", "), "\n")
cat("shared in top-10:", length(intersect(tv_c$Lipid, tv_l$Lipid)), "\n")
