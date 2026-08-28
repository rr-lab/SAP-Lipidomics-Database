# ==============================================================================
# Figure 4 -- GO terms enriched among GWAS candidate genes, by lipid class.
#
# Built from the locus table written by 10_SuppTableS16_GO_collapse_and_rank.R, which has
# already done three things: kept only terms passing the LD-aware permutation
# (q_LD < 0.05), merged GO parent/child terms sharing an identical candidate
# gene set, and merged those gene sets across lipid class, trait layer and
# ontology so that one locus is one row rather than ten.
#
# Selection here is on how SPECIFIC the GO term is, not on fold enrichment and
# not on interval count. A term drawn from a background of thousands of genes
# cannot be interpreted whatever its statistics: zinc ion binding (1015 genes),
# transmembrane transport (526), nucleic acid binding (522) and mRNA binding
# (432) all pass the permutation test and none of them says anything about
# lipids. Requiring a background of at most MAX_BACKGROUND genes removes them
# and keeps the terms that name a process.
#
# Two further terms are dropped by name. Detection of stimulus is vague despite
# a small background, and DNA topological change is topoisomerase activity with
# no defensible link to lipid biology.
#
# Rows are grouped by the lipid class whose GWAS candidates produced them. That
# grouping is the point of the figure: SQDG carries both nutrient terms, nitrate
# and phosphate starvation, which is the sulfolipid-for-phospholipid response to
# phosphate limitation.
#
# Input : table/go_enrichment/Table_GO_enrichment_loci.tsv
# Output: fig/main/Figure5_GO_BP_LDaware.png
#         (filename kept because that is the path main.tex includes; the float
#          prints as Figure 4)
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({ library(stringr); library(forcats); library(tibble) })

loci_file <- Sys.getenv("GO_LOCI",
  file.path(REPO, "table/go_enrichment/Table_GO_enrichment_loci.tsv"))
out_name  <- "Figure5_GO_BP_LDaware.png"
stopifnot(file.exists(loci_file))

MAX_BACKGROUND <- 30
DROP_TERMS <- c("detection of stimulus", "DNA topological change")

go <- vroom(loci_file, delim = "\t", show_col_types = FALSE)

sel <- go %>%
  filter(Background <= MAX_BACKGROUND, !Term %in% DROP_TERMS) %>%
  mutate(
    Condition = factor(Condition, c("CTL", "LIN")),
    # tidy the merged labels for display
    Label = Term %>%
      str_replace(" \\(\\+\\d+ related terms\\)$", "") %>%
      str_replace("nicotianamine aminotransferase activity and L-tyrosine-2-oxoglutarate transaminase activity",
                  "nicotianamine aminotransferase / L-tyrosine transaminase") %>%
      str_replace("quercetin 7-O-glucosyltransferase activity and 3-O-glucosyltransferase activity",
                  "quercetin O-glucosyltransferase activity") %>%
      str_replace("purine nucleoside transmembrane transport",
                  "purine nucleoside / nucleobase transport"),
    Group = ifelse(N_classes > 1, paste0(Classes, "\n(shared)"), Classes)
  )

# Single-class groups first, ordered by how many terms they carry, then the
# gene sets shared across classes.
grp_order <- sel %>%
  group_by(Group) %>%
  summarise(shared = any(N_classes > 1), n = dplyr::n(),
            top = max(Fold_max), .groups = "drop") %>%
  arrange(shared, desc(n), desc(top)) %>% pull(Group)

sel <- sel %>%
  mutate(Group = factor(Group, levels = grp_order),
         Label = str_wrap(Label, 44)) %>%
  arrange(Group, Fold_max) %>%
  mutate(row = row_number(), Label = fct_reorder(Label, row))

fig <- ggplot(sel, aes(Fold_max, Label, colour = Condition)) +
  geom_segment(aes(x = 1, xend = Fold_max, yend = Label),
               colour = "grey80", linewidth = .5) +
  geom_point(size = 4.5) +
  geom_text(aes(label = sprintf("%.0f", Fold_max)),
            hjust = -0.55, size = 4, fontface = "bold", show.legend = FALSE) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = condition_colors) +
  scale_x_log10(breaks = c(1, 3, 10, 30, 60),
                expand = expansion(mult = c(0.02, 0.16))) +
  labs(x = "Fold enrichment among GWAS candidate genes (log scale)", y = NULL) +
  plot_theme +
  # The legend keeps plot_theme's inside top-right position. Bottom-right sits
  # exactly where the 63-fold purine row lands and hides its point.
  theme(axis.text.y        = element_text(size = 12),
        panel.grid.major.y = element_blank(),
        strip.placement    = "outside",
        strip.text.y.left  = element_text(angle = 0, face = "bold", size = 13,
                                          margin = margin(r = 8)),
        panel.spacing.y    = unit(4, "pt"))

save_fig(fig, out_name, width = 13, height = 11)

cat("\n-- terms shown ------------------------------------------------------\n")
print(as.data.frame(sel %>% arrange(Group, desc(Fold_max)) %>%
        mutate(Term = substr(Term, 1, 52)) %>%
        dplyr::select(Group, Condition, Term, Genes, Background, Fold_max, q_LD)))
cat(sprintf("\n  %d loci shown of %d in the table (background <= %d genes, 2 dropped by name)\n",
            nrow(sel), nrow(go), MAX_BACKGROUND))
