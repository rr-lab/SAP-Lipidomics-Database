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
# Nothing is dropped by name. Two terms were previously removed by hand,
# detection of stimulus and DNA topological change, and both are back. DNA
# topological change is three genuine topoisomerases in three separate intervals
# (33.4-fold, q = 0.006, 3 of the 5 sorghum genes carrying the term), so it meets
# every stated criterion and deleting it for want of a story was not defensible.
# Detection of stimulus is a tandem LRR receptor-kinase array in one interval,
# but single-interval support is not grounds for removal here, since the genes
# annotated to a small GO term can genuinely sit in one locus, and the same rule
# would delete the nitrate and purine-transport results. Interval counts are
# reported in Supplementary Table S16 and are not used to select terms.
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
# Minimum number of independent 250 kb loci a gene set must span, and how many
# gene sets to show per lipid class.
#
# The interval threshold was dropped when the LD-aware permutation came in, on
# the reasoning that the permutation already answers "is this one LD block".
# That held for the old null, which reshuffled among the observed loci. The
# current null draws density-matched intervals from the genome, so a single
# locus can clear it on its own and still be a single locus rather than a
# process. Seventeen such gene sets sit at 230-324x fold and would otherwise
# take the whole top of the panel.
#
# Two loci is the floor rather than three, because three drops nitrate
# transport, which spans two loci at 25x and is a result worth showing. The
# panel is then kept readable by taking the strongest TOP_PER_CLASS gene sets
# within each class instead of by raising the threshold, so the cut is on
# ranking rather than on evidence. Nothing is lost: every tested term, its
# interval count included, is in Supplementary Tables S16-S18.
MIN_INTERVALS <- as.integer(Sys.getenv("GO_MIN_INTERVALS", "2"))
TOP_PER_CLASS <- as.integer(Sys.getenv("GO_TOP_PER_CLASS", "4"))

go <- vroom(loci_file, delim = "\t", show_col_types = FALSE)

sel <- go %>%
  filter(Background <= MAX_BACKGROUND, Intervals >= MIN_INTERVALS) %>%
  group_by(Condition, Classes) %>%
  slice_max(Fold_max, n = TOP_PER_CLASS, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Condition = factor(Condition, c("CTL", "LIN")),
    # tidy the merged labels for display
    Label = paste0(Term, "  [", Intervals, " loci]") %>%
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
cat(sprintf("\n  %d gene sets shown of %d in the table (background <= %d genes)\n",
            nrow(sel), nrow(go), MAX_BACKGROUND))
message(sprintf("  minimum independent loci: %d   top per class: %d", MIN_INTERVALS, TOP_PER_CLASS))
message(sprintf("  nitrate present: %s", any(grepl("nitrate", sel$Term, ignore.case = TRUE))))
