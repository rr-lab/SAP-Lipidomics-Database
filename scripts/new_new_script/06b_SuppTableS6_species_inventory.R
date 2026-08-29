# ==============================================================================
# Supplementary Tables S6a, S6b, S6c -- the lipid-species inventory, and the
# frozen species set the inventory is built from.
#
#   Rscript scripts/new_new_script/06b_SuppTableS6_species_inventory.R
#
# THE COUNTING RULE, stated once so it never has to be rediscovered.
#
# A feature counts as a "lipid species" if its annotated name contains a
# parenthesis. That is the rule the original analysis used, and it is what makes
# the paper's 243 / 152 / 49 / 42 reproduce exactly. It keeps everything named in
# lipid shorthand -- PC(16:0/18:2), TG(...), AEG(o-16:2/16:0), Cer(d18:2/20:1) --
# and drops features carried under a trivial chemical name, which in this data
# set are the carotenoids, sterols, tocopherols, quinones and free sphingoid
# bases (23 features in CTL, 28 in LIN).
#
# It is a naming convention, not a biological criterion, and it is applied ONLY
# to the species inventory. The composition analyses in 06_Fig2 use every
# annotated feature, which is why the terpenoid superclass is ~4% of %TIC even
# though the carotenoids are not counted as species here. Methods says so.
#
# Two artefact features, Phytosphingosine and SM(d18:1/17:0), were removed
# upstream and are already absent from the fitted matrices; nothing here re-drops
# them. See _legacy_pipeline/22_lipidome_class_composition.R for why.
#
# Class labels use normalize_lipid_name from _common.R, so PC(18:0/0:0) is
# counted as LPC here and not as PC. The total is unaffected -- the species keeps
# its parenthesis either way -- but the per-class counts differ from the legacy
# table, which counted those three CTL species as PC.
#
# Inputs
#   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
#   data/lipid_class/final_lipid_classes.csv
#
# Outputs
#   data/final_species_set/species_inventory.csv   <- the frozen set, one row per species
#   table/supp/SuppTable_S6a_Species_Summary.csv
#   table/supp/SuppTable_S6b_Species_by_Class.csv
#   table/supp/SuppTable_S6c_Species_by_SuperClass.csv
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({ library(tidyr) })

class_csv <- Sys.getenv("LIPID_CLASS_CSV",
  file.path(DATA_ROOT, "lipid_class/final_lipid_classes.csv"))
stopifnot(file.exists(class_csv))

SET_DIR <- file.path(DATA_ROOT, "final_species_set")
dir.create(SET_DIR, recursive = TRUE, showWarnings = FALSE)

is_species <- function(x) grepl("\\(", x)          # <- the rule

feature_names <- function(path) {
  n <- names(read_trial(path))[-1]
  n[is_species(n)]
}

ctl <- feature_names(CTL_CSV)
lin <- feature_names(LIN_CSV)

all_sp <- sort(unique(c(ctl, lin)))
ann <- vroom(class_csv, show_col_types = FALSE) %>%
  transmute(key = tolower(normalize_lipid_name(Lipids)), SuperClass = Class) %>%
  distinct(key, .keep_all = TRUE)

inv <- tibble::tibble(
    Species    = all_sp,
    Normalized = normalize_lipid_name(all_sp),
    Class      = lipid_class(all_sp),
    key        = tolower(normalize_lipid_name(all_sp)),
    In_CTL     = all_sp %in% ctl,
    In_LIN     = all_sp %in% lin) %>%
  left_join(ann, by = "key") %>%
  mutate(SuperClass = ifelse(is.na(SuperClass), "Unclassified", SuperClass),
         Status = dplyr::case_when(In_CTL & In_LIN ~ "Common",
                                   In_CTL          ~ "CTL only",
                                   TRUE            ~ "LIN only")) %>%
  select(Species, Normalized, Class, SuperClass, In_CTL, In_LIN, Status)

write.csv(inv, file.path(SET_DIR, "species_inventory.csv"), row.names = FALSE)
message("Saved: ", file.path(SET_DIR, "species_inventory.csv"), "  (", nrow(inv), " species)")

# ---- S6a ---------------------------------------------------------------------
s6a <- tibble::tibble(
  Metric = c("Control Total", "LowInput Total", "Common",
             "Control Only", "LowInput Only", "All Unique"),
  Count  = c(length(ctl), length(lin), sum(inv$Status == "Common"),
             sum(inv$Status == "CTL only"), sum(inv$Status == "LIN only"), nrow(inv)))
save_table(s6a, "SuppTable_S6a_Species_Summary.csv")

# ---- S6b / S6c ---------------------------------------------------------------
tally <- function(col) {
  inv %>% group_by(Group = .data[[col]]) %>%
    summarise(Control = sum(In_CTL), LowInput = sum(In_LIN), .groups = "drop") %>%
    arrange(desc(Control + LowInput)) %>%
    rename(!!col := Group)
}
save_table(tally("Class"),      "SuppTable_S6b_Species_by_Class.csv")
save_table(tally("SuperClass"), "SuppTable_S6c_Species_by_SuperClass.csv")

cat("\n-- species inventory --\n"); print(as.data.frame(s6a))
cat("\nfeatures excluded by the rule (no parenthesis in the name)\n")
cat("  CTL ", sum(!is_species(names(read_trial(CTL_CSV))[-1])),
    " | LIN ", sum(!is_species(names(read_trial(LIN_CSV))[-1])), "\n", sep = "")
