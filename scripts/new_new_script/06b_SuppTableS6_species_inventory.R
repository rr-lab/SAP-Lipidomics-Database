# ==============================================================================
# Supplementary Tables S6a, S6b, S6c -- the lipid-species inventory, and the
# frozen species set the inventory is built from.
#
#   Rscript scripts/new_new_script/06b_SuppTableS6_species_inventory.R
#
# WHAT COUNTS AS A SPECIES. Every annotated feature in a trial matrix counts.
# read_trial already drops the four non-lipid columns (LineRaw, PlotID, row,
# col), so what remains is the annotated lipidome as measured.
#
# THIS REPLACES AN EARLIER RULE, 2026-09-16. Until now a feature counted as a
# species only if its name contained a parenthesis, which kept the lipid
# shorthand -- PC(16:0/18:2), TG(...), AEG(o-16:2/16:0) -- and silently dropped
# everything carried under a trivial chemical name: the sterols, carotenoids,
# tocopherols, prenylquinones and a few oxylipins and hormones, 20 features in
# CTL and 26 in LIN. That rule was never stated in the manuscript and it did not
# match the rest of the pipeline. The composition analyses have always used every
# annotated feature, which is why the terpenoid superclass carries about 4% of
# %TIC while the carotenoids were absent from the species count, and the LION
# input in 06a also uses every feature, which is why it reports 164 shared where
# this script used to report 146. One definition now, applied everywhere.
#
# NAMES ARE NORMALISED BEFORE ANYTHING IS COUNTED. This fixes a real error in
# the previous version, which intersected the two trials on raw names. CTL writes
# a lyso species as CLASS(x:y/0:0) and LIN writes it LCLASS(x:y), so LPC(18:2)
# appeared as PC(18:2/0:0) under CTL and LPC(18:2) under LIN and was counted as
# CTL-only AND LIN-only rather than once as common. _common.R says in its own
# comment to normalise before doing anything else; this script now does.
#
# The old rule and the missing normalisation together give 194 / 190 / 146 /
# 48 / 44 / 238. Counting every feature and normalising first gives
# 214 / 216 / 164 / 50 / 52 / 266.
#
# Two artefact features, Phytosphingosine and SM(d18:1/17:0), were removed
# upstream and are already absent from the fitted matrices; nothing here re-drops
# them. See _legacy_pipeline/22_lipidome_class_composition.R for why.
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

# Every annotated feature, normalised before it is compared to anything.
feature_names <- function(path) {
  unique(normalize_lipid_name(names(read_trial(path))[-1]))
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
cat("\nspecies not named in lipid shorthand, which the old rule discarded\n")
non_shorthand <- sort(unique(inv$Species[!grepl("\\(", inv$Species)]))
cat("  ", length(non_shorthand), " of ", nrow(inv), ": ",
    paste(non_shorthand, collapse = ", "), "\n", sep = "")
