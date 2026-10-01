# ==============================================================================
# The frozen lipid species set, as one supplementary table.
#
#   Rscript scripts/new_new_script/79_build_frozen_species_table.R
#
# WHY THIS EXISTS. main.tex cites one table for two different things
#
#   "Classes follow the frozen species set"                 the class of each species
#   "The frozen species set is given in Supplementary ..."  the set itself
#   "...following the LIPID MAPS classification"            category, class, subclass
#
# and until now all three pointed at data/metadata/final_lipid_classes.csv, which
# is the class LOOKUP, not the species set. It carries 316 rows because it still
# holds every annotation the pipeline ever saw, including the 50 that the
# 2026-09-16 curation merged away or dropped -- 13-Keto-9Z,11E-octadecadienoic
# acid among them. A reader counting rows in that table gets 316 where the paper
# says 266, and finds names the paper says are not in the set.
#
# This joins the curated set (data/final_species_set/species_inventory.csv, the
# 266 species, 214 under CTL and 216 under LIN) to the subclass tiers of the
# lookup, so one table answers all three citations and its row count is the
# number the Results quote.
#
# The lookup's own Class column is the LIPID MAPS category, not the shorthand
# class, so it is dropped in favour of the inventory's two explicit columns:
# Class is the shorthand carried in the species name (TG, PC, LPC, ... , Other),
# Category is the LIPID MAPS category.
#
# Input   data/final_species_set/species_inventory.csv
#         data/metadata/final_lipid_classes.csv
# Output  table/supp/SuppTable_S32_Frozen_species_set.csv
# ==============================================================================
suppressPackageStartupMessages(library(data.table))

REPO <- Sys.getenv("SOLD_REPO", ".")
inv  <- fread(file.path(REPO, "data/final_species_set/species_inventory.csv"))
luk  <- fread(file.path(REPO, "data/metadata/final_lipid_classes.csv"))

stopifnot(nrow(inv) == 266, all(inv$Species == inv$Normalized))
miss <- setdiff(inv$Species, luk$Lipids)
if (length(miss)) stop(length(miss), " curated species absent from the class lookup: ",
                       paste(head(miss, 5), collapse = ", "))

out <- merge(inv[, .(Species, Class, Category = SuperClass, In_CTL, In_LIN, Status)],
             luk[, .(Species = Lipids, CommonName, Category_Code, SubClass, Sub_subclass)],
             by = "Species", all.x = TRUE)
setcolorder(out, c("Species", "CommonName", "Class", "Category", "Category_Code",
                   "SubClass", "Sub_subclass", "In_CTL", "In_LIN", "Status"))
setorder(out, Category, Class, Species)

fwrite(out, file.path(REPO, "table/supp/SuppTable_S32_Frozen_species_set.csv"))
message("Saved: table/supp/SuppTable_S32_Frozen_species_set.csv  (", nrow(out), " species)")
cat("\n  detected under CTL:", sum(out$In_CTL), "   under LIN:", sum(out$In_LIN),
    "   both:", sum(out$In_CTL & out$In_LIN), "\n\n")
print(out[, .N, by = .(Category, Class)][order(Category, -N)])
