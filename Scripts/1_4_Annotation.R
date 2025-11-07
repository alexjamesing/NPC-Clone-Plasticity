library(Seurat)
library(SeuratObject)
library(pheatmap)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(Matrix)
library(R6)
library(dplyr)
library(tidyverse)
library(purrr)
library(harmony)
library(plyr)
library(ggridges)
#library(treedataverse)
library(ComplexHeatmap)
library(circlize)

## ------------------------------------------------------------------
## 0. Define I/O paths and load both meta-cell and single-cell data
## ------------------------------------------------------------------
base_folder <- "/omics/odcf/analysis/OE0574_projects/brainbreaks/single_cell_BB_E17-5_Alex/result_all_102025/R_files/Results"
ref <- readRDS(file.path(base_folder, "1_4_Annotation", "rds", "Linnarsson_ref_data_rd.rds"))
query <- readRDS(file.path(base_folder, "1_3_Metacells", "rds", "meta_cell.rds"))
single_cell_path <- file.path(base_folder, "1_3_Metacells", "rds", "sample_clusters.rds")
meta_map_path <- file.path(base_folder, "1_3_Metacells", "rds", "add_meta.rds")
out <- path.expand(file.path(base_folder, "new_annotation"))

if (!file.exists(single_cell_path)) {
  stop("Single-cell object not found at: ", single_cell_path)
}
if (!file.exists(meta_map_path)) {
  stop("Meta-group lookup not found at: ", meta_map_path)
}

meta_map <- readRDS(meta_map_path)
stopifnot("meta_group" %in% colnames(meta_map))

# Ensure that downstream export folders exist so annotation artefacts never overwrite existing results.
dir.create(file.path(out, "files"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "rds"), recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------------
## 1. Prepare the meta-cell object for label transfer
## ------------------------------------------------------------------
query <- UpdateSeuratObject(query)
# Drop stale reductions/graphs so the cell order check passes.
query <- DietSeurat(
  query,
  assays = "RNA",
  dimreducs = NULL,
  graphs = NULL
)

query[["percent.mt"]] <- PercentageFeatureSet(query, pattern = "^mt[-\\.]")
query <- NormalizeData(query)

hv_features <- intersect(VariableFeatures(ref), rownames(query))
anchors <- FindTransferAnchors(
  reference = ref,
  query = query,
  dims = 1:30,
  npcs = 30,
  reduction = "cca",
  features = hv_features
)
predictions_region <- TransferData(
  anchorset = anchors,
  refdata = ref$Region,
  dims = 1:30,
  weight.reduction = "cca"
)
readr::write_csv(
  predictions_region,
  file.path(out, "files", "Linnarsson_prediction_region.csv")
)
query$Linnarsson_region_transfer <- predictions_region$predicted.id
query$Linnarsson_region_scores <- predictions_region$prediction.score.max

saveRDS(query,file.path(out,"rds","query_Linnarsson_region.rds"))

## ------------------------------------------------------------------
## 2. Transfer subclass labels from the La Manno reference
## ------------------------------------------------------------------
predictions_subclass <- TransferData(
  anchorset = anchors,
  refdata = ref$Subclass,
  dims = 1:30,
  weight.reduction = "cca"
)
readr::write_csv(
  predictions_subclass,
  file.path(out, "files", "Predictions_class.csv")
)

## ------------------------------------------------------------------
## 3. Assemble region + subclass predictions for automated annotation
## ------------------------------------------------------------------
# Rely on atlas transfer scores (instead of manual lookups) while keeping
# confidence thresholds explicit for downstream QC.
prediction_summary <- predictions_subclass %>%
  tibble::rownames_to_column("cell_id") %>%                        # keep explicit cell ids
  dplyr::transmute(                                                # retain only relevant subclass info
    cell_id,
    subclass_pred = predicted.id,
    subclass_score = prediction.score.max
  ) %>%
  dplyr::left_join(
    predictions_region %>%
      tibble::rownames_to_column("cell_id") %>%                    # align on identical barcodes
      dplyr::transmute(
        cell_id,
        region_pred = predicted.id,
        region_score = prediction.score.max
      ),
    by = "cell_id"
  )
subclass_threshold <- 0.5
region_threshold <- 0.5
prediction_summary <- prediction_summary %>%
  dplyr::mutate(
    subclass_final = ifelse(
      subclass_score >= subclass_threshold,
      subclass_pred,
      "Unassigned_Subclass"
    ),
    region_final = ifelse(
      region_score >= region_threshold,
      region_pred,
      "Unassigned_Region"
    )
  )

## ------------------------------------------------------------------
## 4. Score lineage marker panels to complement atlas predictions
## ------------------------------------------------------------------
marker_panels <- list(
  Endothelia = c("Cldn5", "Igfbp7"),
  Microglia = c("Tmem119", "P2ry12", "Aif1", "C1qb", "Cx3cr1", "Cd68"),
  Pericyte = c("Pdgfrb", "Cspg4", "Rgs5", "Col3a1"),
  Radial_glia = c("Sox2", "Pax6", "Tnc", "Ptprz1", "Hopx", "Ascl1", "Egfr", "Olig2", "Mki67", "Top2a"),
  Neuroblast = c("Npa3", "Adgrb3", "Fgf15", "Wnt7b", "Tnc", "Hes6"),
  Glioblast = c("Aldh1l1", "Efna5", "Rspo1", "Rspo2", "Rspo3", "Wnt1", "Wnt9a", "Hes1"),
  Oligodendrocyte_precursor_cell = c("Olig2", "Pdgfra"),
  Choroid_plexus = c("Htr2c", "Enpp2"),
  Intermediate_progenitor = c("Eomes", "Neurog2", "Sox2"),
  Neuronal_general = c("Mki67", "Nes", "Dcx", "Neurod2", "Neurod6", "Map2", "Tubb3"),
  GABAergic = c("Isl1", "Slc32a1", "Gad2", "Gad1"),
  Glutamatergic = c("Neurod6", "Slc17a7", "Slc17a6", "Slc17a5"),
  Glycinergic = c("Slc6a5"),
  Cholinergic = c("Chat", "Slc18a3"),
  Serotonergic = c("Slc6a4", "Tph2"),
  Dopaminergic = c("Slc6a3", "Drd2", "Th"),
  Noradrenergic = c("Slc6a2", "Dbh")
)

marker_panels <- lapply(marker_panels, function(genes) {
  genes <- unique(genes)
  genes[genes %in% rownames(query)]
})
marker_panels <- marker_panels[lengths(marker_panels) > 0]

if (length(marker_panels) == 0) {
  stop("No marker genes were found in the query object; cannot perform marker-based annotation.")
}

score_prefix <- "marker_auto_"
query <- AddModuleScore(
  object = query,
  features = marker_panels,
  name = score_prefix,
  assay = DefaultAssay(query)
)

score_cols <- paste0(score_prefix, seq_along(marker_panels))
clean_label <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
renamed_scores <- paste0("score_", clean_label(names(marker_panels)))
for (i in seq_along(score_cols)) {
  query@meta.data[[renamed_scores[i]]] <- query@meta.data[[score_cols[i]]]
}
query@meta.data[score_cols] <- NULL

score_matrix <- as.matrix(query@meta.data[, renamed_scores, drop = FALSE])
max_idx <- max.col(score_matrix, ties.method = "first")
max_score <- score_matrix[cbind(seq_len(nrow(score_matrix)), max_idx)]
second_score <- apply(score_matrix, 1, function(x) {
  if (length(x) < 2) {
    return(NA_real_)
  }
  sort(x, decreasing = TRUE)[2]
})

marker_threshold <- 0
marker_labels <- names(marker_panels)[max_idx]
marker_labels[is.na(marker_labels)] <- "Unassigned_Marker"
marker_labels[max_score < marker_threshold | is.na(max_score)] <- "Unassigned_Marker"

query$marker_celltype <- marker_labels
query$marker_score_max <- max_score
query$marker_score_margin <- max_score - second_score

## ------------------------------------------------------------------
## 5. Blend marker + region/subclass calls into a final label
## ------------------------------------------------------------------
prediction_summary <- prediction_summary %>%
  tibble::column_to_rownames("cell_id")
prediction_summary <- prediction_summary[colnames(query), , drop = FALSE]

query$predicted_region <- prediction_summary$region_final
query$predicted_region_score <- prediction_summary$region_score
query$predicted_subclass <- prediction_summary$subclass_final
query$predicted_subclass_score <- prediction_summary$subclass_score

region_clean <- ifelse(
  is.na(query$predicted_region),
  "Unassigned_Region",
  query$predicted_region
)
region_clean[query$predicted_region_score < region_threshold | is.na(query$predicted_region_score)] <- "Unassigned_Region"

celltype_clean <- ifelse(
  is.na(query$marker_celltype),
  "Unassigned_Marker",
  query$marker_celltype
)

query$auto_annotation <- dplyr::case_when(
  celltype_clean == "Unassigned_Marker" & region_clean == "Unassigned_Region" ~ "Unassigned",
  celltype_clean == "Unassigned_Marker" ~ paste("Unassigned", region_clean, sep = " | "),
  region_clean == "Unassigned_Region" ~ paste(celltype_clean, "Unassigned", sep = " | "),
  TRUE ~ paste(celltype_clean, region_clean, sep = " | ")
)

Idents(query) <- factor(query$auto_annotation)

annotation_counts <- query@meta.data %>%
  tibble::as_tibble(rownames = "cell_id") %>%
  dplyr::count(auto_annotation, sort = TRUE)
readr::write_csv(annotation_counts, file.path(out, "files", "auto_annotation_counts.csv"))

marker_scores_long <- query@meta.data %>%
  tibble::as_tibble(rownames = "cell_id") %>%
  dplyr::select(cell_id, dplyr::starts_with("score_"), marker_celltype, marker_score_max, marker_score_margin)
readr::write_csv(marker_scores_long, file.path(out, "files", "marker_scores.csv"))

saveRDS(query, file.path(out, "rds", "query_annotations_auto.rds"))

## ------------------------------------------------------------------
## 6. Propagate meta-cell annotations back to every barcode
## ------------------------------------------------------------------
meta_lookup <- query@meta.data %>%
  tibble::as_tibble(rownames = "meta_group") %>%
  dplyr::select(
    meta_group,
    marker_celltype,
    auto_annotation,
    predicted_region,
    predicted_region_score,
    predicted_subclass,
    predicted_subclass_score,
    marker_score_max,
    marker_score_margin
  )

single_cells <- readRDS(single_cell_path)
single_cells <- AddMetaData(single_cells, meta_map)

missing_groups <- setdiff(unique(single_cells$meta_group), meta_lookup$meta_group)
if (length(missing_groups) > 0) {
  warning(
    "Some meta_group identifiers were not found in the annotated meta-cell object: ",
    paste(missing_groups, collapse = ", ")
  )
}

single_cells@meta.data <- single_cells@meta.data %>%
  tibble::as_tibble(rownames = "cell_id") %>%
  dplyr::left_join(meta_lookup, by = "meta_group") %>%
  tibble::column_to_rownames("cell_id")

if ("auto_annotation" %in% colnames(single_cells@meta.data)) {
  Idents(single_cells) <- factor(single_cells$auto_annotation)
}

propagated_meta <- single_cells@meta.data %>%
  tibble::as_tibble(rownames = "cell_id") %>%
  dplyr::select(
    cell_id,
    meta_group,
    marker_celltype,
    auto_annotation,
    predicted_region,
    predicted_region_score,
    predicted_subclass,
    predicted_subclass_score,
    marker_score_max,
    marker_score_margin
  )

saveRDS(single_cells, file.path(out, "rds", "single_cell_annotations.rds"))
readr::write_csv(
  propagated_meta,
  file.path(out, "files", "single_cell_annotations.csv")
)

# The meta-cell object now carries atlas-derived labels, and the same labels are
# mirrored on every barcode so downstream scripts (DE, proportions, pathways) can
# keep operating on the single-cell object without modification.
