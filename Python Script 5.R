
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

#---------------------------------------------------------------------------------
excel_file = "/Users/nishitha/Dropbox/My Mac (Nishitha MacBook Pro)/Desktop/CLIENT_DIR_NEC/Heatmaps/NEC_HEATMAP_EXCEL_FILE_SHORT.xlsx"
ann_path = "/Users/nishitha/Dropbox/My Mac (Nishitha MacBook Pro)/Desktop/CLIENT_DIR_NEC/Heatmaps/NEC_Annotation_Data.tsv"
setwd(dirname(excel_file))
dir.create("Heatmap_Plot_Selected_Ant")
setwd("Heatmap_Plot_Selected_Ant")

# Read all sheets into a list of data frames
sheet_names <- excel_sheets(excel_file)
sheet_list <- lapply(sheet_names, function(sheet) {
  data.frame(read_excel(excel_file, sheet = sheet), row.names = 1)
})

# reading the annotation data 
ann_data_1 = read.delim(ann_path)

# Access each data frame along with the sheet name
for (i in seq_along(sheet_list)) {
  ant = sheet_names[i]
  df = sheet_list[[i]]
  print(ant)
  print(nrow(ann_data_1))
  print(nrow(df))
  ann_data <- ann_data_1[row.names(ann_data_1) %in% row.names(df), ]
  print(nrow(ann_data))
  if (nrow(df) == nrow(ann_data)){
    #adding row annotation
    row_ann = rowAnnotation(Sex = ann_data$Sex,
                            Diabetes = ann_data$Diabetes,
                            Hypertension = ann_data$Hyertension,
                            Subject.Outcome = ann_data$Subject_Outcome,
                            Change.in.Antibiotic = ann_data$Change_in_Antibiotic,
                            Sequence.Type = ann_data$Sequence_Type,
                            Phylogroup = ann_data$Phylogroup)

    main_mat_1 = as.matrix(df)
    main_heatmap <-  Heatmap(main_mat_1, name = ant,
                             border = FALSE,
                             na_col = "red",
                             cluster_rows = TRUE, show_row_dend = TRUE, show_row_names = TRUE,
                             row_names_gp = gpar(fontsize = 10), row_names_side = "left",
                             col = circlize::colorRamp2(c(0, 1), colors = c("grey", "black")),
                             cell_fun = function(j, i, x, y, width, height, fill){
                               grid.rect(gp = gpar(fill = NA, col = "white", lwd = 1),
                                         x = x, y = y, width = width, height = height)},
                             right_annotation = row_ann)
    img = paste0(ant,"_NEC_Heatmap_Selected.png")
    png(img, width = 22, height = 8, units = "in", res = 300)
    print(main_heatmap)
    dev.off()
  } else {
    print("NOT SUFFICIANT")
  }
}



