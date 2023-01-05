library(tidyverse)
DB <- read.csv("/Desktop/PanglaoDB_markers_27_Mar_2020.tsv", sep = "\t")
DBfil <- DB %>% 
  select(species, official.gene.symbol, cell.type) %>% 
  filter(species %in% c("Hs", "Mm Hs")) %>% 
  filter(cell.type %in% c("Adipocytes", "B cells", "B cells memory", "B cells naive", "Basal cells",
                          "Basophils", "Chondrocytes", "Dendritic cells", "Endothelial cells", "Eosinophils",
                          "Epithelial cells", "Fibroblasts", "Keratinocytes", "Luminal epithelial cells", 
                          "Macrophages", " Mammary epithelial cells", "Mast cells", "Melanocytes", "Monocytes",
                          "Neutrophils", "NK cells", "Pericytes", "Smooth muscle cells", "T cells"))
table(DBfil$cell.type)
DBcom <- DBfil[,2:3]

celltypes <- list()
i = "Adipocytes"
for (i in unique(DBcom$cell.type)) {
  print(i)
  celltype <- DBcom %>% filter(cell.type == i)
  celltypes[[i]] <- celltype$official.gene.symbol
}

saveRDS(celltypes, "Desktop/pagodaMouseMicroMarkers2.rds")
