library(tidyverse)
pagoda <- readRDS("Desktop/pagodaMouseMicroMarkers.rds")
DB <- read.csv("Desktop/PanglaoDB_markers_27_Mar_2020.tsv", sep = "\t")
DBfil <- DB %>% 
  select(species, official.gene.symbol, cell.type) %>% 
  filter(species %in% c("Mm", "Mm Hs")) %>% 
  filter(cell.type %in% c("Adipocytes", "B cells", "B cells memory", "B cells naive", "Basal cells",
                          "Basophils", "Chondrocytes", "Dendritic cells", "Endothelial cells", "Eosinophils",
                          "Epithelial cells", "Fibroblasts", "Keratinocytes", "Luminal epithelial cells", 
                          "Macrophages", " Mammary epithelial cells", "Mast cells", "Melanocytes", "Monocytes",
                          "Neutrophils", "NK cells", "Pericytes", "Smooth muscle cells", "T cells"))
table(DBfil$cell.type)

# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host = "https://apr2018.archive.ensembl.org")
  mouse = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", host = "https://apr2018.archive.ensembl.org")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

DBcom <- DBfil[,2:3]

celltypes <- list()
#i = "Adipocytes"
for (i in unique(DBcom$cell.type)) {
  print(i)
  celltype <- DBcom %>% filter(cell.type == i)
  genes <- convertHumanGeneList(celltype$official.gene.symbol)
  celltypes[[i]] <-  genes
}

saveRDS(celltypes, "Desktop/pagodaMouseMicroMarkers2.rds")
