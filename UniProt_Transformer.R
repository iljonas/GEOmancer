library(dplyr)
library(tidyr)

geneLocSep <- function(gList){
     gList <- list(gList)
     geneLocs <- as.integer(sapply(gList, function(x) {grep("^[a-z]",x)}))
     genes <- unlist(gList)[geneLocs]
     loci <- ifelse(is.na(geneLocs), unlist(gList), unlist(gList)[-geneLocs])
     return(list(genes,loci))
}

user <- as.character(Sys.info()["user"])

groupedData <- read.csv("C:\\Users\\iljonas\\Downloads\\Capstone Files\\Source Files\\_Locus Names\\My_List_E-modified.txt", sep = '\t') %>%
     mutate(Sep.names = strsplit(as.character(Gene.names), ';')) %>%
     unnest(Sep.names) %>%
     mutate(Sep.groups = strsplit(Sep.names, ' '))

split.names <- sapply(groupedData$Sep.groups, geneLocSep) %>%
     t %>%
     as.data.frame

names(split.names) <- c("Genes", "Loci")

newData <- data.frame(groupedData, split.names) %>%
     unnest(Loci, .drop = FALSE) %>%
     unnest(Genes, .drop = FALSE) %>%
     select(Search.Term:Protein.names, Loci, Genes, Organism:Protein.families, -Gene.names...primary..)