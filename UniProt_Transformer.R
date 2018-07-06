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

groupedData <- read.csv("C:\\Users\\iljonas\\Downloads\\Capstone Files\\Source Files\\_Locus Names\\My_List_E-modified.txt", 
                        sep = '\t', stringsAsFactors = FALSE) %>%
     mutate(Gene.names = gsub('; ', ';', Gene.names)) %>%
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
     select(Search.Term:Protein.names, Loci, Genes, Organism:Protein.families, -Gene.names...primary..) %>%
     mutate(Loci = gsub(Loci, '(\\..|[a-z])$', '')) %>%
     mutate(Gene = gsub(Gene, '_', ''))

#master <- read.csv(paste0("C:\\Users\\", Sys.info()["user"], 
#                "\\Documents\\Capstone Files\\GEO-Antimicrobial-Adjunct-Project\\Output_Files\\Exports\\output.tsv.gz"), 
#                sep = '\t') %>%

master <- read.csv("C:\\Users\\iljonas\\Downloads\\Capstone Files\\Source Files\\_Averaged Files\\Combined_Results-Temp.txt", 
                    sep = '\t', stringsAsFactors = FALSE)

uni.master <- merge(newData, master, by.x = 'Search.Term', by.y = 'ï..Name', all.y = TRUE) %>%
     unique %>% 
     mutate(Search.Term = if_else(!is.na(Loci), Loci, if_else(Source == 'Protein' & !is.na(Entry) & Search.Term != Entry
                                                              , Entry, Search.Term)))