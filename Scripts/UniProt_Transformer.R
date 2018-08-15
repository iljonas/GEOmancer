suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

geneLocSep <- function(gList){
     gList <- list(gList)
     geneLocs <- as.integer(sapply(gList, function(x) {grep("^[a-z]",x)}))
     genes <- unlist(gList)[geneLocs]
     loci <- ifelse(is.na(geneLocs), unlist(gList), unlist(gList)[-geneLocs])
     return(list(genes,loci))
}

user <- as.character(Sys.info()["user"])
user.input <- commandArgs(trailingOnly = TRUE)

primary.pathway <- file.path('C:', 'Users', user, 'Documents', 'GEOmancer', 'Output_Files')
uniprot.pathway <- file.path(primary.pathway, 'UniProt_Downloads', 'uniProt-') %>%
     paste0(user.input[1], '.tsv.gz')
master.pathway <- file.path(primary.pathway, 'Exports', user.input[1]) %>%
     paste0('.tsv.gz')

groupedData <- read.csv(gzfile(uniprot.pathway), sep = '\t', stringsAsFactors = FALSE) %>%
     unnest(Expression.ID = strsplit(as.character(Expression.ID), ',')) %>%
     mutate(Gene.names = gsub('; ', ';', Gene.names)) %>%
     unnest(Sep.names = strsplit(as.character(Gene.names), ';')) %>%
     mutate(Sep.groups = strsplit(Sep.names, ' '))

split.names <- sapply(groupedData$Sep.groups, geneLocSep) %>%
     t %>%
     as.data.frame
names(split.names) <- c("Genes", "Loci")

newData <- data.frame(groupedData, split.names) %>%
     unnest(Loci, .drop = FALSE) %>%
     unnest(Genes, .drop = FALSE) %>%
     select(Expression.ID, Entry:Protein.names, Loci, Genes, Protein.families:Sequence) %>%
     mutate(Loci = gsub('(\\..|[a-z])$', '', Loci)) %>%
     mutate(Gene = gsub('_', '', Genes))

master <- read.csv(gzfile(master.pathway), sep = '\t', stringsAsFactors = FALSE)

uni.master <- merge(newData, master, by = 'Expression.ID', all.y = TRUE) %>%
     unique %>% 
     mutate(Expression.ID = if_else(!is.na(Loci), Loci, if_else(Source == 'Protein' & !is.na(Entry) & Expression.ID != Entry
                                                              , Entry, Expression.ID)))

write.table(uni.master, file = gzfile(paste0(file.path(primary.pathway, 'Exports', 'master-'), user.input[1], '.tsv.gz')),
            sep = '\t', na = 'NA', row.names = FALSE)