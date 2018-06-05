library(reshape2)
library(dplyr)
library(tidyr)
library(rvest)

setMeta.study <- function(source){
     temp.meta <- read.table(source, sep = "\t", 
                             comment.char = "", fill = TRUE, stringsAsFactors = FALSE)
     colnames(temp.meta) <- c("Field", "Description")
     temp.meta <- subset(temp.meta, grepl("#", Field))
     temp.meta$Field <- gsub("#","",temp.meta$Field)
     return(temp.meta)
}

setMeta.samples <- function(source, m.study){
     temp.meta <- read.csv(source, comment.char = "#", sep = '\t', na.strings = 'NA', stringsAsFactors = FALSE)
     row.names(temp.meta) <- subset(m.study, Field == "Col_Meaning", Description, drop = TRUE) %>%
          strsplit(",") %>%
          unlist
     return(temp.meta)
}

setSeries <- function(seriesID, m.samples){
     series.folder <- "Series_Source"
     dir.create(series.folder, showWarnings = FALSE)
     local.path <- paste0(series.folder ,"\\" , seriesID, ".txt.gz")
     
     if(!file.exists(local.path)){
          ftp.path <- paste("ftp://ftp.ncbi.nlm.nih.gov/geo/series", 
                            gsub("...$", "nnn", seriesID), seriesID, "matrix",
                            paste0(seriesID, "_series_matrix.txt.gz"), 
                            sep = "/")
          download.file(ftp.path, destfile = local.path)
     }
     
     temp.series <- read.csv(gzfile(local.path), comment.char = "!", sep = '\t', na.strings = c("", "null"), stringsAsFactors = FALSE)
     colnames(temp.series)[-1] <- sapply(m.samples,  
                                         FUN = function(x){paste(x[!is.na(x)], collapse = "|")})
     colnames(temp.series) <- gsub(" ", "_", colnames(temp.series))
     return(temp.series)
}

setPlatform <- function(platformID){
     pURL <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=", platformID)
     pText <- read_html(pURL) %>%
          html_nodes("pre") %>%
          html_text
     read.csv(text = pText, comment.char = "#", sep = "\t", quote = "", na.strings = "")
}

newControls <- function(c.source, m){
     c.study <- setMeta.study(c.source)
     c.samples <- setMeta.samples(c.source, c.study)
     
     c_series.id <- subset(c.study, Field == "Study", Description, drop = TRUE)
     c_series <- setSeries(series.id, c.samples)
     
     c_series.controls <- select(sec_series, Expression.ID, starts_with("Test"))
     merge(m, c_series.control, by = "Expression.ID")
}

findMatches <- function(dfname, df){
     if(grepl("^Test|Neither", dfname)) {
          sub.df <- df[,names(df) == "Expression.ID" | gsub("^[^|]*\\|", "", dfname) == gsub("^[^|]*\\|", "", names(df))]
          
          if(TRUE %in% grepl("^Test|Neither", names(sub.df))){ #change test to control
               return(sub.df)
          }
     }
}

combine.groups <- function(g){
     combined.group <- select(g, Expression.ID)
     all.group <- as.data.frame(c(Expression.ID = combined.group, select(g, contains("160_min"))))
     control.group <- select(g, contains("120_min"))
     test.group <- select(g, contains("|60_min"))
     counter <- ncol(all.group)
     for(t in test.group){
          for(c in control.group){
               counter <- counter + 1
               all.group[, counter] <- t - c
          }
     }
     
     all.combined <- apply(all.group[, -1], 1, FUN = mean, na.rm = TRUE)
     combined.group$placeholder <- all.combined
     #new.name <- unique(gsub("^[^|]*\\|", "", names(g)[-1]))
     names(combined.group)[2] <- unique(gsub("\\|Dapt.*", "", names(g)[-1]))
     
     return(all.combined)
}

meta.source <- "C:\\Users\\Isaac\\Documents\\Capstone Files\\Temp_Source.txt"

meta.study <- setMeta.study(meta.source)
meta.samples <- setMeta.samples(meta.source, meta.study)

series.id <- subset(meta.study, Field == "Study", Description, drop = TRUE)
series <- setSeries(series.id, meta.samples)

platform.id <- subset(meta.study, Field == "Platform", Description, drop = TRUE)
platform <- setPlatform(platform.id)

subdelim <-  subset(meta.study, Field == "Sub_Delim", Description, drop = TRUE)
if(subdelim == "") {subdelim <- ","}

new.platform <- platform %>%
     select(which(grepl("^ID$|^Strain|^ORF|^SPOT_ID|^PT_ACC$", names(platform)))) %>%
     melt(id.vars = "ID", na.rm = TRUE) %>%
     rowwise %>%
     mutate(value = gsub("::(?(?<=E).|[^+-])*[+-]", "", value, perl = TRUE)) %>%
     mutate(Expression.ID = strsplit(as.character(value), subdelim)) %>%
     unnest(Expression.ID) %>%
     select(-c(value, variable)) %>%
     unique

master <- merge(new.platform, series, by.x = "ID", by.y = "ID_REF", all = TRUE) %>%
     select(-(ID))

outside.controls <- subset(meta.study, Field == "Study_Controls", Description, drop = TRUE)
if(outside.controls != ""){
     master <- newControls(outside.controls, master)
}

sep.master<- sapply(colnames(master), findMatches, df = master)

new.master <- sapply(sep.master, FUN = combine.groups)