library(reshape2)
library(dplyr)
library(tidyr)
library(rvest)

trim <- function(t){gsub("^\\s+|\\s+$", "", t)}

setMeta.study <- function(source){
     temp.meta <- read.table(source, sep = "\t", 
                             comment.char = "", fill = TRUE, stringsAsFactors = FALSE)
     colnames(temp.meta) <- c("Field", "Description")
     temp.meta <- subset(temp.meta, grepl("#", Field))
     temp.meta$Field <- gsub("#","",temp.meta$Field)
     return(temp.meta)
}

setMeta.samples <- function(source){
     temp.meta <- read.table(source, comment.char = "#", sep = '\t', na.strings = 'NA', stringsAsFactors = FALSE)
     row.names(temp.meta) <- temp.meta$V1
     return(temp.meta[-1])
}

setSeries <- function(seriesID, m.samples){
     series.folder <- "Series_Source"
     dir.create(series.folder, showWarnings = FALSE)
     local.path <- paste0(series.folder ,"\\" , seriesID, ".txt.gz")
     
     if(!file.exists(local.path)){
          if(!(grepl("^ftp:\\\\", seriesID))) {
               ftp.path <- paste("ftp://ftp.ncbi.nlm.nih.gov/geo/series", 
                                 gsub("...$", "nnn", seriesID), seriesID, "matrix",
                                 paste0(seriesID, "_series_matrix.txt.gz"), 
                                 sep = "/")
          } else {
               ftp.path = seriesID
          }
          download.file(ftp.path, destfile = local.path)
     }
     
     temp.series <- read.csv(gzfile(local.path), comment.char = "!", sep = '\t', 
                             na.strings = c("", "null"), stringsAsFactors = FALSE)
     temp.series <- temp.series[, order(names(temp.series))]
     temp.series <- select(temp.series, ID_REF, everything())
     m.samples <- m.samples[order(row.names(m.samples)), ]
     if(!(FALSE %in%
          (sapply(names(temp.series)[-1], function(x){trim(x)}) 
               == sapply(row.names(m.samples), function(x){trim(x)})))) {
          names(temp.series)[-1] <- apply(m.samples, 1, FUN = function(x){paste(trim(x), collapse = "|")})
          names(temp.series) <- gsub(" ", "_", names(temp.series))
     }
     else{
          stop("The samples names in the form don't match the sample names in the series data.
               Please check over the entries and rerun the program.")
     }
     return(temp.series)
}

setPlatform <- function(platformID){
     if(grepl("[^.:\\\\]", platformID)) {
          pURL <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=", platformID)
          pText <- read_html(pURL) %>%
               html_nodes("pre") %>%
               html_text
          read.csv(text = pText, comment.char = "#", sep = "\t", quote = "", na.strings = "", stringsAsFactors = FALSE)
     } else {
          read.csv(platformID, comment.char = "#", sep = "\t", quote = "", na.strings = "", stringsAsFactors = FALSE)
     }
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

combine.groups <- function(g, m.study){
     combined.group <- select(g, Expression.ID)
     all.group <- as.data.frame(c(Expression.ID = combined.group, select(g, contains("160_min"))))
     control.group <- select(g, contains("120_min"))
     test.group <- select(g, contains("|60_min"))
     diff.calc <- subset(meta.study, Field == "Diff_Calc", Description, drop = TRUE)
     comp.calc <- subset(meta.study, Field == "Comparison_Calc", Description, drop = TRUE)
     counter <- ncol(all.group)
     for(t in test.group){
          for(c in control.group){
               counter <- counter + 1
               all.group[, counter] <- if_else(diff.calc == "Subtract", t - c, 
                                               if_else(diff.calc == "Log2", log(t / c, 2), 
                                               stop("Error in difference calculation; check JavaScript")))
          }
     }
     
     all.combined <- apply(all.group[, -1], 1, if_else(comp.calc == "Mean", FUN = mean,
                                                      if_else(comp.calc == "Sum", FUN = sum,
                                                              stop("Error in combination calculation; check JavaScript"))) 
                           , na.rm = TRUE)
     combined.group$placeholder <- all.combined
     #new.name <- unique(gsub("^[^|]*\\|", "", names(g)[-1]))
     names(combined.group)[2] <- unique(gsub("\\|Dapt.*", "", names(g)[-1]))
     
     return(all.combined)
}

meta.source <- "C:\\Users\\Isaac\\Documents\\Capstone Files\\GEO-Antimicrobial-Adjunct-Project\\Temp_Source-Rotated.txt"

meta.study <- setMeta.study(meta.source)
meta.samples <- setMeta.samples(meta.source)

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

new.master <- sapply(sep.master, FUN = combine.groups, m.study = meta.study)