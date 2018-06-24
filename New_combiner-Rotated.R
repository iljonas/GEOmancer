library(dplyr)
library(tidyr)
library(rvest)

trim <- function(t){gsub("^\\s+|\\s+$", "", t)}

logDiv <- function(x, y){log(x / y, 2)}

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
     local.path <- paste0(series.folder ,"\\" , toupper(seriesID), ".txt.gz")
     
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
          names(temp.series)[-1] <- apply(m.samples, 1, FUN = function(x){paste(trim(toupper(x)), collapse = "|")})
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
          pURL <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=", toupper(platformID))
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
     #If the column name starts with Test or Neither, assemble all related columns into a data frame (controls are ignored in this operation)
     #To find matches, first make a temp.df to be to hold any non-control columns that match dfname (as well as dfname itself)
     #Then remove all Test/Neither columns from the original data frame, leaving only the controls
     #and remove the group identifier from the beginning of the remaining names
     #, keep only the Expression.ID and any controls (with their group id removed) that match df name when the treatment, dosage, 
     #and exposure time columns are removed
     if(grepl("^(Test|Neither)", dfname)) {
          test.df <- select(df, which(gsub("^[^|]*\\|", "", names(df)) == gsub("^[^|]*\\|", "", dfname)))
          cont.df <- df %>% 
               select(-(grep("^(Test|Neither)", names(df)))) %>%
               gsub("^[^|]*\\|", "", names(df)) %>%
               select(which(names(df) == "Expression.ID" | 
                                 gsub("([^|]*\\|){2}[^|]*$", "", dfname) == gsub("([^|]*\\|){2}[^|]*$", "", names(df))))
          
          #If there are no Test columns in the data frame, you don't need the controls. So these groups are only merged if a Test frame is present
          if(TRUE %in% (grepl("^Test", names(test.df)))) {
               return(merge(test.df, cont.df, by = "Expression.ID"))
          } else {
               return(test.df)
          }
     }
}

combine.groups <- function(g, m.study){
     #initiates combined group with just the expression id. This will be the output data frame
     combined.group <- select(g, Expression.ID)
     
     #initiates the all group, which will contain all individual comparisons
     #starts with the expression id and "Neither" groups, as they don't need to be compared
     all.group <- as.data.frame(c(Expression.ID = combined.group, select(g, startswith("Neither"))))
     
     #separate the control and test groups
     control.group <- select(g, startswith("Control"))
     test.group <- select(g, startswith("Test"))
     
     #get user selected method for calculating difference between each test and control 
     #and method for combining all groups     
     diff.calc <- subset(meta.study, Field == "Diff_Calc", Description, drop = TRUE)
     comp.calc <- meta.study %>% 
          subset(Field == "Comparison_Calc", Description, drop = TRUE) %>%
          tolower
     
     #counter for adding new columns to the all group (after any current columns
     counter <- ncol(all.group) 
     
     #determine diff.calc function based on user input
     if(diff.calc == "Subtraction") {
          diff.calc <- "-"
     } else if(diff.calc == "Log2") {
          diff.calc <- "logDiv"
     }
     
     #loop through each test and compare to each control (via the user selected method)
     for(t in test.group){
          for(c in control.group){
               counter <- counter + 1
               all.group[, counter] <- do.call(diff.calc, x, y)
          }
     }
     
     #combine all groups (except for expression id) via the user selected method
     all.combined <- apply(all.group[, -1], 1, FUN = do.call(comp.calc), na.rm = TRUE)
     
     #normalize the data and add to the combined group
     combined.group$placeholder <- all.combined / max(all.combined)
     
     #with the control groups and expression id removed, all names should be the same once the group type is removed
     #use this common name for the data title
     names(combined.group)[2] <- names(g)[-1] %>%
          select(-(startswith("Control"))) %>%
          gsub("^[^|]*\\|", "") %>%
          unique
     
     return(all.combined)
}

meta.source <- "C:\\Users\\iljonas\\Documents\\Capstone\\GEO-Antimicrobial-Adjunct-Project\\Temp_Source-Rotated.txt"

meta.study <- setMeta.study(meta.source)
meta.samples <- setMeta.samples(meta.source)

series.id <- subset(meta.study, Field == "Study", Description, drop = TRUE)
series <- setSeries(series.id, meta.samples)

platform.id <- subset(meta.study, Field == "Platform", Description, drop = TRUE)
platform <- setPlatform(platform.id)

subdelim <-  subset(meta.study, Field == "Sub_Delim", Description, drop = TRUE)
if(subdelim == "") {subdelim <- ","}

new.platform <- platform %>%
     select(grep("^ID$|^Strain|^ORF|^SPOT_ID|^PT_ACC$", names(platform))) %>%
     gather(key, value, -ID, na.rm = TRUE) %>%
     rowwise %>%
     mutate(value = gsub("::(?(?<=E).|[^+-])*[+-]", "", value, perl = TRUE)) %>%
     mutate(Expression.ID = strsplit(as.character(value), subdelim)) %>%
     unnest(Expression.ID) %>%
     mutate(Source = if_else(key == "PT_ACC", "Protein", "Gene")) %>%
     select(-c(value, key)) %>%
     unique

master <- merge(new.platform, series, by.x = "ID", by.y = "ID_REF", all = TRUE) %>%
     select(-ID)

outside.controls <- subset(meta.study, Field == "Study_Controls", Description, drop = TRUE)
if(outside.controls != ""){
     master <- newControls(outside.controls, master)
}

#
i <- 2
sep.master <- vector()
repeat{
     group <- findMatches(names(master)[i], master)
     group.tests <- select(group, grep("^(Test|Neither)", names(group)))
     master <- select(master, -(names(group.tests)))
     if(!(TRUE %in% (grepl("^Control", master[-1])))){
          break
     }
     i <- i + 1
}
#sep.master<- sapply(colnames(master), findMatches, df = master)

comb.master <- sapply(sep.master, FUN = combine.groups, m.study = meta.study)