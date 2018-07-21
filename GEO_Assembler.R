library(dplyr)
library(tidyr)
library(rvest)

trim <- function(t){gsub("^\\s+|\\s+$", "", t)}

logDiv <- function(x, y){log(x / y, 2)}

folderPath <- function(){
     userID <- Sys.info()["user"]
     folder <- paste0("C:\\Users\\", userID, 
                           "\\Documents\\Capstone Files\\GEO-Antimicrobial-Adjunct-Project\\GEO_Source_Files\\")
}

comparisonCalc <- function(calc.source){
     calc <- calc.source %>% 
          subset(Field == "Comparison Calculation", Description, drop = TRUE) %>%
          tolower
}

clipEnds <- function(cName, exp1, exp2){
     cName <- gsub(exp1, '', cName)
     cName <- gsub(exp2, '', cName)
     return(cName)
}

setMeta.study <- function(source){
     temp.meta <- read.table(source, sep = "\t", 
                             comment.char = "", fill = TRUE, stringsAsFactors = FALSE)
     colnames(temp.meta) <- c("Field", "Description")
     temp.meta <- subset(temp.meta, grepl("#", Field))
     temp.meta$Field <- gsub("#","",temp.meta$Field)
     return(temp.meta)
}

setMeta.samples <- function(source){
     temp.meta <- read.csv(source, comment.char = "#", sep = '\t', na.strings = c('', 'NA'), stringsAsFactors = FALSE)
     row.names(temp.meta) <- temp.meta$Sample_ID
     return(temp.meta[-1])
}

setSeries <- function(seriesID, m.samples){
     series.folder <- paste0(folderPath(), "\\GEO_Series")
     series.path <- paste0(series.folder ,"\\" , toupper(seriesID), ".txt.gz")
     
     if(!file.exists(series.path)){
          ftp.path <- paste("ftp://ftp.ncbi.nlm.nih.gov/geo/series", 
                            gsub("...$", "nnn", seriesID), seriesID, "matrix",
                            paste0(seriesID, "_series_matrix.txt.gz"),
                            sep = "/")
          download.file(ftp.path, destfile = series.path)
     }
     
     temp.series <- read.csv(gzfile(series.path), comment.char = "!", sep = '\t', 
                             na.strings = c("", "null"), stringsAsFactors = FALSE)
     temp.series <- temp.series[, order(names(temp.series))]
     temp.series <- select(temp.series, ID_REF, everything())
     m.samples <- m.samples[order(row.names(m.samples)), ]

     if(!(FALSE %in%
          (sapply(names(temp.series)[-1], function(x){trim(x)}) 
           == sapply(row.names(m.samples), function(x){trim(x)})))) {
          names(temp.series)[-1] <- apply(m.samples, 1, FUN = function(x){paste(trim(toupper(x)), collapse = "|")})
          
          sep.num <- ncol(m.samples) - 1
          for (i in 2:ncol(temp.series)){
               if(nchar(names(temp.series)[i]) - nchar(gsub('\\|', '', names(temp.series)[i])) == sep.num){
                    cols <- which(names(temp.series)[i] == names(temp.series))
                    names(temp.series)[cols] <- paste(names(temp.series)[i], seq_along(cols), sep = '|')
               }
          }
          names(temp.series) <- gsub(" ", "_", names(temp.series))
     }
     else{
          stop("The samples names in the form don't match the sample names in the series data.
               Please check over the entries and rerun the program.")
     }
     
     return(temp.series)
}

setPlatform <- function(platformID){
     platform.folder <- paste0(folderPath(), "\\GEO_Platforms")
     platform.path <- paste0(platform.folder ,"\\" , toupper(platformID), ".txt.gz")
     
     if(!file.exists(platform.path)) {
          pURL <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=", toupper(platformID))
          pText <- read_html(pURL) %>%
               html_nodes("pre") %>%
               html_text
          file.create(platform.path)
          writeLines(pText, con = gzfile(platform.path))
          closeAllConnections()
     }
     
     read.csv(gzfile(platform.path), comment.char = "#", sep = "\t", quote = "", na.strings = "", stringsAsFactors = FALSE)
}

adjustPlatform <- function(platDF, colnums, delim){
     base.plat <- eval(parse(text = paste0("select(platDF, 1, ", colnums, ")")))
     new.plat <- base.plat %>%
          gather(key, value, -ID, na.rm = TRUE) %>%
          rowwise %>%
          mutate(value = gsub("::(?(?<=E).|[^+-])*[+-]", "", value, perl = TRUE)) %>%
          mutate(Expression.ID = strsplit(as.character(value), delim)) %>%
          unnest(Expression.ID) %>%
          mutate(Source = if_else(key == "PT_ACC", "Protein", "Gene")) %>%
          select(-c(value, key)) %>%
          unique
}

assembleSeries <- function(source, is.control = FALSE, ignore.field = ''){
     meta.source <- paste0(folderPath(), "GEO_Forms\\",
                           source, ".tsv")
     
     meta.study <- setMeta.study(meta.source)
     meta.samples <- select(setMeta.samples(meta.source), ifelse(ignore.field != '', -one_of(ignore.field), everything()))
     
     series.id <- subset(meta.study, Field == "Study ID", Description, drop = TRUE)

     series <- setSeries(series.id, meta.samples)
     
     if(is.control == TRUE){
          series <- select(series, ID_REF, starts_with("CONTROL"))
     }
     
     platform.id <- subset(meta.study, Field == "Platform ID", Description, drop = TRUE)
     platform <- setPlatform(platform.id)
     
     col.numbers <- subset(meta.study, Field == "Platform Data Column Numbers", Description, drop = TRUE)
     
     subdelim <- subset(meta.study, Field == "In-column Delimiter", Description, drop = TRUE)
     if(subdelim == "") {subdelim <- ","}
     
     adj.platform <- adjustPlatform(platform, col.numbers, subdelim)

     comp.calc <- comparisonCalc(meta.study)
     adj.series <- merge(adj.platform, series, by.x = "ID", by.y = "ID_REF", all = TRUE) %>%
          select(-ID) %>%
          group_by(Expression.ID, Source) %>%
          summarize_all(eval(parse(text = comp.calc)), na.rm = TRUE)
     adj.series[is.na(adj.series)] <- NA
     
     return(list(meta.study, adj.series))
}

findMatches <- function(dfname, df){
     #If the column name starts with Test or Neither, assemble all related columns into a data frame (controls are ignored in this operation)
     #To find matches, first make a temp.df to be to hold any non-control columns that match dfname (as well as dfname itself)
     #Then remove all Test/Neither columns from the original data frame, leaving only the controls
     #and remove the group identifier from the beginning of the remaining names
     #, keep only the Expression.ID and any controls (with their group id removed) that match df name when the treatment, dosage, 
     #and exposure time columns are removed
     if(grepl("^(TEST|NEITHER)", dfname)) {
          test.df <- select(df, Expression.ID, Source, which(unlist(lapply(names(df), FUN = clipEnds, exp1 = "^[^|]*\\|", exp2 = "\\|[^|]*$"))
                                                     == clipEnds(dfname, exp1 = "^[^|]*\\|", exp2 = "\\|[^|]*$")))
                                                     #gsub("^[^|]*\\|", "", names(df)) == gsub("^[^|]*\\|", "", dfname)))
          cont.df <- df %>% 
               select(-(grep("^(TEST|NEITHER)", names(df)))) %>%
               #setNames(gsub("^[^|]*\\|", "", names(.))) %>%
               select(Expression.ID, Source, 
                      which(unlist(lapply(names(.), FUN = clipEnds, exp1 = "^[^|]*\\|", exp2 = "([^|]*\\|){3}[^|]*$"))
                                 == clipEnds(dfname, "^[^|]*\\|", "([^|]*\\|){3}[^|]*$")))
                                 #gsub("([^|]*\\|){3}[^|]*$", "", dfname) == gsub("([^|]*\\|){3}[^|]*$", "", names(.))))

          #If there are no Test columns in the data frame, you don't need the controls. So these groups are only merged if a Test frame is present
          if(TRUE %in% (grepl("^TEST", names(test.df)))) {
               return(merge(test.df, cont.df, by = c("Expression.ID", "Source")))
          } else {
               return(test.df)
          }
     }
}

combine.groups <- function(g, m.study){
     #initiates combined group with just the expression id and source. This will be the output data frame
     combined.group <- select(g, Expression.ID, Source)

     #initiates the all group, which will contain all individual comparisons
     #starts with the expression id and "Neither" groups, as they don't need to be compared
     all.group <- select(g, Expression.ID, Source, starts_with("NEITHER"))

     #separate the control and test groups
     control.group <- select(g, starts_with("CONTROL"))
     test.group <- select(g, starts_with("TEST"))

     #get user selected method for calculating difference between each test and control 
     #and method for combining all groups     
     diff.calc <- subset(m.study, Field == "Difference Calculation", Description, drop = TRUE)
     comp.calc <- comparisonCalc(m.study)

     #determine diff.calc function based on user input
     if(diff.calc == "Subtraction") {
          diff.calc <- "-"
     } else if(diff.calc == "Log2") {
          diff.calc <- "logDiv"
     }

     #counter for adding new columns to the all group (after any current columns
     counter <- ncol(all.group)
     
     #loop through each test and compare to each control (via the user selected method)
     for(t in test.group){
          for(c in control.group){
               counter <- counter + 1
               all.group[, counter] <- do.call(diff.calc, list(t, c))
          }
     }

     #combine all groups (except for expression id) via the user selected method
     all.combined <- eval(parse(text = paste0('apply(all.group, 1,', comp.calc, ', na.rm = TRUE)')))

     #normalize the data and add to the combined group
     combined.group$placeholder <- all.combined / max(all.combined, na.rm = TRUE)

     #with the control groups and expression id removed, all names should be the same once the group type is removed
     #use this common name for the data title
     names(combined.group)[3] <- names(g)[-(1:2)] %>%
          .[-grep("^CONTROL", .)] %>%
          clipEnds("^[^|]*\\|", "\\|[^|]*$") %>%
          unique

     return(combined.group)
}

user.input <- commandArgs(trailingOnly = TRUE)

source.name <- user.input[1]
exclude <- user.input[2]

master <- assembleSeries(source.name, ignore.field = exclude)

outside.controls <- subset(master[[1]], Field == "Outside study controls", Description, drop = TRUE)
if(outside.controls != ""){
     master[[2]] <- merge(master[[2]], assembleSeries(outside.controls, TRUE, exclude)[[2]])
}

#
i <- 3
comb.master <- select(master[[2]], Expression.ID, Source)

cat("Combined columns created:", '\n')
repeat{
     if(!grepl("^CONTROL", names(master[[2]])[i])){
          sep.master <- findMatches(names(master[[2]])[i], master[[2]])
          master[[2]] <- select(master[[2]], -one_of(names(sep.master)[grep("^(TEST|NEITHER)", names(sep.master))]))
          i <- i - 1
          comb.master <- merge(comb.master, combine.groups(sep.master, master[[1]]), by = c("Expression.ID", "Source"), all = TRUE)
          cat(names(comb.master)[ncol(comb.master)], '\n')
     }

     if(!(FALSE %in% (grepl("^CONTROL", names(master[[2]])[-(1:2)])))
        | i == ncol(master[[2]])){
          break
     }
     i <- i + 1
}

write.table(comb.master, 
          file = gzfile(paste0("C:\\Users\\", Sys.info()["user"], 
                 "\\Documents\\Capstone Files\\GEO-Antimicrobial-Adjunct-Project\\Output_Files\\Exports\\",
                 "output", ".tsv.gz")),
          sep = '\t', na = 'NA', row.names = FALSE)