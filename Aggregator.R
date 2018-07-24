#Isaac Jonas
#Analyzer-Comb_File.R
#This program exports analysis views of all prompted datasets. Each user provided dataset is separated analyzed as provided (analyses used
#described below) and analyzed with only its bacteriostatic columns, bactericidal columns, MRSA columns, and MSSA columns. Each of these is
#exported to its own generated file

library(dplyr)

folderPath <- function(){
     userID <- Sys.info()["user"]
     folder <- paste0("C:\\Users\\", userID, 
                      "\\Documents\\Capstone Files\\GEO-Antimicrobial-Adjunct-Project\\GEO_Source_Files\\")
}

#get user input for import file and read data
getCSV <- function(file.name){
     fullPath <- paste(file.path(folderPath(), file.name, fsep = "\\"), ".tsv.gz", sep = "")
     return(read.csv(gzfile(fullPath), header = TRUE, sep = "\t", fill = TRUE, comment.char = "#"))
}

#creates a separate dataset from base of only positive values then analyzes base dataset for its used of columns, number of null columns, and
#percent of columns that are positive. Both sets are analyzed for their average, quartiles (these are returned as a vector, so they are be
#exported one at a time), and standard deviation. The score is then added to the dataset based on the calculated values
setCSV_analyze <- function(combReduc, fullPath){
     combReduc_pos <- combReduc
     combReduc_pos[combReduc_pos < 0] <- NA
     
     #row number is removed from ID before being exported by removing everything after the ';', which separates the ID name and the row number
     combAgg <- data.frame(Id = sapply(strsplit(row.names(combReduc),";"), '[',1), col_num = apply(!is.na(combReduc), 1, FUN=sum), 
                           col_null = apply(is.na(combReduc), 1, FUN=sum), 
                           perc_pos = apply(!is.na(combReduc_pos), 1, na.rm = TRUE, FUN=sum) / 
                                apply(!is.na(combReduc), 1, na.rm = TRUE, FUN=sum), 
                           avg_all = apply(combReduc, 1, na.rm = TRUE, FUN=mean), 
                           avg_pos = apply(combReduc_pos, 1, na.rm = TRUE, FUN=mean), 
                           quar_0 = apply(combReduc, 1, na.rm = TRUE, FUN=quantile)[1,], 
                           quar_25 = apply(combReduc, 1, na.rm = TRUE, FUN=quantile)[2,], 
                           quar_50 = apply(combReduc, 1, na.rm = TRUE, FUN=quantile)[3,], 
                           quar_75 = apply(combReduc, 1, na.rm = TRUE, FUN=quantile)[4,], 
                           quar_100 = apply(combReduc, 1, na.rm = TRUE, FUN=quantile)[5,], 
                           pos_quar_0 = apply(combReduc_pos, 1, na.rm = TRUE, FUN=quantile)[1,], 
                           pos_quar_25 = apply(combReduc_pos, 1, na.rm = TRUE, FUN=quantile)[2,], 
                           pos_quar_50 = apply(combReduc_pos, 1, na.rm = TRUE, FUN=quantile)[3,], 
                           pos_quar_75 = apply(combReduc_pos, 1, na.rm = TRUE, FUN=quantile)[4,], 
                           pos_quar_100 = apply(combReduc_pos, 1, na.rm = TRUE, FUN=quantile)[5,], 
                           stddev_all = apply(!is.na(combReduc), 1, FUN=sd), 
                           stddev_pos = apply(!is.na(combReduc_pos), 1, FUN=sd))
     
     combAgg_scored <- data.frame(combAgg, score = (combAgg$col_num * combAgg$perc_pos * combAgg$avg_pos) / combAgg$stddev_pos)
     
     return(write.table(gzfile(combAgg_scored), fullPath, quote = FALSE, sep = "\t", row.names = FALSE))
}

#prompt user for export file name and separates data columns for various views of analyzes (methods of action and resistance types)
setCSV <- function(cFile){
     fileName <- readline(prompt = "Enter filename: ")
     #some row genes are not unique by R's non-case sensitive standards and cannot be used as rows. This combines genes with the row number so they become unique
     row.names(cFile) <- paste(cFile[,1],row.names(cFile),sep=";")
     #first row, containing IDs, cannot be included because it cannot be aggregated. Not included with any of these analyses and is saved as row name
     #send full dataset and write to base desired file name
     setCSV_analyze(cFile[,-1], paste(file.path(".\\", fileName), ".tsv.gz", sep = ""))
}

#run program
main <- function(){
     user.input <- commandArgs(trailingOnly = TRUE)
     #read the tsv file based on file name provided by user, then select the grouping column and data columns
     source.file <- getCSV(user.input[1]) %>%
          select(user.input[2], user.input[3])
     agg.file <- eval(parse(text = paste0('aggregate(source.file[2:ncol(source.file)], by = list(IDs = source.file$Expression.ID), FUN = ',
                        user.input[4], ', na.rm=TRUE)')))
     setCSV(agg.file)
}

main()