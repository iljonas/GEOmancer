suppressMessages(library(dplyr))

read.files <- function(pathway, min.value){
     all.results <- read.csv(pathway, sep = '\t')
     all.results$score.perc <- all.results$score / max(all.results$score)
     trans.value <- (100 - min.value) / 100
     unlist(lapply(all.results, is.character))
     reduced.results <- all.results %>%
          filter(score.perc >= trans.value) %>%
          select(-(col_num:stddev_pos), -(score.perc))
}

user.input <- c('Raw_Master_File3', '13:94', 'mean', 'output', '5')

folder.path <- file.path('C:', 'Users', 'Isaac', 'Documents', 'Capstone Files', 
                         'GEO-Antimicrobial-Adjunct-Project')

#upload the master file and remove all of the study columns
master <- file.path(folder.path, user.input[1]) %>%
     paste0('.tsv.gz') %>%
     read.csv(sep = '\t', stringsAsFactors = FALSE) %>%
     select(grep('[^|]', colnames))

folder.files <- list.files(folder.path, full.names = TRUE)
#all.tables <- lapply(list.files, function(x) {read.files(x, 5)})
top.results <- data.frame(Expression.ID = character())

repeat{
     end.of.file.name <- tail(strsplit(folder.path[1], '-', fixed = TRUE))
     sub.path <- folder.files[grep(end.of.file.name, folder.files)]
     data.list <- list(master,
          unlist(lapply(sub.path, 
                        function(x){read.files(x, as.numeric(user.input[5]))})))
     combined.data <- Reduce(function(x,y) {merge(x,y)}, data.list)
     top.results <- merge(top.results, combined.data, all = TRUE)
     
     folder.files <- folder.files[!(grep(end.of.file.name, folder.files))]
     
     if (length(folder.files) == 0){
          break
     }
}
