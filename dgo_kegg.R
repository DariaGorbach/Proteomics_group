#KEGG ANNOTATION
############## KeGG_ID ###################

library(KEGGREST)
library(dplyr)
#dgo_test <- pea_table
#ara2_all <- dgo_test

test_blast$PATHWAY <- NA
test_blast$BRITE <- NA

test_blast$MODULE <- NA
test_blast$DEFINITION <- NA

#ara2_all <- subset(ara2_all, !(ara2_all$kegg_id == '')) #subset all proteins which have KEGG-ID, because keggGet function doesn't work with NA
#later we will combine this kegg contating dataset with the dataset without kegg ids

for (i in (1:NROW(pea_table))){
  test <- as.vector(str_split_fixed(test_blast$KO[i],";",Inf))
  
  pathways <- vector()
  for (j in (1:length(test))){
    val <- test[j]
    df <- keggGet(val)[[1]]
    
    if (is.null(df$PATHWAY)==F) {
      test2 <- paste(as.vector(df$PATHWAY),collapse=";")
      pathways <- c(pathways,test2)
    } else {
      next
    }
  }
  
  pathways <- unique(as.vector(str_split_fixed(pathways,";",Inf)))
  test_blast$PATHWAY[i] <- as.character(paste(as.vector(pathways),collapse=";"))
  
  
  brites <- vector()
  for (j in (1:length(test))){
    val <- test[j]
    df <- keggGet(val)[[1]]
    
    if (is.null(df$BRITE)==F) {
      test3 <- paste(as.vector(df$BRITE),collapse=";")
      brites <- c(brites,test3)
    } else {
      next
    }
  }
  
  brites <- unique(as.vector(str_split_fixed(brites,";",Inf)))
  test_blast$BRITE[i] <- as.character(paste(as.vector(brites),collapse=";"))
  
  
  modules <- vector()
  for (j in (1:length(test))){
    val <- test[j]
    df <- keggGet(val)[[1]]
    
    if (is.null(df$MODULE)==F) {
      test4 <- paste(as.vector(df$MODULE),collapse=";")
      modules <- c(modules,test4)
    } else {
      next
    }
  }
  modules <- unique(as.vector(str_split_fixed(modules,";",Inf)))
  test_blast$MODULE[i] <- as.character(paste(as.vector(modules),collapse=";"))
  
  definition <- vector()
  for (j in (1:length(test))){
    val <- test[j]
    df <- keggGet(val)[[1]]
    
    if (is.null(df$DEFINITION)==F) {
      test5 <- paste(as.vector(df$DEFINITION),collapse=";")
      definition <- c(definition,test5)
    } else {
      next
    }
  }
  definition <- unique(as.vector(str_split_fixed(definition,";",Inf)))
  test_blast$DEFINITION[i] <- as.character(paste(as.vector(definition),collapse=";"))
  
  
  print(i)
}



###################################
library(pbapply)
test_blast$BRITE <- pbapply::pbsapply(test_blast$BRITE,function(brite_full){
  
  patt <- 
    c(seq(09101,09112,1) %>% gsub('^','0',.),
      seq(09121,09124,1) %>% gsub('^','0',.),
      seq(09131,09133,1) %>% gsub('^','0',.),
      seq(09141,09145,1) %>% gsub('^','0',.)
    ) %>% as.character(.) %>% paste(.,collapse = "|")
  
  vec <- str_split_fixed(string = brite_full,pattern = ';',n = Inf) %>% 
    table(.) %>% 
    as.data.frame(.) %>% 
    dplyr::rename('brite_d'=1) %>% 
    dplyr::select(1)
  
  vec <- subset(vec,grepl(patt,vec$brite_d))
  vec <- paste(vec$brite_d %>% str_trim(., "left"),collapse = ';')
  return(vec)
})
test_blast$BRITE %>% gsub('^$','not_assigned',.) -> test_blast$BRITE
