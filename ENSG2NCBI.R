library(org.Hs.eg.db)

formatToEntrez <- function(inFile,outFile){
  
  ens2egs <- as.list(org.Hs.egENSEMBL2EG)
  df <- read.table("salmon.merged.gene_counts.tsv",head=T,sep="\t",quote="")
  ens <- df[,"gene_id"]
  print(paste("number of ensgs:",length(ens)))
  print(paste("number of unique ensgs:",length(unique(ens))))
  
  # mapping ensg to eg(s) (can be more than one)
  eg <- sapply(ens,function(x){r <- ens2egs[[x]];if(is.null(r)){r <- "NA"};r})
  # below, most ls are 1 with a few greater than 1
  ls <- lengths(eg)
  df <- cbind(ls,df)
  # expand each row i ls[i] times
  df <- df[rep(row.names(df), df$ls),]
  print(paste("number of rows after expansion:",length(df$ls)))
  # assigns 1 eg per row
  eg <- unlist(eg)
  df <- cbind(eg,df)
  # finds rows (ensgs) which are not mapped to entrez 
  toremove <- which(eg == "NA")
  # removes unmapped rows
  df <- df[-toremove,]
  print(paste("number of rows after removing unmapped ensgs:",length(df$eg)))
  print(paste("number of unique egs:",length(unique(df$eg))))
  toremove <- which(colnames(df) == "gene_id" | colnames(df) == "ls")
  df <- df[,-toremove]
  df <- aggregate(df[,2:length(df)],by=list(df$eg),FUN=mean)
  colnames(df)[1] <- "entrez"
  # saving
  write.table(df,outFile,sep="\t",quote=F,row.names=F)
}

formatToEntrez("salmon.merged.gene_counts.tsv","EntrezData.txt")