library("isobar")
library("plyr")

#Write a CSV table for each unique spectra file found in the PSM report
writeCSVTables <- function(files, psms){
  
  for(fname in files){
    
    name <- tools::file_path_sans_ext(fname)
    
    write.table(psms[grep(name, psms$file), c("accession", "sequence", "peptide",
                                              "modif", "fixed", "file", "spectrum")],
                paste(name,".corr.csv", sep =""), sep ="\t", row.names = FALSE, quote = FALSE)
  }
}

#Process the PSM report and split it into tables for each MGF file
processPSMReport <- function(fileName){
  
  #Load file
  psms <- read.csv2(fileName, sep = "\t", fill = TRUE)
  
  #Rename columns for more intuitive naming
  names(psms)[c(2:3, 7:11)] <- c("accession", "sequence", "peptide", "modif", 
                                 "fixed", "file", "spectrum")
  
  #Make corrections so isobar can process
  psms$spectrum <- unlist(sub("File:","File:\"", psms$spectrum))
  psms$spectrum <- unlist(sub("NativeID:", "NativeID:\"", psms$spectrum))
  psms$spectrum <- unlist(sub(",", "\",", psms$spectrum))
  psms$spectrum <- paste(psms$spectrum, "\"", sep="")
  
  #Filter for most confident protein identifications
  psms <- psms[psms$Validation == "Confident",]
  
  
  files<- unique(psms$file)
  writeCSVTables(files, psms)
  
  results <- vector("list", length(files))
  
  #Separate the report into csv tables belonging to each MGF file
  for(i in seq(length(files)))
  {
   fname <- tools::file_path_sans_ext(files[i])
   ib <- readIBSpectra(type = "TMT6plexSpectra",id.file = paste(fname,".corr.csv", sep = ""), 
                       peaklist.file = paste("spectra_files/",fname,".mgf",sep =""),decode.titles=TRUE)
   
   results[i] <- ib #More efficient to build result list like this instead of concatenating with it's previous value
  }
  
  return (results)
}

#Correct isotope impurities, normalize data, eliminate any aditional noise
#This is important for correcting any errors that might have happened during the pipeline process
preprocessSpectras <- function(spectras){
  
    result <- vector("list", length(spectras))

    for(i in seq(length(spectras))){
    s <- spectras[[i]]
    
    s <- correctIsotopeImpurities(s)
    s <- normalize(s)
    s <- subtractAdditiveNoise(s)
    
    result[i] <- s  # Throws warning of deprecated implicit embedding
    }
    
    return (result)
}

#Calculate protein ratios for each IBSpectra object
calculateProteinRatios <- function(spectras){
  
  classLabels <- c("126", "127", "128", "129", "130", "131")
  cMethod <- "versus.channel"
  
  remove(result)
  result <- vector("list", length(spectras))
  for(i in seq(length(spectras))){
    
    noiseModel <- NoiseModel(s)
    pr <- proteinRatios(spectras[[i]], noiseModel, cl = classLabels, combn.method = cMethod)
    result[[i]] <- pr
  }
  
  return (result)
}


#Save Protein ratios tables so don't have to rerun the comands again and again
writeProteinRatiosCSV <- function(proteinRatios){
  for(i in seq(length(proteinRatios))){
    obj <- proteinRatios[i]
    write.table(x = proteinRatios[[i]], file = paste("Protein Ratios", i,".csv",sep = ""),sep = "\t",
                col.names = TRUE, row.names =TRUE )
  }
}

#File name
#ident_file <- "Malaria_Natali_TMT6_MGMS_1_Default_PSM_Report.txt"

#Process the PSM report and generate the IBSpectra objects
#spectras <- processPSMReport(ident_file)

#Correct spectras
#correctedSpectras <- preprocessSpectras(spectras)

#Calculate protein quantitative values
#proteinRatios <- calculateProteinRatios(correctedSpectras)

#Save the protein ratios in files so don't have to always run the PSM process methods
#writeProteinRatiosCSV(proteinRatios)

pr1 <- read.csv2("Protein Ratios1.csv", sep= "\t")
pr2 <- read.csv2("Protein Ratios2.csv", sep= "\t")


#Merge the two data samples on the amino acid sequences and the tag the 126 is compared to
merged<- merge(pr1,pr2,by = c("ac", "r2"))

#Statisticall analysis

## Create the table with the data of interest
table<- data.frame(merged$ac,merged$r2 , merged$lratio.x, merged$lratio.y)

## Convert the two ratio columns from factors to numeric types for t.tests
table$merged.lratio.x<- as.numeric(as.character(table$merged.lratio.x))
table$merged.lratio.y<- as.numeric(as.character(table$merged.lratio.y))

#Filter out any missing values
table<- table[!is.na(table$merged.lratio.x),]
table<- table[!is.na(table$merged.lratio.y),]

#For each row in the table, apply function

#t.test results
test <- data.frame(matrix (ncol=4, nrow=0))
colnames(test) <- c("Amino acids", "TValue", "df", "PValue")                   

#Perform paired t.test for all the tag ratios for each protein between the first sample and the second
for(peptide in unique(table$merged.ac)){
  x<- table[which(table$merged.ac == peptide),3]
  y<- table[which(table$merged.ac == peptide),4]
  
  t <- t.test( x, y,  paired = TRUE)
  
  df <- data.frame(peptide, sprintf("%.3f", t$statistic), t$parameter , sprintf("%.3f", t$p.value))
  colnames(df) <- c("Amino acids", "TValue", "df", "PValue")
                   
  test <- rbind( test,df)
  }

# Clustering

#Go term analysis