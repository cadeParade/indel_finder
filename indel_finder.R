## This program is meant to take in a set of FASTA sequences and 
## find out whether a certain set of flanking sequences is present,
## how many bases are in between them, and what the sequence is.

## There are no functions yet because without functions I can access the 
## variables in the console and play with them. At the end there will be 
## functions!

############## SETUP ################################
library("seqinr")
library("Biostrings")

#gria search terms
# search.terms <- DNAStringSet(c("first_binding_site" = "TCGTCCAATAGCTTCT",
#                                "target_site"="CAGTCACGCACGCCTgt", 
#                                "second_binding_site"="gagtttctgctcttta", 
#                                "primer1"="GCGCAGGAATTTCAAAAACACTA", 
#                                "primer2"="ACATTTGCCTGAATGGAGGAGT", 
#                                "upstream_site"="AAGTGGGTGGA"))

#nrxn search terms
search.terms <- DNAStringSet(c(	"first_binding_site"="ATCTTCAGC", 
           				          	"target_site"="CATAAAA", 
           				          	"second_binding_site"="GATGAGGTT", 
           				          	"primer1"="GCGCAGGAATTTCAAAAACACTA", 
                     					"primer2"="ACATTTGCCTGAATGGAGGAGT", 
                     					"upstream_site"="AAGTGGGTGGA"))	
#search.terms.dict <- PDict(search.terms)

bp.cutoff.value <- 250
percent.n.cutoff.value <- 50
filename <- "nrxn1.seq"
table.length <- NA
input.search.terms<- function(){}
read.file <- function(){}
############## Find reverse complement of search terms #######################

find.reverse.complement.DNAStringSet <- function(set.to.reverse){
  #finds reverse complement of search terms
  rev.comp <- DNAStringSet()
  for (i in seq_along(set.to.reverse)){
    rev.comp[i] <- reverseComplement(set.to.reverse[i])
  }
  
  #names reverse complement items in search term
  names(rev.comp)<-paste(names(set.to.reverse),"rev",sep="_")
  
  forward.and.reverse.set <- c(set.to.reverse, rev.comp)
  return(forward.and.reverse.set)
}

all.search.terms <- find.reverse.complement.DNAStringSet(search.terms)


############## Divides raw text into table ################################
read.fasta.file <- function(filename){
  #reads in file
  sequences <- read.fasta(filename, 
                          seqtype = "DNA", 
                          as.string = TRUE, 
                          set.attributes = FALSE, 
                          forceDNAtolower = FALSE)
}

raw.sequence.file <- read.fasta.file(filename)

parse.sequence.name.and.divide.into.table <- function(raw.sequence.file){

  #declares vectors for loop below
  plate.row<-vector(mode="character")
  plate.column<- vector(mode="numeric")
  well.number<-vector(mode="character")
  sequence.column <- vector(mode="character")
  
  # puts well number and sequence in different columns
  # splices well number
  for(i in seq_along(raw.sequence.file)){
  	current<-names(raw.sequence.file[i])
  	split.well.number <- regmatches(current[1],	regexec("([A-H])(\\d\\d)",current[1]))
  	split.well.number<-unlist(split.well.number)
  	
  	well.number[i]<-split.well.number[1]
  	plate.row[i]<- split.well.number[2]
  	plate.column[i]<-split.well.number[3]
  	
  	sequence.column[i] <- as.character(raw.sequence.file[i])                              
  }
  
  sequence.table <- as.data.frame(cbind("well_no"=well.number, 
  						                          "plate_row"=plate.row, 
                                        "plate_column"=plate.column,
  						                          "sequence"=sequence.column), 
                                    stringsAsFactors= FALSE)
  
  #trims sequence string to remove some leading and trailing Ns
  sequence.table$trimmed_sequences <- substring(sequence.table$sequence,25, 
                                                nchar(sequence.table$sequence)-50)
  return(sequence.table)
}

sequence.table <- parse.sequence.name.and.divide.into.table(raw.sequence.file)

#defines length variable to be used in loops
table.length <- nrow(sequence.table)

############## Find percentage of N's in each sequence ###########
#converts sequence column into set of DNAString class
DNA.sequences<- DNAStringSet(sequence.table$sequence)
DNA.sequences.trimmed <- DNAStringSet(sequence.table$trimmed_sequences)

find.sequence.basic.stats <- function(DNA.set){

  sequence.length.vector <- width(DNA.set)
  num.ns.vector <- vcountPattern("N", DNA.set)
  
  n.percent.vector <- format((num.ns.vector/sequence.length.vector)*100, digits=2)
  stats.table <- cbind("well_no"=sequence.table$well_no, 
                   "total_bps"=sequence.length.vector, 
                   "num_Ns"=num.ns.vector, 
                   "percent_ns"=n.percent.vector)
  return(stats.table)
}

n.table <- find.sequence.basic.stats(DNA.sequences)

############## Create  big table to store stats ####################
master.table <- merge(sequence.table, n.table, by="well_no")
#master.table$sequence <- as.character(master.table$sequence)

master.table$total_bps <-as.numeric(as.character(master.table$total_bps))
master.table$num_Ns <-as.numeric(as.character(master.table$num_Ns))
master.table$percent_ns <-as.numeric(as.character(master.table$percent_ns))

############## Screens sequence char stats ####################################
master.table$initial_screen <- "OK"

#annotates master.table
for(i in 1:table.length){
  if (master.table$total_bps[i] < bp.cutoff.value && 
      master.table$percent_ns[i] > percent.n.cutoff.value){
        master.table$initial_screen[i] <- "Omitted: Too many N's + low bp count"
   } else if (master.table$total_bps[i] < bp.cutoff.value){
      master.table$initial_screen[i] <- "Omitted: Low bp count"
   } else if( master.table$percent_ns[i] > percent.n.cutoff.value){
      master.table$initial_screen[i] <- "Omitted: Too many 'N's "
   } 
}


 

############## Perform matchLRPattern ########################################

results <- list(vector(length=table.length))
master.table$sequence_is_reverse <- NA
master.table$LRmatch_length <- 0

for(i in 1:table.length){
  if(master.table$initial_screen[i] == "OK"){
    results[i] <- matchLRPatterns(all.search.terms$first_binding_site, 
                                  all.search.terms$second_binding_site, 
                                  20, 
                                  DNA.sequences.trimmed[[i]], 
                                  Lfixed=F, 
                                  Rfixed=F)  
    if ( length(width(results[[i]])) > 0){
     master.table$sequence_is_reverse[i] <- "forward"
     #master.table$LRmatch_length[i] <- width(results[[i]])
     print(width(results[[i]]))
    }
   else { 
     results[i] <- matchLRPatterns(all.search.terms$second_binding_site_rev, 
                                   all.search.terms$first_binding_site_rev, 
                                   20, 
                                   DNA.sequences.trimmed[[i]], 
                                   Lfixed=F, 
                                   Rfixed=F)
     
      if (length(width(results[[i]])) == 1){
         master.table$sequence_is_reverse[i] <- "reverse"
         master.table$LRmatch_length[i] <- width(results[[i]])
      }
    }
  } 
}



# ############# Gets reverse complement of reversed sequences ##################
# #sequences.to.reverse <- master.table[master.table$reversed,]
# master.table$reversed_sequence <- "not reversed"
# 
# #reverses sequences that are set to be analyzed and marked as reverse
# for (i in 1:table.length){
#   if (master.table$analyze[i] == TRUE && master.table$direction[i] == "reverse"){
#     master.table$reversed_sequence[i]<-toupper(c2s(rev(comp(s2c(master.table$sequence[i]), ambiguous = T))))
#   }
# }


# ############## Input search terms #############################################
# #input.search.terms <- function(){
# 	#print(head(sequence.table))
# 	#print("Input searchable sequences")
# 	#print("1: sequence before first primer")
# 	#print("2: first binding site")
# 	#print("3: target site")
# 	#print("4: second binding site")
# 	#print("5: first primer")
# 	#print("6: second primer")
# 	#print("7: plasmid splice site")
# 	#prompts user in console for searchable strings
# 	#sequence.markers <-scan("",what="character", nmax=7)
# #}
# 	
# 
# ############## Write to files #################################################
# 
 write.table(master.table, file = "summary.all.csv", sep = ",", col.names = NA,
                          qmethod = "double")
# 
# #align.table <- as.data.frame(matrix(NA, nrow=nrow(working.table), ncol=2) )
# 
# #align.table[,1] <- working.table$well_no
# #align.table[,2] <- working.table$sequences_to_align
# 
# #print(head(align.table))
# 
# #write.table(align.table, file="align.csv", row.names=FALSE, col.names=FALSE, sep=",")
# #print(head(working.table))
# #print(nrow(working.table))
# 
# FASTA <- file("FASTA", open="w")
# 
# for (i in 1:nrow(working.table)){
#   writeLines(c(paste(">", working.table$well_no[i]), working.table$sequences_to_align[i]), FASTA)    
# }
# close(FASTA)