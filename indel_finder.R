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
percent.n.cutoff.value <- 32
filename <- "nrxn1.seq"
table.length <- NA


############## FUNCTION DEFINITIONS ##########################################


input.search.terms<- function(){}
read.fasta.file <- function(filename){
  #reads in file
  sequences <- read.fasta(filename, 
                          seqtype = "DNA", 
                          as.string = TRUE, 
                          set.attributes = FALSE, 
                          forceDNAtolower = FALSE)
}
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
find.pct.ns <- function(DNAString.object){
  if(length(DNAString.object) > 1){
    pct.n.vector <- vector(length=length(DNAString.object))
    for(i in 1:length(DNAString.object)){
      total.bps <- width(DNAString.object[i])
      num.ns <- vcountPattern("N", DNAString.object[i])
      pct.ns <- (num.ns/total.bps) * 100
      pct.n.vector[i] <- pct.ns
    }
    return(pct.n.vector)
  }  
  else {
    total.bps <- width(DNAString.object)
    num.ns <- vcountPattern("N", DNAString.object)
    pct.ns <- (num.ns/total.bps) * 100
    pct.n <- pct.ns
    return(pct.n)
  }
}
find.sequence.basic.stats <- function(DNA.set){
  
  n.percent.vector <- find.pct.ns(DNA.set)
  sequence.length.vector <- width(DNA.set)
  num.ns.vector <- vcountPattern("N", DNA.set)
  stats.table <- cbind("well_no"=sequence.table$well_no, 
                       "total_bps"=sequence.length.vector, 
                       "num_Ns"=num.ns.vector, 
                       "percent_ns"=n.percent.vector)
  return(stats.table)
}
findLRmatch <- function(L.search.term, R.search.term, sequences.to.be.searched, max.gap.length = 20, l.and.r.are.fixed = F ){
  matches <- list(vector(length=length(sequences.to.be.searched)))
  for(i in 1:length(sequences.to.be.searched)){
    matches[i] <- matchLRPatterns(L.search.term, 
                                  R.search.term,
                                  max.gap.length, 
                                  sequences.to.be.searched[[i]], 
                                  Lfixed=l.and.r.are.fixed, 
                                  Rfixed=l.and.r.are.fixed)  
  }
  return(matches)
} 
convert.to.numeric <- function(vector){
  converted <- as.numeric(as.character(vector))
  return(converted)
}

############## THINGS ARE HAPPENING ####################

all.search.terms <- find.reverse.complement.DNAStringSet(search.terms)

raw.sequence.file <- read.fasta.file(filename)

sequence.table <- parse.sequence.name.and.divide.into.table(raw.sequence.file)

table.length <- nrow(sequence.table)

#converts sequence column into set of DNAString class
DNA.sequences<- DNAStringSet(sequence.table$sequence)
names(DNA.sequences) <- sequence.table$well_no
DNA.sequences.trimmed <- DNAStringSet(sequence.table$trimmed_sequences)
names(DNA.sequences.trimmed) <- names(DNA.sequences)
n.table <- find.sequence.basic.stats(DNA.sequences)

############## Create  big table to store stats ####################
master.table <- merge(sequence.table, n.table, by="well_no")

master.table$total_bps <-convert.to.numeric(master.table$total_bps)
master.table$num_Ns <-convert.to.numeric(master.table$num_Ns)
master.table$percent_ns <-convert.to.numeric(master.table$percent_ns)

############## Screens sequence char stats ####################################
master.table$initial_screen <- "OK"
#annotates master.table
for(i in 1:table.length){
  if (master.table$total_bps[i] < bp.cutoff.value && 
      master.table$percent_ns[i] > percent.n.cutoff.value){
        master.table$initial_screen[i] <- "Omitted: Too many N's + low bp count"
   } 
  else if (master.table$total_bps[i] < bp.cutoff.value){
      master.table$initial_screen[i] <- "Omitted: Low bp count"
   } 
  else if( master.table$percent_ns[i] > percent.n.cutoff.value){
      master.table$initial_screen[i] <- "Omitted: Too many 'N's "
   } 
}

############## Perform matchLRPattern ########################################

#results <- list(vector(length=table.length))
master.table$direction <- "omitted"
#master.table$results <- NA

forward.matches <- list(vector(length=table.length))
reverse.matches <- list(vector(length=table.length))


for(i in 1:table.length){
  if(master.table$initial_screen[i] == "OK"){
    
    reverse.matches[i] <- findLRmatch(all.search.terms$second_binding_site_rev,
                                   all.search.terms$first_binding_site_rev,
                                   DNA.sequences.trimmed[i])
    
    forward.matches[i] <- findLRmatch(all.search.terms$first_binding_site,
                                   all.search.terms$second_binding_site,
                                   DNA.sequences.trimmed[i])
    
 
   if (length(reverse.matches[[i]]) > 0 && length(forward.matches[[i]]) >0){
      master.table$direction[i] <- "both directions"
    }
    else if(length(reverse.matches[[i]]) > 0){
      master.table$direction[i] <- "reverse"
    }
    else if (length(forward.matches[[i]]) > 0){
      master.table$direction[i] <- "forward"
    }
    else { master.table$direction[i] <- "not found"
    }
  }
}


for (i in 1:table.length){
  master.table$for.matches[i] <- length(forward.matches[[i]])
  master.table$rev.matches[i] <- length(reverse.matches[[i]])
}


for (i in 1:table.length){
  
  if(length(forward.matches[[i]]) > 0 && length(reverse.matches[[i]]) > 0){

    f.num.matches <- length(forward.matches[[i]])
    r.num.matches <- length(reverse.matches[[i]])
    
    f.current <- forward.matches[[i]]
    r.current <- reverse.matches[[i]]
    
    f.min.pct <- 101
    r.min.pct <- 101
    f.best.index <- NA
    r.best.index <- NA
    
    for (j in 1:f.num.matches){
      #find percentage of Ns in each.
      current.width <- width(f.current[j])
      current.Ns <- countPattern("N", f.current[j])
      pct.ns <- (current.Ns/current.width)*100
      if(pct.ns < f.min.pct) {
        f.best.index <- j
        f.min.pct <- pct.ns
      }
    }
     for (k in 1:r.num.matches){
       #find percentage of Ns in each.
       current.width <- width(r.current[k])
       current.Ns <- countPattern("N", r.current[k])
       pct.ns <- (current.Ns/current.width)*100
      if(pct.ns < r.min.pct) {
        r.best.index <- k
        r.min.pct <- pct.ns
      }
     }
    
    if (f.min.pct < r.min.pct){
      master.table$direction[i] <- "forward"
      forward.matches[[i]] <- forward.matches[[i]][f.best.index]
      reverse.matches[[i]] <- NA
    }
    else if (r.min.pct < f.min.pct){
      master.table$direction[i] <- "reverse"
      reverse.matches[[i]] <- reverse.matches[[i]][r.best.index]
      forward.matches[[i]] <- NA
    }
  }
}

master.table$match_length <- 0

for(i in 1:table.length){
  if(master.table$direction[i] == "forward"){
    master.table$match_length[i] <- width(forward.matches[[i]])
  }
  else if (master.table$direction[i] == "reverse"){
    master.table$match_length[i] <- width(reverse.matches[[i]])
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



