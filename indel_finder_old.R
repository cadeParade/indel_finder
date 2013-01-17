## This program is meant to take in a set of FASTA sequences and 
## find out whether a certain set of flanking sequences is present,
## how many bases are in between them, and what the sequence is.

## There are no functions yet because without functions I can access the 
## variables in the console and play with them. At the end there wil be 
## functions!

############## SETUP ################################
library("seqinr")
library("Biostrings")

search.terms <- DNAStringSet(c(	"first_binding_site"="ATCTTCAGC", 
            				          	"second_binding_site"="CATAAAA", 
            				          	"target_site"="GATGAGGTT", 
            				          	"primer1"="GCGCAGGAATTTCAAAAACACTA", 
                      					"primer2"="ACATTTGCCTGAATGGAGGAGT", 
                      					"upstream_site"="AAGTGGGTGGA"))	
#search.terms.dict <- PDict(search.terms)

bp.cutoff.value <- 250
percent.n.cutoff.value <- 50
filename <- "nrxn1.seq"
table.length <- NA
input.search.terms<- function(){}

############## Divides raw text into table ################################

#reads in file
raw.sequence.file <- readLines(filename)
#removes line breaks
collapsed.sequence.file <- paste(raw.sequence.file, collapse="")
#splits by ">" (where the beginning of each sequence is)
split.sequence.file <- strsplit(collapsed.sequence.file, ">")
#converts list to vector with ~96 objects
vector.sequence.file <- unlist(split.sequence.file)
#splits again to separate title string
split.title.and.sequence <- strsplit(vector.sequence.file, "seq")
#remove first item
split.title.and.sequence<-split.title.and.sequence[2:length(split.title.and.sequence)]

#declares vectors for loop below
plate.row<-vector(mode="character")
plate.column<- vector(mode="numeric")
well.number<-vector(mode="character")
sequence.column <- vector(mode="character")

# puts well number and sequence in different columns
#splices well number
for(i in seq_along(split.title.and.sequence)){
	current<-split.title.and.sequence[[i]]
	split.well.number <- regmatches(current[1],	regexec("([A-H])(\\d\\d)",current[1]))
	split.well.number<-unlist(split.well.number)
	
	well.number[i]<-split.well.number[1]
	plate.row[i]<- split.well.number[2]
	plate.column[i]<-split.well.number[3]
	
	sequence.column[i] <- current[2]                              
}

sequence.table <- as.data.frame(cbind("well_no"=well.number, 
						                          "plate_row"=plate.row, 
                                      "plate_column"=plate.column,
						                          "sequence"=sequence.column), 
                                stringsAsFactors= FALSE)



#defines length variable to be used in loops
table.length <- nrow(sequence.table)

############## Find percentage of N's in each sequence ###########
sequence.table$trimmed_sequences <- substring(sequence.table$sequence,
                                              25, nchar(sequence.table$sequence)-50)

#converts sequence column into set of DNAString class
DNA.sequences <- DNAStringSet(sequence.table$sequence)
#trims N's off begin and end of each sequence
DNA.sequences <- trimLRPatterns(Lpattern="NNNNNNNNNNNNNNNNNNNNNNNNN", Rpattern="NNNNNNNNNNNNNNNNNNNNNNNNN", 
                                DNA.sequences, Lfixed = F, Rfixed=F)



sequence.length.vector <- width(DNA.sequences)
num.ns.vector <- vcountPattern("N", DNA.sequences)

n.percent.vector <- format((num.ns.vector/sequence.length.vector)*100, digits=2)
n.table <- cbind("well_no"=sequence.table$well_no, 
                 "total_bps"=sequence.length.vector, 
                 "num_Ns"=num.ns.vector, 
                 "percent_ns"=n.percent.vector)

############## Create  big table to store stats ####################
master.table <- merge(sequence.table, n.table, by="well_no")
#master.table$sequence <- as.character(master.table$sequence)

master.table$total_bps <-as.numeric(as.character(master.table$total_bps))
master.table$num_Ns <-as.numeric(as.character(master.table$num_Ns))
master.table$percent_ns <-as.numeric(as.character(master.table$percent_ns))

 
############## Find reverse complement of search terms #######################

#finds reverse complement of search terms
rev.comp.search.terms <- DNAStringSet()
for (i in seq_along(search.terms)){
  rev.comp.search.terms[i] <- reverseComplement(search.terms[i])
}

#names reverse complement items in search term
names(rev.comp.search.terms)<-paste(names(search.terms),"rev",sep="_")

all.search.terms <- c(search.terms, rev.comp.search.terms)

############## Perform matchLRPattern ########################################



 lrmatch.vec <- DNA.sequences
 target1 <- search.terms[[1]]
 target2 <- search.terms[[2]]
 

for(i in 1:5){
   current.sequence <- DNA.sequences[[i]]
   print(matchLRPatterns(target1, target2, 20, current.sequence
                         , Lfixed=F, Rfixed=F)) 
                         #max.Lmismatch=3, max.Rmismatch=3))
 }
 

############## Decides whether sequence is forward or reverse   ############################
# 
# #converts sequence column into set of DNAString class
# DNA.sequences <- DNAStringSet(sequence.table$sequence)
# 
# #sets up logical table to fill with T/F whether each search term is found in a sequence
# search.terms.found <- data.frame(matrix(NA, nrow = table.length,
#                                         ncol = length(all.search.terms)))
# names(search.terms.found) <- names(all.search.terms)
# t#emp.vector <- vector(length=table.length)
# 
# 
# for(i in 1:length(all.search.terms)){
# #  print("hi")
#   temp.vector[i] <- vcountPattern(as.character(all.search.terms[i]),
#                                DNA.sequences[i],
#                                #max.mismatch=round(0.1*width(all.search.terms[i])), 
#                                #fixed=FALSE)  
#                                             )
#   search.terms.found[,i] <- temp.vector
# }
#print(head(search.terms.found))

#for(i in 1:12)

#for each search term, enters T/F if found in that sequence
#for (i in 1:length(all.search.terms)){
#  search.terms.found[,i] <- grepl(all.search.terms[i], master.table$sequence, 
#                                  ignore.case = TRUE)
#}

#splits search.terms.found into forward and reverse search terms
#sequence.test.forward <- search.terms.found[,1:6]
#sequence.test.reverse <- search.terms.found[,7:12]

#Adds direction column to master table. Fills with "not found"
#master.table$direction <- "not_found"


#Tests each row to see if any forward sequences are found, any reverse sequences
#are found
#for (i in 1:table.length){
#  if(any(sequence.test.forward[i,] == TRUE)){
#    master.table$direction[i] <- "forward"
#   } 
#  if(any(sequence.test.reverse[i,] == TRUE)){
#    master.table$direction[i] <- "reverse"
#   }
#}

# ############## Screens sequence char stats ####################################
# master.table$initial_screen <- "OK"
# 
# #annotates master.table
# for(i in 1:table.length){
#   if( master.table$total_bps[i] < bp.cutoff.value &&
#         master.table$percent_ns[i] > percent.n.cutoff.value){
#     master.table$initial_screen[i] <- "Omitted: Too many N's + low bp count"
#   } else if (master.table$total_bps[i] < bp.cutoff.value){
#     master.table$initial_screen[i] <- "Omitted: Low bp count"
#   } else if( master.table$percent_ns[i] > percent.n.cutoff.value){
#     master.table$initial_screen[i] <- "Omitted: Too many 'N's "
#   } #else (it.will.analyze.this.row[i] <- TRUE)
# }
# 
# ############# Marks bad rows ##################################################
# 
# master.table$analyze <- FALSE
# 
# for(i in 1:table.length){
#   if(master.table$initial_screen[i] == "OK" 
#      && master.table$direction[i] != "not_found"){
#         master.table$analyze[i] <- TRUE
#   }
# }
# 
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
# 
# ############# Makes master search column ##############################
# 
# master.table$sequences_to_align <- "-"
# for(i in 1:table.length){
#   if (master.table$direction[i] == "forward"){ 
#     master.table$sequences_to_align[i] <- master.table$sequence[i] 
#   }
#   if (master.table$direction[i] == "reverse"){
#     master.table$sequences_to_align[i] <- master.table$reversed_sequence[i]
#   }
# }
# 
# #this is for printing to the FASTA file. 
# working.table <- subset(master.table, master.table$analyze)
# #working.table <- master.table[,master.table$analyze == TRUE]
# 
# ############# Finds matches of target sequences ##########################
# #master.table$match_length <- "No match"
# 
# #match <- matchPattern(WT, C3, max.mismatch=5,  with.indels = TRUE)
# 
# #match.length <- width(match)
# 
# #for (i in 1:table.length){
# #  master.table$match_length[i]<- matchPattern(search.terms[3], master.table$sequences_to_align[i])
# #}
# #print(head(master.table$match_length))
# 
# 
# ############## Find search terms in sequences #################################
# 
# #matches <- list()
# 
# #for(i in seq_along(all.search.terms)){
# #  current.search.term <- all.search.terms[i]
# #  matches[[i]] <- grep(current.search.term, 
# #                  master.table$sequence, 
# #                  ignore.case = TRUE)  
# #}
# 
# #names(matches) <- names(all.search.terms)
# #print (matches)
# 
# #binding1 <- unique(c(matches[[1]], matches[[7]]))
# #binding2 <- unique(c(matches[[2]], matches[[8]]))
# #target <- unique(c(matches[[3]], matches[[9]]))
# #primer1 <- unique(c(matches[[4]], matches[[10]]))
# #primer2 <- unique(c(matches[[5]], matches[[11]]))
# #upstream <- unique(c(matches[[6]], matches[[12]]))
# 
# #print(binding1)
# #print("Rows that match binding site 2: ")
# #print(binding2)
# #print("Rows that match target site: ")
# #print(target)
# 
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
# write.table(master.table, file = "summary.all.csv", sep = ",", col.names = NA,
#                          qmethod = "double")
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