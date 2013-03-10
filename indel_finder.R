## This program is meant to take in a set of FASTA sequences and 
## find out whether a certain set of flanking sequences is present,
## how many bases are in between them, and what the sequence is.

############## SETUP ################################
library("seqinr")
library("Biostrings")

#ptpmt1.3 search terms
search.terms <- DNAStringSet(c("first_binding_site" = "GAGTGGCAGTCTGTT",
                               "second_binding_site" = "TGGACACTGTGGACCT",
                               "target_site" = "GGCGTTGAGCAGATCAGAT"))

filename <- "ptpmt1.3.txt"   
wt.sequence <- "gccaccgttgaatgaataaataaaaaaaacaagcaaacaaacaaatgtaaatattagaaaggaataacaatattttagctttggggtttattttttttgtcttttgccctggaaatgattttaatctgtgacttatttacatgttgcacaacaaagcatttcagaattagatttttttaaaaaaaaggtaaaactacttcagattaaaatggttgaatatttctttttgttccactgtagGTGCAGCTCAGAGCGGCGTAAGGAGAAATCTCGTGATGCCGCGCGCTGCAGACGGAGTAAAGAGACAGAGGTGTTTTATGAACTGGCTCATCATCTTCCCCTTCCACACAGCATCAGCTCACATTTGGATAAAGCGTCCATCATGAGACTGGCTATCAGCTTCCTGCGGACACGCAAACTCGTCAACTCAGgtacacagtcagtatatgacaattattaattcaaaccagctttattatattgaacaagaaggtcacataaactgcaatg"


# #epas search terms
# search.terms <- DNAStringSet(c("first_binding_site" = "TACAATACTCCCACTGAA",
#                                "second_binding_site" = "ACTCATGGACAGTTGGTA",
#                                "target_site" = "ATGACAGATGCAGACAG"))
# filename <- "epas1b.seq"   
# wt.sequence <- "cactgttgttaggagggttcagagtagcaggatgaagttgctgttgtttattttggatgtgagccaagggcttgagagctagaacaagactagtatagtgtgcacacacactaacttgcattctaaaactcttgtgtttgtgctgtattgcagGCTACAATACTCCCACTGAAATGACAGATGCAGACAGACTCATGGACAGTTGGTATCTGAAGTCACTCGGTGGCTTTATTACAGTGGTAACATCAGATGGAGACATGATCTTCTTATCGGAGAACATCAACAAtagtaacgcacactgtatcaacacatgaatcga"

# #gria search terms
# search.terms <- DNAStringSet(c("first_binding_site" = "TCGTCCAATAGCTTCT",
#                                "second_binding_site"="gagtttctgctcttta",
#                                "target_site" = "CAGTCACGCACGCCTgt"))
# filename <- "gria_sequence.seq"  
# wt.sequence <- "cactactactgctgtccttcactcgaacaagtctgaagtgaagtgtatgttcttaacccctctgatctcgcatgcagGTGGTCTGTTCATGCGCTCCACGGTCCAGGAGCACAGCGCGTTCCGCTTCGCCGTTCAGCTCTACAACACCAATCAGAACATCACTGAGAAACCTTTCCATCTCAATTACAACGTGGACAATCTGGAGTCGTCCAATAGCTTCTCAGTCACGCACGCCTgtgagtttctgctctttatctttcccccacactcaaaaaaataccaccaccaccatcaaaaataaaacaaatctcaccacagtctctatattgcttgaagactcattgcggtattgaaaagaaaagctttttattctcaccaaagacgtgggtttgatttgcaacttgttagggacatattgtgaggaggcaggagacaaccctcggattctaatttt"


#nrxn search terms
# search.terms <- DNAStringSet(c(	"first_binding_site"="ATCTTCAGC", 
#            				          	  "second_binding_site"="GATGAGGTT",
#                                 "target_site" = "CATAAAA"))
# filename <- "nrxn1.seq"
# wt.sequence <- "GCGCAGGAATTTCAAAAACACTACTTTAGTTGTGGACGAAGAAATCAAGTGGGTGGAGGTAAAGTCGAAACGGAGGGACATGACGGTCTTCAGCCATTTATTCTTAGGGGGGATACCTCCTGAACTGCGATCTGTAGCATTACGCCTCACATCTTCAGCCATAAAAGATGAGGTTCCCTACAAAGGATGGATAACCAACCTGAGAGTGAACGGCTCGGAGCCGGTGCTTATCGGTAGCGATGGAGTCAACAGCGACATTTGCGAAGCCGACCACATTTGCCTGAATGGAGGAGT"

bp.cutoff.value <- 275
percent.n.cutoff.value <- 32
table.length <- NA
max.mismatch <- 2


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
  if(class(DNAString.object) == "XStringViews" ||
     class(DNAString.object) == "DNAString" ||
     class(DNAString.object) == "DNAStringSet"){
    if(length(DNAString.object) > 1){
      pct.n.vector <- vector(length=length(DNAString.object))
      for(i in 1:length(DNAString.object)){
        total.bps <- width(DNAString.object[i])
        num.ns <- vcountPattern("N", DNAString.object[i])
        pct.ns <- (num.ns/total.bps) * 100
        pct.n.vector[i] <- signif(pct.ns, digits=3)
      }
      return(pct.n.vector)
    }  
    else if (length(DNAString.object) == 1
             && class(DNAString.object) != "logical"){
      total.bps <- width(DNAString.object)
      num.ns <- countPattern("N", DNAString.object)
      pct.ns <- (num.ns/total.bps) * 100
      pct.n <- signif(pct.ns, digits = 3)
      return(pct.n)
    }
  }
  else( return("wrong class") )
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
screen.for.low.bps.and.ns <- function(){
  initial_screen <- vector(length=table.length)
  initial_screen[1:nrow(master.table)] <- "OK"
  for(i in 1:table.length){
    if (master.table$total_bps[i] < bp.cutoff.value && 
          master.table$percent_ns[i] > percent.n.cutoff.value){
      initial_screen[i] <- "Omitted: Too many N's + low bp count"
    } 
    else if (master.table$total_bps[i] < bp.cutoff.value){
      initial_screen[i] <- "Omitted: Low bp count"
    } 
    else if( master.table$percent_ns[i] > percent.n.cutoff.value){
      initial_screen[i] <- "Omitted: Too many 'N's "
    } 
  }  
  return(initial_screen)
}
findLRmatch <- function(L.search.term, R.search.term, 
                        sequences.to.be.searched, 
                        max.L.mismatch, max.R.mismatch, 
                        max.gap.length = 50, l.and.r.are.fixed = F ){
  matches <- vector("list",length(sequences.to.be.searched))
 
  for(i in 1:length(sequences.to.be.searched)){
    matches[[i]] <- matchLRPatterns(L.search.term, 
                                  R.search.term,
                                  max.gap.length, 
                                  sequences.to.be.searched[[i]], 
                                  Lfixed=l.and.r.are.fixed, 
                                  Rfixed=l.and.r.are.fixed,
                                  max.Lmismatch=max.L.mismatch,
                                  max.Rmismatch=max.R.mismatch)  
  }
  return(matches)
} 
convert.to.numeric <- function(vector){
  converted <- as.numeric(as.character(vector))
  return(converted)
}
find.best.LR.match <- function(match.list,i){
  current <- match.list[[1]]
  table.column <- vector(length=length(match.list))
  #if more than one match, finds one with least Ns and gets rid of others
  if(length(current) > 0){
    
    num.matches <- length(current)
    current <- trim(current)
    min.pct <- 101
    best.index <- NA
    
    for (j in 1:num.matches){
      current.width <- width(current[j])
      current.Ns <- countPattern("N", current[j])
      pct.ns <- (current.Ns/current.width)*100
      if(pct.ns < min.pct) {
        best.index <- j
        min.pct <- pct.ns
      }
    }
    match.list[i] <- current[best.index] 
    return(match.list[i])
  }  
  else{
    empty.view <- Views(DNAString("NNN"), start=1,end=1)
    match.list[i] <-empty.view
    return(match.list[i])
  }

}
find.f.and.r.matches.for.each.sequence <- function(){
  matches <- list("f.mismatch.0" = 0, "f.mismatch.1" = 0,"f.mismatch.2"= 0,
                  "r.mismatch.0"= 0,"r.mismatch.1"= 0,"r.mismatch.2"= 0)
  f.mismatch.0 <- vector("list", table.length)
  f.mismatch.1 <- vector("list", table.length)
  f.mismatch.2 <- vector("list", table.length)
  
  r.mismatch.0 <- vector("list", table.length)
  r.mismatch.1 <- vector("list", table.length)
  r.mismatch.2 <- vector("list", table.length)
  
  matches$f.mismatch.0 <- f.mismatch.0
  matches$f.mismatch.1 <- f.mismatch.1
  matches$f.mismatch.2 <- f.mismatch.2
  
  matches$r.mismatch.0 <- r.mismatch.0
  matches$r.mismatch.1 <- r.mismatch.1
  matches$r.mismatch.2 <- r.mismatch.2
  
  
#   f.match.mismatch.0  <- vector("list", table.length)
#   print(paste("length f.match.mismatch.0 empty)", length(f.match.mismatch.0)))
#   f.match.mismatch.1  <- vector("list", table.length)
#   f.match.mismatch.2  <- vector("list", table.length)
#   
#   r.match.mismatch.0  <- vector("list", table.length)
#   r.match.mismatch.1  <- vector("list", table.length)
#   r.match.mismatch.2  <- vector("list", table.length)
  
  forward.matches <- vector("list", table.length)
  reverse.matches <- vector("list", table.length)
#  print(paste("length of forward.matches initially", length(forward.matches)))
#  print(paste("length of reverse.matches initially", length(reverse.matches)))
  
  
  #print(paste("matches",matches))
  for(i in 1:table.length){
  #for(i in 0:table.length){
    if(master.table$initial_screen[i] == "OK"){
      for(mismatch.num in 0:max.mismatch){
        forward.matches[i] <- findLRmatch(all.search.terms$first_binding_site,
                                          all.search.terms$second_binding_site,
                                          DNA.sequences.trimmed[i], mismatch.num,mismatch.num)
        
        reverse.matches[i] <- findLRmatch(all.search.terms$second_binding_site_rev,
                                          all.search.terms$first_binding_site_rev,
                                          DNA.sequences.trimmed[i], mismatch.num, mismatch.num)
        
        if(mismatch.num == 0){
#          print(paste("length matches[[1]] inside if 1 BEFORE", length(matches[[1]]), i))
          matches$f.mismatch.0[i] <- find.best.LR.match(forward.matches[i],i)
          matches$r.mismatch.0[i] <- find.best.LR.match(reverse.matches[i],i)
#          print(paste("length matches[[1]] inside if 1", length(matches[[1]]), i))
        }
        else if (mismatch.num == 1){
          matches$f.mismatch.1[i] <- find.best.LR.match(forward.matches[i],i)
          matches$r.mismatch.1[i] <- find.best.LR.match(reverse.matches[i],i)  
#          print(paste("length matches[[1]] inside if 2", length(matches[[1]]), i))
          
        }
        else if (mismatch.num == 2){
          matches$f.mismatch.2[i] <- find.best.LR.match(forward.matches[i],i)
          matches$r.mismatch.2[i] <- find.best.LR.match(reverse.matches[i],i)
#          print(paste("length matches[[1]] inside if 3", length(matches[[1]]), i))
          
        } 
      }
    }
  }
#  print(paste("length matches[[1]] at end of find.f.and.r.matches.for.each.sequence()", length(matches[[1]])))
  return(matches)
}
find.best.match.in.a.direction <- function(mismatch.0, mismatch.1, mismatch.2){
  best.matches <- vector("list", table.length)
  #best.matches[1:length(best.matches)]= NULL
  for(i in 1:table.length){
    if(is.na(mismatch.0[i]) == FALSE &&
       is.null(mismatch.0[[i]]) == FALSE &&
       width(mismatch.0[[i]]) > 1 && 
       class(mismatch.0[[i]]) == "XStringViews"){
      best.matches[i] <- mismatch.0[[i]]
    }
  }
  
  for(i in 1:table.length){
    if(is.na(mismatch.1[i]) == FALSE &&
       is.null(mismatch.1[[i]]) == FALSE &&
       class(best.matches[[i]]) != "XStringViews" &&
       width(mismatch.1[[i]]) > 1){
      best.matches[i] <- mismatch.1[[i]]
    }
  }
  
  for(i in 1:table.length){
    if(is.na(mismatch.2[i]) == FALSE &&
       is.null(mismatch.2[[i]]) == FALSE &&
       class(best.matches[[i]]) != "XStringViews" &&
       width(mismatch.2[[i]]) > 1){ 
      best.matches[i] <- mismatch.2[[i]]
    }
  }
  return(best.matches)
}
find.pct.ns.of.match <- function(best.matches){
  pct.ns <- vector(length=table.length)
  for(i in 1:table.length){
    if(#width(best.matches[[i]]) == 1 #||
       is.null(best.matches[[i]]) == TRUE){
    #ß){
      pct.ns[i] <- 1000
    }
    else{
      pct.ns[i] <- find.pct.ns(best.matches[[i]])
    }
  }
  return(pct.ns)
}
find.average.start.of.matches <- function(direction.string){
  add <- 0
  starting.position.sum <- 0
  total.number.of.sequences <- sum(master.table$direction == direction.string)
  for(i in 1:table.length){
    if(master.table$direction[i] == direction.string){
      add <- start(best.matches.final[[i]])
      starting.position.sum <- starting.position.sum + add
    }
  }
  average <- starting.position.sum/total.number.of.sequences
  return(average)
}
make.all.sequences.forward <- function(){
  sequences.for.alignment <- vector(length=table.length)
  sequences.for.alignment[1:length(sequences.for.alignment)] <- "probably not it"
  for(i in 1:table.length){
   
    if(master.table$direction[i] == "reverse"){
      end.of.match.index <- end(best.matches.final[[i]])
      sequence.for.trimming <- as.character(DNA.sequences.trimmed[[i]])
      sequences.for.alignment[i] <- substr(sequence.for.trimming, 1, 
                                           end.of.match.index+average.f.start)
      sequences.for.alignment[i] <- as.character(reverseComplement(
                                                 DNAString(sequences.for.alignment[i])))
    }
    else if(master.table$direction[i] == "forward" ){
      sequences.for.alignment[i] <- as.character(DNA.sequences.trimmed[[i]])
    }
  }
  return(sequences.for.alignment)
}
create.vector.of.match.lengths <- function(){
  match_length <- vector(length=table.length)
  acceptable.pct.n.for.length.of.binding.site <- 20
  for(i in 1:table.length){
    if(master.table$direction[i] != "probably not it" && 
         master.table$initial_screen[i] == "OK" &&
         best.pct.ns[i] < acceptable.pct.n.for.length.of.binding.site){
      match_length[i] <- width(best.matches.final[[i]])
    }
    else{
      match_length[i] <- "no good match"
    }
  }
  return(match_length)
}
select.only.notable.matches <- function(){
  WT.match.length <- length(c(search.terms$first_binding_site, 
                              search.terms$second_binding_site, 
                              search.terms$target_site))
  align.table <- cbind("well_no"= master.table$well_no, "sequences"=sequences.for.alignment)
  for(i in table.length:1){
    if(align.table[i,2] == "probably not it" || 
       align.table[i,2] == FALSE || 
       master.table$match_length[i] == "no good match" || 
       master.table$match_length[i] == WT.match.length){
      align.table <- as.data.frame(align.table[-i,])
    }
  }
  return(align.table)
}
label.sequence.direction <- function(forward.pct.ns, reverse.pct.ns,forward.best.matches, reverse.best.matches){
  best.matches.final <- vector("list", table.length)
  best.pct.ns <- vector(length=table.length)
  direction <- vector(length = table.length)
  for(i in 1:table.length){
    if(reverse.pct.ns[i] < forward.pct.ns[i]){
      best.matches.final[[i]]<- reverse.best.matches[[i]]
      direction[i] <- "reverse"
      #master.table$match_length[i] <- width(best.matches.final[[i]])
      best.pct.ns[i] <- reverse.pct.ns[i]
    }
    if(forward.pct.ns[i] < reverse.pct.ns[i]){
      best.matches.final[[i]] <- forward.best.matches[[i]]
      direction[i] <- "forward"
      #master.table$match_length[i] <- width(best.matches.final[[i]])
      best.pct.ns[i] <- forward.pct.ns[i]
    }
    if(forward.pct.ns[i] == reverse.pct.ns[i]){
      best.matches.final[[i]] <- "probably not it"
      direction[i] <-"probably not it"
      #master.table$match_length[i] <- "probably not it"
      #best.pct.ns[i] <- "probably not it"
      best.pct.ns[i] <- forward.pct.ns[i]
    }
  }
  results <- list("best.matches" = best.matches.final, 
                  "best.pct.ns" = best.pct.ns, 
                  "direction" = direction)
  return(results)
}

############## THINGS ARE HAPPENING ####################

############## setting up table with sequences
all.search.terms <- find.reverse.complement.DNAStringSet(search.terms)
raw.sequence.file <- read.fasta.file(filename)
sequence.table <- parse.sequence.name.and.divide.into.table(raw.sequence.file)
#print(paste("length sequence.table", nrow(sequence.table)))
table.length <- nrow(sequence.table)

#converts sequence column into set of DNAString class
DNA.sequences<- DNAStringSet(sequence.table$sequence)
names(DNA.sequences) <- sequence.table$well_no
DNA.sequences.trimmed <- DNAStringSet(sequence.table$trimmed_sequences)
names(DNA.sequences.trimmed) <- names(DNA.sequences)

############## logs basic info and puts them in a table
n.table <- find.sequence.basic.stats(DNA.sequences)

############## Create  big table to store stats 
master.table <- merge(sequence.table, n.table, by="well_no")

############## Converts columns to numeric
master.table$total_bps <-convert.to.numeric(master.table$total_bps)
master.table$num_Ns <-convert.to.numeric(master.table$num_Ns)
master.table$percent_ns <-convert.to.numeric(master.table$percent_ns)

############## Screens sequence char stats 
master.table$initial_screen <- screen.for.low.bps.and.ns()

############## Perform matchLRPattern
matches <- find.f.and.r.matches.for.each.sequence()
f.match.mismatch.0 <- matches$f.mismatch.0
f.match.mismatch.1 <- matches$f.mismatch.1
f.match.mismatch.2 <- matches$f.mismatch.2
r.match.mismatch.0 <- matches$r.mismatch.0
r.match.mismatch.1 <- matches$r.mismatch.1
r.match.mismatch.2 <- matches$r.mismatch.2
#print(paste("length f.mismatch.0",length(f.match.mismatch.0)))


forward.best.matches <- find.best.match.in.a.direction(f.match.mismatch.0,
                                                       f.match.mismatch.1,
                                                       f.match.mismatch.2) 
reverse.best.matches <- find.best.match.in.a.direction(r.match.mismatch.0,
                                                       r.match.mismatch.1,
                                                       r.match.mismatch.2) 

forward.pct.ns <- find.pct.ns.of.match(forward.best.matches)
reverse.pct.ns <- find.pct.ns.of.match(reverse.best.matches)


############### decides if forward or reverse 
matches.and.pct.ns <- label.sequence.direction(forward.pct.ns, reverse.pct.ns, 
                                               forward.best.matches,reverse.best.matches)
### splits return from function into 3 pieces
best.matches.final<- matches.and.pct.ns$best.matches
best.pct.ns <- matches.and.pct.ns$best.pct.ns
master.table$direction <- matches.and.pct.ns$direction

############### finds average forward and reverse starting position
average.f.start <- find.average.start.of.matches("forward")
average.r.start <- find.average.start.of.matches("reverse")

############## Makes all sequences forward 
sequences.for.alignment <- make.all.sequences.forward()
master.table$match_length <- create.vector.of.match.lengths()

############## makes table for FASTA ALIGNMENT
align.table <- select.only.notable.matches()

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

write.table(master.table, file = "summary.all.csv", sep = ",", col.names = NA,
                        qmethod = "double")


FASTA <- file("FASTA", open="w")
writeLines(c(paste(">", "WT"), wt.sequence), FASTA)
for (i in 1:nrow(align.table)){
  writeLines(c(paste(">", align.table[i,1]), 
               as.character(align.table[i,2])), 
             FASTA)    
}
close(FASTA)
