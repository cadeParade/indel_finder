This program is written specifically to solve one problem in one lab. In the course of the scientists' studies, they must go through a large number of DNA sequences (typically 96 at a time) and analyze each one to see whether there is an insertion or a deletion (indel) in a specified area. This would be easy if the bounding sequences of the indel region were always clear, but often the sequence file will come back with "N"s instead of bases or sometimes there are polymorphisms in the binding site which makes a simple find function inadequate. 

This program finds the indel bounding region, including Ns or polymorphisms and returns how long the match is. If the match is the expected length (the length of the two bounding regions plus the wildtype sequence between them) then there are no indels. If the match length is different than this, then there is a very good chance an indel has been created. 

The program uses primarily the matchLRPatterns function from the R Biostrings package. 

matchLRPatterns takes two input sequences, in this case the left and right binding sites, and finds within a target sequence any region that looks like it's flanked by the left and right sequences. 
We can then find how many bases are between the matched pairs in order to see if there are indels in the target sequence. 

I will explain the full process the program goes through so you can trust it as far as you wish. 

1. The raw sequence data is read from a FASTA formatted text file. 
2. The sequences are put into a column of their own. The beginning 25 bases and the ending 50 bases are trimmed from each sequence because they are almost always full of Ns. Without these deleted, the matches become overwhelming.
3. Each sequence is has some initial characteristics recorded -- the total number of bases and the percentage of N's each sequence has. 
4. Sequences are omitted from analysis at this stage if they have too few bases or too high percentage of N's. The threshold for these cutoff values can be changed to alter the analysis parameters. 
5. The matchLRPatterns function is run on each sequence with both forward and reverse binding site sequences. The matchLRPattern function has several parameters you can alter to make it more or less specific. It has been set to allow matching with "N"s as sometimes the sequence is messy but still readable as a match. It has also been set to allow up to two mismatches in both the left and right binding site. This is to account for polymorphisms that sometimes happen there. 
6. The matches are recorded separately for forward and reverse and since there are often still N's at the beginning and ends of the sequences, there are often many matches. To find the best one, the match with the least percent N's in each set of matches is kept while the rest are discarded. 
7. Then, for each well, the percent N's of the forward match and reverse match (if there are both) are compared and the lowest wins. If there are no matches in either direction, the sequence is labeled as "not found".
8. Each well is then labeled as forward or reverse, based on which direction binding site the best match came from. 
9. The match length column is the most important. This tells how long the segment including the left binding site, target site (altered or not), and right binding site is. If the length is different than the expected length calculated by adding the length of the two binding sites and the target site, then there is a good chance of an index in the target site. 
 
The method is not perfect yet and so there may be some false positives. This method also will not recognize point mutations in the target site. 

The program creates a csv file with a table listing all of the results listed above. 
