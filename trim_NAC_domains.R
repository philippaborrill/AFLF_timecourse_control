# script to trim NAC domain alignment
# Philippa Borrill # 02-05-2017

setwd("Y:/PB_AFLF/control_timecourse/TF_analysis/phylogenetics")
library("seqinr")

NAC_aln <- read.alignment(file = "wheat_barley_rice_NAC_msa_NAC_domain_no_numbers.fa", format = "fasta") # read in alignment of NACs

NAC_aln$seq

printMultipleAlignment <- function(alignment, chunksize=60)
{
  # this function requires the Biostrings package
  require("Biostrings")
  # find the number of sequences in the alignment
  numseqs <- alignment$nb
  # find the length of the alignment
  alignmentlen <- nchar(alignment$seq[[1]])
  starts <- seq(1, alignmentlen, by=chunksize)
  n <- length(starts)
  # get the alignment for each of the sequences:
  aln <- vector()
  lettersprinted <- vector()
  for (j in 1:numseqs)
  {
    alignmentj <- alignment$seq[[j]]
    aln[j] <- alignmentj
    lettersprinted[j] <- 0
  }
  # print out the alignment in blocks of 'chunksize' columns:
  for (i in 1:n) { # for each of n chunks
    for (j in 1:numseqs)
    {
      alnj <- aln[j]
      chunkseqjaln <- substring(alnj, starts[i], starts[i]+chunksize-1)
      chunkseqjaln <- toupper(chunkseqjaln)
      # Find out how many gaps there are in chunkseqjaln:
      gapsj <- countPattern("-",chunkseqjaln) # countPattern() is from Biostrings package
      # Calculate how many residues of the first sequence we have printed so far in the alignment:
      lettersprinted[j] <- lettersprinted[j] + chunksize - gapsj
      print(paste(chunkseqjaln,lettersprinted[j]))
    }
    print(paste(' '))
  }
}

printMultipleAlignment(NAC_aln, 60)




cleanAlignment <- function(alignment, minpcnongap, minpcid)
{
  # make a copy of the alignment to store the new alignment in:
  newalignment <- alignment
  # find the number of sequences in the alignment
  numseqs <- alignment$nb
  # empty the alignment in "newalignment")
  for (j in 1:numseqs) { newalignment$seq[[j]] <- "" }
  # find the length of the alignment
  alignmentlen <- nchar(alignment$seq[[1]])
  # look at each column of the alignment in turn:
  for (i in 1:alignmentlen)
  {
    # see what percent of the letters in this column are non-gaps:
    nongap <- 0
    for (j in 1:numseqs)
    {
      seqj <- alignment$seq[[j]]
      letterij <- substr(seqj,i,i)
      if (letterij != "-") { nongap <- nongap + 1}
    }
    pcnongap <- (nongap*100)/numseqs
    # Only consider this column if at least minpcnongap % of the letters are not gaps:
    if (pcnongap >= minpcnongap)
    {
      # see what percent of the pairs of letters in this column are identical:
      numpairs <- 0; numid <- 0
      # find the letters in all of the sequences in this column:
      for (j in 1:(numseqs-1))
      {
        seqj <- alignment$seq[[j]]
        letterij <- substr(seqj,i,i)
        for (k in (j+1):numseqs)
        {
          seqk <- alignment$seq[[k]]
          letterkj <- substr(seqk,i,i)
          if (letterij != "-" && letterkj != "-")
          {
            numpairs <- numpairs + 1
            if (letterij == letterkj) { numid <- numid + 1}
          }
        }
      }
      pcid <- (numid*100)/(numpairs)
      # Only consider this column if at least %minpcid of the pairs of letters are identical:
      if (pcid >= minpcid)
      {
        for (j in 1:numseqs)
        {
          seqj <- alignment$seq[[j]]
          letterij <- substr(seqj,i,i)
          newalignmentj <- newalignment$seq[[j]]
          newalignmentj <- paste(newalignmentj,letterij,sep="")
          newalignment$seq[[j]] <- newalignmentj
        }
      }
    }
  }
  return(newalignment)
}


cleaned_NAC_aln <- cleanAlignment(NAC_aln, 10, 0)

printMultipleAlignment(cleaned_NAC_aln, 60)

write.fasta(cleaned_NAC_aln$seq, cleaned_NAC_aln$nam, "wheat_barley_rice_NAC_msa_NAC_domain_no_numbers_10perc.fa", open="w" )

cleaned_NAC_aln50 <- cleanAlignment(NAC_aln, 10, 0)
write.fasta(cleaned_NAC_aln50$seq, cleaned_NAC_aln50$nam, "wheat_barley_rice_NAC_msa_NAC_domain_no_numbers_50perc.fa", open="w" )


# now do the trimming for wheat only NAC domains
wheat_NAC_aln <- read.alignment(file = "wheat_NAC_msa_NAC_domain_no_numbers_shortID.fas", format = "fasta") # read in alignment of NACs
wheat_NAC_aln$seq

cleaned_wheat_NAC_aln10 <- cleanAlignment(wheat_NAC_aln, 10, 0)
write.fasta(cleaned_wheat_NAC_aln10$seq, cleaned_wheat_NAC_aln10$nam, "wheat_NAC_msa_NAC_domain_no_numbers_10perc_shortID.fa", open="w" )


