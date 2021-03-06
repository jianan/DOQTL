################################################################################
# Main function for running the haplotype HMM for DO mice.
# You must supply either a set of founders with your data or a set of 
# founder means and variances.  This means that either is.founder.F1 must have
# true values for at least one sample from each founder or init.means and 
# init.covars must be populated with state mean and covariance estimates for all
# SNPs.
# Arguments: data: list, containing geno, sex and gen.
#            chr: character vector, with chromosomes to run.  Must match
#                 the chromosome IDs in the snps table.  NULL means
#                 run all.
#            founders: list, containing genotype calls for the founders
#                      and F1s.
#            snps: data.frame, with 4 columns containing SNPs to use. 
#                  SNP IDs, chr, Mb position, cM position in columns 1:4.
#            output.dir: character, with the file directory for the final
#                        theta & rho mean and variance files.
#            trans.prob: list, containing transition probabilities between
#                        each SNP. Only required for custom sample types.
#            plot: boolean that is true if the user would like to plot a sample
#                  chromosome as the model progresses.
calc.genoprob.alleles = function(data, chr, founders, snps, output.dir = ".",
                        trans.prob.fxn = do.trans.probs, plot = FALSE,
                        write.gp36=FALSE, write.b=FALSE) {
  # No NAs in geno.
  if(any(is.na(data$geno))) {
    data$geno[is.na(data$geno)] = "N"
    message("Changing NA values in genotypes to 'N'.")
  } # if(any(is.na(data$geno))
  # Select the chromosomes to run.  If chr.to.run = "all", then run all chrs.
  if(chr[1] == "all") {
    chr = names(snps)
  } # else
  # Extract the sample and founder temporary data file names.
  tmpfiles = list(dg = data$geno, fg = founders$geno)
  # Loop through each chromosome.
  for(curr.chr in chr) {
    print(paste("CHR", curr.chr))
    # SNPs on the current chromosome.
    cur.snps = snps[[which(names(snps) == curr.chr)]]
    # Read in the sample and founder data.
    chr.pattern = paste("chr", curr.chr, "\\.", sep = "")
	d = NULL # d is the variable that is loaded in the load() statements below.
    load(tmpfiles$dg[grep(chr.pattern, tmpfiles$dg)])
    data$geno = d
    load(tmpfiles$fg[grep(chr.pattern, tmpfiles$fg)])
    founders$geno = d
    rm(d)
    gc()
    # If this is the X chromosome, split the samples by sex, using 36 states
    # for the females and 8 states for the males.
    if(curr.chr == "X") {
      # Run the females first.
      females = which(data$sex == "F")
      female.prsmth = NULL
      female.b = NULL
      # Only run if there are samples that are female.
      if(length(females) > 0) {
        print("Females")
        cur.data = list(geno = data$geno[females,], 
                   sex = data$sex[females], gen = data$gen[females])
        attributes(cur.data) = attributes(data)
        founder.subset = which(founders$sex == "F")
        cur.founders = list(geno = founders$geno[founder.subset,],
                       sex = founders$sex[founder.subset],
                       code = founders$code[founder.subset],
                       states = founders$states$X$F)
        tmp = hmm.allele(data = cur.data, founders = cur.founders,
              sex = "F", snps = cur.snps, chr = curr.chr, 
              trans.prob.fxn = trans.prob.fxn)
        female.prsmth = tmp$prsmth
        female.b = tmp$b
        rm(tmp)
      } # if(length(females) > 0)
      # Run the males, which only have founder states, not F1s because the 
      # males are hemizygous.
      males = which(data$sex == "M")
      male.prsmth = NULL
      male.b = NULL
      # Only run if there are samples that are male.  If only founders are
      # male, there's no point in running this.
      if(length(males) > 0) {
        print("Males")
        cur.data = list(geno = data$geno[males,], 
                   sex = data$sex[males], gen = data$gen[males])
        attributes(cur.data) = attributes(data)
        founder.subset = which(founders$sex == "M")
        cur.founders = list(geno = founders$geno[founder.subset,],
                       sex = founders$sex[founder.subset],
                       code = founders$code[founder.subset],
                       states = founders$states$X$M)
        cur.founders = keep.homozygotes(cur.founders)
        tmp = hmm.allele(data = cur.data, founders = cur.founders,
              sex = "M", snps = cur.snps, chr = curr.chr, 
              trans.prob.fxn = trans.prob.fxn)
        male.prsmth = tmp$prsmth
        male.b = tmp$b
        rm(tmp)
      } # if(length(males) > 0)
	  
      # Combine the male and female prsmth data.
      if(length(females) > 0) {
        if(length(males) > 0) {
          prsmth = array(-.Machine$double.xmax, c(length(founders$states$X$F), 
                   nrow(data$geno), nrow(cur.snps)), dimnames = 
                   list(founders$states$X$F, rownames(data$geno), 
                   cur.snps[,1]))
          m = match(dimnames(female.prsmth)[[2]], dimnames(prsmth)[[2]])
          prsmth[,m,] = female.prsmth
          m = match(dimnames(male.prsmth)[[2]], dimnames(prsmth)[[2]])
          m2 = match(paste(dimnames(male.prsmth)[[1]], dimnames(male.prsmth)[[1]], sep = ""),
               dimnames(prsmth)[[1]])
          prsmth[m2,m,] = male.prsmth
          rm(female.prsmth, male.prsmth)
        } else {
         prsmth = female.prsmth
         rm(female.prsmth)
        } # else
      } else if(length(males) > 0) {
         prsmth = male.prsmth
         rm(male.prsmth)
      } # else if(length(males) > 0)
      
      # Write out the smoothed probabilities and the emission probabilities.
      write.results(prsmth = prsmth, b = female.b, output.dir = output.dir, chr = curr.chr,
                    all.chr = chr, sex = "F", write.gp36=write.gp36, write.b=write.b)
      write.results(b = male.b, output.dir = output.dir, chr = curr.chr,
                    all.chr = chr, sex = "M", write.gp36=write.gp36, write.b=write.b)
					
    } else if (curr.chr == "Y") {
      stop("Chr Y not implemented yet.")
    } else if (curr.chr == "M") {
      stop("Chr M not implemented yet.")
    } else {
      # Autosomes.
      data$geno[data$geno == "-"] = "N"
      cur.founders = founders
      cur.founders$states = cur.founders$states$auto
      tmp = hmm.allele(data = data, founders = cur.founders, sex = "F", 
            snps = cur.snps, chr = curr.chr, trans.prob.fxn = trans.prob.fxn)
      prsmth = tmp$prob.mats$prsmth
      b = tmp$b
      # Write out the smoothed probabilities and the founder state meand and 
      # variances.
      write.results(prsmth = tmp$prsmth, b = tmp$b, output.dir = output.dir, 
	                chr = curr.chr, all.chr = chr, write.gp36=write.gp36,
                        write.b=write.b)
      rm(tmp)
    } # else
    gc()
	
  } # for(chr)
} # calc.genoprob.alleles()
