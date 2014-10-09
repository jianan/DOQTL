################################################################################
# Write out the smoothed probabilities for each sample and the founder state
# means and variances. This will write out a file for all samples on the 
# autosomes.  On the X Chr, you must specify which sex is being written out.
# 
write.results = function(prsmth, theta.rho.means, theta.rho.covars, b, output.dir,
                         chr, all.chr, sex, write.gp36=FALSE, write.b=FALSE) {
  if(!missing(prsmth)) {
    dimnames(prsmth)[[2]] = make.names(dimnames(prsmth)[[2]])
    prsmth = exp(prsmth) ## in *.genotype.probs.txt, values are log-ed.
    if(write.gp36){
      save(prsmth, file=paste0(output.dir, "Chr", chr, ".founder.probs.36.Rdata"))
    }
    
    # Create a matrix that we can use to multiply the 36 state probabilities.
    spl = strsplit(dimnames(prsmth)[[1]], split = "")
    names(spl) = dimnames(prsmth)[[1]]
    spl = lapply(spl, factor, levels = sort(unique(unlist(spl))))
    spl = lapply(spl, table)
    mat = matrix(unlist(spl), length(spl[[1]]), length(spl), dimnames = 
        list(names(spl[[1]]), names(spl)))
    mat = mat * 0.5
    
    # Place the founder probabilities in the matrix.
    model.probs <- array(NA, dim=c(dim(prsmth)[2], 8, dim(prsmth)[3]))
    for(i in 1:dim(prsmth)[2]){
      model.probs[i,,] = mat %*% prsmth[,i,]
    }
    dimnames(model.probs)[[2]] <- LETTERS[1:8]
    dimnames(model.probs)[c(1,3)] <- dimnames(prsmth)[c(2,3)]
    save(model.probs, file=paste0(output.dir, "Chr", chr, ".founder.probs.Rdata"))

  } # if(!missing(prsmth))
  # Write out the genotype mean and variance estimates.
  if(!missing(theta.rho.means)) {  
    if(chr != "X") {
      # Write out the genotype mean and variance estimates.
      save(theta.rho.means, file = paste(output.dir, "chr", chr,
           ".final.means.Rdata", sep = ""))
      save(theta.rho.covars, file = paste(output.dir, "chr", chr,
           ".final.covars.Rdata", sep = ""))
    } else {
    
      if(is.null(sex)) {
        stop("write.results: sex cannot be null in write.results if chr = X.")
      } # if(is.null(sex))
    
      if(sex == "M") {
        sex = "male"
      } else {
        sex = "female"
      } # else
    
      save(theta.rho.means, file = paste(output.dir, "chr", chr, ".", sex,
           ".final.means.Rdata", sep = ""))
      save(theta.rho.covars, file = paste(output.dir, "chr", chr, ".", sex,
           ".final.covars.Rdata", sep = ""))
    } # else
  } else if(!missing(b) & write.b) {
    if(chr != "X") {
      # Write out the genotype mean and variance estimates.
      save(b, file = paste(output.dir, "chr", chr, ".emission.probs.Rdata", 
	       sep = ""))
    } else {
    
      if(is.null(sex)) {
        stop("write.results: sex cannot be null in write.results if chr = X.")
      } # if(is.null(sex))
    
      if(sex == "M") {
        sex = "male"
      } else {
        sex = "female"
      } # else
    
      save(b, file = paste(output.dir, "chr", chr, ".", sex, 
	       ".emission.probs.Rdata", sep = ""))
    } # else    
  } # else if(!missing(b))
  
} # write.results()
