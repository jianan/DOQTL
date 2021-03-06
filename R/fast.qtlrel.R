################################################################################
# A faster impelmentation of the QTLRel algorithm for additive covariates with
# a kinship matrix.  For interactive covariates, use the main QTLRel function.
# This is only about twice as fast as QTLRel.
# Code adapted from Riyan Cheng's QTLRel package on CRAN.
# NOTE: No NAs are allowed anywhere! Filter your data before calling this.
# Daniel Gatti
# Dan.Gatti@jax.org
# Aug, 2, 2013
################################################################################
# Arguments: pheno: numeric vector with no NAs containing the phenotype.
#            probs: 3D numeric array of haplotype probabilities.
#            K: numeric matrix.
fast.qtlrel = function(pheno, probs, K, addcovar, snps) {

  # If the SNPs are in bp, rescale them to Mb.
  if(max(snps[,3], na.rm = TRUE) > 200) {
    snps[,3] = snps[,3] * 1e-6
  } # if(max(snps[,3]) < 200)

  prdat = list(pr = probs, chr = snps[,2], dist = snps[,3], snp = snps[,1])
  class(prdat) = c(class(prdat), "addEff")
  err.cov = NULL

  if(!missing(K)) {

    K = as.matrix(K)
    vTmp = list(AA = 2 * K, DD = NULL, HH = NULL, AD = NULL, MH = NULL,
                EE = diag(nrow(K)))
    vc = NULL
    if(missing(addcovar)) {
      vc = estVC(y = pheno, v = vTmp)
    } else {
      vc = estVC(y = pheno, x = addcovar, v = vTmp)
    } # else

    err.cov = matrix(0, nrow(K), ncol(K))
    for(j in which(vc$nnl)) {
      err.cov = err.cov + vTmp[[j]] * vc$par[names(vTmp)[j]]
    } # for(j)

    rm(vTmp)

    # Invert the covariance matrix.
    eW = eigen(err.cov, symmetric = TRUE)
    if (min(eW$values) < 0 && abs(min(eW$values)) > sqrt(.Machine$double.eps)) {
      stop("fast.qtlrel: W is not positive definite")
    } else {
      eW$values[eW$values <= 0] = Inf
    } # else
    err.cov = eW$vector %*% diag(eW$values^-0.5) %*% t(eW$vector)
    rm(eW)

  } # if(!missing(K))

  # Add the intercept.
  if(missing(addcovar)) {
    addcovar = matrix(1, nrow = length(pheno), ncol = 1, dimnames =
               list(rownames(pheno), "Intercept"))
  } else {
    addcovar = as.matrix(cbind(rep(1, length(pheno)), addcovar))
    colnames(addcovar)[1] = "Intercept"
  } # else

  # Remove A as the basis.
  probs = probs[,-1,]

  # Null model.
  ss.null = 0
  if(!is.null(err.cov)) {
    ytmp = err.cov %*% pheno
    xtmp = err.cov %*% addcovar
    qr.null = qr(xtmp)
    ss.null = sum(qr.resid(qr.null, ytmp)^2)
  } else {
    qr.null = qr(addcovar)
    ss.null = sum(qr.resid(qr.null, pheno)^2)
  } # else

  # Additive model for all SNPs.
  addx = cbind(addcovar, probs[,,1])
  ss = rep(0, nrow(snps))
  coef = matrix(0, nrow(snps), ncol(addx), dimnames = list(snps[,1],
         colnames(addx)))
  rng = (ncol(addcovar)+1):ncol(addx)

  perc.var = 0
  lrs = 0
  lod = 0

  if(!is.null(err.cov)) {

    for(s in 1:nrow(snps)) {
      addx[,rng] = probs[,,s]
      xtmp = err.cov %*% addx
      qr.add = qr(xtmp)
      ss[s] = sum(qr.resid(qr.add, ytmp)^2)
      coef[s,] = qr.coef(qr.add, ytmp)
    } # for(s)

  } else {

    for(s in 1:nrow(snps)) {
      addx[,rng] = probs[,,s]
      qr.add = qr(addx)
      ss[s] = sum(qr.resid(qr.add, pheno)^2)
      coef[s,] = qr.coef(qr.add, pheno)
    } # for(s)

  } # else

  perc.var = 100 * (1.0 - (ss / ss.null))
  lrs = -length(pheno) * log(ss / ss.null)
  lod = lrs / (2 * log(10))

  # Get the p-value from the LRS.
  p = pchisq(q = -length(pheno) * log(ss / ss.null), 
      df = ncol(addx) - ncol(addcovar), lower.tail = FALSE)

  return(list(lod = cbind(snps[,1:4], perc.var  = perc.var, lrs = lrs, 
         lod = lod, p = p, neg.log10.p = -log(p, 10)), coef = coef))

} # fast.qtlrel()

