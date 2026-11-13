r <- readRDS("sens_spec_test/results.rds"); res <- r$res
cn <- c("orig99","orig95","fwd99","fwd95","fwd95_d0")
g <- paste(substr(r$groups$MCH,1,4), r$groups$group)
tab <- function(par, col="median", f=1) { m <- sapply(res, function(x) {s <- x$summary; s[[col]][startsWith(s$variable, paste0(par,"["))]*f}); colnames(m) <- cn; rownames(m) <- g; m }
cat("\nIFR median per 10,000\n"); print(round(tab("p", f=1e4),2))
cat("\np_HEV_d median (orig = apparent, fwd = true)\n"); print(round(tab("p_HEV_d"),3))
cat("\nexpected community HEV deaths total\n"); print(round(colSums(tab("deaths_exp_phevd","mean")),1))
k <- c("alpha_mu","alpha_sigma","alpha_mu_d","alpha_sigma_d")
m <- sapply(res, function(x){s<-x$summary; i <- match(k, s$variable); sprintf("%.2f (%.2f, %.2f)", s$median[i], s$q2.5[i], s$q97.5[i])}); rownames(m) <- k; colnames(m) <- cn; print(m, quote=FALSE)
