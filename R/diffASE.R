
diffASE <- function(a1_counts, tot_counts, estimates, min_counts = 0, min_cells = 5, metadata, cores = NULL){
    
    require(doParallel)
  
    lbetabin <- function(df, glob_params, mu, theta){
      
      min_theta=1e-06
      theta <- pmax(theta, min_theta)
      
      y <- df[,1]
      n <- df[,2]
      
      alpha <- mu/theta
      beta <- (1-mu)/theta
      
      sum(lchoose(n, y) + lgamma(alpha+beta) - lgamma(n+alpha+beta) + 
            lgamma(y+alpha) - lgamma(alpha) + lgamma(n-y+beta) - lgamma(beta))
      
    }
    
    lbetabin_alphabeta <- function(df, inits){
      
      y <- df[,1]
      n <- df[,2]
      
      #mu = inits[1]
      #theta = inits[2]
      
      alpha = inits[1]
      beta = inits[2]
      
      #min_theta=1e-06
      #theta <- pmax(theta, min_theta)
      
      #alpha = mu/theta
      #beta = (1-mu)/theta
      
      sum(lchoose(n, y) + lgamma(alpha+beta) - lgamma(n+alpha+beta) + 
            lgamma(y+alpha) - lgamma(alpha) + lgamma(n-y+beta) - lgamma(beta))
      
    }
    
    
    # A function to rbind the data frames in the list, leaving out one at a time
    rbind_consec <- function(df_list, idx) {
      df_to_rbind <- df_list[-idx]  # Leave out the data frame at the specified index
      do.call(rbind, df_to_rbind)  # Use do.call to rbind all data frames in the list
    }
    
    groups <- split(metadata, f = metadata$cell_type)
    
    res <- list()
    
  
    if (!is.null(cores)) {
      
      cl <- parallel::makePSOCKcluster(cores)
      registerDoParallel(cl)
    
  
      tmp <-  foreach::foreach(k = 1:nrow(a1_counts), .combine = rbind, .multicombine = F) %dopar%  {   
      #for  (k in 1:nrow(a1_counts)) {
      
      y <-  a1_counts[k,]
      n <- tot_counts[k,]
      a2 <- n - y
      
      df <- data.frame(y = y, n = n)
      df <- na.omit(df)
      #df <- df[df$n >= min_counts,]
      df_subset <- df[df$n >= min_counts,]
      
      #log2FC[k] = log2(mean(y)) - log2(mean(a2))
      #AR[k] = mean(y/n, na.rm = T)
      #N[k] = dim(df[df$n >= 5,])[1]
      #tot_gene_mean[k] = mean(df$n)
      
      if (dim(df)[1] >= min_cells){
        
        #obtaining null likelihood - allelic imbalance is the same across the cell types
        #without dispersion correction
        mu <- estimates[k, "mean_reestim"]
        theta <- estimates[k, "theta_reestim"]
        nul.lik_orig <- lbetabin(df, mu = mu, theta = theta)
        
        #optim.nul <- tryCatch(optim(par = c(mu, theta), lbetabin,
        #                            hessian=T, df = df, method = "L-BFGS-B",
        #                            lower = c(1e-8), upper=c(1e6),
        #                            control = list( fnscale=-1 )), error=function(e) e)
        
        #with dispersion correction
        theta_adj <- estimates[k, "thetaCorrected"]
        nul.lik_adj <- lbetabin(df, mu = mu, theta = theta_adj)
        
        #this dataset contains counts for each individual cell type
        df_split <- lapply(groups, function(q) df_subset[rownames(df_subset) %in% rownames(q),])
        
        #this dataset contains counts for the cell type excluding one of interest
        # Loop through the indices of the data frames in the list
        df_oneout <- lapply(seq_along(df_split), function(i) rbind_consec(df_split, i))
        
        #fitting glm to obtain initial values
        celltype.binom <- lapply(df_split, function(q) tryCatch(
                                                        glm(y/n ~ 1, family = "binomial", 
                                                            weights = n, data = q), error=function(e) e))
        celltype.inits <- lapply(celltype.binom, function(q) c(q$fitted.values[1],
                                                               1-q$fitted.values[1]))
        
        #for each cell type obtaining parameters alpha and beta
        celltype.optim <- mapply(function(p, q) tryCatch(optim(par = p, lbetabin_alphabeta,
                                                         hessian=T, df = q, method = "L-BFGS-B",
                                                         lower = c(1e-8), upper=c(1e6),
                                                         control = list( fnscale=-1 )), error=function(e) NA),
                                 celltype.inits, df_split, SIMPLIFY = F)
        #celltype.optim <- lapply(df_split, function(q) tryCatch(optim(par = c(mu,theta), lbetabin_alphabeta,
        #                                                           hessian=T, df = q, method = "L-BFGS-B",
        #                                                           lower = c(1e-8), upper=c(1e6),
        #                                                           control = list( fnscale=-1 )), error=function(e) e))
        
        #they are used to calculate mean AR
        celltype.pars <- lapply(celltype.optim, function(q) q[[1]])
        celltype.means <- lapply(celltype.pars, function(q) q[1]/(q[1] + q[2]))
        #using mean AR and dispersion estimated across all cells obtain likelihood values
        celltype.lik_orig <- mapply(function(p,q) lbetabin(p, mu = q, theta = theta),
                                    df_split, celltype.means, SIMPLIFY = F)
        #adjusted likelihood values are obtained with shrunk dispersion 
        celltype.lik_adj <- mapply(function(p, q) lbetabin(p, mu = q, theta = theta_adj),
                                   df_split, celltype.means, SIMPLIFY = F)
        
        #fitting glm to obtain initial values for all other cell types
        rest.binom <- lapply(df_oneout, function(q) tryCatch(
          glm(y/n ~ 1, family = "binomial", 
              weights = n, data = q), error=function(e) e))
        rest.inits <- lapply(rest.binom, function(q) c(q$fitted.values[1],
                                                      1-q$fitted.values[1]))
        #for each cell type obtaining parameters alpha and beta
        rest.optim <- mapply(function(p, q) tryCatch(optim(par = p, lbetabin_alphabeta,
                                                           hessian=T, df = q, method = "L-BFGS-B",
                                                           lower = c(1e-8), upper=c(1e6),
                                                           control = list( fnscale=-1 )), error=function(e) NA),
                                 rest.inits, df_oneout, SIMPLIFY = F)
        
        #repeating the same for other cell types
        #rest.optim <- lapply(df_oneout, function(q) tryCatch(optim(par = c(mu,theta), lbetabin_alphabeta,
        #                                                              hessian=T, df = q, method = "L-BFGS-B",
        #                                                              lower = c(1e-8), upper=c(1e6),
        #                                                              control = list( fnscale=-1 )), error=function(e) e))
        rest.pars <- tryCatch(lapply(rest.optim, function(q) q[[1]]), error=function(e) NA)
        rest.means <- tryCatch(lapply(rest.pars, function(q) q[1]/(q[1] + q[2])), error=function(e) NA)
        
        rest.lik_orig <- mapply(function(p,q) lbetabin(p, mu = q, theta = theta),
                                df_oneout, rest.means, SIMPLIFY = F)
        rest.lik_adj <- mapply(function(p,q) lbetabin(p, mu = q, theta = theta_adj),
                               df_oneout, rest.means, SIMPLIFY = F)
        
        #rest.lik_orig <- lapply(df_oneout, function(q) lbetabin(q, mu = mu, theta = theta))
        #rest.lik_adj <- lapply(df_oneout, function(q) lbetabin(q, mu = mu, theta = theta_adj))
                               
        #obtaining alternative likelihood - allelic imbalance is not the same between the cell types
        #summing the values
        alt.lik_orig <- Map(sum, celltype.lik_orig, rest.lik_orig)
        alt.lik_adj <- Map(sum, celltype.lik_adj, rest.lik_adj)
        
        llr_orig <- mapply(function(p,q) p - q, nul.lik_orig, alt.lik_orig, SIMPLIFY = F)
        pval_orig <- lapply(llr_orig, function(q)
                            pchisq(-2*(q), df = 1, lower.tail = FALSE))
        
        llr_adj <- mapply(function(p,q) p - q, nul.lik_adj, alt.lik_adj, SIMPLIFY = F)
        pval_adj <- lapply(llr_adj, function(q)
          pchisq(-2*(q), df = 1, lower.tail = FALSE))
        
        master_list <- Map(cbind, celltype.means, rest.means, nul.lik_orig, alt.lik_orig, llr_orig, pval_orig,
                   nul.lik_adj, alt.lik_adj, llr_adj, pval_adj)
        
        master_list <- lapply(master_list, function(p) {lapply(p, as.data.frame)})  
        final <- do.call(Map, c(f = rbind, master_list))
        
        final <- lapply(final, function(q) {colnames(q) <- c("mean_ct", "mean_other", "null_lik_orig", "alt_lik_orig",
                                                             "orig_llr", "orig_pval", "null_lik_adj", "alt_lik_adj",
                                                             "adj_llr", "adj_pval");
        return(q)})
        
        final <- lapply(final, function(q) {rownames(q) <- rownames(a1_counts); return(q)})
        
        tmp <- final
        #for (i in 1:length(res)){
        #  for (j in 1:length(res[[i]]))
        #  res[[i]][[j]] <- as.data.frame(res[[i]][[j]])
        #}
        
      } else {
        
        tmp <- NULL
      
      }
        
      }
      
      return(tmp)
      stopCluster(cl)   
  }   
      

}
