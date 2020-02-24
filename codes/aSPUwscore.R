#' Doubley Adaptive Sum of Powered Score tests (daSPU) test for single trait with multiple weights for each genetic marker
#'
#' It returns p-values and test statistics
#'
#' @param U Score vector for the marker set we are interested in. (N by 1 matrix)
#'
#' @param V Corresponding covariance matrix for the score vector. (N by N matrix)
#'
#' @param weight Multiple weights for each genetic markers. (N by M matrix)
#'
#' @param pow1 SNP or single CpG sites specific power(gamma values) used in daSPU test.
#'
#' @param pow2 Specific power(gamma values) used for different weight.
#'
#' @param n.perm number of permutations.
#'
#' @export
#' @return P-values for SPUMpath tests and aSPUMpath test.
#'
#' @author Chong Wu and Wei Pan

waSPU <- function(U,V,weight, pow = c(1:8,Inf),n.perm = 1000,version = "V2") {
    
    pow1 = pow
    weight = as.matrix(weight)
    U = as.matrix(U)
    # remove SNPs corresponding to zero weight
    weight.tmp = abs(weight)
    index = rowSums(weight.tmp) > 0
    U = U[index,]
    V = V[index,]
    V = V[,index]
    weight = weight[index,]
    weight = as.matrix(weight)
    U = as.matrix(U)
    V = as.matrix(V)
    
    if(dim(weight)[2] == 1) {
        sim.aSPUw3(U,V,weight,pow1,n.perm)
    } else {
        if(version == "V1") {
            sim.aSPUw4(U,V,weight,pow1,n.perm)
        } else if(version =="V2") {
            sim.aSPU_pan(U,V,weight,pow1,n.perm)
        } else {
            sim.aSPUwV3(U,V,weight,pow1,n.perm)
        }
    }
}

# New, Dr. Pan wants to check if the TWAS-onmibus is right or not.

combined.score <- function(U,V,weight){
    weight = as.matrix(weight)
    cov = t(weight) %*% V %*% weight
    
    df = qr(cov)$rank
    T = t(t(weight) %*% U ) %*% solve(cov) %*% t(weight) %*% U
    
    pvalue = 1 - pchisq(T, df = df)
    return(list(T = T, pvs = pvalue))
}




sim.aSPUw3 <- function(U,V,weight, pow=c(1:8,Inf),n.perm = 1000){
    
    Ts <- rep(0, length(pow))
    for(j in 1:length(pow)){
        if (pow[j] < Inf) {
            Ts[j] <- sum((U * weight)^pow[j])
        } else {
            Ts[j] <- max(abs(weight * U))
        }
    }
    
    eV <- eigen(V)
    eV$values[eV$values<0] = 0
    
    CovSsqrt <- t(eV$vectors %*% (t(eV$vectors) * sqrt(eV$values)))
    pow[pow==Inf] <- 0 # pass 0 as infitiy
    n.pow = length(pow)
    T0s <- big.matrix(n.perm,n.pow,type = "double")
    
    Ts.abs <- abs(Ts)
    
    Res = calcT0Wsim3(as.matrix(CovSsqrt), as.matrix(weight),as.matrix(pow),as.matrix(Ts.abs),T0s@address,n.perm)
    
    minp0 = Res$minp0
    pPerm0 = Res$pPerm0
    
    Paspu <- (sum(minp0 <= min(pPerm0)) + 1) / (n.perm+1)
    pvs <- c(pPerm0, Paspu)
    
    Ts <- c(Ts, min(pPerm0))
    pow[pow==0] <- Inf
    names(Ts) <- c(paste("SPU(", pow,")", sep=""), "aSPU")
    names(pvs) <- names(Ts)
    
    return(pvs)
}

sim.aSPUw4 <- function(U,V,weight, pow=c(1:8,Inf),n.perm = 1000){
    
    npow = length(pow)
    nweight = dim(weight)[2]
    Ts <- rep(0,npow * nweight)
    
    for(i in 1:nweight) {
        for(j in 1:npow) {
            if (pow[j] < Inf) {
                Ts[(j-1) * nweight + i] <- sum((U * weight[,i])^pow[j])
            } else {
                Ts[(j-1) * nweight + i] <- max(abs(U * weight[,i]))
            }
        }
    }
    
    eV <- eigen(V)
    eV$values[ eV$values < 0 ] = 0
    
    CovSsqrt <- t(eV$vectors %*% (t(eV$vectors) * sqrt(eV$values)))
    pow[pow==Inf] <- 0 # pass 0 as infitiy
    
    T0s.size = npow * nweight
    T0s <- big.matrix(n.perm,T0s.size,type = "double")
    T1s <- big.matrix(n.perm,npow,type = "double")
    
    
    minp0_sign <- big.matrix(n.perm,nweight,type = "double")
    
    cov.res = calcT0Wsim4(as.matrix(CovSsqrt), as.matrix(weight),as.matrix(pow),T0s@address,minp0_sign@address, n.perm)
    
    
    #may change the first part into analytical form
    cov.1 = t(weight) %*% V %*% weight
    cov.res[1:nweight,1:nweight] = cov.1
    
    # calculate the final test statistics and its corresponding distribution
    calc_test_ch(cov.res,as.matrix(weight),as.matrix(pow),T0s@address,T1s@address,n.perm)
    
    Ts.chong = calc_test_ts(cov.res, as.matrix(weight), as.matrix(pow), as.matrix(Ts))
    
    Res = calc_p_ch(as.matrix(pow),as.matrix(Ts.chong),T1s@address, n.perm)
    
    minp0 = Res$minp0
    pPerm0 = Res$pPerm0
    
    Paspu <- (sum(minp0 <= min(pPerm0)) + 1) / (n.perm+1)
    pvs <- c(pPerm0, Paspu)
    
    pow[pow==0] <- Inf
    names(pvs) <- c(paste("SPU(", pow,")", sep=""), "aSPU")
    
    return(pvs)
}





sim.aSPUwV3 <- function(U,V,weight, pow=c(1:8,Inf),n.perm = 1000){
    
    npow = length(pow)
    nweight = dim(weight)[2]
    Ts <- rep(0,npow * nweight)
    
    for(i in 1:nweight) {
        for(j in 1:npow) {
            if (pow[j] < Inf) {
                Ts[(j-1) * nweight + i] <- sum((U * weight[,i])^pow[j])
            } else {
                Ts[(j-1) * nweight + i] <- max(abs(U * weight[,i]))
            }
        }
    }
    
    eV <- eigen(V)
    eV$values[ eV$values < 0 ] = 0
    
    CovSsqrt <- t(eV$vectors %*% (t(eV$vectors) * sqrt(eV$values)))
    pow[pow==Inf] <- 0 # pass 0 as infitiy
    
    T0s.size = npow * nweight
    T0s <- big.matrix(n.perm,T0s.size,type = "double")
    T1s <- big.matrix(n.perm,npow,type = "double")
    
    
    minp0_sign <- big.matrix(n.perm,nweight,type = "double")
    
    res.tmp = calcT0WsimV3(as.matrix(CovSsqrt), as.matrix(weight),as.matrix(pow),T0s@address,minp0_sign@address, n.perm)
    
    
    cov.res = res.tmp$cov
    
    mean.res = res.tmp$mean
    
    #may change the first part into analytical form
    cov.1 = t(weight) %*% V %*% weight
    cov.res[1:nweight,1:nweight] = cov.1
    
    # calculate the final test statistics and its corresponding distribution
    calc_test_ch_V3(cov.res,as.vector(mean.res), as.matrix(weight),as.matrix(pow),T0s@address,T1s@address,n.perm)
    
    Ts.chong = calc_test_ts_V3(cov.res, as.vector(mean.res),as.matrix(weight), as.matrix(pow), as.matrix(Ts))
    
    Res = calc_p_ch(as.matrix(pow),as.matrix(Ts.chong),T1s@address, n.perm)
    
    minp0 = Res$minp0
    pPerm0 = Res$pPerm0
    
    Paspu <- (sum(minp0 <= min(pPerm0)) + 1) / (n.perm+1)
    pvs <- c(pPerm0, Paspu)
    
    pow[pow==0] <- Inf
    names(pvs) <- c(paste("SPU(", pow,")", sep=""), "aSPU")
    
    return(pvs)
}



sim.aSPU_pan <- function(U,V,weight, pow=c(1:8,Inf),n.perm = 1000){
    
    npow = length(pow)
    nweight = dim(weight)[2]
    Ts <- rep(0,npow * nweight)
    
    for(i in 1:nweight) {
        for(j in 1:npow) {
            if (pow[j] < Inf) {
                Ts[(j-1) * nweight + i] <- sum((U * weight[,i])^pow[j])
            } else {
                Ts[(j-1) * nweight + i] <- max(abs(U * weight[,i]))
            }
        }
    }
    
    eV <- eigen(V)
    eV$values[ eV$values < 0 ] = 0
    
    CovSsqrt <- t(eV$vectors %*% (t(eV$vectors) * sqrt(eV$values)))
    pow[pow==Inf] <- 0 # pass 0 as infitiy
    
    T0s.size = npow * nweight
    T0s <- big.matrix(n.perm,T0s.size,type = "double")
    
    T1s <- big.matrix(n.perm,npow,type = "double")
    
    
    minp0_sign <- big.matrix(n.perm,nweight,type = "double")
    
    cov.res = calcT0Wsim4(as.matrix(CovSsqrt), as.matrix(weight),as.matrix(pow),T0s@address,minp0_sign@address, n.perm)
    
    ## Dr.Pan idea
    # calculate test statistics and its corresponding distribution
    Ts.abs <- abs(Ts)
    
    
    T2s <- big.matrix(n.perm,T0s.size,type = "double")
    
    Res = calc_p_pan(as.matrix(Ts.abs),T0s@address,T2s@address, minp0_sign@address,T0s.size,nweight, n.perm)
    
    Ts2 = Res$pPerm0
    cov.res = Res$cov
    Ts2 = qnorm(1- Ts2/2) * rep(sign(Ts[1:nweight]),npow)
    
    
    rm(T0s)
    calc_test_pan(cov.res,as.matrix(weight),as.matrix(pow),T1s@address,T2s@address, n.perm)
    
    Ts.chong = calc_test_ts(cov.res, as.matrix(weight), as.matrix(pow), as.matrix(Ts2))
    
    Res = calc_p_ch_pan(as.matrix(pow),as.matrix(Ts.chong),T1s@address, n.perm)
    
    minp0 = Res$minp0
    pPerm0 = Res$pPerm0
    
    Paspu <- (sum(minp0 <= min(pPerm0)) + 1) / (n.perm+1)
    pvs <- c(pPerm0, Paspu)
    pow[pow==0] <- Inf
    
    names(pvs) <- c(paste("SPU(", pow,")", sep=""), "aSPU")
    
    return(pvs)
}
