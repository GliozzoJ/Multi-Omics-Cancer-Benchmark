# Functions to compute:
# - log-rank test
# - enrichment analysis on clinical variables (chi-squared test for 
# categorical variables, kruskal-wallis test for continuous)
# All p-values are estimated using permutation tests (i.e., we permuted the 
# cluster labels between samples and used the test statistic to obtain an 
# empirical p-value). 
# DISCLAIMER: the code is copied (with slight modifications to run it outside
# the benchmark workflow) from 
# https://github.com/Shamir-Lab/Multi-Omics-Cancer-Benchmark
# Further information about the code can be find in the reference paper and
# related supplementary:
# Nimrod Rappoport , Ron Shamir, Multi-omic and multi-view clustering algorithms: 
# review and cancer benchmark, Nucleic Acids Research, Volume 46, Issue 20, 
# 16 November 2018, Pages 10546–10562, https://doi.org/10.1093/nar/gky889
# August 2023


# Load libraries
library("survival")
library("parallel")


#' Test survival curve differences using survival package
#' 
#' Nothe that NA values for "Survival" and "Death" variables are set to zero.
#'
#' @param groups vector. Named vector with obtained clusterings.
#' @param survival.file.path string. Path with survival data (accepted table or
#' .rds format containing a dataframe with rownames). Data need to be samples x features.
#' You can also pass directly a dataframe to this variable.
#' @param keep.NA boolean. Do you want to keep NAs in survival time and event? (def. T).
#' If FALSE, NAs are substituted with 0. 
#' @param death.name string. Name of the column in data that has the survival event
#' as 1 = Death, 0 = Alive.
#' @param surv.name string. Name of the column in data that has the survival time.
#'
#' @return List with various elements, including:
#' - chisq: the chisquare statistic for a test of equality.
#' - var: the variance matrix of the test.
#' - pvalue: the p-value corresponding to the Chisquare statistic
#' @export
#'
#' @examples
check.survival <- function(groups, survival.file.path, keep.NA=T, death.name="cdr.os", 
                           surv.name="cdr.os.time") {
    
    # Read RDS or table
    if(is.character(survival.file.path)){
        
        ext <- strsplit(basename(survival.file.path), ".", fixed=T)[[1]][-1]
        if(ext == "rds"){
            survival.data = readRDS(survival.file.path)
            patient.names.in.file = rownames(survival.data)
        } else{
            survival.data = read.table(survival.file.path, header = TRUE)
            patient.names.in.file = as.character(survival.data[, 1])
        }
        
    } else{
        
        survival.data = survival.file.path
        patient.names.in.file = rownames(survival.data)
    }
    
    patient.names = names(groups) #clustering names
    patient.names.in.file = toupper(substring(patient.names.in.file, 1, 12)) # data names
    
    stopifnot(all(patient.names %in% patient.names.in.file))
    
    # reorder data to match clustering
    indices = match(patient.names, patient.names.in.file)
    ordered.survival.data = survival.data[indices,]
    
    # Add column with clustering to survival data
    ordered.survival.data["cluster"] <- groups
    
    # Retrieve index of column for survival time and event
    idx.death <- which(colnames(ordered.survival.data) == death.name)
    idx.surv <- which(colnames(ordered.survival.data) == surv.name)
    
    # In survival event and time, substitute NAs with 0
    if(!keep.NA){
        ordered.survival.data[,idx.surv][ is.na(ordered.survival.data[,idx.surv]) ] = 0
        ordered.survival.data[,idx.death][ is.na(ordered.survival.data[,idx.death]) ] = 0
    }
    
    # Compute survival analysis
    x = "cluster"
    y = paste("Surv(", surv.name,",", death.name,")" )
    form = as.formula(paste(y, "~", x))
    surv <- survdiff(form, data=ordered.survival.data)
    
    return(surv)
    
}


#' Compute empirical survival using a permutation test
#'
#' @param clustering vector. Named vector with obtained clusterings.
#' @param seed integer. Set seed for reproducibility. 
#' @param survival.file.path string. Path with survival data (accepted table or
#' .rds format containing a dataframe with rownames). Data need to be samples x features.
#' You can also pass directly a dataframe to this variable.
#' @param keep.NA boolean. Do you want to keep NAs in survival time and event? (def. T).
#' If FALSE, NAs are substituted with 0. 
#' @param death.name string. Name of the column in data that has the survival event
#' as 1 = Death, 0 = Alive.
#' @param surv.name string. Name of the column in data that has the survival time.
#' @param mc.cores integer. Number of cores to use (def. half of available cores).
#' @param verbose boolean. Do you want to print results?
#'
#' @return List with p-value, confidence interval and number of permutations
#' @export
get.empirical.surv <- function(clustering, seed=42, 
                               survival.file.path, keep.NA=T, 
                               death.name="cdr.os", 
                               surv.name="cdr.os.time", mc.cores=detectCores()/2,
                               verbose=F) {
    
    set.seed(seed, "L'Ecuyer-CMRG")
    surv.ret = check.survival(clustering, survival.file.path=survival.file.path, 
                              keep.NA=keep.NA, death.name=death.name, 
                              surv.name=surv.name) # Compute survival diff
    orig.chisq = surv.ret$chisq # Extract chi2
    orig.pvalue = get.logrank.pvalue(surv.ret) # compute log-rank p-value
    
    # The initial number of permutations to run
    num.perms = round(min(max(10 / orig.pvalue, 1000), 1e6))
    should.continue = T
    
    total.num.perms = 0
    total.num.extreme.chisq = 0
    
    while (should.continue) {
        if(verbose == T){
            print('Another iteration in empirical survival calculation')
            print(num.perms)
        }
        # Compute the chi2 statistic for each permutation
        perm.chisq = as.numeric(mclapply(1:num.perms, function(i) {
            cur.clustering = sample(clustering)
            names(cur.clustering) = names(clustering)
            cur.chisq = check.survival(cur.clustering, 
                                       survival.file.path=survival.file.path, 
                                       keep.NA=keep.NA, death.name=death.name, 
                                       surv.name=surv.name)$chisq
            return(cur.chisq)
        }, mc.cores=mc.cores))
        
        total.num.perms = total.num.perms + num.perms #update number of perms computed
        total.num.extreme.chisq = total.num.extreme.chisq + sum(perm.chisq >= orig.chisq) #number of extreme chisq
        
        binom.ret = binom.test(total.num.extreme.chisq, total.num.perms)  #binomial test
        cur.pvalue = binom.ret$estimate # the estimated probability of success
        cur.conf.int = binom.ret$conf.int # confidence interval for the probability of success
        
        if(verbose == T){
            print(c(total.num.extreme.chisq, total.num.perms))
            print(cur.pvalue)
            print(cur.conf.int)
        }
        
        sig.threshold = 0.05
        is.conf.small = ((cur.conf.int[2] - cur.pvalue) < min(cur.pvalue / 10, 0.01)) & ((cur.pvalue - cur.conf.int[1]) < min(cur.pvalue / 10, 0.01))
        is.threshold.in.conf = cur.conf.int[1] < sig.threshold & cur.conf.int[2] > sig.threshold
        #if ((is.conf.small & !is.threshold.in.conf) | (total.num.perms > 2e7)) {
        if ((!is.threshold.in.conf) | (total.num.perms > 5e6)) {
            should.continue = F
        } else {
            num.perms = 1e5 #1e5
        }
    }
    
    return(list(pvalue = cur.pvalue, conf.int = cur.conf.int, total.num.perms=total.num.perms, 
                total.num.extreme.chisq=total.num.extreme.chisq))
}


#' Enrichment analysis on clinical variables (discrete and continuous).
#'
#' @param clustering vector. Named vector with obtained clusterings.
#' @param clinical.data.path string. Path to clinical data (table, rds) or
#' directly a dataframe.
#' @param clinical.metadata list. A list indicating for each clinical variable
#' to evaluate if it is 'DISCRETE' or 'NUMERIC'. Following an example:
#'     clinical.metadata = list(gender='DISCRETE', age_at_initial_pathologic_diagnosis='NUMERIC', 
#'     pathologic_M='DISCRETE', pathologic_N='DISCRETE', pathologic_T='DISCRETE', 
#'     pathologic_stage='DISCRETE')
#' @param seed integer. Set seed for reproducibility.
#' @param mc.cores integer. Number of cores to use (def. half of available cores).
#'
#' @return List with p-values, confidence intervals and number of permutations
#' for each tested clinical variable.
#' @export
check.clinical.enrichment <- function(clustering, clinical.data.path, 
                                      clinical.metadata, seed=42, 
                                      mc.cores=detectCores()/2) {
    
    # Read RDS or table
    if(is.character(clinical.data.path)){
    
        ext <- strsplit(basename(clinical.data.path), ".", fixed=T)[[1]][-1]
        if(ext == "rds"){
            clinical.params = readRDS(clinical.data.path)
            
        } else {
            clinical.params = read.table(clinical.data.path, 
                                         header=T, stringsAsFactors = F)
        }
        
    } else {
        
        clinical.params = clinical.data.path
    }
    
    
    pvalues = c()
    total.num.iters = c()
    total.num.extreme = c()
    conf.int = list()
    
    params.being.tested = c()
    
    c=1 #counter
    for (clinical.param in names(clinical.metadata)) { # iter across clinical variables
        
        print(paste("Analysing clinical variable:", clinical.param))
        
        # Skip clinical variable that are not in the clinical.metadata given as
        # input
        if (!(clinical.param %in% colnames(clinical.params))) {
            print(paste0('WARNING: ', clinical.param, ' skipped from analysis (not in provided clinical.metadata.)'))
            next
        }
        
        # Extract clinical value with patients in the same order of clustering
        clinical.values = clinical.params[names(clustering), clinical.param]
        
        # Check if variable is discrete or numeric
        is.discrete.param = clinical.metadata[clinical.param] == 'DISCRETE' #boolean for discrete
        is.numeric.param = clinical.metadata[clinical.param] == 'NUMERIC' #boolean for numeric
        stopifnot(is.discrete.param | is.numeric.param)
        
        # skip parameter if many missing values
        # if a clinical value has too many NAs is skipped. It is skipped also
        # if all values are identical.
        if (is.numeric.param) {
            numeric.entries = !is.na(as.numeric(clinical.values))
            if (2 * sum(numeric.entries) < length(clinical.values)) {
                print(paste0('WARNING: skipping on ', clinical.param, ' because too many NAs: ', 
                             (length(clinical.values) - numeric.entries)))
                next
            }
        } else {
            not.na.entries = !is.na(clinical.values)
            should.skip = F
            if (2 * sum(not.na.entries) < length(clinical.values)) {
                should.skip = T
            } else if (length(table(clinical.values[not.na.entries])) == 1) {
                should.skip = T
            }
            
            if (should.skip) {
                print(paste0('WARNING: skipping on ', clinical.param, ' because too many NAs or only one value'))
                next
            }
        }
        
        params.being.tested = c(params.being.tested, clinical.param) # collect params tested
        
        if (is.discrete.param) { #chi-squared test
            #clustering.with.clinical = cbind(clustering, clinical.values)
            #tbl = table(as.data.frame(clustering.with.clinical[!is.na(clinical.values),]))
            #test.res = chisq.test(tbl)
            #pvalue = test.res$p.value
            res = get.empirical.clinical(clustering[!is.na(clinical.values)], 
                                            clinical.values[!is.na(clinical.values)], T, 
                                            seed=seed, mc.cores=mc.cores)
            
        } else if (is.numeric.param) { #kruskal wallis
            #test.res = kruskal.test(as.numeric(clinical.values[numeric.entries]),
            #				clustering[numeric.entries])
            #pvalue = test.res$p.value
            if(length(table(clustering[numeric.entries])) > 1){
                
                res = get.empirical.clinical(clustering[numeric.entries], 
                                             as.numeric(clinical.values[numeric.entries]), F, 
                                             seed=seed, mc.cores=mc.cores)
                
            } else {
                
                res = list("pvalue"=1, "conf.int"=c(NA,NA), "total.num.iters"=NA, 
                           "total.num.extreme"=NA)
            
            }
        }
        
        pvalues = c(pvalues, res$pvalue)
        total.num.iters = c(total.num.iters, res$total.num.iters)
        total.num.extreme = c(total.num.extreme, res$total.num.extreme)
        conf.int[[c]] = res$conf.int
        
        c=c+1 #update counter
    }
    names(pvalues) = params.being.tested
    names(total.num.iters) = params.being.tested
    names(total.num.extreme) = params.being.tested
    names(conf.int) = params.being.tested
    
    return(list(pvalues=pvalues, conf.int=conf.int, 
           total.num.iters=total.num.iters, total.num.extreme=total.num.extreme))
}


#' Compute empirical chi2/kruskal-wallis using a permutation test
#'
#' @param clustering vector. Clustering for each sample.
#' @param clinical.values vector. Vector with values for a clinical variable.
#' @param is.chisq boolean. True to compute chi2-test on discrete variables, False
#' to compute Kruskal-Wallis Test on numerical variables. 
#' @param seed integer. Set seed for reproducibility.
#' @param mc.cores integer. Number of cores to use (def. half of available cores).
#'
#' @return P-value for considered test, confidence interval and number of
#' iterations.
#' @export
get.empirical.clinical <- function(clustering, clinical.values, is.chisq, 
                                   seed=42, mc.cores=detectCores()/2) {
    
    set.seed(seed, "L'Ecuyer-CMRG")
    
    if (is.chisq) {
        clustering.with.clinical = cbind(clustering, clinical.values)
        tbl = table(as.data.frame(clustering.with.clinical))
        test.res = chisq.test(tbl)
    } else {
        test.res = kruskal.test(as.numeric(clinical.values), clustering)
    }
    orig.pvalue = test.res$p.value
    num.iter = 1000
    total.num.iters = 0
    total.num.extreme = 0
    should.continue = T
    
    while (should.continue) {
        #print('another iteration in empirical clinical')
        perm.pvalues = as.numeric(mclapply(1:num.iter, function(i) {
            cur.clustering = sample(clustering)
            names(cur.clustering) = names(clustering)
            
            if (is.chisq) {
                clustering.with.clinical = cbind(cur.clustering, clinical.values)
                tbl = table(as.data.frame(clustering.with.clinical))
                test.res = chisq.test(tbl)
            } else {
                test.res = kruskal.test(as.numeric(clinical.values), cur.clustering)
            }
            cur.pvalue = test.res$p.value
            return(cur.pvalue)
        }, mc.cores=mc.cores))
        total.num.iters = total.num.iters + num.iter
        total.num.extreme = total.num.extreme + sum(perm.pvalues <= orig.pvalue)
        
        binom.ret = binom.test(total.num.extreme, total.num.iters)
        cur.pvalue = binom.ret$estimate
        cur.conf.int = binom.ret$conf.int
        
        sig.threshold = 0.05
        is.threshold.in.conf = cur.conf.int[1] < sig.threshold & cur.conf.int[2] > sig.threshold
        if (!is.threshold.in.conf | total.num.iters > 1e5) {
            should.continue = F
        }
    }
    return(list(pvalue=cur.pvalue, conf.int=cur.conf.int, 
                total.num.iters=total.num.iters, total.num.extreme=total.num.extreme))
}

###################
# Other functions #
###################

#' Compute log-rank p-value
#'
#' @param survdiff.res list. Result of Test Survival Curve Differences from
#' "survival" package.
#'
#' @return Log-rank p-value
#' @export
get.logrank.pvalue <- function(survdiff.res) {
    1 - pchisq(survdiff.res$chisq, length(survdiff.res$n) - 1)  
}



