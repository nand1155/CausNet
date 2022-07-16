#' Compute index r for some set T, which is a subset of a set with n elements
#'
#' @param n Number of elements in Superset of T.
#' @param T Subset of the Superset.
#' @return index r for some set T.
subsetr = function(n, T){
  r = 0L
  for(i in seq_len(n)){
    if( is.element(i, T)) r = r + 2L ^ (n - i)
  }
  return(r)
}

#' Find a subset set T of a Superset with n elements, using index r
#'
#' @param n Number of elements in Superset of T.
#' @param r Index of the Subset of the Superset.
#' @return a subset set T that has index r as a subset of Superset with n elements.
subsetur = function(n, r){
  T = NULL
  for(i in n:1){
    if( r %% 2 == 1 ) T = c(i,T)
    r = r %/% 2
  }
  return( T )
}


#' Find possible parents for each variable based on FDR or correlation curoff.
#' @importFrom stats cor.test
#' @param mydata Data frame of data.
#' @param alpha Cutoff value of FDR or correlation.
#' @param n.var Number of elements in data.
#' @param fdr True/False - whether FDR to be used.
#' @return a list of possible parents for each variable.

mypp = function(mydata, alpha, n.var, fdr){
  pp = vector('list', n.var)
  for(v in 1:n.var){
    for(p in 1:n.var){
      if(p != v){
        if (fdr == TRUE){
          p.value = cor.test(mydata[, v], mydata[, p])$p.value
          if(p.value < alpha) pp[[v]] = c(pp[[v]], p)
        }
        else{
          corr = cor(mydata[, v], mydata[, p])
          if(corr > alpha) pp[[v]] = c(pp[[v]], p)
        }

      }
    }
  }
  return(pp)
} # end mypp

#' Find possible parents for each variable based on FDR or correlation
#' cutoff for phenotype driven two-level possible parents search.
#' @importFrom stats cor.test
#' @param mydata Data frame of data.
#' @param alpha1 Cutoff value of FDR or correlation for first-level parents for outcome variable.
#' @param alpha1 Cutoff value of FDR or correlation for second-level parents for outcome variable.
#' @param n.var Number of elements in data.
#' @param y Outcome variable.
#' @param fdr True/False - whether FDR to be used.
#' @return a list of possible parents for each variable.

mypp1 = function(mydata, alpha1, alpha2, n.var, y, fdr){
  pp = vector('list', n.var)
  nSubset1 = NULL
  nSubset = y
    for(p in 1:n.var){
      if(p != y){

        if (fdr == TRUE){
          p.value = cor.test(mydata[, y], mydata[, p])$p.value
          if(p.value < alpha1){
            pp[[y]] = c(pp[[y]], p)
            nSubset1 = c(nSubset1,p)
          }
        else{
          corr = cor(mydata[, y], mydata[, p])
          if(corr > alpha1){
            pp[[y]] = c(pp[[y]], p)
            nSubset1 = c(nSubset1,p)
          }
        }
      }
    }
  }
  for(q in 1:(n.var-1)){# not including outcome variable
    for(r in nSubset1){
      if(q != r){
        if (fdr == TRUE){
          p.value = cor.test(mydata[, q], mydata[, r])$p.value
          corr = cor(mydata[, r], mydata[, q])
          if(p.value < alpha2){
            pp[[r]] = c(pp[[r]], q)

          }
        }
        else{
          corr = cor(mydata[, r], mydata[, q])
          if(corr > alpha2){
            pp[[r]] = c(pp[[r]], q)

          }
        }

      }
    }
  }
  return(pp)
} # end mypp


#' Find possible parents for each variable based on FDR or correlation
#' cutoff for phenotype driven two-level possible parents of survival outcome variable.
#' @importFrom stats cor.test
#' @param mydata Data frame of data with last two colums of survival outcome and time.
#' @param alpha1 Cutoff value of FDR or correlation for first-level parents for outcome variable.
#' @param alpha1 Cutoff value of FDR or correlation for second-level parents for outcome variable.
#' @param Covariates The column names of Covariates in mydata.
#' @param y Outcome variable.
#' @return a list of possible parents for each variable.

mypp_surv = function(mydata, Covariates, alpha1, alpha2){
  n.var = ncol(mydata)
  n.var = n.var -1 # last col survival time
  pval = NULL
  i=0

  for (g in Covariates){
    i=i+1
    f = coxph(Surv(as.numeric(unlist(mydata[n.var+1])),
                   as.numeric(unlist((mydata[n.var]))))
              ~  as.matrix(mydata[g]),mydata)
    an = anova(f)
    pval = c(pval,an$`Pr(>|Chi|)`[4])
  }
  pval1=p.adjust(pval, "BH")
  ppDis = NULL
  i=0
  for (k in pval1){ # fixed pval1
    i = i+1
    if (k < alpha1) ppDis = c(ppDis,i)
  }
  sGenes = Covariates[ppDis]
  parentsDis = ppDis # col numbers in Covariates

  mydata1 = mydata[Covariates]
  mydata=cbind(mydata1,surdata)
  cnames = colnames(mydata)
  n.var = ncol(mydata)
  surv = 0
  if(!is.null(surdata)){
    n.var = n.var -1 # last col survival time
    surv = 1
  }

  pp = vector('list', n.var)
  nSubset1 = parentsDis
  y = n.var
  pp[[y]] = nSubset1

  for(v in 1:(n.var-1)){
    for(p in 1:(n.var-1)){
      if(p != v){
        #p.value = cor.test(mydata[, v], mydata[, p])$p.value
        p.value = cor(mydata[, v], mydata[, p])
        if(p.value > alpha2) pp[[v]] = c(pp[[v]], p)

      }
    }
  }

  return( pp )

}

#' Reformat a list of possible parents, pp, into a list of possible offspring, po
#' @param pp List of possible parents.
#' @return a list of possible offsprings for each variable.


pofun = function(pp){
  po = vector('list', length(pp))
  for(i in 1:length(pp)){
    for(j in 1:length(pp)){
      if( is.element(i, pp[[ j ]]) ) po[[i]] = c( po[[i]], j )
    }
  }
  return( po )
}

#' Create a list of all possible subsets of a given Set upto specified cardinality.
#' @param vec Vector of the Superset.
#' @param maxSubsetSize The maximum cardinality of the subsets of a given Superset to be found.
#' @return A list of all possible subsets of a given Set upto specified cardinality.
#' @importFrom utils combn

comb1 = function(vec, maxSubsetSize) {
  n = length(vec)
  out = vector('list', min( maxSubsetSize, n))
  for (j in 1:min( maxSubsetSize, n)) { # only upto 2 elements

    if( length(vec) > 1 ) out[[ j ]] = combn(vec, j) else out[[ j ]] = matrix(vec, nrow=1, ncol=1)
  }
  return(out)
}

#' Create all subsets of possible parents subsets for each node.
#' @param pp List of possible parents.
#' @return A list of all subsets of possible parents subsets for each node.
#' @importFrom utils combn

pp.sets = function(pp){
  n = length(pp)
  pps = vector('list', n)
  for(j in 1:n){
    if (!is.null(pp[[ j ]]))
    pps[[ j ]] = comb1( pp[[ j ]],n )
  }
  return(pps)
}


#' Compute the BIC score of a node with given set of parents.
#' @param y Child node.
#' @param x List of parents of y.
#' @return BIC score.
#' @importFrom stats BIC lm

score.bic.lm = function(y, x, mydata) {
  y.nm = colnames(mydata)[ y ]
  if( is.element(x[1], 1:ncol(mydata)) ) x.nms = colnames(mydata)[ x ] else x.nms = "1"
  fit = lm(paste0(y.nm, ' ~ ', paste(x.nms, collapse=' + ')), data = mydata)
  bic = -(1/2)*BIC(fit)
  return(bic)
}

#' Compute the BIC score of the survival outcome node with given set of parents.
#' @param y Child node.
#' @param x List of parents of y.
#' @return BIC score.
#' @importFrom stats BIC lm
#' @importFrom stats BIC survival

score.bic.surv = function(y, x, mydata) {
  n = ncol(mydata)
  y.nm = colnames(mydata)[ y ]
  delta = colnames(mydata)[ n ]
  X = NULL
  if(!is.na(x) ) {
    x.nms = colnames(mydata)[ x ]
    for(v in x.nms){
      X = c(X,v)
    }

    l = as.matrix(mydata[X])
    fit = coxph(Surv(as.numeric(unlist(mydata[unlist(y)])), unlist(mydata[delta])) ~ l)

  }

  else{
    x.nms = "1"
    fit = coxph(Surv(as.numeric(unlist(mydata[unlist(y)])), unlist(mydata[delta])) ~ 1)
    }

  bic = -(1/2)*BIC(fit)
  return(bic)
}

#' Compute the scores for all sets of possible parents for each node.
#' @param mydata Data frame of data with last column as output variable in case of binary/continuous outcome;
#' in case of survival data, data without the survival outcome and time.
#' @param pps List of possible parents' subsets.
#' @param surv Whether Survival outcome data, 0 or 1.
#' @param scoreFn Score Function to be used. Values 'bic' or 'bge'.
#' @return List of scores for all sets of possible parents for each node.

pp.sets.s = function(mydata, pps, surv, scoreFn){
  n = length(pps)
  ppss = vector('list', n) # possible parent set scores
  ppss1 = vector('list', n) # possible parent best sets
  ppss2 = vector('list', n) # possible parent best set scores
  for(v in 1:n){
    n.pp = ncol(pps[[v]][[1]])
    if (!is.null(n.pp)){

    for(set.size in 1:n.pp){
      ppss[[v]][[set.size]] = rep(NA, ncol(pps[[v]][[set.size]]) )
      ppss1[[v]][[set.size]] = rep(NA, ncol(pps[[v]][[set.size]]) )
      ppss2[[v]][[set.size]] = rep(NA, ncol(pps[[v]][[set.size]]) )
      for(pset.i in 1:ncol(pps[[v]][[set.size]])){
        v.pset = pps[[v]][[set.size]][ , pset.i ]
        if(v==n & surv ==1){
            ppss[[v]][[set.size]][ pset.i ] = score.bic.surv( v, v.pset, mydata )
            }
          else{

            if(scoreFn=="bic"){
              ppss[[v]][[set.size]][ pset.i ] = score.bic.lm( v, v.pset, mydata )
              }
            else{
              ppss[[v]][[set.size]][ pset.i ] = score_bge( v, v.pset, mydata )

            }

          }

        }
    }
    }
  }
  return(ppss)
}

#' Compute the scores for nodes with no parents for each node.
#' @param mydata Data frame of data with last column as output variable in case of binary/continuous outcome;
#' in case of survival data, data with last two columns the survival outcome and time.
#' @param vset List of Covariates.
#' @param surv Whether Survival outcome data, 0 or 1.
#' @param scoreFn Score Function to be used. Values 'bic' or 'bge'.
#' @return List of scores for nodes with no parents for each node.


mscores = function(vset, mydata, surv = 0, scoreFn){
  n = ncol(mydata)
  mscore = NULL
  for(v in vset){
    if(v==n-1 & surv ==1){
      mscore = c(mscore, score.bic.surv( v, NA, mydata ))
    }
    else{

      if(scoreFn=="bic"){
        mscore = c(mscore, score.bic.lm( v, NA, mydata ))
      }
      else{
        mscore = c(mscore, score_bge( v, NA, mydata ))

      }

    }
  }
  return(mscore)
}

#' Get the score of a node, given its parent set, and the list of parent set scores.
#' @param mydata Data frame of data with last column as output variable in case of binary/continuous outcome;
#' in case of survival data, data with last two columns the survival outcome and time.
#' @param v node variable.
#' @param pset Parents of v.
#' @param pps Possible parents sets' subsets.
#' @param ms List of scores of nodes with no parents.
#' @return The score of a node, given its parent set, and the list of parent set scores.

get.score = function(v, pset, pps, ppss, ms){
  if(is.null(pset)){ # FIX check for null NOT length<1
    myscore = ms[v]
  } else {
    l = length(pset)
    aa = apply(pps[[v]][[l]],2,setequal,y=pset)
    myscore = ppss[[v]][[l]][ aa ]
  }
  return(myscore)
} # end get.score

#' Find the best parent set of a node, for parent set pset, given the
#' list of Possible parents sets' subsets, list of Possible parents sets' subsets
#' scores, and list of scores of all nodes with no parents.
#' @param v node variable.
#' @param pset Parents of v.
#' @param pps List of Possible parents sets' subsets.
#' @param ppss List of Possible parents sets' subsets scores.
#' @param ms List of scores of all nodes with no parents.
#' @param max_parents Maximum in-degree.
#' @param surv Whether Survival outcome data, 0 or 1.
#' @return The score of a node, given its parent set, and the list of parent set scores.

bestps = function(v, pset, pps, ppss, ms, max_parents, surv){
  # all subsets of pset
  best.score = get.score(v, NULL, pps, ppss, ms)
  best.set = NULL
  sb = comb1( pset, max_parents ) # max_parents -maxSubsetSize
  for(j in 1:length(sb)){
    for(q in 1:ncol(sb[[j]])){
      tmp.set = sb[[j]][, q]
      tmp.score = get.score(v, tmp.set, pps, ppss, ms)
      if(surv == 1){
        if(q == 1) {
          best.score = tmp.score
          best.set = tmp.set
        }
      }
      if( tmp.score > best.score ){
        best.set = tmp.set
        best.score  = tmp.score
      }
    }
  }
  outp = vector('list', 2)
  outp[[1]] = best.set
  outp[[2]] = best.score
  return(outp)
} # end function bestps


#' Find the best parent sets and best scores for all sets of possible parents for each node.
#' @param pps List of Possible parents sets' subsets.
#' @param ppss List of Possible parents sets' subsets scores.
#' @param ms List of scores of all nodes with no parents.
#' @param max_parents Maximum in-degree.
#' @param surv Whether Survival outcome data, 0 or 1.
#' @return The list of best parent sets and best scores for all sets of possible parents for each node.


pp.sets.bs = function(pps, ppss, ms, max_parents, surv){
  bps = ppss # best parent sets
  bpss = ppss # best parent set scores
  for(v in 1:length(pps)){ # length(pps) - no. of vars
    # for each possible parent set, compare score w/ all subsets to identify
    # the best possible parent sets
    n.pp = ncol(pps[[v]][[1]]) #n.pp number of pp for v
    #print("n.pp")
    # print(n.pp)
    if (!is.null(n.pp)){ # for each v, n.pp number of pp .. for each of
      # pps for v, need to find all subsets, find their score and then choose
      # the best score and best parents for that pps
      for(set.size in 1:n.pp){ # n.pp number of pp
      #print("#######")
      #for(set.size in 1:min(n.pp,max_parents)){ # max parents NOT here
        n.sets = ncol(pps[[v]][[set.size]])
        # n.sets- no of sets of size n.sets
        #print(n.sets)
        if (!is.null(n.sets)){
          bps[[v]][[set.size]] = vector('list', n.sets)
          for(k in 1:n.sets){
              tmp.set = pps[[v]][[set.size]][, k]
              # tmp.set one set of size set.size at a time

              #print(tmp.set)
              tmp = bestps(v, tmp.set, pps, ppss, ms, max_parents,surv)
              # print("tmp")
              #print(tmp)
              bps[[v]][[set.size]][[k]] = tmp[[1]]
              bpss[[v]][[set.size]][k] = tmp[[2]]
            }
          }
        }
      }
  }
  outp = vector('list', 2)
  outp[[ 1 ]] = bps
  outp[[ 2 ]] = bpss
  return( outp )
} # end function pp.sets.bs

#' Given a node s, and a set w, compute score for s with its best parent set in w.
#' argument bps is generated by pp.sets.bs
#' @param s node variable.
#' @param w Possible Parents subset for node s.
#' @param pps List of Possible parents sets' subsets.
#' @param bps List of best parent subset in each Possible parents sets' subsets.
#' @return The score for s with its best parent set in w.

swscore = function(s, w, pp, pps, bps){
  pset = w[ is.element( w, pp[[s]] ) ] # find possible parents of s in w
  l = length(pset)

  aa = apply(pps[[s]][[l]],2,setequal,y=pset)
  best.ps = bps[[1]][[s]][[l]][aa]
  best.pss = bps[[2]][[s]][[l]][aa]

  outp = vector('list', 2)
  outp[[ 1 ]] = best.ps
  outp[[ 2 ]] = best.pss
  return(outp)
} # end swscore

#' Given a set w of nodes, and v a possible offspring of at least one of the nodes in w,
#' compute network scores for all possible
#' sinks of sets, w1 = {w + v}.
#' @param w A set of nodes.
#' @param w.networkscore The score of the best network for w.
#' @param pp List of possible parents.
#' @param po List of possible offsprings.
#' @param bps List of best parent subset in each Possible parents sets' subsets.
#' @param m Number of variables in the whole data.
#' @return A list of network scores for all possible
#' sinks of setsof type w1 = {w + v}.


wsink.scores = function(w, w.networkscore, pp, po, pps, bps, m){

  # Find possible offspring for w not already in w
  wpo = NULL
  for(v in w){
    wpo = unique( c(wpo, po[[v]][ !is.element( po[[v]], w ) ]) )
  }

  if( length(wpo) > 0 ){
    windx <- k <- sink <- wscore <- numeric(length(wpo))
    # Expand w by one po node, a possible sink, and compute score
    rowno = 1
    for(s in wpo){
      s.score = swscore(s, w, pp, pps, bps)[[2]]
      wscore[rowno] = s.score + w.networkscore
      windx[rowno] = subsetr(m, c(w, s))
      k[rowno] = length(w) + 1
      sink[rowno] = s
      rowno = rowno + 1
    }
  } else {
    w1sinks = NULL    # end if length(wpo)
    return( w1sinks )
  }

  w1sinks <- list(wscore = wscore,
                  windx = windx,
                  k = k,
                  sink = sink)
  return( w1sinks )
} # end wsink.scores

#' Given a list of best sinks, find a best network.
#' @param bsinks A list of best sinks.
#' @param m Number of variables in the whole data.
#' @return A data frame of a best network.


bestnet = function(bsinks, m){
  nms = c("windx", "k", "sink","subset","wscore", "component") # adding subset
  bestnets = as.data.frame(matrix(NA, nrow=0, ncol=length(nms)))
  names(bestnets) = nms

  ## Order best sinks, removing a sink at each step
  # Loop over connected components
  mycomp = 1
  rowno = 1
  while(nrow(bsinks) > 0){

    ks = unique( bsinks$k )
    ks = ks[order(ks, decreasing=TRUE)]
    q = max(ks)
    q1 = min(ks)
    aa = match(q, bsinks$k) #match gets all rows with k=q
    tmp.s = bsinks[ aa, ]   #match gets all rows with k=q
    bestnets[rowno, c("windx", "k", "subset","sink","wscore")] = tmp.s[1, c("windx", "k", "subset","sink","wscore")]
    bestnets[rowno, "component"] = mycomp
    # print("tmp.s[1, windx]")
    # print(tmp.s[1, "windx"])
    myw = subsetur(m, tmp.s[1, "windx"]) # just pick 1st row
    w1 = myw[ !is.element(myw, tmp.s[1, "sink"]) ] # remove sink in the top row
    w1indx = subsetr(m, w1) # myw w/o sink

    rowno = rowno + 1
    wlen = length(w1)
    #print("wlen===")
    #print(wlen)
    if(wlen >1){
      for(d in wlen:q1){ # FIX for partial - replaced 2 by q1
        aa = match(w1indx, bsinks$windx)
        tmp.s = bsinks[ aa, ]
        ######## aa gets only one top row
        ######## so inside this loop you just get the top rows for each subset
        bestnets[rowno, c("windx", "k","subset", "sink","wscore")] = tmp.s[1, c("windx", "k","subset", "sink","wscore")]
        bestnets[rowno, "component"] = mycomp
        w = subsetur(m, tmp.s[1, "windx"])
        w1 = w[ !is.element(w, tmp.s[1, "sink"]) ]
        w1indx = subsetr(m, w1)
        rowno = rowno + 1
      }
    }


    # remove all rows in bsinks with sets that have elements of the largest w in bestnets
    aa = NULL
    for(j in 1:nrow(bsinks)){
      wj = subsetur(m, bsinks[j, "windx"])
      aa = c(aa, sum( is.element(wj, myw) ) == 0)
    }
    bsinks = bsinks[aa, ]
    break
  } # end while loop

  return(bestnets)
} # end bestnet

#' Find edges of a best network connected components from ordering of sinks.
#' @param bnets A list of best orderings of sinks.
#' @param pp List of Possible parents for all nodes.
#' @param pps List of Possible parents sets' subsets.
#' @param bps List of best parent subset in each Possible parents sets' subsets.
#' @return A data frame of edges of a best network.

sink2net = function(bnets, pp, pps, bps){
  m = length(bps[[1]])
  nms = c("node.source", "node.sink", "component")
  mynets = as.data.frame(matrix(NA, nrow=0, ncol=length(nms)))
  names(mynets) = nms
  n.comp = length(unique(bnets$component))
  rowno = 1
  for(c in 1:n.comp){
    tmpn = bnets[ is.element(bnets$component, c), ]
    for(q in max(tmpn$k):min(tmpn$k)){ # REPLCED k by q, 2 by min(tmpn$k)
      bp.set = NULL
      tmp = tmpn[ is.element(tmpn$k, q), ]
      s = tmp[ 1, "sink"]
      w = subsetur(m, tmp[ 1, "windx"])
      bp.set1 = swscore(s, w, pp, pps, bps) ## BUG FIX [[1]]
      bp.set = swscore(s, w, pp, pps, bps)[[1]] # ??????? OK???
      if(!is.null(bp.set[[1]])){
        src = bp.set[[1]]
        #print(length(bp.set[[1]]))
        snk = rep(s, length(bp.set[[1]]))
        cmp = rep(c, length(bp.set[[1]]))
        mynets[ c(rowno:(rowno+length(bp.set[[1]])-1)),] = cbind(src, snk, cmp)

      }
      rowno = rowno+length(bp.set[[1]])
    }


  }
  return(mynets)
} # end sink2net

#' Compute the BGE score of a given child node with given set of parents.
#' @param y Child node.
#' @param x List of parents of y.
#' @return Numeric BGE score
#' @export
score_bge <- function(y, x, mydata) {
  #print("score bge")
  names <- colnames(mydata)

  if (any(is.na(x))) {
    model <- paste0("[",
                    names[y],
                    "]")
    model.net <- bnlearn::model2network(model)
    model_data <- mydata[y]
    bnlearn::modelstring(model.net) <- model
    res <- bnlearn::score(model.net, model_data, type = "bge")

  } else {
    names1 = names[y]
    names2 = names[x]
    names3 <- c(names1,names2)
    namesLeft = names[!(names %in% names3)]
    namesUse = c(namesLeft,names2)
    model <- paste0(paste0("[", paste0(names[x], collapse = "]["), "]"),
                    "[",
                    names[y],
                    "|",
                    paste0(names[x], collapse = ":"),
                    "]")

    model.net <- bnlearn::model2network(model)
    model_data <- mydata[c(y, x)]
    bnlearn::modelstring(model.net) <- model

    model1 <- paste0("[", paste0(names[x], collapse = "]["), "]")
    #print(model1)
    model.net1 <- bnlearn::model2network(model1)
    model_data1 <- mydata[names[x]]
    #print(model.net1)

    bnlearn::modelstring(model.net1) <- model1
    bnlearn::modelstring(model.net) <- model
    #print(model)
    networkScore <- bnlearn::score(model.net, model_data, type = "bge")
    sourceScore <- bnlearn::score(model.net1, model_data1, type = "bge")
    res = networkScore - sourceScore
  }
  #print(res)
  res
}



