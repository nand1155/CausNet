

## create all truncated sets of possible parents for each node
## k levels
pp.sets1 = function(pp, k = NULL){
  n = length(pp)
  pps = vector('list', n)
  for(j in 1:n){
    if (!is.null(pp[[ j ]]))
      pps[[ j ]] = comb11( pp[[ j ]],n , k )
  }
  return(pps)
}

# create list object with all combinations of vec {all} upto card n-k

#' @importFrom utils combn
comb11 = function(vec, N, k) {
  n = length(vec)
  out = vector('list', n)
  for (j in 1 :(N-k)) { # for all nodes

    if( length(vec) > 1 ) out[[ j ]] = combn(vec, j) else out[[ j ]] = matrix(vec, nrow=1, ncol=1)
  }
  return(out)
}

## create object with scores for all sets of possible parents for each node
pp.sets.s.pOrd = function(mydata, pps, surv, score, pOrd){
  n = length(pps)
  ppss = vector('list', n) # possible parent set scores
  ppss1 = vector('list', n) # possible parent best sets
  ppss2 = vector('list', n) # possible parent best set scores
  for(v in 1:n){
    n.pp = ncol(pps[[v]][[1]])
    #print("%%%%%%%% n.pp=")
    #print(n.pp)
    if (!is.null(n.pp)){
      #print("22222222")
      for(set.size in 1:(n-pOrd)){

        ppss[[v]][[set.size]] = rep(NA, ncol(pps[[v]][[set.size]]) )
        ppss1[[v]][[set.size]] = rep(NA, ncol(pps[[v]][[set.size]]) )
        ppss2[[v]][[set.size]] = rep(NA, ncol(pps[[v]][[set.size]]) )
        for(pset.i in 1:ncol(pps[[v]][[set.size]])){
          v.pset = pps[[v]][[set.size]][ , pset.i ]
          if(v==n & surv ==1){
            ppss[[v]][[set.size]][ pset.i ] = score.bic.surv( v, v.pset, mydata )
          }
          else{

            if(score=="bic"){
              #print(v)
              #print("BIC")
              ppss[[v]][[set.size]][ pset.i ] = score.bic.lm( v, v.pset, mydata )
            }
            else{
              #print(v)
              #print("Bge")
              ppss[[v]][[set.size]][ pset.i ] = score_bge( v, v.pset, mydata )

            }

          }
        }

      }# end for(set.size in 1:n.pp)
    }
  }
  return(ppss)
}

