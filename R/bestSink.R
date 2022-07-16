#' Find best sinks for all subsets of variables in data.
#' @importFrom stats cor.test
#' @param mydata Data frame of data.
#' @param ms List of scores of nodes with no parents.
#' @param pp List of possible parent sets for all variables.
#' @param po List of possible offspring sets for all variables.
#' @param pps List of Possible parents sets' subsets.
#' @param ppss List of Possible parents sets' subsets' scores.
#' @param bps List of best parents sets' for each parents subsets'.
#' @param surv Whether Survival outcome data, 0 or 1.
#' @return A data frame of best sinks for each subset with score.

bestSinks = function(pp, ms, po, pps, ppss, bps, mydata, surv){

  m = ncol(mydata) # number of nodes
  if(surv == 1){
    m = m-1
  }
  nms = c("windx", "k", "sink", "wscore") #subset
  sinks.tmp = as.data.frame(matrix(NA, nrow=0, ncol=length(nms)))
  rownames(sinks.tmp) <- NULL
  names(sinks.tmp) = nms

  # best sinks and scores for subnetworks of one node, which is the node itself and its score
  for(s in 1:m){
    sinks.tmp[s, "windx"] = subsetr(m, s)
    sinks.tmp[s, "k"] = 1
    sinks.tmp[s, "sink"] = s
    sinks.tmp[s, "wscore"] = ms[s]
    }
  sinks.tmp = round(sinks.tmp,4)
  mysinks = sinks.tmp
  bsinks = sinks.tmp[0, ] # names row

  for(q in 2:m){
    sinks.tmp1 = list()
    wscore <- windx <- k <- sink <- numeric(m*m)
    index <- 1

    for(j in seq_len(nrow(sinks.tmp))) {
      w = subsetur(m, sinks.tmp[j, "windx"])
      w.networkscore = sinks.tmp[j, "wscore"]
      w1sinks = wsink.scores(w, w.networkscore, pp, po, pps, bps, m)
      index_subset <- seq_along(w1sinks$wscore)-1+index
      wscore[index_subset] <- w1sinks$wscore
      windx[index_subset] <- w1sinks$windx
      k[index_subset] <- w1sinks$k
      sink[index_subset] <- w1sinks$sink
      index <- index + length(index_subset)
    }
    sinks.tmp1  <- data.frame(wscore = wscore[seq_len(index-1)],
                              windx = windx[seq_len(index-1)],
                              k = k[seq_len(index-1)],
                              sink = sink[seq_len(index-1)])
    sinks.tmp1 = round(sinks.tmp1,4)
    #print(sinks.tmp1,quote = TRUE, row.names = FALSE)
    # break q loop if there are no more offspring for any sets
    if( nrow(sinks.tmp1) == 0 ) break

    # for each subset w, find the best sinks
    # keep the row/rows with max score
    myws = unique( sinks.tmp1$windx )
    for(wind in 1:length(myws)){
      myw = myws[ wind ]
      tmp = sinks.tmp1[ is.element( sinks.tmp1$windx, myw ), ]
      tmp1 = tmp[ tmp$wscore >= max(tmp$wscore), ]
      bsinks = rbind( bsinks, tmp1 )
    }
    bsinks = unique(bsinks)
    sinks.tmp = bsinks[ is.element( bsinks$k, q ), ]
    sinks.tmp = unique(sinks.tmp)

  } # end q for loop

  #print(bsinks,quote = TRUE, row.names = FALSE)
  bsinks11 = bsinks
  wsubCol = NULL
  if(nrow(bsinks)>0){
    for(i in 1:nrow(bsinks)){
      if(is.na(bsinks[i,"windx"])) wsubset = NA
      else wsubset <- paste0(subsetur(m,bsinks[i,"windx"]), collapse = ",")
      wsubCol = c(wsubCol,wsubset)
    }

  }

  bsinks["subset"] = wsubCol
  #print(bsinks,quote = TRUE, row.names = FALSE)

  return(bsinks)
} # end bestSinks
