# function to loop through offsprings to find best sinks
bestSinksPartialOrds = function(pp, ms, po, pps, ppss, bps, mydata, surv, max_parents, pOrd){
  print("inside bestSinksPartialOrds")
  m = ncol(mydata) # number of nodes
  if(surv == 1){
    m = m-1
  }
  nms = c("windx", "k", "sink", "wscore") #subset
  sinks.tmp = as.data.frame(matrix(NA, nrow=0, ncol=length(nms)))
  rownames(sinks.tmp) <- NULL
  names(sinks.tmp) = nms

  # best sinks and scores for subnetworks of pOrd nodes and its score
  pOrdsubs = combn(m, (m-pOrd))
  print(pOrdsubs)
  index = 1
  for(s1 in 1:ncol(pOrdsubs)){
    w = pOrdsubs[,s1]
    print("w")
    print(w)
      w.networkscore = 0
      w1sinks = wsink.scores(w, w.networkscore, pp, po, pps, bps, m)
      if(!is.null(w1sinks)){
        print("w1sinks")
        print(w1sinks,quote = TRUE, row.names = FALSE)
        ks = unique( w1sinks$wscore )
        q = max(ks)
        print("q")
        print(q)
        sinks.temp1  <- data.frame(wscore = w1sinks$wscore,
                                  windx = w1sinks$windx,
                                  k = w1sinks$k,
                                  sink = w1sinks$sink)
        #sinks.temp1 = lapply(sinks.temp1,as.numeric)
        d = sinks.temp1[(sinks.temp1$wscore==q),]
        rowsToDelete = list()
        for (j in nrow(d)) {
          # find bps for sink in w
          # find w in pps
          s = d[j,]$sink
          l = length(w)
          print("s")
          print(s)
          print("l")
          print(l)
          if(length(pps[[s]])>= l){
            aa = apply(pps[[s]][[l]],2,setequal,y=w)
            best.ps = bps[[1]][[s]][[l]][aa]
            if(length(best.ps) != 0) # due to problem with gen ordering
               {
                 #print("best.ps===" )
                 #print(best.ps )
                 best.pss = bps[[2]][[s]][[l]][aa]

                 print("best.ps")
                 print(best.ps[[1]])

                 for (b in best.ps[[1]]) {
                   print("b")
                   print(b)
                   d[j,]$wscore = d[j,]$wscore + ms[b]
                 }
                 d[j,]$wscore = d[j,]$wscore + best.pss

            }# end length(pps[[s]])>= l
            else{
              rowsToDelete = c(rowsToDelete, j)
            }
        }# end if(length(pps[[s]][[l]])!= 0)
        else{
          rowsToDelete = c(rowsToDelete, j)
        }
      }
      if(length(rowsToDelete) != 0) {
        ix <- which(rowsToDelete %in% rownames(d))
        d = d[-ix,]
      }
      sinks.tmp = rbind( sinks.tmp, d)
    } # end if(!is.null(w1sinks))
  }

  mysinks = sinks.tmp
  print("mysinks")
  #sinks.tmp = lapply(sinks.tmp,as.numeric)
  #sinks.tmp = round(sinks.tmp,3)
  print(mysinks,quote = TRUE, row.names = FALSE)

  bsinks = sinks.tmp[0, ] # names row
  bsinks = mysinks


  print("bsinks")
  #print(bsinks,quote = TRUE, row.names = FALSE)
  # best sinks and scores for subnetworks of size 2:m
  for(q in (m-pOrd+1):m){
    print(q)
    sinks.tmp1 = list()
    wscore <- windx <- k <- sink <- numeric(m*m)
    index <- 1

    for(j in seq_len(nrow(sinks.tmp))) {
      w = subsetur(m, sinks.tmp[j, "windx"])
      w.networkscore = sinks.tmp[j, "wscore"]
      w1sinks = wsink.scores(w, w.networkscore, pp, po, pps, bps, m)
      print("w1sinks")
      print(w1sinks)
      index_subset <- seq_along(w1sinks$wscore)-1+index
      print("index_subset")
      print(index_subset)
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
    print(sinks.tmp1)
    #sinks.tmp1 = round(sinks.tmp1,3)
    print(sinks.tmp1,quote = TRUE, row.names = FALSE)
    # break q loop if there are no more offspring for any sets
    if( nrow(sinks.tmp1) == 0 ) break

    # for each subset w, find the best sinks
    # keep the row/rows with max score
    myws = unique( sinks.tmp1$windx )
    for(wind in 1:length(myws)){
      print("wind")
      print(wind)
      myw = myws[ wind ]
      tmp = sinks.tmp1[ is.element( sinks.tmp1$windx, myw ), ]
      print(tmp)
      scrs = unique( tmp$wscore )
      sc = max(scrs)
      tmp1 = tmp[(tmp$wscore==sc),]
      print(tmp1)
      bsinks = rbind( bsinks, tmp1 )
    }
    print("bsinks11111111111111111111")
    print(bsinks,quote = TRUE, row.names = FALSE)
    bsinks = unique(bsinks)
    sinks.tmp = bsinks[ is.element( bsinks$k, q ), ]
    # NEW Sep 29 -- remove duplicates in sinks.tmp - sink
    #sinks.tmp = sinks.tmp[ , -which(names(sinks.tmp) %in% c("sink"))]
    # print("bsinks22222222222222222")
    # print(sinks.tmp,quote = TRUE, row.names = FALSE)
    sinks.tmp = unique(sinks.tmp)
    # NEW Sep 29 --keep only the row/rows with max score
    # for each subset with card q in sinks.tmp
  } # end q for loop

  print(bsinks,quote = TRUE, row.names = FALSE)
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
  print(bsinks,quote = TRUE, row.names = FALSE)

  return(bsinks)
} # end bestSinks
