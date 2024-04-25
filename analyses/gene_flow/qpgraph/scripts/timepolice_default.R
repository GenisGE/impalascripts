
getAnc <- function(x,df,rev=FALSE){

    anc <-  df[df$to==tail(x,1),"from"]
    if(length(anc)>2){
        cat("2> anc. To",tail(x,1),"\n")
        print(anc)
    }
    if(length(anc)==2){
        if(rev)
            return( c(getAnc(c(x,anc[1]),df=df),getAnc(anc[2],df=df)))
        else
            return( c(getAnc(c(x,anc[2]),df=df),getAnc(anc[1],df=df)))
    }
    if(length(anc)==0)
        return(x)
  
    getAnc(c(x,anc),df=df)  
}

# returns the node with time fuckup

timepoliceOne <- function(df){

    nodes <- unique(df$to)
    ancestors <-lapply(nodes,getAnc,df=df)
    ancestors2 <-lapply(nodes,getAnc,df=df,rev=T)

    ##is the first direct ancestor node also part of the path of the second ancestror
    timeFuck <- sapply(ancestors,function(x) x[2]%in%x[-2] )
    ## same but swich ancestors
    timeFuck2 <- sapply(ancestors2,function(x) x[2]%in%x[-2] )

    nodes[timeFuck | timeFuck2]
}

timepoliceTwo <- function(df){
    ## remove scenarious with two simultationus admixture events
    #    w <- which(igraph::degree(g, mode='in') == 2)
    w <- df$to%in%df[duplicated(df$to),2]
    if(any(duplicated(df$from[w])))
        return(TRUE)

    return(FALSE) ## no simultanous admixture
}


timepolice <- function(graph) {


  # admixed node with whos two ancestors are ancestors of each other.   e.g.
  # A,B -> X (X mix of A and B), A->Z -> Y->B (B has A as an ancestor)
  df <- as.data.frame(igraph::as_edgelist(graph), stringsAsFactors = FALSE)
  names(df) <- c("from", "to")
  if (length(timepoliceOne(df)) > 0) {
    return(FALSE)
  }

  # node with has to admixed direct descendands
  # A,B -> X (X mix of A and B), A,C -> Y (A has two admixed direct descendants)
  if (timepoliceTwo(df)) {
    return(FALSE)
  }



  return(TRUE)
}
