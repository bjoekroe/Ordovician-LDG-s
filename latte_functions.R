############ functions for latte_final

###############
### basic counts

bcounts <- function(occs) 
{counts_ordo <- matrix(,3,7)
for (i in 2:(nrow(xts)-1)) 
{
  cns <- as.character(xts[i,1])
  ox <- subset(occs, grepl(cns, occs$ts))
  counts_ordo[3,i-1] <- length(unique((ox$geoplate)))
  counts_ordo[2,i-1] <- length(unique((ox$collection_no)))
  counts_ordo[1,i-1] <- length(unique((ox$accepted_name)))
}
return(counts_ordo)
}

######################
### rarefaction trends

rd.fun <- function(ooclat, quota) {
  
  xocc.2a <- xtabs(~ooclat$ts +ooclat$accepted_name,drop.unused.levels = TRUE);
  xocc.2a  <- t(xocc.2a )
  xocc.2a  <- as.data.frame.matrix(xocc.2a)
  
  md2a <- matrix(,(nrow(xts)-2),6)
  bin.name <- array()
  rf.s <- array()
  
  for (i in 2:(nrow(xts)-1)) {
    cns <- as.character(xts[i,1])
    bin.name.a <- cns
    bin.name <- rbind(bin.name, bin.name.a)
    if (cns %in% colnames(xocc.2a ) == TRUE) {
      rf.s <- xocc.2a [,cns]
      occi <- subset(ooclat, ooclat$ts==cns)
      md2a[i-1, 4] <- diversity(rf.s)/log(specnumber(rf.s))
      rarx <- rarefy(rf.s, quota, se=TRUE) 
      md2a[i-1,3] <- rarx[1]
      md2a[i-1,1] <- length(unique(occi$collection_no))
      md2a[i-1,2] <- specnumber(rf.s)
      md2a[i-1,5] <- rarx[1]+rarx[2]
      md2a[i-1,6] <- rarx[1]-rarx[2]
    }
  }
  bin.name <- bin.name[-1]
  row.names(md2a) <- bin.name
  colnames(md2a) <- c("collections", "taxa", "rarediv", "evenness", "rar+err", "rar-err")
  md2a[is.na(md2a)] <- 0

  md2bins <- as.data.frame(md2a)
  return(md2bins)	
}

#####################
### rarefaction LDG's

rarplus <- function(occs, lt, mmcol)
{
  ### create arrays
  rd.a <- array(,dim =ncol(occs)) # rarefied diversity
  se_rd.a <- array(, dim =ncol(occs))
  ss <- array(, dim =ncol(occs)) # number of genera in bin
  so <- array(, dim =ncol(occs)) # number of occurrences in bin
  
  ### (1) calculates S raw (ss) and number of collections (so) per bin
  for (i in 1:ncol(occs)) 
  {
    rfs <- occs[,i]
    ss[i] <- length(which(rfs>0)) 
    so[i] <- sum(rfs)
    rar.x <- rarefy(rfs, mmcol, se=TRUE)
    se_rd.a[i] <- rar.x[2]
    rd.a[i] <- rar.x[1]
  }
  
  ### (2) does permutated resampling and calculates min collections of individual resamples per time bin
  tot.res <- rbind(ss, so, rd.a, se_rd.a)
  tot.res <- round(tot.res, digits=2)
  colnames(tot.res) <- c(1:ncol(occs))
  return(tot.res)
}

####################
### alpha, beta, gamma trends per Latzone

bplusa_tot <- function(occs)
{ 
  occs <- of
  rp <- array(, c(9,7,4))  # array 1: ts; 2: zone, 3: values per zone
  for (i in 1:4) ### ### loop for each latitudinal zone
  {
    
    occi <- subset(occs, occs$paleolat==i)
    for (k in 2:(nrow(xts)-1)) ## loop for each time slice
    {
      cns <- as.character(xts[k,1])
      ots <- subset(occi, grepl(cns, occi$ts))
      otl <- droplevels(ots)
      otl <- na.omit(otl[,c("geoplate","accepted_name")])
      ox <- table(otl$geoplate, otl$accepted_name)
      oalpha <- table(otl$accepted_name)
      beta.d <- d(ox, lev="beta", boot="TRUE")
      gamma.d <- d(ox, lev="gamma", boot="TRUE")
      alpha.d <- d(ox, lev="alpha", boot="TRUE")
      
      rp[1,k-1,i] <- sum(oalpha) # nr of occurrences
      rp[2,k-1,i] <- nrow(oalpha) # nr of genera
      rp[3,k-1,i] <- nrow(ox) # nr of plates
      rp[4,k-1,i] <- as.numeric(as.character(gamma.d[1]))
      rp[5,k-1,i] <- as.numeric(as.character(gamma.d[2])) # SE
      rp[6,k-1,i] <- as.numeric(as.character(beta.d[1]))
      rp[7,k-1,i] <- as.numeric(as.character(beta.d[2])) # SE
      rp[8,k-1,i] <- as.numeric(as.character(alpha.d[1])) # SE
      rp[9,k-1,i] <- as.numeric(as.character(alpha.d[2])) # SE
    }
  }
  return(rp)
}

#####################
### alpha, beta, gamma LDG's

betaplusgamma <- function(occs)
{ 
  rk <- array(, c(9,4,7))  # array 1: ts; 2: zone, 3: values per zone
  for (i in 2:(nrow(xts)-1)) ### loop for each time slice !! avoid i=5 (nrow(xts)-1)
  {
    cns <- as.character(xts[i,1])
    
    for (k in 1:4) ### loop for each latitudinal zone
    {
      ots <- subset(occs, grepl(cns, occs$ts))
      otl <- ots[which(ots$paleolat==k),]
      otl <- droplevels(otl)
      otl <- na.omit(otl[,c("geoplate","accepted_name")])
      ox <- table(otl$geoplate, otl$accepted_name)
      
      oalpha <- table(otl$accepted_name)
      beta.d <- d(ox, lev="beta", boot="TRUE")
      gamma.d <- d(ox, lev="gamma", boot="TRUE")
      alpha.d <- d(ox, lev="alpha", boot="TRUE")
      
      rk[1,k,i-1] <- sum(oalpha) # nr of occurrences
      rk[2,k,i-1] <- nrow(oalpha) # nr of genera
      rk[3,k,i-1] <- nrow(ox) # nr of plates
      rk[4,k,i-1] <- as.numeric(as.character(gamma.d[1]))
      rk[5,k,i-1] <- as.numeric(as.character(gamma.d[2])) # SE
      rk[6,k,i-1] <- as.numeric(as.character(beta.d[1]))
      rk[7,k,i-1] <- as.numeric(as.character(beta.d[2])) # SE
      rk[8,k,i-1] <- as.numeric(as.character(alpha.d[1]))
      rk[9,k,i-1] <- as.numeric(as.character(alpha.d[2])) # SE
    }
  }
  return(rk)
}


############################
######## function calcs erates 
eratescalc <- function(of)
{ 
  e.rates <- matrix(0,7,4)
  re <- array(, c(6,4,7))  # array 1: values; 2: zone, 3: ts
  for (i in 2:(nrow(xts)-1)) 
  {
    cns <- as.character(xts[i,1])
    cns3 <- as.character(xts[i+1,1])
    cns1 <- as.character(xts[i-1,1])
    ots <- subset(of, grepl(cns, of$ts))
    ots3 <- subset(of, grepl(cns3, of$ts))
    ots1 <- subset(of, grepl(cns1, of$ts))
    ots.u <- unique(ots[,2:4])
    ots1.u <- unique(ots1[,2:4])
    ots3.u <- unique(ots3[,2:4])
    me.immis <- merge(ots1.u, ots.u, by='accepted_name') #  = bc, all bottom crossers
    me.tc <- merge(x=ots.u, y=ots3.u, by="accepted_name") # top crossers
    me.throughnull <- merge(me.immis, ots3.u, by = 'accepted_name')  # rangethrough
    me.onull <- merge(x=ots.u, y=ots1, by="accepted_name", all = TRUE) # all lower
    me.singlnull <- merge(x=me.onull, y=ots3, by="accepted_name", all = TRUE)
    singlnull <- me.singlnull[is.na(me.singlnull$ts.y)&is.na(me.singlnull$ts),] # singletons
    immis.null <-  me.immis[which(!as.numeric(me.immis[,5])-as.numeric(me.immis[,3])==0),]
    extirps.null <-  me.tc[which(!as.numeric(me.tc[,5])-as.numeric(me.tc[,3])==0),]
    
    genus.count <- length(unique(ots$accepted_name))
    sing.n <- length(unique(singlnull$accepted_name))
    tc.n <- length(unique(me.tc$accepted_name))
    bc.n <- length(unique(me.immis$accepted_name))
    through.n <- length(unique(me.throughnull$accepted_name))
    msd <- ((through.n*2)+tc.n+bc.n)/2
    ore <- (-log(through.n/(tc.n+through.n)))/d.iff[i-1]
    exe <- (-log(through.n/(bc.n+through.n)))/d.iff[i-1]
    e.rates[i-1,1] <- genus.count
    e.rates[i-1,2] <- msd
    e.rates[i-1,3] <- ore
    e.rates[i-1,4] <- exe
    
    for (k in 1:4) {
      
      extirps.l <- extirps.null[extirps.null$paleolat.x==k,]
      extirps.r <-  mean(extirps.l$paleolat.y-extirps.l$paleolat.x)# as relation
      
      immis.l <- immis.null[immis.null$paleolat.y==k,]
      immis.r <- mean(immis.l$paleolat.y-immis.l$paleolat.x) # # as relation
      
      tc.lnx <- me.tc[me.tc$paleolat.x==k,]
      tc.ln <- length(unique(tc.lnx$accepted_name)) 
      tc.l <- tc.ln-(length(unique(extirps.l$accepted_name)))
      
      bc.lnx <- me.immis[me.immis$paleolat.y==k,]
      bc.ln <- length(unique(bc.lnx$accepted_name))
      bc.l <- bc.ln-(length(unique(immis.l$accepted_name)))
      
      sing.lnx <- singlnull[singlnull$paleolat.y==k,]
      sing.l <- length(unique(sing.lnx$accepted_name)) 
      
      through.lnx <- me.throughnull[me.throughnull$paleolat.y==k,]
      through.l <- length(unique(through.lnx$accepted_name))
      
      #extirps.r <- length(unique(extirps.l$accepted_name))/genus.count # as relation
      #immis.r <- length(unique(immis.l$accepted_name))/genus.count # # as relation
      
      msd.l <- ((through.l*2)+tc.l+bc.l)/2
      ore.l <- -log(through.l/(tc.l+through.l))/d.iff[i-1]
      exe.l <- -log(through.l/(bc.l+through.l))/d.iff[i-1]
      
      re[1,k,i-1] <- genus.count
      re[2,k,i-1] <- msd.l
      re[3,k,i-1] <- ore.l
      re[4,k,i-1] <- exe.l
      re[5,k,i-1] <- extirps.r
      re[6,k,i-1] <- immis.r
    }
  }
  return(re)
  return(e.rates)
}
