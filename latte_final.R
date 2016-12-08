###############################################################################################################################
###### This file calculates diversity measures based on results from pb_rnames_match.r ###################
###############################################################################################################################

#### results: r1-4, rarified diversity trends for latzones; rl, rarified diversity LDG's, rk, gamma LDG's; rt, gamma total; 
#### results: rp, gamma- latitudinal zones; re, evolutionary rates trend

library(vegan)
library(vegetarian)
source("/Users/bkroger/Documents/r/rnames/rname_finals/rname_functions.r");
setwd("/Users/bkroger/Documents/r/rnames/rname_finals/");
source("latte_functions.r");

##############
#####LOAD DATA 

ordo_res <-read.table(file="ordo_match_chrono.csv", sep=',', header=T)
cts <- cbind(Chrono_ts,row.names(Chrono_ts))
xts <- cts
d.iff <- c.diff

#http2<-paste("http://paleobiodb.org/data1.2/colls/list.txt?interval=Ordovician&show=ref,loc,stratext,paleoloc")
#coll_tot_1 <- read.table(http2, sep=',', header=T)
#http2<-paste("http://paleobiodb.org/data1.2/colls/list.txt?interval=Silurian&show=ref,loc,stratext,paleoloc")
#coll_tot_2 <- read.table(http2, sep=',', header=T)
#http2<-paste("http://paleobiodb.org/data1.2/colls/list.txt?interval=Cambrian&show=ref,loc,stratext,paleoloc")
#coll_tot_3 <- read.table(http2, sep=',', header=T)
#coll_tots <- rbind(coll_tot_1, coll_tot_2, coll_tot_3)
#coll_tots <- coll_tots[,c("collection_no", "paleolat", "geoplate")]
# write.csv(coll_tots, file = "coll_tots_erates_lat.csv")
coll_tots <- read.table(file = "coll_tots_erates_lat.csv", sep=',', header=T)
coll_tots <- coll_tots[,2:4]

#gen_C <- read.table("http://paleobiodb.org/data1.2/occs/list.txt?interval=Cambrian&show=loc,class", sep=',', header=T)
#gen_C <- subset(gen_C, gen_C$early_interval == "Furongian" | gen_C$late_interval == "Furongian")
#gen_C <- subset(gen_C, gen_C$accepted_rank == "genus" | gen_C$accepted_rank == "species")
#gen_C <- as.data.frame(gen_C,stringsAsFactors=FALSE)
# write.csv(gen_C, file = "gen_C_erates_lat.csv")
gen_C <- read.table(file = "gen_C_erates_lat.csv", sep=',', header=T)

#gen_S <- read.table("http://paleobiodb.org/data1.2/occs/list.txt?interval=Silurian&show=stratext,geo,loc", sep=',', header=T)
#gen_S <- subset(gen_S, gen_S$early_interval == "Rhuddanian")
#gen_S <- subset(gen_S, gen_S$accepted_rank == "genus" | gen_S$accepted_rank == "species")
# write.csv(gen_S, file = "gen_S_erates_lat.csv")
gen_S <- read.table(file = "gen_S_erates_lat.csv", sep=',', header=T)

# gen_O <- subset(ordo_res, ordo_res$accepted_rank == "genus" | ordo_res$accepted_rank == "species")
# gen_O <- subset(gen_O, gen_O$ts_count==1) # subset of accuracy == 1 time bin
# gen_O <- subset(gen_O, ordo_res$identified_no>0)
# write.csv(gen_O, file = "gen_O_erates_lat.csv")
gen_O <- read.table(file = "gen_O_erates_lat.csv", sep=',', header=T)

##########################
##########################
###preparation of data tables

g.Ca <- matrix("Camb",nrow(gen_C), 1 )
g.Cb <- matrix(,nrow(gen_C), 1 )
g.Cc <- as.character(gen_C$collection_no)
g.c <-cbind(g.Cb, g.Ca, g.Cc)
for (i in 1:nrow(gen_C)) {
  if(grepl(" ",gen_C$accepted_name[i])){g.c[i,1] <- strsplit(as.character(gen_C$accepted_name[i])," ")[[1]][1]}
  else {g.c[i,1] <- as.character(gen_C$accepted_name[i])}
}

g.sa <- matrix("Sil",nrow(gen_S), 1 )
g.sb <- matrix(,nrow(gen_S), 1 )
g.sc <- as.character(gen_S$collection_no)
g.s <-cbind(g.sb, g.sa, g.sc)
for (i in 1:nrow(gen_S)) {
  if(grepl(" ",gen_S$accepted_name[i])){g.s[i,1] <- strsplit(as.character(gen_S$accepted_name[i])," ")[[1]][1]}
  else {g.s[i,1] <- as.character(gen_C$accepted_name[i])}
}

g.oa <- as.character(gen_O$oldest)
g.ob <- matrix(,nrow(gen_O), 1 )
g.oc <- as.character(gen_O$collection_no)
g.o <-cbind(g.ob, g.oa, g.oc)
for (i in 1:nrow(gen_O)) {
  if(grepl(" ",gen_O$accepted_name[i])){ g.o[i,1] <- strsplit(as.character(gen_O$accepted_name[i])," ")[[1]][1]}
  else {g.o[i,1] <- as.character(gen_O$accepted_name[i])}
}

g.c <-as.data.frame(g.c)
g.s <-as.data.frame(g.s)
g.o <-as.data.frame(g.o)
colnames(g.c) <- c("accepted_name", "ts", "collection_no")
names(g.s) <- c("accepted_name", "ts", "collection_no")
names(g.o) <- c("accepted_name", "ts", "collection_no")

occs2 <- rbind(unique(g.c), unique(g.s), unique(g.o))
occs2 <- na.omit(occs2)
occs2 <- drop.levels(occs2)

ordo_resl <-  merge(occs2, coll_tots, by='collection_no')
ordo_resl$paleolat <- as.numeric(as.character(ordo_resl$paleolat))

ooclat1 <- subset(ordo_resl, abs(ordo_resl$paleolat)<=15 )
ooclat2 <- subset(ordo_resl, abs(ordo_resl$paleolat)>15&abs(ordo_resl$paleolat)<=30)
ooclat3 <- subset(ordo_resl, abs(ordo_resl$paleolat)>30&abs(ordo_resl$paleolat)<=45)
ooclat4 <- subset(ordo_resl, abs(ordo_resl$paleolat)>45&abs(ordo_resl$paleolat)<=60)
ooclat1 <- drop.levels(ooclat1)
ooclat2 <- drop.levels(ooclat2)
ooclat3 <- drop.levels(ooclat3)
ooclat4 <- drop.levels(ooclat4)

of <- array()
of <- within(ordo_resl, paleolat[abs(ordo_resl$paleolat)<=15] <- 1)
of <- within(of, paleolat[abs(of$paleolat)>15&abs(of$paleolat)<=30] <- 2)
of <- within(of, paleolat[abs(of$paleolat)>30&abs(of$paleolat)<=45] <- 3)
of <- within(of, paleolat[abs(of$paleolat)>45&abs(of$paleolat)<=60] <- 4)
of <- of[which(of$paleolat>=1 & of$paleolat<=4),]

ots <- matrix(0,1,ncol(of))
colnames(ots) <- colnames(of)
bin.name <- array(, dim =nrow(xts)-2);
ol <- list()
lt <- c(1:4); # four latidunal zones
min_col <- array(,nrow(xts)-2,)

##########################
##########################
########## counts

######### counts general numbers
bc1 <- bcounts(ooclat1)
bc2 <- bcounts(ooclat2)
bc3 <- bcounts(ooclat3)
bc4 <- bcounts(ooclat4)
bca <- rbind(bc1, bc2, bc3, bc4)
bcc <- bcounts(ordo_resl)

############ counts of minimal genus counts in Latitudinal zones
for (i in 2:(nrow(xts)-1)) 
{
  cns <- as.character(xts[i,1])
  bin.name[i-1] <- cns
  ots <- subset(of, grepl(cns, of$ts))
  ox <- table(ots$accepted_name, ots$paleolat)
  ol[[i-1]] <- ox
}

for (i in 1:(nrow(xts)-2)) 
{
  rx <- ol[[i]]
  rx <- replace(rx, rx>0, 1)
  min_col[i] <-min(colSums(rx))
}
mmcol <- min(min_col)

#######################
# calculation of minimal value for bootstrapping

ri.min <- numeric()  
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
  me.12 <- merge(ots1.u, ots.u, by='accepted_name')
  ri.min[i-1] <- nrow(me.12)
}
ri.min <- min(ri.min)

##########################
##########################
########## calculations

### rarefied diversity trends
r.rt <- rd.fun(ordo_resl, 500)

r1 <- rd.fun(ooclat1, mmcol)
r2 <- rd.fun(ooclat2, mmcol)
r3 <- rd.fun(ooclat3, mmcol)
r4 <- rd.fun(ooclat4, mmcol)
rares <- rbind(r1,r2, r3, r4)


### rarefied diversity LDG's

rl <- array(, c(nrow(xts)-2,length(lt),4)); # array 1: ts; 2: zone, 3: values per zone (1, gen; 2, occs, 3, rare, 4rareci)
for (i in 1:(nrow(xts)-2)) 
{
  om <- ol[[i]]
  rl2 <- rarplus(om, lt, mmcol)
  for (k in 1:ncol(rl2)) ### durch zonen
  {
    rl[i,k,] <- rl2[,k]
  }
}

### alpha, beta, gamma trends

r.gam <- matrix(0,9,7) ### total trend
for (k in 2:(nrow(xts)-1)) 
{
  cns <- as.character(xts[k,1])
  ots <- subset(of, grepl(cns, of$ts))
  otl <- droplevels(ots)
  otl <- na.omit(otl[,c("geoplate","accepted_name")])
  ox <- table(otl$geoplate, otl$accepted_name)
  
  oalpha <- table(otl$accepted_name)
  beta.d <- d(ox, lev="beta", boot="TRUE")
  gamma.d <- d(ox, lev="gamma", boot="TRUE")
  alpha.d <- d(ox, lev="alpha", boot="TRUE")
  
  r.gam[1,k-1] <- sum(oalpha) # nr of occurrences
  r.gam[2,k-1] <- nrow(oalpha) # nr of genera
  r.gam[3,k-1] <- nrow(ox) # nr of plates
  r.gam[4,k-1] <- as.numeric(as.character(gamma.d[1]))
  r.gam[5,k-1] <- 1.95*as.numeric(as.character(gamma.d[2])) # CI
  r.gam[6,k-1] <- as.numeric(as.character(beta.d[1]))
  r.gam[7,k-1] <- 1.95*as.numeric(as.character(beta.d[2])) # CI
  r.gam[8,k-1] <- as.numeric(as.character(alpha.d[1]))
  r.gam[9,k-1] <- 1.95*as.numeric(as.character(alpha.d[2])) # CI
}

rp <- bplusa_tot(of) ## per latzone

### alpha, beta, gamma LDG's

rk <- betaplusgamma(of) # rk[value,LatZone,TS]; value: : 1 occs, 2 gen; 3 plates; 4,5 gamma; 6,7 beta; 8,9 alpha

#######################
# immigrations bootstrapped

ri.m <- numeric()  
ri.ce <- numeric()
ri.e <- numeric()  
ri.ee <- numeric()

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
  me.12 <- merge(ots1.u, ots.u, by='accepted_name')
  me.23 <- merge(ots.u, ots3.u, by='accepted_name')
  me.i <- as.numeric(me.12[,5])-as.numeric(me.12[,3])
  me.i3 <- as.numeric(me.23[,3])-as.numeric(me.23[,5])
  for (k in 1:1000) {
    me.mean <- sample(me.i, size=length(me.i), replace=TRUE)
    me.mean3 <- sample(me.i3, size=length(me.i), replace=TRUE)
    me.q <- qt(0.95,df=length(me.mean)-1)*sd(me.mean)/sqrt(length(me.mean))
    me.q3 <- qt(0.95,df=length(me.mean3)-1)*sd(me.mean3)/sqrt(length(me.mean3))
  }
  ri.m[i-1] <- mean(as.numeric(me.12[,5])-as.numeric(me.12[,3]))
  ri.ce[i-1] <- me.q
  ri.e[i-1] <- mean(as.numeric(me.23[,3])-as.numeric(me.23[,5]))
  ri.ee[i-1] <- me.q3
}

#######################
# evolutionary rates complete set

re <- eratescalc(of) # re[,1,]: gencont, 2 mds, 3 ori, 4 exe, 5 extirps, 6 immis


############################
############################
##########################
###tests

########### spearmans's rank after generalised differencing (McKinney, 1990) 

trend.fun <- function(ipua, ipub) {
  me <- data.frame(matrix(,7,2))
  colnames(me) <- c("a", "b")
  me[,1] <- ipua
  me[,2] <- ipub
  
  me.plus <- rbind(me[1,], me[])
  me.minus <- rbind(me[], me[1,])
  
  qmax.a <- 1/(max(as.numeric(as.character(me$a))))
  qmax.b <- 1/(max(as.numeric(as.character(me$b))))
  malpha.a <- na.omit(qmax.a*(as.numeric(as.character(me.plus[,1]))-as.numeric(as.character(me.minus[,1]))))
  malpha.b <- na.omit(qmax.b*(as.numeric(as.character(me.plus[,2]))-as.numeric(as.character(me.minus[,2]))))
  spear.c <- cor(cbind(malpha.a, malpha.b), method="spearman")
  spear.t <- cor.test(malpha.a, malpha.b, method="spearman")
  resis <- rbind(spear.c, spear.c)
  neList <- list("rho"=spear.c, "p"=spear.t)
  return(neList)
}

#beta/gamma
(trend.fun(rp[6,,1], rp[4,,1]) + trend.fun(rp[6,,2], rp[4,,2]) 
+ trend.fun(rp[6,,3], rp[4,,3]) + abs(trend.fun(rp[6,,4], rp[4,,4]))+ trend.fun(r.gam[8,], r.gam[4,]))/5
#beta/rarefied
(trend.fun(rp[6,,1], r1[,3]) + trend.fun(rp[6,,2], r2[,3]) 
+ trend.fun(rp[6,,3], r3[,3]) + abs(trend.fun(rp[6,,4], r4[,3]))+ trend.fun(r.gam[8,], r.rt[,3]))/5
#gamma/rarefied
(trend.fun(rp[4,,1], r1[,3]) + trend.fun(rp[4,,2], r2[,3]) 
+ trend.fun(rp[4,,3], r3[,3]) + abs(trend.fun(rp[4,,4], r4[,3]))+trend.fun(r.gam[4,], r.rt[,3]))/5

#total
#beta/gamma
trend.fun(r.gam[8,], r.gam[4,])
#beta/rarefied
trend.fun(r.gam[8,], r.rt[,3])
#gamma/rarefied
trend.fun(r.gam[4,], r.rt[,3])
