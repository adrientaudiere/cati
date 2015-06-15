###############################################################################################################################################
# FDchange: R function to compute change in functional diversity following disturbance                                                        # 
#                                                                                                                                             #  
# INPUTS:                                                                                                                                     #
#                                                                                                                                             #
#       - "abundances" = matrix (2 x S) with the abundances of the S species "before" and "after" disturbance                                 #
#                                                                                                                                             #
#       - "coord" = matrix (S x T) with coordinates of the S species present in the "abundances" matrix                                       #
#                           on the T axes of the functional space (NA not allowed)                                                            #
#                                                                                                                                             # 
#        NB: number of species before and after disturbance must be strictly higher than T                                                    #                                
#                                                                                                                                             #
#    R Libraries 'ape', 'geometry' and 'rcdd' need to be downloaded  before running this function                                             #
#                                                                                                                                             #
#                                                                                                                                             #
# OUTPUTS: a list with                                                                                                                        #
#               $ NbSp: a vector (3) with number of species before and after disturbance and change between the 2 periods (after-before)      #
#                                                                                                                                             #
#               $ sp_ab: a matrix (S*3), with species relative abundances (%) before and after disturbance,                                   #
#                         and their status given their changes in abundance : "extirpated", "loser", "unchanged", "winner", "introduced"      #
#                                                                                                                                             #
#               $ FId: a matrix (4*T) with  Community Weighted Mean (CWM) values before and after disturbance for each functional axis        #
#                       delta= change in CWM (after-before), and delta_st= delta CWM standardized by trait range in the species pool          #
#                                                                                                                                             #
#               $ FD: a matrix (3*7) with indices values before and after disturbance and their temporal change (delta=after-before)          #
#                   - FRic: functional richness (expressed as a percentage of the functional space filled by the S species)                   # 
#                   - FEve: functional evenness                                                                                               #
#                   - FDiv: functional divergence                                                                                             #
#                   - FDis: functional dispersion (expressed as a percentage of the maximal pairwise distance observed in the species pool)   #
#                   - FEnt: functional entropy (Rao's quadratic entropy expressed as an equivalent number of species)                         #
#                   - FSpe: functional specialization (expressed as a percentage of the maximal specialization observed in the species pool)  #
#                   - FOri: functional originality (expressed as a percentage of the maximal originality observed in the species pool)        #
#                                                                                                                                             #
#               $ FShift: a single value vector with functional shift expressed as a percentage of non-overlap                                #
#                                       between convex hulls shaping communities before and after disturbance respectively                    #
#                                                                                                                                             #
#               $ details: a list with details about FD indices                                                                               #
#                     - vertices_before and vertices_after: vectors with identities of species being vertices before and after disturbance    #
#                     - mstvect_before and mstvect_after: matrix describing the topology of Minimum Spanning Trees for each period            #
#                     - B_before and B_after: vector(length T) with coordinates of the center of gravity (B) of the vertices for each period  #
#                     - meandDB_before and mean_dB_after: mean distance to B for each period                                                  #
#                     - O: vector (length T) with coordinates of the center of gravity of the species pool                                    #
#                     - speS: vector (length S) with functional specialization of each of the S species                                       #
#                     - NN: matrix 0/1 (S*S) with pairs of nearest neighbours (focal species=rows)                                            #
#                     - oriS: vector (length S) with functional originality of each of the S species                                          #
#                                                                                                                                             #
#  An example of how using this function is provided at this end of the script (cf figure in Box2)                                            #
#                                                                                                                                             #
###############################################################################################################################################


FDchange<-function(coord,abundances)  {

#loading required libraries
require(ape)
require(geometry)
require(rcdd)  
############################################################################################

# identities of species present before and after disturbance
sp_before<-which(abundances[1,]!=0)
sp_after<-which(abundances[2,]!=0)

# number of species

NbSp<-c( length(sp_before), length(sp_after), length(sp_after)-length(sp_before)  )
names(NbSp)<-c("before","after","delta")

# check coord matrix
if (is.numeric(coord)==F) stop ("coord values must be numeric")
if(is.na(sum(coord))) stop ("NA are not allowed in 'coord' matrix")
T<-ncol(coord) # T = number of functional axes
if (dim(abundances)[2]!=dim(coord)[1]) stop(" error : different number of species in 'coord' and 'abundances' matrices ")
if (NbSp["before"]<=T)  stop(paste("error : community 'before' must contain at least ",T+1, " species",sep=""))
if (NbSp["after"]<=T)  stop(paste("error : community 'after' must contain at least ",T+1, " species",sep=""))

############################################################################################

# computing relative abundances
relab<-abundances/apply(abundances,1,sum,na.rm=T)

# status given changes in abundances
status<-rep(NA,nrow(coord)) 
status[which(relab["Before",] == relab["After",])]<-"unchanged"
status[which(relab["Before",] > relab["After",])]<-"loser"
status[which(relab["Before",] <  relab["After",])]<-"winner"
status[which(relab["Before",]!=0 & relab["After",]==0)]<-"extirpated"
status[which(relab["Before",]==0 & relab["After",]!=0)]<-"introduced"

sp_ab<-data.frame(relab_before=round(relab[1,]*100),relab_after=round(relab[2,]*100),status=as.factor(status) )
row.names(sp_ab)<-row.names(coord)


# definition of matrix FId
FId<-matrix(0,4,T,dimnames=list(c("before","after","delta","delta_st"),colnames(coord)) )

# computing changes in community weighted mean trait values
FId["before",]<-relab["Before",sp_before] %*% coord[sp_before,] 
FId["after",]<-relab["After",sp_after] %*% coord[sp_after,]
FId["delta",]<-FId["after",]-FId["before",]
FId["delta_st",]<-FId["delta",] / ( apply(coord,2,max) - apply(coord,2,min) )

###################################################################################################
# function to compute functional richness, evenness and divergence
# ab: relative abundances of species (no 0) ; tr: matrix of species coord
FRicFEveFDiv<-function(ab,tr)  {

# S = number of species
S<-length(ab)
    
# FRIC
FRic<-round(convhulln(tr,"FA")$vol,6)

# identity of vertices
vert0<-convhulln(tr,"Fx TO 'vert.txt'")
vert1<-scan("vert.txt",quiet=T)
vertices<-(vert1+1)[-1]

# FEve
# computation of inter-species euclidian distances
distT<-dist(tr, method="euclidian")

# computation of Minimum Spanning Tree and conversion of the 'mst' matrix into 'dist' class
linkmst<-mst(distT) 
mstvect<-as.dist(linkmst)

# computation of the pairwise cumulative relative abundances and conversion into 'dist' class
ab2<-matrix(0,nrow=S,ncol=S)
for (q in 1:S)
for (r in 1:S)
ab2[q,r]<-ab[q]+ab[r] # end of q,r 
ab2vect<-as.dist(ab2)

# computation of EW for the (S-1) segments
EW<-rep(0,S-1)
flag<-1
for (m in 1:((S-1)*S/2))
{if (mstvect[m]!=0) {EW[flag]<-distT[m]/(ab2vect[m]) ; flag<-flag+1}}  # end of m

# computation of the PEW and comparison with 1/S-1, finally computation of FEve
minPEW<-rep(0,S-1)  ;  OdSmO<-1/(S-1)
for (l in 1:(S-1))
minPEW[l]<-min( (EW[l]/sum(EW)) , OdSmO )  # end of l
            
FEve<-round( ( (sum(minPEW))- OdSmO) / (1-OdSmO ) ,6)


# FDiv
# coord values of vertices
trvertices<-tr[vertices,]

# coordinates of the center of gravity of the vertices (B)
B<-apply(trvertices,2,mean)

# computing euclidian dstances to B (dB)
dB<-apply(tr, 1, function(x) { (sum((x-B)^2) )^0.5} )

# mean of dB values, deviations to mean and relative abundances-weighted mean deviation
meandB<-mean(dB)
devdB<-dB-meandB
abdev<-ab*devdB
ababsdev<-ab*abs(devdB)

# computation of FDiv
FDiv<-round( (sum(abdev)+meandB) / (sum(ababsdev)+meandB) ,6)
            
RED<-c(FRic=FRic, FEve=FEve, FDiv=FDiv, details=list(vertices=vertices, mst=linkmst, B=B, meandB=meandB) ) 
invisible(RED)              
} # end of function FRicFEveFDiv
############################################################################################

# definition of matrix FD
FD<-matrix(0,7,3,dimnames=list(c("FRic","FEve","FDiv","FDis","FEnt","FSpe","FOri"),c("before","after","delta")) )


# computing convex hull volume of the total pool of species
CHV_pool<-as.numeric( FRicFEveFDiv(ab=rep(1,nrow(coord)), tr=coord)["FRic"] )


# computing FRic FEve and FDiv before and after disturbance
F_RED_before<-FRicFEveFDiv( relab["Before",sp_before] , coord[sp_before,]  )
F_RED_after<-FRicFEveFDiv( relab["After",sp_after] , coord[sp_after,]  )
FD["FRic",c("before","after")]<-c(F_RED_before$FRic , F_RED_after$FRic) / CHV_pool
FD["FEve",c("before","after")]<-c(F_RED_before$FEve, F_RED_after$FEve)
FD["FDiv",c("before","after")]<-c(F_RED_before$FDiv , F_RED_after$FDiv)
vertices.before<-row.names(coord[sp_before,])[F_RED_before$details.vertices]
vertices.after<-row.names(coord[sp_after,])[F_RED_after$details.vertices]

# FDis: abundance-weighted mean distance to abundance-weighted centroid 
FDis_before<-relab["Before",sp_before] %*% apply(coord[sp_before,], 1, function(x) { (sum((x-FId["before",])^2) )^0.5} ) 
FDis_after<-relab["After",sp_after] %*% apply(coord[sp_after,], 1, function(x) { (sum((x-FId["after",])^2) )^0.5} ) 
FD["FDis",c("before","after")]<-c(FDis_before , FDis_after)/ max(dist(coord))
 
 
# standardized distance between species
dist2<- as.matrix(dist(coord)) / max(dist(coord))

# Rao's quadratic entropy (Q): sum of abundance-weighted distances before and after disturbance
d2ab_before<-relab["Before",sp_before] %*% dist2[sp_before,sp_before] %*% relab["Before",sp_before]
d2ab_after<-relab["After",sp_after] %*% dist2[sp_after,sp_after] %*% relab["After",sp_after]

# Q expressed as an equivalent number of species 
Q_before<- 1/(1-d2ab_before)
Q_after<- 1/(1-d2ab_after)
FD["FEnt",c("before","after")]<-c(Q_before , Q_after)
 
 
# computing functional specialization of all species : distance to hypothetical mean species
O<-apply(coord, 2, mean)
speS<-apply(coord, 1, function(x) { (sum((x-O)^2) )^0.5} ) 

# computing functional specialization before and after disturbance
# values are standardized by dividing by the maximal value of specialization
FD["FSpe","before"]<-speS[sp_before] %*% relab["Before",sp_before] / max(speS)
FD["FSpe","after"]<-speS[sp_after] %*% relab["After",sp_after]  / max(speS)


# computing functional originality of all species and identity of nearest neighbour
dist_F<-as.matrix(dist(coord,method="euclidean")) ; dist_F[which(dist_F==0)]<-NA
oriS<-apply(dist_F, 1, min, na.rm=T )
NN<-dist_F ; NN<-NN-apply(NN,1,min,na.rm=T) ; NN[which(NN!=0)]<-NA   ; NN[which(NN==0)]<-1


# computing functional specialization and originality before and after disturbance
# values are standardized by dividing by the maximal value of originality
FD["FOri","before"]<-oriS[sp_before] %*% relab["Before",sp_before] / max(oriS)
FD["FOri","after"]<-oriS[sp_after] %*% relab["After",sp_after]  / max(oriS)

FD[,"delta"]<-FD[,"after"]-FD[,"before"]
############################################################################################
# compute functional shift between periods
set1<-coord[sp_before,]
set2<-coord[sp_after,]

# tranforming species coordinates in the multidimensional space in true rational number written as character string
# reduce set of points to vertices only using redundant function
# changing polytope representation: vertices to inequality constraints
set1rep <- d2q(cbind(0, cbind(1, set1)))
polytope1 <- redundant(set1rep, representation = "V")$output
H_chset1 <- scdd(polytope1, representation = "V")$output

set2rep <- d2q(cbind(0, cbind(1, set2)))
polytope2 <- redundant(set2rep, representation = "V")$output
H_chset2 <- scdd(polytope2, representation = "V")$output

# intersection between the two polytopes
H_inter <- rbind(H_chset1, H_chset2)
V_inter <- scdd(H_inter, representation = "H")$output

# extracting coordinates of vertices
vert_1n2 <- q2d(V_inter[ , - c(1, 2)])

# computing convex hull volume of the intersection (if it exists)
vol_inter<-0
if (is.matrix(vert_1n2)) # vector if one vertex in common
    if( nrow(vert_1n2)>ncol(vert_1n2) ) vol_inter<-convhulln(vert_1n2,"FA")$vol

FShift<- 1- (vol_inter / CHV_pool )
############################################################################################

# grouping results
res<-list(NbSp=NbSp, sp_ab=sp_ab, FId=FId, FD=FD, FShift=FShift, 
          detailsFD=list(vertices_before=vertices.before, vertices_after=vertices.after,
                       mst_before=F_RED_before$details.mst, mst_after=F_RED_after$details.mst,
                       B_before=F_RED_before$details.B, B_after=F_RED_after$details.B,
                       meandB_before=F_RED_before$details.meandB, meandB_after=F_RED_after$details.meandB,
                       O=O, speS=speS, NN=NN, oriS=oriS)  
         )
invisible(res)

}# end of function FDchange


