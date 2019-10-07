#### Small-scale fishing modelling using Bayesian networks (a.k.a. belief propagation networks)

## Author: Phillip P.A. Staniczenko  pstaniczenko@brooklyn.cuny.edu
## Please cite: 
## Social ties explain catch portfolios of small-scale fishers in the Caribbean
## Alexander, S.M., Staniczenko, P.P.A. & Bodin, O.
## Fish and Fisheries, vol. X, ppXX--XX (2019)

## Date: 7th October 2019

#### Notes
## Before running this R script, you must first compile for your specific computer the C source file
## fNMLofDAG.c
## to generate the executable file
## fNMLGet
##
## To run this R script, ensure all required files are in your working directory, then in R type
## source("SSF_BNAnalysis_Staniczenko_distribute.R")

#### Outputs
## Total lengths for model-data combinations;
## prints results to screen and also saves results as a csv file in the working directory

#### Required program files
## ZipData.py
## fNMLGet executable (must be compiled for you computer from source file, fNMLofDAG.c)

#### Required data files
## manID.csv
## manTrait_metadata.csv (heirarchy is column 3)
## fishID.csv
## fishTrait_metadata.csv
## catchportfolio_manfishedgelist.csv
## socialnetwork_edgelist.csv

# Clear workspace
rm(list=ls(all=TRUE))

## User settings

myfishgroup_BN <- "fishonegroup"
# "fishonegroup", "sizeclass", "pricecateg", "fisheachgroup", "sizeclass-pricecateg"
# pricecateg is column 2 (e.g., 1=cheap, 2=expensive)
# sizeclass is column 3 (e.g., 1=small, 2=medium, 3=large)
# fisheachgroup is column 1

mymangroup_BN <- "manonegroup"
# "manonegroup", "geartype", "landingsite", "geartype-landingsite", "maneachgroup"
# geartype is column 4 (e.g., 1,2,3)
# landingsite is column 5 (e.g., 1,2)
# maneachgroup is column 1

BN_model <- 2  # 0: AND
			   # 1: OR
			   # 2: SUM
			   # 3: FULL
# For BNs, the DAG in a FULL model should have a maximum of 32 columns

# Number of randomised BNs to try
numrandomisations <- 10

#### Libraries
require(gsl)
library(igraph)
					 
#### Data processing

# Load data for fish and man
manID <- read.csv("manID.csv", header=FALSE, sep=",", colClasses=c("numeric","character"))
fishID <- read.csv("fishID.csv", header=FALSE, sep=",", colClasses=c("numeric","character"))
manfishedgelist <- read.csv("catchportfolio_manfishedgelist.csv", header=FALSE, sep=",", colClasses=c("character","character"))
#
manTrait <- read.csv("manTrait_metadata.csv", header=TRUE, sep=",",colClasses="character")
fishTrait <- read.csv("fishTrait_metadata.csv", header=TRUE, sep=",",colClasses="character")

REORDER <- as.numeric(manTrait[,3])
REORDER <- sort(REORDER, decreasing=TRUE, index.return=TRUE)
REORDER <- REORDER$ix

# Get bipartite network from edge list
mybipartitematrix <- matrix(rep(0, dim(fishID)[1] * dim(manID)[1]), nrow=dim(fishID)[1], ncol=dim(manID)[1])
# fill matrix
for (myrow in 1:dim(manfishedgelist)[1]){
  mybipartitematrix[fishID[which(fishID[,2]==manfishedgelist[myrow, 2], arr.ind=TRUE), 1], manID[which(manID[,2]==manfishedgelist[myrow, 1], arr.ind=TRUE), 1]] <- 1
}

# REORDER everything
mybipartitematrix <- mybipartitematrix[, REORDER]
manID <- manID[REORDER, ]
manTrait <- manTrait[REORDER, ]

numrows <- dim(mybipartitematrix)[1]
numcols <- dim(mybipartitematrix)[2]
							  
write.table(mybipartitematrix, file="mybipartitematrix.txt", row.names=FALSE, col.names=FALSE, sep=" ")

# load social network data
socialnetworkedgelist <- read.csv("socialnetwork_edgelist.csv", header=FALSE, sep=",", colClasses=c("character","character"))

# Get social network adjacency matrix from edge list
socialnetwork <- matrix(rep(0, dim(mybipartitematrix)[2] * dim(mybipartitematrix)[2]), nrow=dim(mybipartitematrix)[2], ncol=dim(mybipartitematrix)[2])
# fill matrix
for (myrow in 1:dim(socialnetworkedgelist)[1]){
  socialnetwork[manID[which(manID[,2]==socialnetworkedgelist[myrow, 1], arr.ind=TRUE), 1], manID[which(manID[,2]==socialnetworkedgelist[myrow, 2], arr.ind=TRUE), 1]] <- 1
  socialnetwork[manID[which(manID[,2]==socialnetworkedgelist[myrow, 2], arr.ind=TRUE), 1], manID[which(manID[,2]==socialnetworkedgelist[myrow, 1], arr.ind=TRUE), 1]] <- 1
}

# REORDER
socialnetwork <- socialnetwork[REORDER, REORDER]

# No social interactions
socialnetwork_empty <- matrix(rep(0, dim(mybipartitematrix)[2] * dim(mybipartitematrix)[2]), nrow=dim(mybipartitematrix)[2], ncol=dim(mybipartitematrix)[2])

bipartitetoadjacency <- function(mybipartitematrix){
	numrows <- dim(mybipartitematrix)[1]
	numcols <- dim(mybipartitematrix)[2]
	adj1 <- cbind(matrix(rep(0, numrows * numrows), numrows, numrows), mybipartitematrix)
	adj2 <- cbind(t(mybipartitematrix), matrix(rep(0, numcols * numcols), numcols, numcols))
	adj <- rbind(adj1, adj2)
	return(adj)
}
							  
myadjacencymatrix <- bipartitetoadjacency(mybipartitematrix)

#### Functions

# Compute complexity of random graph of size X
Complexity <- function(X){
  if (X == 0){
    return (0.0)
  }
  # for small X, use the exact value
  XX = X
  if (XX < 100.){
    return (log2(exp(XX + log(gamma_inc(XX,XX)) - (XX - 1.0) * log(XX)) + 1.0))
  }
  # for large X, use the approximation
  return (log2(1. + sqrt(pi * XX / 2.) - 1./3. + sqrt(2. * pi) / (24. * sqrt(XX)) 
               - 4. / (135. * XX) +  sqrt(2. * pi) / (576. * sqrt(XX * XX * XX)) +
                 8. / (2835. * XX * XX)))
}

# Compute NML (really, total length) for a slice containing Ones ones and Zeros zeros
NMLSlice <- function(Ones, Zeros){
  NML = Complexity(Ones + Zeros)
  p = 0.0
  # if either Ones or Zeros = 0, the log_likelihood is 0
  if (Ones * Zeros > 0){
    # compute maximum log likelihood
    p = Ones / (Ones + Zeros)
    NML = NML - Ones * log2(p) - Zeros * log2(1.0 - p)
  }
  return (NML)
}

# Compute NML (really, total length) for a matrix-SBM combination
getNML <- function(mybipartitematrix, rowgroups, colgroups){
  # iterate over groups
  NMLmatrix <- 0.0
  for (i in 1:length(rowgroups)){
  	for (j in 1:length(colgroups)){
  	  Ones <- 0
  	  mycounter <- 0
  	  # iterate over elements within a rowgroup-colgroup pair
  	  for (k in 1:length(rowgroups[[i]])){
  	  	for (l in 1:length(colgroups[[j]])){
  	  	  mycounter <- mycounter + 1
  	  	  Ones <- Ones + mybipartitematrix[rowgroups[[i]][k],colgroups[[j]][l]]
  	  	}
  	  }
  	  NMLmatrix <- NMLmatrix + NMLSlice(Ones, mycounter - Ones)
  	}
  }
  return(NMLmatrix)
}

# Check if graph is DAG
isDAG <- function(myDAG){
  M <- myDAG
  M[is.na(M)] <- 0
  M[M == -1] = 1
  graph <- graph.adjacency(M, mode="directed",)
  return(is.dag(graph))
}

# Get DAG in format for fNMLGet executable
getreducedDAGform <- function(myDAG, myDAGorder){
  numnodes <- length(myDAGorder)
  reducedDAGform <- numnodes - myDAGorder
  errorcheck <- c(0,0)
  for (i in numnodes:2){
  	reducedDAGform <- c(reducedDAGform, rev(myDAG[1:(i - 1), i]))
  	errorcheck <- c(errorcheck, rev(myDAG[1:(i - 1), i]))
  }
  #print(sum(errorcheck))
  write.table(reducedDAGform, file="myDAGreduced.txt", row.names=FALSE, col.names=FALSE, sep=" ")
  return(0)
}

gettotallengthforblockDAGcombination <- function(mybipartitematrix, rownodes, colnodes, socialnetworkbasis, DAGorderingbasis){
	if (length(colnodes)==1){  # only one fisherman
		TotalLength <- NMLSlice(sum(mybipartitematrix[rownodes, colnodes]), length(mybipartitematrix[rownodes, colnodes]) - sum(mybipartitematrix[rownodes, colnodes]))
	} else {
		mybipartitematrix_reduced <- mybipartitematrix[rownodes, colnodes]
		DAGorderingbasis_reduced <- DAGorderingbasis[colnodes,]
		socialnetworkbasis_reduced <- socialnetworkbasis[colnodes, colnodes]
		# build DAG
		if (sum(socialnetworkbasis_reduced) > 0){
			for (k in 1:dim(socialnetworkbasis_reduced)[1]){
				for (l in 1:dim(socialnetworkbasis_reduced)[2]){
					if (socialnetworkbasis_reduced[k, l] > 0){
						if (as.numeric(DAGorderingbasis_reduced[k, 2]) >= as.numeric(DAGorderingbasis_reduced[l, 2])){
							socialnetworkbasis_reduced[k, l] <- 0
						}
					}
				}
			}
		}
		myDAGcols <- socialnetworkbasis_reduced
		if (length(rownodes==1)){  # only one fish
			myDAGorder <- 1:length(colnodes)
		} else {
			myDAGorder <- 1:dim(mybipartitematrix_reduced)[2]
		}
		
		if (is.null(dim(mybipartitematrix_reduced))==TRUE){
			mybipartitematrix_reduced <- matrix(c(mybipartitematrix_reduced, 1), nrow=1, ncol=length(mybipartitematrix_reduced)+1)
			write.table(mybipartitematrix_reduced, file="mybipartitematrix_reduced.txt", row.names=FALSE, col.names=FALSE, sep=" ")
			# write DAG in reduced form
    		getreducedDAGform(myDAGcols, myDAGorder)
    		# call fNMLGet executable
    		zipmatrixnumrows <- "1"
    		zipmatrixnumcols <- toString(length(mybipartitematrix_reduced) - 1)
    		myexpression <- paste("./fNMLGet ", zipmatrixnumrows, " ", zipmatrixnumcols, " ", "mybipartitematrix_reduced.txt", " myDAGreduced.txt ", toString(BN_model), " 123 100 100", sep="")
    		system(myexpression, ignore.stdout=TRUE, ignore.stderr=TRUE)
    		mypattern <- "mybipartitematrix_reduced.txt"
    		myfiles <- list.files(path=".", pattern=mypattern, full.names=FALSE)
    		BN_model_results <- unlist(regmatches(myfiles[2], gregexpr('\\(?[0-9,.]+', myfiles[2])))
    		TotalLength <- as.numeric(BN_model_results[2])
		} else {
			write.table(mybipartitematrix_reduced, file="mybipartitematrix_reduced.txt", row.names=FALSE, col.names=FALSE, sep=" ")
			# write bipartite matrix in zipped format
    		system("python ZipData.py mybipartitematrix_reduced.txt")
    		# write DAG in reduced form
    		getreducedDAGform(myDAGcols, myDAGorder)
    		# call fNMLGet executable
    		mypattern <- "mybipartitematrix_reduced.txt-z-"
    		myfiles <- list.files(path=".", pattern=mypattern, full.names=FALSE)
    		zipmatrixnumrows <- unlist(regmatches(myfiles[1], gregexpr('\\(?[0-9,]+', myfiles[1])))[1]
    		zipmatrixnumcols <- unlist(regmatches(myfiles[1], gregexpr('\\(?[0-9,]+', myfiles[1])))[2]
    		myexpression <- paste("./fNMLGet ", zipmatrixnumrows, " ", zipmatrixnumcols, " ", myfiles[1], " myDAGreduced.txt ", toString(BN_model), " 123 100 100", sep="")
    		system(myexpression, ignore.stdout=TRUE, ignore.stderr=TRUE)
    		myfiles <- list.files(path=".", pattern=mypattern, full.names=FALSE)
    		BN_model_results <- unlist(regmatches(myfiles[2], gregexpr('\\(?[0-9,.]+', myfiles[2])))
    		TotalLength <- as.numeric(BN_model_results[4])
    	}
	
		# Clean up files
  		system("rm mybipartitematrix_reduced.txt*")
  		system("rm myDAGreduced.txt")
  	}
  	return(TotalLength)
}

gettotallengthforblockDAGcombination_plusrandomBNs <- function(mybipartitematrix, rownodes, colnodes, socialnetworkbasis, DAGorderingbasis, numrandomisations){
	TLforrandomBNs <- -99
	if (length(colnodes)==1){  # only one fisherman
		TotalLength <- NMLSlice(sum(mybipartitematrix[rownodes, colnodes]), length(mybipartitematrix[rownodes, colnodes]) - sum(mybipartitematrix[rownodes, colnodes]))
	} else {
		mybipartitematrix_reduced <- mybipartitematrix[rownodes, colnodes]
		DAGorderingbasis_reduced <- DAGorderingbasis[colnodes,]
		socialnetworkbasis_reduced <- socialnetworkbasis[colnodes, colnodes]
		# build DAG
		if (sum(socialnetworkbasis_reduced) > 0){
			for (k in 1:dim(socialnetworkbasis_reduced)[1]){
				for (l in 1:dim(socialnetworkbasis_reduced)[2]){
					if (socialnetworkbasis_reduced[k, l] > 0){
						if (as.numeric(DAGorderingbasis_reduced[k, 2]) <= as.numeric(DAGorderingbasis_reduced[l, 2])){
							socialnetworkbasis_reduced[k, l] <- 0
						}
					}
				}
			}
		}
		myDAGcols <- socialnetworkbasis_reduced
		if (length(rownodes==1)){  # only one fish
			myDAGorder <- 1:length(colnodes)
		} else {
			myDAGorder <- 1:dim(mybipartitematrix_reduced)[2]
		}
		
		if (is.null(dim(mybipartitematrix_reduced))==TRUE){
			mybipartitematrix_reduced <- matrix(c(mybipartitematrix_reduced, 1), nrow=1, ncol=length(mybipartitematrix_reduced)+1)
			write.table(mybipartitematrix_reduced, file="mybipartitematrix_reduced.txt", row.names=FALSE, col.names=FALSE, sep=" ")
			# write DAG in reduced form
    		getreducedDAGform(myDAGcols, myDAGorder)
    		# call fNMLGet executable
    		zipmatrixnumrows <- "1"
    		zipmatrixnumcols <- toString(length(mybipartitematrix_reduced) - 1)
    		myexpression <- paste("./fNMLGet ", zipmatrixnumrows, " ", zipmatrixnumcols, " ", "mybipartitematrix_reduced.txt", " myDAGreduced.txt ", toString(BN_model), " 123 100 100", sep="")
    		system(myexpression, ignore.stdout=TRUE, ignore.stderr=TRUE)
    		mypattern <- "mybipartitematrix_reduced.txt"
    		myfiles <- list.files(path=".", pattern=mypattern, full.names=FALSE)
    		BN_model_results <- unlist(regmatches(myfiles[2], gregexpr('\\(?[0-9,.]+', myfiles[2])))
    		TotalLength <- as.numeric(BN_model_results[2])
		} else {
			write.table(mybipartitematrix_reduced, file="mybipartitematrix_reduced.txt", row.names=FALSE, col.names=FALSE, sep=" ")
			# write bipartite matrix in zipped format
    		system("python ZipData.py mybipartitematrix_reduced.txt")
    		# write DAG in reduced form
    		getreducedDAGform(myDAGcols, myDAGorder)
    		# call fNMLGet executable
    		mypattern <- "mybipartitematrix_reduced.txt-z-"
    		myfiles <- list.files(path=".", pattern=mypattern, full.names=FALSE)
    		zipmatrixnumrows <- unlist(regmatches(myfiles[1], gregexpr('\\(?[0-9,]+', myfiles[1])))[1]
    		zipmatrixnumcols <- unlist(regmatches(myfiles[1], gregexpr('\\(?[0-9,]+', myfiles[1])))[2]
    		myexpression <- paste("./fNMLGet ", zipmatrixnumrows, " ", zipmatrixnumcols, " ", myfiles[1], " myDAGreduced.txt ", toString(BN_model), " 123 100 100", sep="")
    		system(myexpression, ignore.stdout=TRUE, ignore.stderr=TRUE)
    		myfiles <- list.files(path=".", pattern=mypattern, full.names=FALSE)
    		BN_model_results <- unlist(regmatches(myfiles[2], gregexpr('\\(?[0-9,.]+', myfiles[2])))
    		TotalLength <- as.numeric(BN_model_results[4])
    	}
	
		# Clean up files
  		system("rm mybipartitematrix_reduced.txt*")
  		system("rm myDAGreduced.txt")
  		
  		# Now get total lengths for randomised BNs
  		if (sum(socialnetworkbasis_reduced)>0){
  			TLforrandomBNs <- rep(0, numrandomisations)
  			numDAGedges <- sum(socialnetworkbasis_reduced)
  			numuppertri <- as.integer(0.5*dim(socialnetworkbasis_reduced)[1]*(dim(socialnetworkbasis_reduced)[1]-1))
  			for (myrandomisation in 1:numrandomisations){
  				socialnetworkbasis_reduced_random <- socialnetworkbasis_reduced * 0
  				uppertrivalues <- c(rep(1, numDAGedges), rep(0, numuppertri-numDAGedges))
  				uppertrivalues <- sample(uppertrivalues, length(uppertrivalues))
  				socialnetworkbasis_reduced_random[upper.tri(socialnetworkbasis_reduced_random)] <- uppertrivalues
  				myDAGcols <- socialnetworkbasis_reduced_random
  				if (length(rownodes==1)){  # only one fish
					myDAGorder <- 1:length(colnodes)
				} else {
					myDAGorder <- 1:dim(mybipartitematrix_reduced)[2]
				}
  				# now have randomised BN
  				if (dim(mybipartitematrix_reduced)[1]==1){
  					write.table(mybipartitematrix_reduced, file="mybipartitematrix_reduced.txt", row.names=FALSE, col.names=FALSE, sep=" ")
  					# write DAG in reduced form
    				getreducedDAGform(myDAGcols, myDAGorder)
    				# call fNMLGet executable
    				zipmatrixnumrows <- "1"
    				zipmatrixnumcols <- toString(length(mybipartitematrix_reduced) - 1)
    				myexpression <- paste("./fNMLGet ", zipmatrixnumrows, " ", zipmatrixnumcols, " ", "mybipartitematrix_reduced.txt", " myDAGreduced.txt ", toString(BN_model), " 123 100 100", sep="")
    				system(myexpression, ignore.stdout=TRUE, ignore.stderr=TRUE)
    				mypattern <- "mybipartitematrix_reduced.txt"
    				myfiles <- list.files(path=".", pattern=mypattern, full.names=FALSE)
    				BN_model_results <- unlist(regmatches(myfiles[2], gregexpr('\\(?[0-9,.]+', myfiles[2])))
    				TotalLength_randomBN <- as.numeric(BN_model_results[2])
    				system("rm mybipartitematrix_reduced.txt*")
  					system("rm myDAGreduced.txt")
  					TLforrandomBNs[myrandomisation] <- TotalLength_randomBN
  				} else {
  					write.table(mybipartitematrix_reduced, file="mybipartitematrix_reduced.txt", row.names=FALSE, col.names=FALSE, sep=" ")
					# write bipartite matrix in zipped format
    				system("python ZipData.py mybipartitematrix_reduced.txt")
    				# write DAG in reduced form
    				getreducedDAGform(myDAGcols, myDAGorder)
    				# call fNMLGet executable
    				mypattern <- "mybipartitematrix_reduced.txt-z-"
    				myfiles <- list.files(path=".", pattern=mypattern, full.names=FALSE)
    				zipmatrixnumrows <- unlist(regmatches(myfiles[1], gregexpr('\\(?[0-9,]+', myfiles[1])))[1]
    				zipmatrixnumcols <- unlist(regmatches(myfiles[1], gregexpr('\\(?[0-9,]+', myfiles[1])))[2]
    				myexpression <- paste("./fNMLGet ", zipmatrixnumrows, " ", zipmatrixnumcols, " ", myfiles[1], " myDAGreduced.txt ", toString(BN_model), " 123 100 100", sep="")
    				system(myexpression, ignore.stdout=TRUE, ignore.stderr=TRUE)
    				myfiles <- list.files(path=".", pattern=mypattern, full.names=FALSE)
    				BN_model_results <- unlist(regmatches(myfiles[2], gregexpr('\\(?[0-9,.]+', myfiles[2])))
    				TotalLength_randomBN <- as.numeric(BN_model_results[4])
    				system("rm mybipartitematrix_reduced.txt*")
  					system("rm myDAGreduced.txt")
  					TLforrandomBNs[myrandomisation] <- TotalLength_randomBN
  				}
  			}
  		}	
  	}
  	if (sum(TLforrandomBNs)==-99){
  		TLforrandomBNs <- TotalLength
  	}
  	el <- new.env()
  	assign("TotalLength", TotalLength, envir=el)
  	assign("TotalLength_randomBN", TLforrandomBNs, envir=el)
  	return(as.list(el))
}

#### Main

print("MODELLING SMALL-SCALE FISHING")
beginTime <- Sys.time()

# Assign nodes in rows into groups -- fish
if (myfishgroup_BN=="fishonegroup"){
	groupbytrait <- "0"
	numrowgroups <- 1
  	rowgroups <- vector("list", numrowgroups)
   	rowgroups[[1]] <- 1:dim(mybipartitematrix)[1]
} else if (myfishgroup_BN=="pricecateg"){
    groupingtrait <- fishTrait[, 2]
  	groupbytrait <- unique(groupingtrait)
    numrowgroups <- length(groupbytrait)
    rowgroups <- vector("list", numrowgroups)
    for (i in 1:numrowgroups){
    	rowgroups[[i]] <- which(groupingtrait==groupbytrait[i], arr.ind=TRUE)
    }
} else if (myfishgroup_BN=="sizeclass"){
   	groupingtrait <- fishTrait[, 3]
	groupbytrait <- unique(groupingtrait)
   	numrowgroups <- length(groupbytrait)
   	rowgroups <- vector("list", numrowgroups)
   	for (i in 1:numrowgroups){
   		rowgroups[[i]] <- which(groupingtrait==groupbytrait[i], arr.ind=TRUE)
   	}
} else if (myfishgroup_BN=="sizeclass-pricecateg"){
    groupingtrait1 <- fishTrait[, 3]
    groupingtrait2 <- fishTrait[, 2]
    groupingtraitcombo <- paste(groupingtrait1, "+", groupingtrait2, sep="")
    groupbytrait <- unique(groupingtraitcombo)
    numrowgroups <- length(groupbytrait)
    rowgroups <- vector("list", numrowgroups)
    for (i in 1:numrowgroups){
    	rowgroups[[i]] <- which(groupingtraitcombo==groupbytrait[i], arr.ind=TRUE)
    }
    print ("fish sizeclass-pricecateg groups (1=small, 2=medium, 3=large; + ; 1=cheap, 2=expensive)):")
    print (groupbytrait)
} else if (myfishgroup_BN=="fisheachgroup"){
   	groupingtrait <- fishTrait[, 1]
	groupbytrait <- unique(groupingtrait)
   	numrowgroups <- length(groupbytrait)
   	rowgroups <- vector("list", numrowgroups)
   	for (i in 1:numrowgroups){
   		rowgroups[[i]] <- which(groupingtrait==groupbytrait[i], arr.ind=TRUE)
   	}
}
groupbytrait_row <- groupbytrait
# Assign nodes in columns into groups -- man
if (mymangroup_BN=="manonegroup"){
	groupbytrait <- "0"
	numcolgroups <- 1
   	colgroups <- vector("list", numcolgroups)
   	colgroups[[1]] <- 1:dim(mybipartitematrix)[2]
} else if (mymangroup_BN=="geartype"){
	groupingtrait <- manTrait[, 4]
	groupbytrait <- unique(groupingtrait)
  	numcolgroups <- length(groupbytrait)
   	colgroups <- vector("list", numcolgroups)
   	for (i in 1:numcolgroups){
   		colgroups[[i]] <- which(groupingtrait==groupbytrait[i], arr.ind=TRUE)
   	}
} else if (mymangroup_BN=="landingsite"){
   	groupingtrait <- manTrait[, 5]
	groupbytrait <- unique(groupingtrait)
   	numcolgroups <- length(groupbytrait)
   	colgroups <- vector("list", numcolgroups)
   	for (i in 1:numcolgroups){
   		colgroups[[i]] <- which(groupingtrait==groupbytrait[i], arr.ind=TRUE)
   	}
} else if (mymangroup_BN=="geartype-landingsite"){
   	groupingtrait1 <- manTrait[, 4]
   	groupingtrait2 <- manTrait[, 5]
   	groupingtraitcombo <- paste(groupingtrait1, "+", groupingtrait2, sep="")
   	groupbytrait <- unique(groupingtraitcombo)
   	numcolgroups <- length(groupbytrait)
   	colgroups <- vector("list", numcolgroups)
   	for (i in 1:numcolgroups){
   		colgroups[[i]] <- which(groupingtraitcombo==groupbytrait[i], arr.ind=TRUE)
   	}
   	print ("geartype-landingsite groups:")
   	print (groupbytrait)
} else if (mymangroup_BN=="maneachgroup"){
   	groupingtrait <- manTrait[, 1]
	groupbytrait <- unique(groupingtrait)
   	numcolgroups <- length(groupbytrait)
   	colgroups <- vector("list", numcolgroups)
   	for (i in 1:numcolgroups){
   		colgroups[[i]] <- which(groupingtrait==groupbytrait[i], arr.ind=TRUE)
   	}
}
groupbytrait_col <- groupbytrait

socialnetworkbasis <- socialnetwork

resultsmatrix_BN <- matrix(rep(0, 9 * length(rowgroups) * length(colgroups)), nrow=length(rowgroups) * length(colgroups), ncol=9)
colnames(resultsmatrix_BN) <- c("TL_EmptyBN", "TL_BN", "TL_RandBN", "+/-", "p(TL_RandBN<TL_BN)", "D(BN-EmptyBN)", "D(BN-RandBN)", "fishgroup", "mangroup")
resultsmatrix_BN_partitions <- matrix(rep(0, 2 * length(rowgroups) * length(colgroups)), nrow=length(rowgroups) * length(colgroups), ncol=2)
colnames(resultsmatrix_BN_partitions) <- c("fishgroup", "mangroup")

myrow <- 1
for (i in 1:length(rowgroups)){
	for (j in 1:length(colgroups)){
		resultsmatrix_BN_partitions[myrow, 1] <- groupbytrait_row[i]
		resultsmatrix_BN_partitions[myrow, 2] <- groupbytrait_col[j]
		TLoutput <- gettotallengthforblockDAGcombination_plusrandomBNs(mybipartitematrix, rowgroups[[i]], colgroups[[j]], socialnetworkbasis, manTrait[, c(2,3)], numrandomisations)
		resultsmatrix_BN[myrow, 2] <- TLoutput$TotalLength
		resultsmatrix_BN[myrow, 3] <- mean(TLoutput$TotalLength_randomBN)
		resultsmatrix_BN[myrow, 4] <- sd(TLoutput$TotalLength_randomBN)
		resultsmatrix_BN[myrow, 1] <- gettotallengthforblockDAGcombination(mybipartitematrix, rowgroups[[i]], colgroups[[j]], socialnetwork_empty, manTrait[, c(2,3)])
		resultsmatrix_BN[myrow, 5] <- sum(TLoutput$TotalLength_randomBN<TLoutput$TotalLength)/numrandomisations
		if (myfishgroup_BN == "sizeclass-pricecateg"){
			resultsmatrix_BN[myrow, 8] <- groupbytrait_row[i]
			resultsmatrix_BN[myrow, 9] <- groupbytrait_col[j]
			resultsmatrix_BN[myrow, 6] <- as.numeric(resultsmatrix_BN[myrow, 2]) - as.numeric(resultsmatrix_BN[myrow, 1])
			resultsmatrix_BN[myrow, 7] <- as.numeric(resultsmatrix_BN[myrow, 2]) - as.numeric(resultsmatrix_BN[myrow, 3])
		} else if (mymangroup_BN == "geartype-landingsite"){
			resultsmatrix_BN[myrow, 8] <- groupbytrait_row[i]
			resultsmatrix_BN[myrow, 9] <- groupbytrait_col[j]
			resultsmatrix_BN[myrow, 6] <- as.numeric(resultsmatrix_BN[myrow, 2]) - as.numeric(resultsmatrix_BN[myrow, 1])
			resultsmatrix_BN[myrow, 7] <- as.numeric(resultsmatrix_BN[myrow, 2]) - as.numeric(resultsmatrix_BN[myrow, 3])
		} else {
			resultsmatrix_BN[myrow, 8] <- as.numeric(groupbytrait_row[i])
			resultsmatrix_BN[myrow, 9] <- as.numeric(groupbytrait_col[j])
			resultsmatrix_BN[myrow, 6] <- resultsmatrix_BN[myrow, 2] - resultsmatrix_BN[myrow, 1]
			resultsmatrix_BN[myrow, 7] <- resultsmatrix_BN[myrow, 2] - resultsmatrix_BN[myrow, 3]
		}
		myrow <- myrow + 1
	}
}

print("Fish partition:")
print(myfishgroup_BN)
print("Man partition:")
print(mymangroup_BN)
print("Results (also saved as a csv file in working directory):")
print(resultsmatrix_BN)

savefilename <- paste("resultsmatrix_BN_", myfishgroup_BN,"_", mymangroup_BN,".csv", sep="")
write.table(resultsmatrix_BN, file=savefilename, row.names=FALSE, col.names=TRUE, sep=",")

# Clean up files
system("rm mybipartitematrix.txt")

runTime <- Sys.time() - beginTime
print("Program completed successfully")
print("Run time:")
print(runTime)