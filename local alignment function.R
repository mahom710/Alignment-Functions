########################
## # Local Alignment ##
########################
library(Biostrings)


## Define input sequences
x <- "FIPFSAGPRNCIGQK"
y <- "PFGFGKRSCMGRRLA"

## Vectorize input sequences
x <- substring(x, 1:nchar(x), 1:nchar(x))
y <- substring(y, 1:nchar(y), 1:nchar(y))

## Define gap penality and sub matrix
gp <- 8 # Gap penalty
data("BLOSUM50"); subm <- BLOSUM50 


## Create dynamic programming matrix based on x, y and gp
ma <- matrix(NA, length(y)+1, length(x)+1, dimnames=list(c("gp", y), c("gp", x)))
ma[1,] <- 0
ma[,1] <- 0

#Create second identical matrix that will store the coordinates of where each data point points
indexMa <- matrix(list(NA), length(y)+1, length(x)+1, dimnames=list(c("gp", y), c("gp", x)))
#indexMa[1, 1] <- list(c(1, 2))  # this is an example of how to store a cord
# note that the first col and row should be empty so that the indices of the ma nad indexMa match

#Solve of each index
for (i in 2:(length(x)+1)){
  for (j in 2:(length(y)+1)) {
    
    #substitution penalty             
    aminoX <- colnames(ma)[j] 
    index_position_blosum_i <- grep(aminoX,row.names(subm))
    aminoY <- row.names(ma)[i]
    index_position_blosum_j <- grep(aminoY,colnames(subm))
    sub_pen <- subm[index_position_blosum_i, index_position_blosum_j]
    im1jm1 <- ma[i-1,j-1] + sub_pen
    
    #gap penalties
    im1j <- ma[i-1,j] - gp
    ijm1 <- ma[i,j-1] - gp
    
    #places the score in the matrix
    maxscore <- max(im1jm1,im1j,ijm1,0) 
    ma[i,j] <- maxscore
    
    #stores the cords of where the max score came from
    if(maxscore == im1jm1){
      indexMa[i,j] <- list(c(i-1,j-1))
    } else if (maxscore == im1j){
      indexMa[i,j] <- list(c(i-1,j))
    } else if (maxscore == ijm1){
      indexMa[i,j] <- list(c(i,j-1))
    } 
  }
  #I could improve this part by storing the direction in the index matrix instead of the coordinates.
  #For example, the im1jm1 in the above "if else" loop could just store the direction (diag) and we could
  #cut out the entire "checking to see where the current index points" part into just looking for the direction.
  
}

#fix the first row and column points for the index matrix
indexMa[1,1] <- list(c(0,0))
for (i in 1:length(x)){
  indexMa[1,i+1] = list(c(1,i))
}
for (i in 1:length(y)){
  indexMa[i+1] = list(c(i,1))
}

# Alignment

#find the largest value after converting matrix to linear vector
index_largest <- which.max(ma)

# Convert the index to matrix indices
index_matrix <- arrayInd(index_largest, dim(ma))

curindex <- c(index_matrix[1], index_matrix[2]) 
  
xalign <- c()
yalign <- c()

#tracing back through the matrix until it reaches the 1,1
while (ma[curindex[1],curindex[2]] > 0) {
  
  #checking to see where the current index points
  up <- unlist(indexMa[curindex[1], curindex[2]]) == c(curindex[1]-1,curindex[2]) #tells if the arrow points up
  left <- unlist(indexMa[curindex[1], curindex[2]]) == c(curindex[1],curindex[2]-1) #tells if the arrow points left
  diag <- unlist(indexMa[curindex[1], curindex[2]]) == c(curindex[1]-1,curindex[2]-1) #tells if the arrow points diag
  
  #create a string for the alignments and update the current index
  if(all(up)){
    yalign <- paste0(row.names(ma)[curindex[1]],yalign)
    xalign <- paste0("-", xalign)
    curindex <- c(curindex[1]-1,curindex[2])
  } else if(all(left)){
    yalign <- paste0("-", yalign)
    xalign <- paste0(colnames(ma)[curindex[2]], xalign)
    curindex <- c(curindex[1],curindex[2]-1)
  } else if(all(diag)){
    yalign <- paste0(row.names(ma)[curindex[1]], yalign)
    xalign <- paste0(colnames(ma)[curindex[2]], xalign)
    curindex <- c(curindex[1]-1,curindex[2]-1)
  } 
  
  
}

print(xalign)
print(yalign)


writeLines(c(xalign, yalign), "LocalAlignment.txt")
write.table(ma, file="LocalAlignment.xls", quote=FALSE, na = "", col.names = NA, sep="\t")

seq






















