########################
## # Global Alignment ##
########################

GlobalAlign <- function(FirstSequence, SecondSequence, GapPenalty, LevelOfSimilarity, returnMatrices){
  
  #can add a feature of input for sequences (like fastq, string)
library(Biostrings)

# add some timer to benchmarking
start_time <- Sys.time()

# Possible Errors
if (GapPenalty < 0) {
  stop("Error: Gap penalty cannot be negative.")
}
only_allowed_characters <- function(x, allowed_characters) {
  return(!grepl(paste0("[^ARNDCQEGHILKMFPSTWYVBX*]"), x))
}
if (!only_allowed_characters(x)){
  stop("Error: Ensure first sequence is uppercase and proteins.")
}
if (!only_allowed_characters(y)){
  stop("Error: Ensure second sequence is uppercase and proteins.")
}
if (x == ""){
  stop("Error: First sequence is empty.")
}
if (y == ""){
  stop("Error: Second sequence is empty.")
}

## Define input sequences (by axis in matrix)
x <- FirstSequence; y <- SecondSequence #x <- "HEAGAWGHEE"; y <- "PAWHEAE"

# Vectorize input sequences
x <- unlist(strsplit(x, split = ""))
y <- unlist(strsplit(y, split = ""))

## Define gap penality and sub matrix
gp <- GapPenalty 

if(LevelOfSimilarity == "low"){
  data("BLOSUM50"); subm <- BLOSUM50 
} else if (LevelOfSimilarity == "moderate"){
  data("BLOSUM62"); subm <- BLOSUM62
} else if (LevelOfSimilarity == "high"){
  data("BLOSUM80"); subm <- BLOSUM80
}

## Create dynamic programming matrix based on x, y and gp
ma <- matrix(NA, length(y)+1, length(x)+1, dimnames=list(c("gp", y), c("gp", x)))
ma[1,] <- seq(0, -(length(ma[1,])-1) * gp, -gp)
ma[,1] <- seq(0, -(length(ma[,1])-1) * gp, -gp)

#Create second identical matrix that will store the coordinates of where each data point points. (note) that the first col and row should be empty so that the indices of the ma nad indexMa match
indexMa <- matrix(list(NA), length(y)+1, length(x)+1, dimnames=list(c("gp", y), c("gp", x)))

#Solve of each index 
for (i in 2:(length(y)+1)){
  for (j in 2:(length(x)+1)) {
    
    #substitution penalty             
    sub_pen <- as.numeric(subm[y[i-1],x[j-1]])
    im1jm1 <- ma[i-1,j-1] + sub_pen
    
    #gap penalties
    im1j <- ma[i-1,j] - gp
    ijm1 <- ma[i,j-1] - gp
    
    #places the score in the matrix
    maxscore <- max(im1jm1,im1j,ijm1) 
    ma[i,j] <- maxscore
    
    #stores the cords of where the max score came from
    if(maxscore == im1jm1){
      indexMa[i,j] <- "diag"
    } else if (maxscore == im1j){
      indexMa[i,j] <- "up"
    } else {
      indexMa[i,j] <- "left"
    }
  }
}

#fix the first row and column points for the index matrix (fix by pointers)
indexMa[1,1] <- "END"
for (i in 1:length(x)){
  indexMa[1,i+1] = "left"
}
for (i in 1:length(y)){
  indexMa[i+1] = "up"
}

# Alignment
curindex <- c(length(y)+1, length(x)+1)
finalscore <- ma[curindex[1],curindex[2]]
xalign <- c()
yalign <- c()

#tracing back through the matrix until it reaches the 1,1
while (indexMa[curindex[1],curindex[2]] != "END") {
  
  pointer <- unlist(indexMa[curindex[1],curindex[2]])
  #create a string for the alignments and update the current index
  if(pointer == "up"){
    yalign <- paste0(row.names(ma)[curindex[1]],yalign)
    xalign <- paste0("-", xalign)
    curindex <- c(curindex[1]-1,curindex[2])
  } else if(pointer == "left"){
    yalign <- paste0("-", yalign)
    xalign <- paste0(colnames(ma)[curindex[2]], xalign)
    curindex <- c(curindex[1],curindex[2]-1)
  } else if(pointer == "diag"){
    yalign <- paste0(row.names(ma)[curindex[1]], yalign)
    xalign <- paste0(colnames(ma)[curindex[2]], xalign)
    curindex <- c(curindex[1]-1,curindex[2]-1)
  } 
  
}

# end benchmark
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)

print(paste0("Seq1:",xalign))
print(paste0("Seq2:",yalign))
print(paste0("Final Score:", finalscore))

if (returnMatrices == TRUE){
  return(list(ma,indexMa))
}

}





















