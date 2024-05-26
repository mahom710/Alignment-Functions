CheckForInputErrors <- function(x,y,GapPenalty){
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
}

CreateSubstitutionMatrix <- function(gp,LevelOfSimilarity){
  if(LevelOfSimilarity == "low"){
    data("BLOSUM50"); subm <- BLOSUM50 
  } else if (LevelOfSimilarity == "moderate"){
    data("BLOSUM62"); subm <- BLOSUM62
  } else if (LevelOfSimilarity == "high"){
    data("BLOSUM80"); subm <- BLOSUM80
  }
  return(subm)
}

CreateScoringMatrix <- function(x,y,gp){
  ma <- matrix(NA, length(y)+1, length(x)+1, dimnames=list(c("gp", y), c("gp", x)))
  ma[1,] <- seq(0, -(length(ma[1,])-1) * gp, -gp)
  ma[,1] <- seq(0, -(length(ma[,1])-1) * gp, -gp)
  return(ma)
}

CreateTraceBackMatrix <- function(x,y){
  indexMa <- matrix(list(NA), length(y)+1, length(x)+1, dimnames=list(c("gp", y), c("gp", x)))
  return(indexMa)
}

SolveScoringMatrixAndTraceback <- function(x,y,subm,ma,indexMa,gp){
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
  
  # fix the first col and row
  indexMa[1,1] <- "END"
  for (i in 1:length(x)){
    indexMa[1,i+1] = "left"
  }
  for (i in 1:length(y)){
    indexMa[i+1] = "up"
  }
  
  return(list(scoringMatrix=ma, tracebackMatrix=indexMa))
}

GlobalAlignmentFromTraceBack <- function(x,y,subm,ma,indexMa,gp){
  # setup
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
  return(list(finalScore=finalscore,alignment1=xalign,alignment2=yalign))
  
}



























