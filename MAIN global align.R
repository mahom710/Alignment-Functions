GlobalAlign <- function(FirstSequence, SecondSequence, GapPenalty, LevelOfSimilarity, returnMatrices){
  # Load in functions
  source("./Functions.R")
  
  # Start Timer
  start_time <- Sys.time()
  
  # Define input sequences (by axis in matrix)
  x <- FirstSequence; y <- SecondSequence; 
  
  # Check for input errors
  CheckForInputErrors(FirstSequence,SecondSequence,GapPenalty)
  
  # Vectorize input sequences
  x <- unlist(strsplit(x, split = ""))
  y <- unlist(strsplit(y, split = ""))
  
  # Create the substitution matrix (BLOSUM)
  substitutionMatrix <- CreateSubstitutionMatrix(GapPenalty,LevelOfSimilarity)
  
  # Create scoring matrix
  scoringMatrix <- CreateScoringMatrix(x,y,GapPenalty)
  
  # Create traceback matrix
  tracebackMatrix <- CreateTraceBackMatrix(x,y)
  
  # Populate both the scoring and traceback matrices
  listMatrices <- SolveScoringMatrixAndTraceback(x,y,substitutionMatrix,scoringMatrix,tracebackMatrix,GapPenalty)
  scoringMatrix <- listMatrices$scoringMatrix; tracebackMatrix <- listMatrices$tracebackMatrix
  
  # Conduct the traceback to find the alignments
  listFinal <- GlobalAlignmentFromTraceBack(x,y,substitutionMatrix,scoringMatrix,tracebackMatrix,GapPenalty)
  finalScore <- listFinal$finalScore; alignment1 <- listFinal$alignment1; alignment2 <- listFinal$alignment2
  
  # End Timer
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  print(time_taken)
  
  # Print Results
  print(paste0("Seq1:",alignment1))
  print(paste0("Seq2:",alignment2))
  print(paste0("Final Score:", finalScore))
  
  if (returnMatrices == TRUE){
    return(list(scoringMatrix=scoringMatrix,tracebackMatrix=tracebackMatrix))
  }
}
  


