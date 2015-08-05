wig2jahmm <- function (directory) {
   fnames <- list.files(directory, pattern="\\.wig(\\.gz)?$")
   # Look for "input" or "control" in file names to build the
   # control profile.
   inputs <- grep("(input|control)", fnames,
                 ignore.case=TRUE, value=TRUE)
   targets <- grep("(input|control)", fnames,
                   ignore.case=TRUE, value=TRUE, invert=TRUE)

   # Check that all files are present.
   if (length(inputs) < 1) stop("You need at least one Input profile.")
   else write(c("Input file(s):", inputs), stderr())
   if (length(targets) < 1) stop("You need at least one ChIP profile.")
   else write(c("Target file(s):", targets), stderr())

   write("Reading files...", stderr())
   # Read chromosome information from the first file.
   dfChrom <- wig2df(fnames[1], chrom_only=TRUE)
   # Read input control files and sum them into a single vector.
   dfInputs <- rowSums(sapply(X=inputs, FUN=wig2df))
   # Read ChIP profiles and make a 'data.frame' out of them.
   dfTargets <- data.frame(lapply(X=targets, FUN=wig2df))
   names(dfTargets) <- seq(length(dfTargets))
   dfAll <- data.frame(chrom=dfChrom, input=dfInputs, dfTargets)

   write("Running jahmm...", stderr())
   return(jahmm(dfAll))
}

wig2df <- function(fname, chrom_only=FALSE) {
   input <- readLines(fname)
   input <- grep("^(#|browser|track)", input, value=TRUE, invert=TRUE)
   fileLength <- length(input)

   decLinesInd <- grep("^fixedStep", input)
   wholeDecLines <- input[decLinesInd]
   decLines <- matrix(unlist(strsplit(wholeDecLines, " ")),
                      nrow=length(wholeDecLines), byrow=TRUE)

   if (chrom_only) {
      # Keep the chromosome information.
      chrom <- matrix(unlist(strsplit(decLines[, 2], "=")),
                      ncol=2, byrow=TRUE)[, 2]
      chromRep <- c(diff(decLinesInd) - 1, fileLength - max(decLinesInd))
      chromCol <- rep(chrom, chromRep)
      return(chromCol)
   } else {
      # Only keep the scores.
      scores <- as.numeric(grep("^fixedStep", input,
                              value=TRUE, invert=TRUE))
      return(scores)
   }
}
