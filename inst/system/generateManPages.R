library(R.utils)
library(RUVnormalize)

## Set default author
author <- "Laurent Jacob"
     
setwd("R")
sourceDirectory(".")
rdocFiles <- list.files(".", full.names=TRUE)

## Compile the Rdoc:s into Rd files (saved in the destPath directory)
destPath <- "../man"
destPath <- Arguments$getWritablePath(destPath)

print(rdocFiles)

dummy <- lapply(rdocFiles, Rdoc$compile, destPath=destPath, verbose=TRUE, debug=TRUE)
     
## List the generated Rd files
rdFiles <- list.files(destPath, full.names=TRUE)
print(rdFiles)
     
## Show one of the files
## file.show(rdFiles[1])
     
## Clean up
setwd("..")
