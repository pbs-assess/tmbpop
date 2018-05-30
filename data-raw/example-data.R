Ct               <- read.table("inst/extdata/Ct.txt") # catc = Ct,
ma               <- as.numeric(unlist(read.table("inst/extdata/ma.txt"))) # prop-at-age mature
wa1              <- as.numeric(unlist(read.table("inst/extdata/wa1.txt"))) # weight females
wa2              <- as.numeric(unlist(read.table("inst/extdata/wa2.txt"))) # weight males
S1t              <- read.table("inst/extdata/S1t.txt") # survey 1: West Coast Vancouver Island
S2t              <- read.table("inst/extdata/S2t.txt") # survey 2: National Marine Fisheries Service
S3t              <- read.table("inst/extdata/S3t.txt") # survey 3: GB Reed
patC1            <- read.table("inst/extdata/patC1.txt") # prop-at-age catch females
patC2            <- read.table("inst/extdata/patC2.txt") # prop-at-age catch males
ntC              <- read.table("inst/extdata/ntC.txt") # number trips prop-at-age catch
patS11           <- read.table("inst/extdata/patS11.txt") # prop-at-age survey 1 females
patS12           <- read.table("inst/extdata/patS12.txt") # prop-at-age survey 1 males
ntS1             <- read.table("inst/extdata/ntS1.txt") # number trips prop-at-age survey 1
MiscFixedParam   <- read.table("inst/extdata/MiscFixedParam.txt", # additional fixed param
  header = TRUE
) # selectivity survey 2 and 3, sigmaR

pop_example <- list(
  Ct = Ct,
  ma = ma,
  wa1 = wa1,
  wa2 = wa2,
  S1t = S1t,
  S2t = S2t,
  S3t = S3t,
  patC1 = patC1,
  patC2 = patC2,
  ntC = ntC,
  patS11 = patS11,
  patS12 = patS12,
  ntS1 = ntS1,
  MiscFixedParam = MiscFixedParam
)

usethis::use_data(pop_example, overwrite = TRUE)
