############################################################
# models.R
#
# Script to read in SISCAH operating model history and
# MP estimates.
#
# Last Modified: Aug 22, 2024
#
############################################################

# first, load model histories
histFolder <- file.path("data","SOG_DDM_omGrid")

fitFolders <- c("fit_parBatSOG_MbhGrid_h.71",
                "fit_parBatSOG_OMgrid_h.70_2",
                "fit_parBatSOG_OMgrid_h.70_3",
                "fit_parBatSOG_OMgrid_h.70_4",
                "fit_parBatSOG_OMgrid_h.70_5")

fitPaths <- here(histFolder, fitFolders, paste0(fitFolders, ".rds"))

histRpts <- lapply(X = fitPaths, FUN = readRDS)

scenarios <- c("SOG_Mb0.532_h0.70",
               "SOG_Mb0.562_h0.65",
               "SOG_Mb0.562_h0.70",
               "SOG_Mb0.562_h0.75",
               "SOG_Mb0.584_h0.70")


wtPosts <- readRDS("data/wtPosts.rds")

names(histRpts) <- scenarios

# order hist reports
histRpts <- histRpts[c(3,2,4,1,5)]

# Now weighted ensemble OMs
simFolder <- file.path("data","SOG_wtdPerf_Ugrid")

# Now read the infoFile in each sim folder
dirList   <- list.dirs(simFolder, recursive = FALSE)
dirList   <- dirList[grepl(x = dirList, pattern = "sim_")]

infoList  <- file.path(dirList,"infoFile.txt")
infoList  <- lapply(X = infoList, FUN = lisread)
infoList  <- lapply(X = infoList, FUN = as.data.frame)
info.df   <- do.call(rbind, infoList)

mps   <- c("Ugrid_0.06","Ugrid_0.14")
mps2  <- c("maxTHR_0.06","maxTHR_0.14")

# Filter info.df for chosen MPs
info.df <- info.df |> filter( mp %in% mps )
mpBlobList <- list()

for(k in 1:nrow(info.df))
{
  simLabel <- paste0("sim_",info.df$simLabel[k])
  path <- file.path(simFolder,simLabel,paste0(simLabel,".RData"))
  # Load blob
  load(path)

  # Save to mpBlobList
  mpBlobList[[k]] <- blob

}
names(mpBlobList) <- mps2

# Pull some model dimensions and year labels
fYear   <- blob$ctlList$opMod$fYear
nS      <- blob$om$nS
nP      <- blob$om$nP
nT      <- blob$om$nT
tMP     <- blob$om$tMP
pT      <- blob$ctlList$opMod$pT

# Labels for automated text generation
species <- blob$om$speciesNames
stock   <- blob$om$stockNames

# Get good reps (should be all there for ensemble)
goodRepIdx <- which(blob$goodRep)

# Draw a random replicate for plotting
set.seed(101)
randReplicate <- sample(goodRepIdx, size = 1)
randTraces    <- sample(goodRepIdx, size = 3)
usePosts      <- blob$ctlList$opMod$posteriorSamples

# Load MP application fit report
fit_maxTHR0.14 <- readRDS("data/fit_maxTHR0.14/fit_maxTHR0.14.rds")

# Load ref pts and ensemble parameters
ensRefPtsTable <- read.csv("data/SOG_ensRefPts.csv") |>
  mutate_if(is.numeric, round, 2)
ensParTable <- read.csv("data/ensOM_meanPars.csv") |>
  mutate_if(is.numeric,round,2) |>
  mutate(PBTGtLRP = ifelse(PBTGtLRP > 0.99, 0.99, .)) # if >.99 set = 0.99 else leave alone
