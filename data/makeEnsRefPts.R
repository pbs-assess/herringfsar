# Render tech report for IPHC survey analysis
# See more info on csasdown at:
# https://github.com/pbs-assess/csasdown
library(dplyr)
# scenRefPts <- read.csv("SOG_allScenRefPtsMean.csv")
# scenRefPts$weights <- c(0.34,rep(0.165,4))

# for(k in 2:11)
#   scenRefPts[,k] <- scenRefPts[,k] * scenRefPts$weights

# ensRefPts <- scenRefPts |> summarise_if(is.numeric,sum)

# write.csv(ensRefPts, file = "SOG_ensRefPts.csv")


postParTable <- read.csv("omGrid_postPars.csv")
ensParTable <- read.csv("ensOM_meanPars.csv") |>
                mutate_if(is.numeric,round,2) |>
                mutate( PBTGtLRP = 0.99 )

parTable <- rbind(postParTable,ensParTable)
parTable$X <- NULL

scenRefPtsTable <- read.csv("SOG_allScenRefPts.csv")
ensRefPtsTable <- read.csv("SOG_ensRefPts.csv") |>
                  mutate_if(is.numeric, round, 2)

refPtsTable <- rbind(scenRefPtsTable,ensRefPtsTable)
refPtsTable$X <- NULL



allEstTable <- parTable |> 
                left_join(refPtsTable, by = "Scenario") |>
                dplyr::select(  h,
                                Mb,
                                M0,
                                m1,
                                Mbar,
                                qs,
                                qd,
                                qb,
                                R0, 
                                B0,                                 
                                BT,
                                DT = DT.x,
                                PBTGtLRP,
                                Bmsy,
                                Umsy,
                                MSY,
                                Uusr,
                                USR,
                                Yusr,
                                Ucrash) |>
                t()

colnames(allEstTable) <- c("OM 1","OM 2","OM 3","OM 4","OM 5","Ensemble")


rownames(allEstTable) <- c( "$h$",
                            "$M_b$",
                            "$M_0$",
                            "$m_1$",
                            "$\\overline{M}$",
                            "$q_s$",
                            "$q_d$",
                            "$q_{blend}$",
                            "$R_0$",
                            "$B_0$",
                            "$B_{2023}$",
                            "$B_{2023}/B_0$",
                            "$P(B_{2023} > 0.3 B_0)$",
                            "$B_{MSY}$",
                            "$U_{MSY}$",
                            "$MSY$",
                            "$U_{USR}$",
                            "$USR$",
                            "$Y_{USR}$",
                            "$U_{Crash}$")

