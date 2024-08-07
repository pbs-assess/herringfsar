# Combine SOG data into a single CSV
# See more info on csasdown at:
# https://github.com/pbs-assess/csasdown
library(dplyr)

C.df <- read.csv("catchData.csv", header = TRUE ) |>
          group_by( Year ) |>
          summarise( Catch = sum(Value))

I.df <- read.csv("SplitIdx/SOG_BlendedIdx.csv", header = TRUE ) |>
          mutate(totSpawn = Dive + Surface ) |>
          dplyr::select( Year, totSpawn ) 

final.df <- I.df |>
            left_join(C.df)


write.csv(final.df, file = "SOGIdxCatData.csv")

