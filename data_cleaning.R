library(tidyverse)
library(readr)

linkco2013us_den <- read_csv("Infant Health Datasets/linkco2013us_den.csv")

cols_to_keep = c("mager41",
                 "cig_rec",
                 "fagecomb",
                 "precare",
                 "wtgain",
                 "rf_ghyp",
                 "rf_gest",
                 "rf_eclam",
                 "rf_ppoutc",
                 "bwtr14",
                 "sex")

df_infant_2013 = linkco2013us_den %>% select(all_of(cols_to_keep))

new_df = matrix(0,nrow=0,ncol=ncol(df_infant_2013))

rows_to_keep = rbernoulli(nrow(df_infant_2013),p=0.01)
rows_to_keep = as.integer(rows_to_keep)

df_infant_2013 = df_infant_2013[rows_to_keep,]


rows_to_keep = c()

for (i in 1:nrow(df_infant_2013)){
  if(!any(is.na(df_infant_2013[i,]))){
    rows_to_keep = c(rows_to_keep,i)
  }
}

df_infant_2013 = df_infant_2013[rows_to_keep,]

df_infant_2013 = df_infant_2013 %>%
                filter(fagecomb != 99) %>%
                filter(cig_rec != "U") %>%
                filter(precare != 99) %>%
                filter(wtgain != 99) %>%
                filter(rf_ghyp != "U") %>%
                filter(rf_gest != "U") %>%
                filter(rf_eclam != "U") %>%
                filter(rf_ppoutc != "U") %>%
                filter(bwtr14 != 14)

write.csv(df_infant_2013, "reduced_infant_data.csv")                                
                                        
                    
