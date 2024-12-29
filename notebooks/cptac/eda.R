library(tidyverse)
library(readxl)

df <- read_excel('data/cptac/metadata.xlsx', sheet = 'clinical_data')
names(df)

# filter IDHwt GBM, Treatment-Naive, Primary GBM
df <- df |> 
    filter(path_diagnosis == "Glioblastoma", 
        idh1_r132h_mut_status != "Positive", 
        treatment == "NA", 
        str_detect(sample_type, "Primary"))


