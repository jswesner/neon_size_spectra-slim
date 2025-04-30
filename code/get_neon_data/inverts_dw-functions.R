# functions for estimating macroinvertebrate dry weights 
# written by Justin Pomeranz
# jfpomeranz@gmail.com
# started June 2020
# modified from analysis in Pomeranz et al. 2022 Global change Biology

# function to add published length weight regression coefficients to macroinvertebrate taxonomic data from NEON sites. Only returns observations which has a corresponding LW coefficient. i.e. output may be smaller than input
LW_coef <- function(x, lw_coef, percent = FALSE){
  # x = inv_taxonomyProcessed table from NEON data product "DP1.20120.001"
  # lw_coef = table of taxon specific length_weight coefficients including higher taxonomic classification. Note that this table must have a column named "taxon" which is used to match up with the different taxonomic levels in x
  # percent = Report percent of rows in x which have a taxon-specific LW equation? default is FALSE
  
  
  # record number of rows in x to calculate percent of observations with a LW equation
  x.nrow <- nrow(x)
  # make indexes for each taxonomic level and extract those rows
  # remove rows from x that are accounted for in taxononmic tables
  
  # genus
  genus_index <- which(x$genus %in% lw_coef$taxon)
  # make table of observations which have a genus-specific LW equation
  genus_table <- x[genus_index,]
  # remove these from global taxonomic data
  x <- x[-genus_index,]
  
  # repeat process for other taxonomic levels
  
  # family
  family_index <- which(x$family %in% lw_coef$taxon)
  family_table <- x[family_index,]
  x <- x[-family_index,]
  
  # order
  order_index <- which(x$order %in% lw_coef$taxon)
  order_table <- x[order_index,]
  x <- x[-order_index,]
  
  # subclass
  subclass_index <- which(x$subclass %in% lw_coef$taxon)
  subclass_table <- x[subclass_index,]
  x <- x[-subclass_index,]
  
  # class
  class_index <- which(x$class %in% lw_coef$taxon)
  class_table <- x[class_index,]
  x <- x[-class_index,]
  
  # combine coefficients to taxonomic tables
  lw_cols <- c("taxon",
               "a",
               "b",
               "L_units",
               "dw_units",
               "formula_type",
               "formula")
  genus_table <- merge(genus_table, lw_coef[,lw_cols], 
                       by.x = "genus",
                       by.y = "taxon", all.x = TRUE)
  family_table <- merge(family_table, lw_coef[,lw_cols], 
                        by.x = "family",
                        by.y = "taxon", all.x = TRUE)
  order_table <- merge(order_table, lw_coef[,lw_cols], 
                       by.x = "order",
                       by.y = "taxon", all.x = TRUE)
  subclass_table <- merge(subclass_table, lw_coef[,lw_cols], 
                          by.x = "subclass",
                          by.y = "taxon", all.x = TRUE)
  class_table <- merge(class_table, lw_coef[,lw_cols], 
                       by.x = "class",
                       by.y = "taxon", all.x = TRUE)
  datout <- rbind(genus_table, family_table, order_table, subclass_table, class_table)
  
  if (percent == TRUE){
    message(paste(nrow(datout)/x.nrow * 100,
                  "percent of input data had Length-Weight equations"))
  }
  return(datout)
}

# estimate dry weight (dw) values from length-weight regression coefficients
est_dw <- function(x, fieldData){
  # x = inv_taxonomyProcessed table from NEON data product "DP1.20120.001" with LW coefficients added using the LW_coeff() function
  # fieldData = inv_fieldData table from NEON data product "DP1.20120.001"
  
  # functions uses tidyverse functions
  require(tidyverse)
  
  # simplify fieldData to three columns
  field = fieldData %>%
    select(siteID, sampleID, benthicArea) %>%
    distinct()
  # add benthicArea column from fieldData to x. This is necessary to calculate number per m2 below
  
  # join by siteID and sampleID
  x <- left_join(x, field, by = c("siteID", "sampleID"))
  
  # correct counts to per meter squared. Round to nearest integer
  x <- x %>%
    mutate(no_m2 = round(estimatedTotalCount / benthicArea))
  
  # calculate dw based on different formula types
  # sizeClass = length of individual in mm
  x <- x %>% 
    mutate(dw = case_when(formula_type == 1 ~ a * sizeClass^b,
                          formula_type == 2 ~ exp(a + b * log(sizeClass))))
  return(x)
}
