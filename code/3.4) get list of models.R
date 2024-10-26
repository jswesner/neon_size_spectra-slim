model_files = list.files(
  path = "models/", 
  pattern = "\\.rds$", 
  full.names = TRUE
)

# Keep only models with updated culling via the Clauset method
model_files = model_files[grepl("clauset", basename(model_files))]

model_list = NULL

# Loop through each file and read the RDS into the list
for (file in model_files) {
  model_list[[file]] <- readRDS(file)
}

saveRDS(model_list, file = "models/model_list.rds")

