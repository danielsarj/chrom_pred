library(rhdf5)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/deeplearn_pred/DIY/contribution_scores')
conditions <- c('NI', 'Flu')
datasets <- c('counts', 'profile')

for (dset in datasets){
for (cond in conditions){
  
  # list all .h5 files in the folder
  h5_files <- list.files(getwd(), 
    pattern='^'%&%cond%&%'.*\\.'%&%dset%&%'_scores\\.h5$', 
    full.names=T)

  # initialize a list to store datasets
  datasets <- list()

  # loop through the .h5 files
  for (file in h5_files) {
    cat('Processing file:', file, '\n')
  
    # list keys (datasets) in the current file
    keys <- h5ls(file)$name
  
    # read each dataset and store it
    for (key in keys) {
      data <- h5read(file, key)
    
      # ensure the dataset is numeric and add to list
      if (is.numeric(data)) {
        if (!is.null(datasets[[key]])) {
          datasets[[key]] <- datasets[[key]] + data  # add data
        } else {
          datasets[[key]] <- data  # initialize dataset
        }
      }
    }
  }

  # average the accumulated datasets
  file_count <- length(h5_files)
  for (key in names(datasets)) {
    datasets[[key]] <- datasets[[key]] / file_count
  }

  # save new file
  for (key in names(datasets)) {
  h5write(datasets[[key]], cond%&%'_avg_'%&%dset%&%'_scores.h5', key)
  }

  cat('Average '%&%dset%&%' for '%&%cond%&%' saved.')
}
}