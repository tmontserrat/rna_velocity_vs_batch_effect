# Function to apply is.nan() to a dataframe

# Ref.: https://stackoverflow.com/questions/18142117/how-to-replace-nan-value-with-zero-in-a-huge-data-frame

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))