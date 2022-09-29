### Function to select the proportion of cells of each population

# Function to select random cells keeping known cell populations proportion
select_populations <- function(data, 
                               proportions, 
                               cell_names_column = "clusters",
                               set.seed = 123) {
  
  set.seed(set.seed)
  # Show the different populations present in the data
  print(unique(data[, cell_names_column]))
  
  # Empty vector to save the barcodes
  cells_b1 <- c()
  
  # Counting variable
  count <- 1
  
  # Iterate the different populations
  for (type in unique(data[, cell_names_column])) {
    
    # Shuffle the rows
    data <- data[sample(1:nrow(data)),]
    
    # Choose the cells for the batch
    cells_b1 <- c(cells_b1, 
                  sample(rownames(data[data[, cell_names_column] == type, ]),
                         round(proportions[count] * length(rownames(data[data[, cell_names_column] == type, ]))),
                         prob = rep(1, nrow(data[data[, cell_names_column] == type, ]))
                  )
    )
    
    count <- count + 1
    
  }
  return(cells_b1)
}