# Function to process TPS file and compute scaled distances
process_data <- function(file_path) {
  lines <- readLines(file_path)
  lm_indices <- which(grepl("LM=", lines))
  scale_indices <- which(grepl("SCALE=", lines))
  id_indices <- which(grepl("ID=", lines))
  
  if (length(lm_indices) == 0 || length(scale_indices) == 0 || length(id_indices) == 0) {
    stop("Missing LM=, SCALE=, or ID= fields in TPS file.")
  }
  
  results <- data.frame(ID = character(), Mean_Distance = numeric(), stringsAsFactors = FALSE)
  
  for (i in seq_along(lm_indices)) {
    num_landmarks <- as.numeric(gsub("LM=", "", lines[lm_indices[i]]))
    scale <- as.numeric(gsub("SCALE=", "", lines[scale_indices[i]]))
    specimen_id <- gsub("ID=", "", lines[id_indices[i]])
    specimen_id <- gsub("\\s+", "", specimen_id)
    
    landmark_lines <- lines[(lm_indices[i] + 1):(lm_indices[i] + num_landmarks)]
    landmarks_coords <- do.call(rbind, strsplit(landmark_lines, "\\s+"))
    landmarks_coords <- as.data.frame(landmarks_coords, stringsAsFactors = FALSE)
    landmarks_coords <- dplyr::mutate_all(landmarks_coords, as.numeric)
    
    if (nrow(landmarks_coords) < 2 || is.na(scale)) {
      next
    }
    
    # Calculate all pairwise distances
    distances <- numeric()
    for (j in 1:(nrow(landmarks_coords) - 1)) {
      for (k in (j + 1):nrow(landmarks_coords)) {
        dist <- calculate_distance(landmarks_coords[j, 1], landmarks_coords[j, 2], 
                                   landmarks_coords[k, 1], landmarks_coords[k, 2])
        scaled_distance <- dist * scale
        distances <- c(distances, scaled_distance)
      }
    }
    
    mean_dist <- mean(distances)
    results <- rbind(results, data.frame(ID = specimen_id, Mean_Distance = mean_dist, stringsAsFactors = FALSE))
  }
  
  return(results)
}
