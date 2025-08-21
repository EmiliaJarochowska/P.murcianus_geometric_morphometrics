# Function to process TPS file and compute element lengths

process_data <- function(file_path) {
  lines <- readLines(file_path)
  lm_indices <- which(grepl("LM=", lines))
  scale_indices <- which(grepl("SCALE=", lines))
  id_indices <- which(grepl("ID=", lines))
  
  if (length(lm_indices) == 0 || length(scale_indices) == 0 || length(id_indices) == 0) {
    stop("Missing LM=, SCALE=, or ID= fields in TPS file.")
  }
  
  results <- data.frame(ID = character(), Length = numeric(), stringsAsFactors = FALSE)
  
  for (i in seq_along(lm_indices)) {
    num_landmarks <- as.numeric(gsub("LM=", "", lines[lm_indices[i]]))
    scale <- as.numeric(gsub("SCALE=", "", lines[scale_indices[i]]))
    specimen_id <- gsub("ID=", "", lines[id_indices[i]])
    
    landmark_lines <- lines[(lm_indices[i] + 1):(lm_indices[i] + 2)]
    landmarks_coords <- do.call(rbind, strsplit(landmark_lines, "\\s+"))
    landmarks_coords <- as.data.frame(landmarks_coords, stringsAsFactors = FALSE)
    landmarks_coords <- dplyr::mutate_all(landmarks_coords, as.numeric)
    
    if (nrow(landmarks_coords) < 2 || is.na(scale)) {
      next
    }
    
    # Calculate distance between landmarks 1 and 2 only
    if (nrow(landmarks_coords) >= 2) {
      dist <- calculate_distance(landmarks_coords[1, 1], landmarks_coords[1, 2], 
                                 landmarks_coords[2, 1], landmarks_coords[2, 2])
      scaled_distance <- dist * scale
      results <- rbind(results, data.frame(ID = specimen_id, Length = scaled_distance, stringsAsFactors = FALSE))
    }
  }
  
  results$ID <- results$ID %>% 
    trimws() %>% 
    toupper()
  
  return(results)
}