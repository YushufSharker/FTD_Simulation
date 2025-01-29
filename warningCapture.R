# Capture warning during simulation

library(stringr)

# Function to extract a word by matching it
get_warning <- function(text, word_match) {
  words <- str_split(text, "\\s+")[[1]] # Split the text into words
  selected_word <- words[str_detect(words, fixed(word_match))]

  if (length(selected_word) > 0) {
    return(selected_word)
  } else {
    return(NA) # Return NA if the word is not found
  }
}

get_phrase <- function(message, phrase) {
  if (any(grepl(phrase, message))) {
    return(phrase)
  } else {
    return(NA)
  }
}


