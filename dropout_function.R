
# function for dorpout handling,
# M = vercot of visit time points
dropout <- function(rate = .1, data = dat, M) {
  # Get unique ids
  ids <- unique(dat$id)
  # Sample
  id_rm <- sample(ids,
                  size = floor(length(ids) * rate))
  dat <- subset(dat, !(id %in% id_rm))
  return(dat)
}
for (v in M) {
  ids <- unique(subset(dat, M >= v)$id)
  # Sample
  id_rm <- sample(ids,
                  size = floor(length(ids) * drop_out))
  dat <- subset(dat, !(M >= v & (id %in% id_rm)))
}
