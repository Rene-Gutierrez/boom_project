### Prior Correlation Matrix

Lambda <- function(M, g, a){
  return(t((1 - g) * a + g * t(M)))
}