matrix.double.center <- function(mat){
  mean_row_partial <- apply(mat, MARGIN = 2, 
                            function(x){
                              x - rowMeans(mat)
                            }
  )
  mean_col_partial <- t(apply(mean_row_partial, MARGIN = 1, 
                              function(x){
                                x - colMeans(mat)
                              }))
  double_center_matrix <- mean_col_partial + mean(as.matrix(mat))
  return(double_center_matrix)
}