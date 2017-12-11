findpiv <- function(A, k, p, tol) {
  shape <- dim(A)
  m <- shape[1]
  n <- shape[2]
  r <- which(abs(A) > tol)
  if (length(r) == 0)
    break
  
  #
  r <- r[1]
  j <- as.integer(as.double(r-1)/m )+ 1
  p <- p + j - 1
  
  k <- k + r - (j - 1) * m - 1
  return(list(k, p))
}
m <- findpiv(matrix(data = c(1,0,0,0,1,0,0,0,1), 3, 3), 1,1, .Machine$double.eps)

m
m
switch2rows <- function(A, m, n, i, j) {
  B <- A
  B[m, i:j] <- A[n, i:j]
  B[n, i:j] <- A[m, i:j]
  return(B)
}

plu <- function(A) {
  shape <- dim(A)
  
  m <- shape[1]
  n <- shape[2]
  P <- as.matrix(diag(m))
  L <- as.matrix(diag(m))
  U <- matrix(0, m, n)
  pivcol <- c()
  tol <- sqrt(.Machine$double.eps)
  
  sign <- 1
  
  p <- 1
  for (k in 1:min(m, n)) {
    xy = findpiv(A[k:m, p:n], k, p, tol)

                                   
    r <- xy[[1]]
    p <- xy[[2]]
    if(length(r)==0)return(list(P=P,L=L,U=U))
    if (r!= k) {
      A <- switch2rows(A, k, r, 1, n)
      print(dim(A))
      if (k > 1) 
        L <- switch2rows(L, k, r, 1, k-1)

      P <- switch2rows(P, k, r, 1, m)
      sign <- -sign
    }
    if (abs(A[k, p]) >= tol) {
      pivcol[length(pivcol) + 1] <- p
      for (i in (k+1):m) {
        L[i, k] <- A[i, p] / A[k, p]
        for (j in (k +1):n) {
          A[i, j] <- A[i, j] - L[i, k]*A[k, j]
        }
      }
    }
    for (j in k:n) {
      U[k, j] <- A[k, j] * (abs(A[k, j]) >= tol)
    }
    if (p < n)
      p <- p + 1
  
  }
  
  
  # 
  return(list(P,L,U))
}

A <- matrix(data = c(1,0,2,3,1,0,0,0,1), 3, 3)

n <- plu(A)
n

