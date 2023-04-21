d_vn.mpfr = function(N, alpha, v) {
  #ok
  #N=180;v=2;alpha=.5
  # if(N<NN){
  #   stop("use d_vn ")
  # }
  NN = 168
  zero <- Rmpfr::mpfr(rep(0,N*(N+1)),2)
  N1 <- as.integer(N)
  N2 <- as.integer(N+1)
  qnj <-  matrix(0, nrow = N, ncol = N+1)
  qnj.mpfr <-  new("mpfrMatrix",zero,Dim=c(N1,N2))


  # colnames(qnj) = paste0("j=", seq(0, 120))
  # rownames(qnj) = paste0("N=", seq(1, 120))
  colnames(qnj) = paste0("j=", seq(0, N))
  rownames(qnj) = paste0("N=", seq(1, N))
  colnames(qnj.mpfr) = paste0("j=", seq(0, N))
  rownames(qnj.mpfr) = paste0("N=", seq(1, N))

  for (n in 1:N) {
    qnj[n, n + 1] <- 1
  }

  for (n in 1:N) {
    qnj.mpfr[n, n + 1] <- 1
  }

  for (n in 2:N) {


    for (j_pos in 2:n) {

      # cat("j=",j_pos,"\n")
      cat("n=",n,"\r")
      #j = j_pos-1
      if(n<=NN){
        qnj[n, j_pos] <-
          qnj[n - 1, j_pos - 1] + (n - 1 - (j_pos - 1) * alpha) * qnj[n - 1, j_pos]
      }

      if(n>NN){
        qnj.mpfr[NN,] <- qnj[NN,]
        qnj.mpfr[n, j_pos] <-
          qnj.mpfr[n - 1, j_pos - 1] + (n - 1 - (j_pos - 1) * alpha) * qnj.mpfr[n - 1, j_pos]
      }
    }
  }

  qnj.complet <- Rmpfr::rbind(qnj[1:NN,],qnj.mpfr[-c(1:NN),])

  f_alpha_v_N_1 = alpha ^ ((0:(N-1))) * v ^ (((0:(N-1))) * alpha - (N-1))
  dvn_1 <- (-1) ^ (N-1) * sum(f_alpha_v_N_1 * qnj.complet[N-1,1:N ]) * exp(-v ^ alpha)

  f_alpha_v_N = alpha ^ ((0:N)) * v ^ (((0:N)) * alpha - N)
  dvn <- (-1) ^ N * sum(f_alpha_v_N * qnj.complet[N, ]) * exp(-v ^ alpha)
  return(list(
    dvn = dvn,
    dvn_1 = dvn_1
  ))
}
