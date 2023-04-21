d_vn = function(N,alpha,v){
  #N=1;v=2;alpha=.5

  if(N==1){
    dvn = alpha*v^(alpha-1)
    qnj=1
    f_alpha_v = alpha^((N))*v^(((N))*alpha-N)
    dvn <- (-1)^N*sum(f_alpha_v*qnj)*exp(-v^alpha)
    return(list(dvn=dvn,qnj=qnj,f_alpha_v=f_alpha_v))
  }

  qnj <-  matrix(0, nrow = N, ncol = N+1)
  qnj[,1]  <- 0

  colnames(qnj) = paste0("j=",seq(0,N))
  rownames(qnj) = paste0("N=",seq(1,N))

  for(n in 1:N){
    qnj[n,n+1] <- 1
  }

  for(n in 2:N){

    for(j_pos in 2:(N)){

      #j = j_pos-1
      qnj[n,j_pos] <- qnj[n-1,j_pos-1] + (n-1-(j_pos-1)*alpha)*qnj[n-1,j_pos]
    }
  }

  f_alpha_v = alpha^((0:N))*v^(((0:N))*alpha-N)
  dvn <- (-1)^N*sum(f_alpha_v*qnj[N,])*exp(-v^alpha)
  return(list(dvn=dvn,qnj=qnj,f_alpha_v=f_alpha_v))
}
