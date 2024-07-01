compute_G_game<-function(G_marginal, A, S){
  m_a = max(A$k)
  m_s = max(S$l)
  L   = max(S$s)
  n_firm <- max(A$i)
  output <- matrix(0,nrow = m_a*m_s, ncol = m_s)
  for (k in 1:m_a){
    for (l in 1:m_s){
      A_k <- A$a[A$k == k]
      S_l <- S$s[S$l == l]
      sa_1 <- c(A_k[1],S_l[1])
      sa_2 <- c(A_k[2],S_l[2])
      sa_3 <- c(A_k[3],S_l[3])
      transition_1 <- G_marginal[(sa_1[2]-1)*2+sa_1[1]+1,]
      transition_2 <- G_marginal[(sa_2[2]-1)*2+sa_2[1]+1,]
      transition_3 <- G_marginal[(sa_3[2]-1)*2+sa_3[1]+1,]
      transition <- expand.grid(transition_1,transition_2,transition_3)
      transition$t <- transition$Var1*transition$Var2*transition$Var3 
      output[(l-1)*m_a+k,] <-  transition$t
    }
  }
  return(output)
}
