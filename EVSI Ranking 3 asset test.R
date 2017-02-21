calc_jp_num <- function(n_sims, mu_a, mu_b, mu_c, tau_a, tau_b, tau_c) {

  a <- rnorm(n_sims, mu_a, tau_a^-.5)
  b <- rnorm(n_sims, mu_b, tau_b^-.5)
  c <- rnorm(n_sims, mu_c, tau_c^-.5)

  abc <- mean(a > b & a > c & b > c)
  acb <- mean(a > c & a > b & c > b) 
  bac <- mean(b > a & b > c & a > c)  
  bca <- mean(b > c & b > a & c > a) 
  cab <- mean(c > a & c > b & a > b) 
  cba <- mean(c > b & c > a & b > a) 
  return(c(abc, acb, bac, bca, cab, cba))
  
}

calc_jp <- function(mu_a, mu_b, mu_c, tau_a, tau_b, tau_c) {
  require(mvtnorm)
  pr_ord <- function(ord) {
  
  mu_ <- unname(unlist(mget(paste('mu', ord, sep='_'), inherits=TRUE)))
  var_ <- 1 / unname(unlist(mget(paste('tau', ord, sep='_'), inherits=TRUE)))
  
  muz <- c(mu_[2] - mu_[1], mu_[3] - mu_[2])
  varz <- c(var_[1] + var_[2], var_[2] + var_[3])
    
  pmvnorm(c(-Inf, -Inf), c(0, 0), 
    mean=muz, 
    sigma=rbind(c(varz[1], -var_[2]), c(-var_[2], varz[2])))
  }
  
  return(
    c(pr_ord(c('a', 'b', 'c')), 
      pr_ord(c('a', 'c', 'b')), 
      pr_ord(c('b', 'a', 'c')), 
      pr_ord(c('b', 'c', 'a')), 
      pr_ord(c('c', 'a', 'b')), 
      pr_ord(c('c', 'b', 'a'))
    )
  )
}

upd_tau <- function(sigma2, n, tau) {
  tau + (n / sigma2) 
}

EVPI <- function(mu_a, mu_b, mu_c, tau_a, tau_b, tau_c) {
  
  jps <- calc_jp(mu_a, mu_b, mu_c, tau_a, tau_b, tau_c)

  u <- c(mu_a, mu_b, mu_c)
  
  E_pri_uA <- sum(jps * c(1, 1, .5, 0, .5, 0))
  E_pri_uB <- sum(jps * c(.5, 0, 1, 1, 0, .5))
  E_pri_uC <- sum(jps * c(0, .5, 0, .5, 1, 1))

  1 - pmax(E_pri_uA, E_pri_uB, E_pri_uC)
  
  
}

EVSI <- function(mu_a, mu_b, mu_c, tau_a, tau_b, tau_c, sigma2, N, p) {
  sapply(N, 
    function(N) {
      jps <- calc_jp(mu_a, mu_b, mu_c, tau_a, tau_b, tau_c)
    
      #E_pri_uA <- sum(jps * c(1, 1, .5, 0, .5, 0))
      E_pri_uA <- sum(jps * mu_a)
      #E_pri_uB <- sum(jps * c(.5, 0, 1, 1, 0, .5))
      E_pri_uB <- sum(jps * mu_b)
      #E_pri_uC <- sum(jps * c(0, .5, 0, .5, 1, 1))
      E_pri_uC <- sum(jps * mu_c)
      
      tau_a <- upd_tau(sigma2, N * p[1], tau_a)
      tau_b <- upd_tau(sigma2, N * p[2], tau_b)
      tau_c <- upd_tau(sigma2, N * (1 - p[1] - p[2]), tau_c)
      
      jps <- calc_jp(mu_a, mu_b, mu_c, tau_a, tau_b, tau_c)
    
      #E_post_uA <- sum(jps * c(1, 1, .5, 0, .5, 0))
      E_post_uA <- sum(jps * mu_a)
      #E_post_uB <- sum(jps * c(.5, 0, 1, 1, 0, .5))
      E_post_uB <- sum(jps * mu_b)
      #E_post_uC <- sum(jps * c(0, .5, 0, .5, 1, 1))
      E_post_uC <- sum(jps * mu_c)
      
      pmax(E_post_uA, E_post_uB, E_post_uC) - 
        pmax(E_pri_uA, E_pri_uB, E_pri_uC)
    }
  )  
}

EVSI_opt <- function(mu_a, mu_b, mu_c, tau_a, tau_b, tau_c, sigma2, N) {
   require(parallel)
  simplify2array(mclapply(N, 
     function(N) {   
       sol <- constrOptim(c(1/3, 1/3),
         function(p) {
           EVSI(mu_a, mu_b, mu_c, tau_a, tau_b, tau_c, sigma2, N, p)
        }, NULL, rbind(c(1, 0), c(0, 1), c(-1, -1)), c(0, 0, -1),
        control=list(fnscale=-1)        
       )
       return(c(sol$value, sol$par[1], sol$par[2], 1 - sol$par[1] - sol$par[2]))
     }, mc.cores=detectCores()
   ))
}

# EVSI_opt <- function(mu_a, mu_b, mu_c, tau_a, tau_b, tau_c, sigma2, N) {
#   require(parallel)
#   
#   n.cores <- detectCores()
#   N_ <<- as.list(as.data.frame(matrix(N[seq_len(length(N) %/% n.cores * 
#     n.cores)], nrow=length(N) %/% n.cores)))
#   if (length(N) != length(unlist(N_))) {
#     N_[[n.cores]] <- c(N_[[n.cores]], setdiff(N, unlist(N_)))
#   }  
#   N <- N_
#   
#   do.call(rbind, mclapply(N, 
#     function(N) {   
#       sol <- matrix(NA, length(N), 4)
#       sol_ <- constrOptim(c(1/3, 1/3),
#         function(p) {
#           EVSI(mu_a, mu_b, mu_c, tau_a, tau_b, tau_c, sigma2, N[1], p)
#         }, NULL, rbind(c(1, 0), c(0, 1), c(-1, -1)), c(0, 0, -1),
#         control=list(fnscale=-1)        
#       )
#       
#       sol[1, ] <- c(sol_$value, sol_$par[1], sol_$par[2], 1 - sol_$par[1] - sol_$par[2])
#       
#       sol_ <- constrOptim(c(sol[1, 2], sol[1, 3]),
#         function(p) {
#           EVSI(mu_a, mu_b, mu_c, tau_a, tau_b, tau_c, sigma2, N[1], p)
#         }, NULL, rbind(c(1, 0), c(0, 1), c(-1, -1)), c(0, 0, -1),
#         control=list(fnscale=-1)        
#       )
#       
#       sol[1, ] <- c(sol_$value, sol_$par[1], sol_$par[2], 1 - sol_$par[1] - sol_$par[2])
#       
#       for (i in seq_along(N)[-1]) {
#         sol_ <- constrOptim(c(sol[i - 1, 2], sol[i - 1, 3]),
#           function(p) {
#             EVSI(mu_a, mu_b, mu_c, tau_a, tau_b, tau_c, sigma2, N[i], p)
#           }, NULL, rbind(c(1, 0), c(0, 1), c(-1, -1)), c(0, 0, -1),
#           control=list(fnscale=-1)        
#         )
#         sol[i, ] <- c(sol_$value, sol_$par[1], sol_$par[2], 1 - sol_$par[1] - sol_$par[2])
#       }  
#       sol
#     }, mc.cores=n.cores
#   ))
# } 

library(manipulate)

manipulate({
  p <- c(1/3, 1/3)
  EVPI_ <<- EVPI(mu_a=mu_a, mu_b=mu_b, mu_c=mu_c, tau_a=tau_a, tau_b=tau_b, tau_c=tau_c)
  EVSI_ <<- EVSI(mu_a=mu_a, mu_b=mu_b, mu_c=mu_c, tau_a=tau_a, tau_b=tau_b, tau_c=tau_c,
    sigma2=sigma2, N=seq(3, N, len=sql), p=p)
  EVSI_opt_ <<- t(EVSI_opt(mu_a=mu_a, mu_b=mu_b, mu_c=mu_c, tau_a=tau_a, tau_b=tau_b, 
    tau_c=tau_c, sigma2=sigma2, N=seq(3, N, len=sql)))
  
  par(mfcol=c(2, 2), mar=c(1, 1, 0, 0), oma=c(4, 3.5, 1, 5))
  
  plot.new()  
  plot.window(xlim=c(0, N), ylim=c(0, 1.05))
  box()
  axis(2, las=1)
  mtext('p', side=2, line=3, las=1, font=3)
  abline(h=1/3, lwd=2, col='grey', lty=2)
  points(seq(3, N, len=sql), EVSI_opt_[, 2], type='l', lwd=2, lty=2)
  points(seq(3, N, len=sql), EVSI_opt_[, 3], type='l', lwd=2, col='red', lty=3)
  points(seq(3, N, len=sql), EVSI_opt_[, 4], type='l', lwd=2, col='blue', lty=4)
  mtext('a ', side=3, line=-1.3, adj=1)
  
  plot.new()  
  plot.window(xlim=c(0, N), ylim=c(0, .57))
  box()
  axis(1)
  axis(2, las=1)
  mtext('N', side=1, line=2.5, font=3)
  mtext('--- p=optimum', adj=0, side=1, line=3.5, col='red')
  mtext(sprintf('--- p=%s', round(p, digits=2)), adj=1, side=1, line=3.5, col='black')
  mtext('EVSI', side=2, line=3)
  abline(h=EVPI_, lwd=2)
  points(seq(3, N, len=sql), EVSI_, type='l', lwd=2, lty=3)
  points(seq(3, N, len=sql), EVSI_opt_[, 1], type='l', lty=2, col='red')
  text(0, EVPI_, sprintf('                         EVPI = %s', round(EVPI_, 3)), pos=3)
  mtext('c ', side=3, line=-1.3, adj=1)
   
  plot.new()  
  plot.window(xlim=c(0, 6), ylim=c(
    min(mu_a - 2.5 * (tau_a^-.5), mu_b - 2.5 * (tau_b^-.5), mu_c - 2.5 * (tau_c^-.5)), 
    max(mu_a + 2.5 * (tau_a^-.5), mu_b + 2.5 * (tau_b^-.5), mu_c + 2.5 * (tau_c^-.5))))
  box()
  axis(4, las=1)
  segments(c(1.5, 3.5, 5.5), c(mu_a - 1.96 * (tau_a^-.5), mu_b - 1.96 * (tau_b^-.5), mu_c - 1.96 * (tau_c^-.5)), 
    y1=c(mu_a + 1.96 * (tau_a^-.5), mu_b + 1.96 * (tau_b^-.5), mu_c + 1.96 * (tau_c^-.5)), , col=c('black', 'red', 'blue'))
  points(c(1.5, 3.5, 5.5), c(mu_a, mu_b, mu_c), pch=19, col=c('black', 'red', 'blue'))
  text(1.5, mu_a, bquote({mu*"'"}[a]==.(round(mu_a, 3))), pos=2, font=3, col='black')
  text(3.5, mu_b, bquote({mu*"'"}[b]==.(round(mu_b, 3))), pos=2, font=3, col='red')
  text(5.5, mu_c, bquote({mu*"'"}[c]==.(round(mu_c, 3))), pos=2, font=3, col='blue')
  text(.5, min(mu_a - 2.5 * (tau_a^-.5), mu_b - 2.5 * (tau_b^-.5), mu_c - 2.5 * (tau_c^-.5)) / 2, 
    bquote(sigma == .(round(sqrt(sigma2), 1))))
  mtext('b ', side=3, line=-1.3, adj=1)
   
  plot.new()  
  plot.window(xlim=c(0, N), ylim=c(0, max(EVSI_opt_[, 1] - EVSI_, .001) * 1.2))
  box()
  axis(1)
  axis(4, las=1)
  mtext('N', side=1, line=2.5, font=3)
  mtext(bquote('EVSI'^'p=optimum' - 'EVSI'^.(paste('p =', round(p, digits=2)))), side=4, line=4)
  optN <<- seq(3, N, len=sql)[which.max(EVSI_opt_[, 1] - EVSI_)]
  if (var(EVSI_opt_[, 1]  - EVSI_) == 0) optN <<- 0
  abline(v=optN, col='red')
  text(optN,
    max(EVSI_opt_[, 1] - EVSI_, .001) * 1.18, 
    round(optN), col='red', pos=4)
  abline(h=max(EVSI_opt_[, 1] - EVSI_), col='blue', lty=4)
  text(0, max(EVSI_opt_[, 1] - EVSI_),  
    bquote('                                                   EVSI'^
    'p=optimum' - 'EVSI'^.(paste('p =', round(p, digits=2))) == 
      .(round(max(EVSI_opt_[, 1] - EVSI_), 3))),
    pos=3, col='blue')
  points(seq(3, N, len=sql), EVSI_opt_[, 1] - EVSI_, type='l', lwd=3) 
  mtext('d ', side=3, line=-1.3, adj=1)
  
},
mu_a=slider(-5, 5, step=.01, initial=1),
mu_b=slider(-5, 5, step=.01, initial=0),
mu_c=slider(-5, 5, step=.01, initial=-1),
tau_a=slider(0.001, 5, initial=.1), 
tau_b=slider(0.001, 5, initial=.1),
tau_c=slider(0.001, 5, initial=.1),
N=slider(1, 50000, initial=50),
sigma2=slider(0.001, 225, initial=225),
sql=slider(3, 100, initial=100)
)


