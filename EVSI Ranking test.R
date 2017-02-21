EVPI <- function(mu_a, tau_a, tau_b) {
  PI <- pnorm(mu_a / (sqrt(2) * ((tau_a^-.5) + (tau_b^-.5))))
  1 - pmax(PI, 1 - PI)  
}

EVSI <- function(mu_a, tau_a, tau_b, sigma2, N, p) {
  PI <- pnorm(mu_a / (sqrt(2) * ((tau_a^-.5) + (tau_b^-.5))))
  
  PI_prime <- pnorm(mu_a / (sqrt(2) * ((((tau_a * sigma2 + N * p) / sigma2)^-.5) + 
    (((tau_b * sigma2 + N * (1 - p)) / sigma2)^-.5)))
    )
  pmax(PI_prime, 1 - PI_prime) - pmax(PI, 1 - PI)
}

EVSI_opt <- function(mu_a, tau_a, tau_b, sigma2, N) {
  sapply(N, 
    function(N) {   
      optimize(
        function(p) {
          EVSI(mu_a, tau_a, tau_b, sigma2, N, p)
        },             
        interval=c(0, 1), maximum=TRUE
      )
    }
  )
} 
  




library(manipulate)

manipulate({
  EVPI_ <- EVPI(mu_a=mu_a, tau_a=tau_a, tau_b=tau_b)
  EVSI_ <- EVSI(mu_a=mu_a, tau_a=tau_a, tau_b=tau_b, sigma2=sigma2, N=0:N, p=p)
  EVSI_opt_ <- EVSI_opt(mu_a=mu_a, tau_a=tau_a, tau_b=tau_b, sigma2=sigma2, N=0:N)
  
  par(mfcol=c(2, 2), mar=c(1, 1, 0, 0), oma=c(4, 3.5, 1, 5))
  
  plot.new()  
  plot.window(xlim=c(0, N), ylim=c(0, 1.05))
  box()
  axis(2, las=1)
  mtext('p', side=2, line=3, las=1, font=3)
  abline(h=p, lwd=2, col='grey', lty=2)
  optP <- unlist(EVSI_opt_[1, ])[-1]
  points(1:N, optP, type='l', lwd=2)
  optP <- round(pmax(optP, 1 - optP), 3)
  maxOptP_N <- which.max((1:N)[optP == max(optP)])
  if (tau_a == tau_b) maxOptP_N <- 0
  abline(v=maxOptP_N, col='red')
  text(maxOptP_N, 1, maxOptP_N, pos=4, col='red')
  mtext('a ', side=3, line=-1.3, adj=1)
  
  plot.new()  
  plot.window(xlim=c(0, N), ylim=c(0, .57))
  box()
  axis(1)
  axis(2, las=1)
  mtext('N', side=1, line=2.5, font=3)
  mtext('--- p=optimum', adj=0, side=1, line=3.5, col='red')
  mtext(sprintf('--- p=%s', p), adj=1, side=1, line=3.5, col='black')
  mtext('EVSI', side=2, line=3)
  abline(h=EVPI_, lwd=2)
  points(0:N, EVSI_, type='l', lwd=2, lty=3)
  points(0:N, unlist(EVSI_opt_[2, ]), type='l', lty=2, col='red')
  text(0, EVPI_, sprintf('                   EVPI = %s', round(EVPI_, 3)), pos=3)
  mtext('c ', side=3, line=-1.3, adj=1)
  
  plot.new()  
  plot.window(xlim=c(0, 4), ylim=c(min(mu_a - 2.5 * (tau_a^-.5), -2.5 * (tau_b^-.5)), 
    max(mu_a + 2.5 * (tau_a^-.5), 2.5 * (tau_b^-.5))))
  box()
  axis(4, las=1)
  segments(c(1.5, 3.5), c(mu_a - 1.96 * (tau_a^-.5), -1.96 * (tau_b^-.5)), 
    y1=c(mu_a + 1.96 * (tau_a^-.5), 1.96 * (tau_b^-.5)))
  points(c(1.5, 3.5), c(mu_a, 0), pch=19)
  text(1.5, mu_a, bquote({mu*"'"}[a]==.(round(mu_a, 3))), pos=2, font=3)
  text(3.5, 0, bquote({mu*"'"}[b]==0), pos=2, font=3)
  text(.5, min(mu_a - 2.5 * (tau_a^-.5), -2.5 * (tau_b^-.5)) / 2, 
    bquote(sigma == .(sqrt(round(sigma2, 1))))) 
  mtext('b ', side=3, line=-1.3, adj=1)
  
  plot.new()  
  plot.window(xlim=c(0, N), ylim=c(0, max(unlist(EVSI_opt_[2, ]) - EVSI_, .001) * 1.2))
  box()
  axis(1)
  axis(4, las=1)
  mtext('N', side=1, line=2.5, font=3)
  mtext(bquote('EVSI'^'p=optimum' - 'EVSI'^.(paste('p =',p))), side=4, line=4)
  optN <- (0:N)[which.max(unlist(EVSI_opt_[2, ]) - EVSI_)]
  if (tau_a == tau_b) optN <- 0
  abline(v=optN, col='red')
  text(optN,
    max(unlist(EVSI_opt_[2, ]) - EVSI_, .001) * 1.18, 
       round(optN), col='red', pos=4)
  abline(h=max(unlist(EVSI_opt_[2, ]) - EVSI_), col='blue', lty=4)
  text(0, max(unlist(EVSI_opt_[2, ]) - EVSI_),  
    bquote('                                                   EVSI'^
      'p=optimum' - 'EVSI'^.(paste('p =', p)) == 
      .(round(max(unlist(EVSI_opt_[2, ]) - EVSI_), 3))),
    pos=3, col='blue')
  points(0:N, unlist(EVSI_opt_[2, ]) - EVSI_, type='l', lwd=3)
  mtext('d ', side=3, line=-1.3, adj=1)
  
    
  },
  mu_a=slider(.01, 10, step=.01, initial=1), 
  tau_a=slider(0.001, 5, initial=.001), 
  tau_b=slider(0.001, 5, initial=5),
  N=slider(1, 2000, initial=2000),
  sigma2=slider(0.001, 100, initial=100),
  p=slider(0, 1, initial=.5, step=.01)  
)
  





  