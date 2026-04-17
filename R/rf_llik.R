rf_estot <- function(pars,qq)
{
   ph <- pars[1];
   ii <- pars[2];
   estot <- -ii * ph * log(1 - qq)/(ph - ii * log(1 - qq));
   return(estot);
}

rf_lesk <- function(pars,qq,k)
{
   logfun <- function(x) rf_lesk_int(x,pars,qq,k);
   lesk <- integral_peak(logfun);
   return(lesk)
}

rf_lesk_int <- function(x,pars,qq,k)
{
   ph <- pars[1];
   ii <- pars[2];
   y <- rep(0,length(x));
   for(cnt in 1:length(x))
   {
      y[cnt] <- sum(log(ii * x[cnt] + (0:(k - 1)))) + ii * x[cnt] * log(1 - qq) - ph * x[cnt];
   }
   y <- y + k * log(qq) - lgamma(k + 1) + 2 * log(ph);
   return(y);
}

rf_lestot_int <- function(x,pars,qq)
{
   ph <- pars[1];
   ii <- pars[2:length(pars)];
   iiqq <- -sum(ii * log(1 - qq));
   y <- rep(0,length(x));
   for(cnt in 1:length(x))
   {
      if(x[cnt] > 0)
      {
         y[cnt] <- logminexpm1approx(iiqq,x[cnt],1) - ph * x[cnt];
      } else
      {
         y[cnt] <- -Inf;
      }
   }
   y <- y + 2 * log(ph);
   return(y);
}
