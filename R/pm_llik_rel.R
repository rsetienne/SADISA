pm_estot_rel <- function(pars,qq)
{
   logfun <- function(x) exp(pm_estot_rel_int(x,pars,jj = qq));
   estot <- stats::integrate(f = logfun,lower = 0, upper = Inf)$value;
   return(estot);
}

pm_estot_rel_int <- function(x,pars,jj)
{
   th <- pars[1];
   ii <- pars[2];
   y <- log(th) - th * x - log(x) + log1p(-exp(ii * x * log(ii/(ii + jj))));
   return(y);
}

pm_lesk_rel <- function(pars,qq,k)
{
   logfun <- function(x) pm_lesk_rel_int(x,pars,jj = qq,p = k);
   lesk <- integral_peak(logfun);
   return(lesk);
}

pm_lesk_rel_int <- function(x,pars,jj,p)
{
   th <- pars[1];
   ii <- pars[2];
   y <- log(th) - th * x - log(x) +
        -ii * p + ii * x * log(ii) +
        -lgamma(ii * x) + (ii * x - 1) * log(p) +
        log1p(-exp(-jj * p));
   return(y);
}
