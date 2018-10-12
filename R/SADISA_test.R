#' @title Tests SADISA for data sets included in the paper by Haegeman & Etienne
#' @description Tests SADISA for data sets included in the paper by Haegeman & Etienne
#' @keywords model species-abundance-distribution
#' @references Haegeman, B. & R.S. Etienne (2017). A general sampling formula for community structure data. Methods in Ecology & Evolution. In press.
#' @export

SADISA_test <- function(tol = 1E-6)
{
   datasets = NULL; rm(datasets);
   fitresults = NULL; rm(fitresults);
   utils::data('datasets', package = 'SADISA');
   utils::data('fitresults', package = 'SADISA');
   cat('\n\nTesting PM+DL model - unconditional (Table 1):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset1.abunvec[[i]];
      po <- fitresults$fit1a.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dl'));
      cat('\nThe difference is:',result - fitresults$fit1a.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit1a.llikopt[[i]],tolerance = tol);
   }
   cat('\n\nTesting PM+DL model - conditional (Table 1):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset1.abunvec[[i]];
      po <- fitresults$fit1b.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pmc','dl'));
      cat('\nThe difference is:',result - fitresults$fit1b.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit1b.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting RF+DL model (Table 2):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset1.abunvec[[i]];
      po <- fitresults$fit2.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('rf','dl'));
      cat('\nThe difference is:',result - fitresults$fit2.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit2.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting MDD+DL model (Table 3):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset1.abunvec[[i]];
      po <- fitresults$fit3.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('dd','dl'))
      cat('\nThe difference is:',result - fitresults$fit3.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit3.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting multiple-samples model (Table 4):\n');
   for(i in 1:11)
   {
      nn <- datasets$dset2.abunvec[[i]];
      po <- fitresults$fit4.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dl'), mult = 'ms');
      cat('\nThe difference is:',result - fitresults$fit4.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit4.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting multiple-guilds model (Table 5):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset3.abunvec[[i]];
      po <- fitresults$fit5.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dl'));
      cat('\nThe difference is:',result - fitresults$fit5.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit5.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting PR+DL model (Table S1):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset1.abunvec[[i]];
      po <- fitresults$fit6.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pr','dl'));
      cat('\nThe difference is:',result - fitresults$fit6.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit6.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting PM+LDD model (Table S2):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset1.abunvec[[i]];
      po <- fitresults$fit7.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dd'));
      cat('\nThe difference is:',result - fitresults$fit7.llikopt[[i]],'  ');
      testthat:: expect_equal(result,fitresults$fit7.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting large data sets (Table S3):\n');
   # large data set
   for(i in 1:6)
   {
      nn <- datasets$dset4a.abunvec[[i]];
      po <- fitresults$fit8a.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dl'));
      cat('\nThe difference is:',result - fitresults$fit8a.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit8a.llikopt[[i]],tolerance = tol);
   }
   # small data set
   for(i in 1:6)
   {
      nn <- datasets$dset4b.abunvec[[i]];
      po <- fitresults$fit8b.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dl'));
      cat('\nThe difference is:',result - fitresults$fit8b.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit8b.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting multiple-samples + multiple-guilds model (new):\n')
   nn1 <- datasets$dset2.abunvec[[1]];
   nn2 <- datasets$dset2.abunvec[[2]];
   po1 <- fitresults$fit4.parsopt[[1]];
   po2 <- fitresults$fit4.parsopt[[2]];
   nn <- list();
   po <- list();
   for(i in 1:3)
   {
      nn[[i]] <- list(nn1[[i]],nn2[[i]]);
      po2[[i]][1] <- po1[[i]][1];
      po[[i]] <- list(po1[[i]],po2[[i]]);
   }
   result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dl'), mult = 'both');
   result2 <- SADISA_loglik(abund = nn1,pars = po1,model = c('pm','dl'), mult = 'ms');
   result2 <- result2 + SADISA_loglik(abund = nn2,pars = po2,model = c('pm','dl'), mult = 'ms');
   cat('\nThe difference is:',result - result2,'  ');
   testthat::expect_equal(result,result2,tolerance = tol);

   cat('\n\nTesting multiple-samples model for protracted speciation (new):\n');
   for(i in 1:11)
   {
      nn <- datasets$dset2.abunvec[[i]];
      po <- fitresults$fit4.parsopt[[i]];
      for(j in 1:3) {po[[j]] <- c(po[[j]][1],1e8,po[[j]][-1])}
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pr','dl'), mult = 'ms');
      cat('\nThe difference is:',result - fitresults$fit4.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit4.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting multiple-samples model for diversity-dependence (new):\n');
   for(i in 1:11)
   {
      nn <- datasets$dset2.abunvec[[i]];
      po <- fitresults$fit4.parsopt[[i]];
      for(j in 1:3) {po[[j]] <- c(po[[j]][1],1e-6,po[[j]][-1])}
      result <- SADISA_loglik(abund = nn,pars = po,model = c('dd','dl'), mult = 'ms');
      cat('\nThe difference is:',result - fitresults$fit4.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit4.llikopt[[i]],tolerance = tol);
   }

}
