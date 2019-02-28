---
title: "Using SADISA"
author: "Rampal S. Etienne"
date: "2019-02-27"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{"Using SADISA"}
  %\VignetteEngine{knitr::knitr}
  %\VignetteEncoding{UTF-8}{inputenc}
---

## SADISA introduction

SADISA contains three user-level functions, SADISA_loglik, SADISA_ML and SADISA_test. SADISA_loglik calculates the loglikelihood (sampling formula) of a set of abundances given a mainland-island community model, i.e. a model for the mainland (metacommunity) and the island (local community) that is connected to it via dispersal, and SADISA_ML finds the parameters that maximizes that loglikelihood. Tests that replicate the analyses in Haegema & Etienne (2016) are provided in SADISA_test. The model makes the Independent-Species (IS) assumption, and therefore differs from the sampling formulae that have been derived within the context of neutral theory, which all assume that individuals of species interact, but all in the same way (zero-sum constraint). Mathematically, the sampling formula is different, but in many cases the loglikelihood is practically indistinguishable from that obtained earlier, and so are the maximum likelihood parameters. The IS approach allows for a much wider set of data and models than the zero-sum constraint.

The package contains a fourth documented function, integral_peak, which is a numerical procedure for evaluating the integrals needed in the IS approach. The function is exported for use in other packages, but does not need to be called independently for the purposes of SADISA.

We represent the data, $\mathcal{D}$, of a single community sample as
abundance frequencies $s_{k}$ (number of species that are observed $k$ times
in the sample). The sampling formula is given by a product of Poisson
samples with parameters $\lambda _{k}$  

$\mathbb{P}(\mathcal{D})=\mathrm{e}^{-\Lambda }\prod_{k\mid s_{k}>0}\frac{\lambda _{k}^{s_{k}}}{s_{k}!}$

with $\lambda_{k}$ the expected number of observed species with abundance $k$,

$\lambda _{k} = \mathbb{E}s_{k}=\int_0^{\infty} \mathbb{P}(k|x)\,\rho (x)\,\mathrm{d}x$

and $\Lambda$ the expected number of observed species,

$\Lambda = \sum_{k>0}\mathbb{E}s_{k}=\int \mathbb{P}(\text{obs}|x)\,\rho (x)\,\mathrm{d}x$

where $\mathbb{P}(\text{obs}|x)$ the probability that a species with relative abundance x in the metacommunity is present in the data,

$\mathbb{P}(\text{obs}|x)=1-\mathbb{P}(0|x)$


## SADISA usage

First clear memory and load the package:


```r
rm(list = ls())
library(SADISA)
```

Then, choose an abundance data set.

For a single sample, you can use:


```r
abund <- c(1717,1681,983,788,755,724,681,644,617,381,379,376,364,346,345,325,322,294,289,
           288,285,264,248,244,236,236,229,218,203,201,188,184,177,167,164,163,156,149,147,
           143,121,118,111,101,100,99,98,98,98,93,92,92,88,87,85,85,82,81,80,78,76,70,68,68,
           67,67,64,63,63,61,58,55,55,55,54,52,52,51,50,49,47,45,45,43,43,41,40,39,39,38,38,
           36,33,33,33,33,33,32,31,30,29,29,28,28,28,27,27,27,26,26,26,26,25,25,25,25,23,23,
           23,23,22,22,22,21,21,21,21,20,19,18,17,16,16,16,15,15,15,14,14,14,13,13,13,13,13,
           12,12,12,12,12,12,12,10,10,10,10,10,10,10,9,9,8,8,8,7,7,7,7,7,7,6,5,5,5,5,5,5,5,5,
           4,4,4,4,4,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
           1,1,1,1,1)
```

This is one of the censuses of the 50 ha Barro Colorado Island forest plot in Panama.

Then, provide the model that you want the loglikelihood for. For instance, if we want a metacommunity that is neutral with point speciation, and a local community that is neutral but with limited dispersal from the metacommunity, the model is:


```r
model <- c('pm','dl')
```

with 'pm' for point mutation and 'dl' for dispersal limitation.

Then choose values for the model parameters. These are $\theta$, the fundamental biodiversity parameter, and $I$, the fundamental dispersal parameter. Let us choose:


```r
pars <- c(4.793930230665386e+01, 2.174825700126286e+03)
```

Let us now calculate the loglikelihood and store it in 'result':


```r
result <- SADISA_loglik(abund = abund,pars = pars,model = model)
result
```

```
## [1] -317.7003
```

Have a look at the test file SADISA_test.R included in the package for further examples of SADISA_loglik.

Let us now try to find the maximum likelihood parameters for this model and this data set. We need to initialize the optimization by setting initial parameter values and indicating which parameters need to be optimized, for example:


```r
initpars <- c(40,2000)
idpars <- c(1,1)
```

Here we are optimizing both parameters, because idpars has a 1 in each element. If we just wanted to optimize one parameter, we would use


```r
initpars <- c(40,2000)
idpars <- c(0,1)
```

So here we only optinmize I and give it an initial value of 2000. We must also specify whether some parameters are set equal to others or not. In this example, it does not make much sense to set the parameters equal, but we still show how to do this. We store this information in the variable 'labelpars'. In this case we use


```r
labelpars <- c(1,2)
```

which means that we are considering two parameter values (out of two) with labels 1 and 2. The 1 corresponds to the first parameter in initpars and idpars and the 2 to the second.

To run the optimization we type:


```r
result <- SADISA_ML(abund = abund,initpars = initpars,idpars = idpars, labelpars = labelpars,model = model)
```

```
## The loglikelihood for the initial parameter values is -320.9138.
## 
## Parameters after likelihood maximization:
## [1]   40.000 4016.133
## 
## Maximum loglikelihood:
##  -319.1055
```

Let us look at an example with multiple samples. We load the available data sets:


```r
data(datasets)
```

dataset$dset2.abunvec contains various data sets with multiple samples. Let us select the first one:


```r
abund <- datasets$dset2.abunvec[[1]]
abund
```

```
## [[1]]
##   [1] 1717 1681  983  788  755  724  681  644  617  381  379  376  364  346
##  [15]  345  325  322  294  289  288  285  264  248  244  236  236  229  218
##  [29]  203  201  188  184    0  177  167  164  163  156  149  147  143    0
##  [43]    0   22  121   20  118  111    0    0  101   21  100   99   98   98
##  [57]   98   10   93   92   92   50   88   87   85   85    0   82   81   80
##  [71]   78   76   70   52   68    0   68   67   67   64   63   63   61   58
##  [85]    0   55   55   55   54   12   52   51   49    0    0   47    0   45
##  [99]   45   43   43    0   41   40   39   39   38   38    0   36   25   19
## [113]   33   12   33   33   33   33   32    0   21   31   13   30   29    0
## [127]   29   28   28   28   27   27   27   26   26   26   26    0   25   25
## [141]   25    1   23   23    0   23   23   22   22   21    2    0   21    5
## [155]   15   18    0   17    1   16    4   16   16    0   15    0    0   15
## [169]    0   14   14   14    0    0   13   13   13    0   13    0   12   12
## [183]   12   12   12    0    0    0   10   10   10   10   10   10    9    9
## [197]    8    5    3    0    8    1    8    7    7    7    7    7    7    0
## [211]    2    3    6    0    0    0    0    5    5    0    5    0    5    5
## [225]    5    0    4    0    4    0    0    0    0    4    4    0    3    3
## [239]    3    0    3    0    0    0    3    3    0    0    3    0    0    2
## [253]    0    0    2    0    2    1    0    2    0    2    0    2    0    0
## [267]    0    2    2    0    2    0    2    2    1    0    1    1    0    1
## [281]    0    0    0    1    1    0    0    0    0    1    0    0    0    0
## [295]    0    0    1    0    0    0    0    0    0    1    1    1    0    1
## [309]    0    1    1    1
## 
## [[2]]
##   [1]  16   0   3   0   0   0   0   0   0  38   0   0   0   0   0   0   0
##  [18]   0   0  15   0   0   0   0   0   0   0  38   0   0  36   0 179   0
##  [35]   0   0   0   1   0  20   0   0   0 122   0   0   3   0 102   0   0
##  [52]   0   0   0   0   0   0   0   4   1   0   0   4   1   0   0   0   0
##  [69]   3   0   0   0   0   0   0   0   3   0   0   0   1   2   0   0  57
##  [86]   0  13   0   1   0   0   0   0   0   0   3   0   1   0   0   1  42
## [103]   0   0  36   0   0   1   0   0   0  34   0   0   2   0   0  31   0
## [120]   0  31   0   0   0   0   0  13  14   1   1   0   5   0   0   0   2
## [137]   0   0   0   6   7   0   0   0   0   0   0   0   0   0   0   0   0
## [154]   1  19   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
## [171]   0   0   0   0   0   0   0  13   0   0   0   0   0   0   0  12   0
## [188]  11   0  10   0   0   1   3   1   1   9   0   3   8   0   0   0   0
## [205]   0   0   0   0   1   6   6   6   0   0   0   5   5   0   2   0   0
## [222]   5   0   1   0   0   0   4   1   0   0   0   0   0   0   0   0   0
## [239]   0   0   0   0   0   3   2   0   0   3   0   0   3   0   0   2   2
## [256]   2   0   2   2   0   0   0   0   0   0   2   2   0   0   2   0   2
## [273]   0   0   1   2   0   0   1   0   1   1   1   0   0   1   0   0   1
## [290]   0   1   1   0   1   1   1   1   1   1   0   0   1   1   0   0   1
## [307]   1   0   1   0   0   0
## 
## [[3]]
##   [1]   1   3   0  54   8   0   0   1  58   0   0   1   0 227   0  36   4
##  [18]  10  77   2   0   6   0  53  28 218   0   0   7   0   0   1   0   0
##  [35]   7   5   0   0   0   0  20 129 126   0   3 120   0   0   0 102   0
##  [52] 101   0   0   0   0   8  97   0   0   0  89  57   6   7   0  82   0
##  [69]   0   0   1   0   1  69   0  68   0   2   0   0   0   0   0   2   0
##  [86]  28   0   0   6  54  29  10   3  48  48   0  46   0   0   0   0   0
## [103]   0  28   0   0   0   0  37   0  35   0   0  33   4   0   0   0   0
## [120]  32   0   0  30   0   0  29   0   6   0   3   0   2   1   0   0   0
## [137]   9  26   0   0   0  24   0   0  23   1   0   0   0   5  21  21   1
## [154]  19   0   0  18   0  16   6  16   2   0  15   1  15  15   0  14   3
## [171]   0   8  14  14   0   0   1   0   0  12   0   0   4   0   0   0  12
## [188]   0   0   0   0   0   0   0   0   0   0   9   8   0   0   8   0   0
## [205]   0   7   0   0   0   0   0   0   0   6   5   1   0   0   0   5   0
## [222]   0   0   0   3   5   0   0   1   4   4   4   4   0   0   4   0   3
## [239]   0   3   0   3   3   0   0   2   3   0   0   3   0   0   2   0   0
## [256]   0   0   0   0   0   2   0   2   0   2   0   0   0   0   0   1   0
## [273]   0   0   2   0   0   0   0   0   0   0   0   0   0   0   1   1   0
## [290]   0   1   0   1   0   0   0   0   0   0   1   1   0   0   0   0   0
## [307]   0   0   0   0   0   0
```

Because this data set contains three samples we need to specify parameters for three samples, which we take from the file fitresult.RData provided with the package:


```r
data(fitresults)
pars <- fitresults$fit4.parsopt[[1]]
```

Note that we have named the parameters. Let us fix $theta$ and optimize the the $I$ parameters assuming they are all equal. This can be coded in labelpars, initpars and idpars as follows:


```r
labelpars <- list(c(1,2),c(1,2),c(1,2))
labelpars
```

```
## [[1]]
## [1] 1 2
## 
## [[2]]
## [1] 1 2
## 
## [[3]]
## [1] 1 2
```

```r
initpars <- c(pars[[1]][1],pars[[1]][2])
initpars
```

```
## [1] 259.42034  44.45759
```

```r
idpars <- c(0,1)
```

Because the parameters to be optimized, I, are all assumed equal, they all receive a 2. The initial value of I will then be in initpars[2]. All thetas are equal (but there is of course only one theta), so they all receive a 1. Its fixed value is in initpars[1].

Now we run the optimization by adding the option mult = 'ms':


```r
result <- SADISA_ML(abund = abund,initpars = initpars,idpars = idpars, labelpars = labelpars,model = model,mult = 'ms')
```

```
## The loglikelihood for the initial parameter values is -1116.12.
## 
## Parameters after likelihood maximization:
## [[1]]
## [1] 259.42034  44.45759
## 
## [[2]]
## [1] 259.42034  44.45759
## 
## [[3]]
## [1] 259.42034  44.45759
## 
## 
## Maximum loglikelihood:
##  -1116.12
```

We observe that the starting values were indeed the optimal values.

Now we show an example where we no longer assume that all the I values are the same. We then use;


```r
labelpars <- list(c(1,2),c(1,3),c(1,4))
labelpars
```

```
## [[1]]
## [1] 1 2
## 
## [[2]]
## [1] 1 3
## 
## [[3]]
## [1] 1 4
```

```r
initpars <- c(pars[[1]][1],pars[[1]][2],pars[[2]][2],pars[[3]][2])
initpars
```

```
## [1] 259.42034  44.45759  44.45759  44.45759
```

```r
idpars <- c(0,1,1,1)
```

We run again the optimization by adding the option mult = 'ms':


```r
result <- SADISA_ML(abund = abund,initpars = initpars,idpars = idpars, labelpars = labelpars,model = model,mult = 'ms')
```

```
## The loglikelihood for the initial parameter values is -1116.12.
## 
## Parameters after likelihood maximization:
## [[1]]
## [1] 259.42034  63.43168
## 
## [[2]]
## [1] 259.4203  31.5572
## 
## [[3]]
## [1] 259.42034  35.50356
## 
## 
## Maximum loglikelihood:
##  -1100.448
```

Assume now that the three samples actually represent three different guilds. Again we want to optimize the I parameters, but we keep theta fixed. We use the same structure:


```r
labelpars <- list(c(1,2),c(1,3),c(1,4))
labelpars
```

```
## [[1]]
## [1] 1 2
## 
## [[2]]
## [1] 1 3
## 
## [[3]]
## [1] 1 4
```

```r
initpars <- c(pars[[1]][1],pars[[1]][2],pars[[2]][2],pars[[3]][2])
initpars
```

```
## [1] 259.42034  44.45759  44.45759  44.45759
```

```r
idpars <- c(0,1,1,1)
```

and run the optimization, but with the option mult = 'mg':


```r
result <- SADISA_ML(abund = abund,initpars = initpars,idpars = idpars, labelpars = labelpars,model = model,mult = 'mg')
```

```
## The loglikelihood for the initial parameter values is -530.6148.
## 
## Parameters after likelihood maximization:
## [[1]]
## [1] 259.4203  61.4082
## 
## [[2]]
## [1] 259.42034  32.55765
## 
## [[3]]
## [1] 259.42034  36.77088
## 
## 
## Maximum loglikelihood:
##  -522.001
```

Finally, for a combination of multiple samples and multiple guilds, we set parameters as follows:


```r
abund1a <- abund[[1]]
abund1b <- abund[[2]]
abund1c <- abund[[3]]
abund2 <- list(list(abund1a,abund1a),list(abund1b,abund1b),list(abund1c,abund1c))
labelpars <- list(list(c(1,2),c(1,3)),list(c(1,2),c(1,3)),list(c(1,2),c(1,3)))
labelpars
```

```
## [[1]]
## [[1]][[1]]
## [1] 1 2
## 
## [[1]][[2]]
## [1] 1 3
## 
## 
## [[2]]
## [[2]][[1]]
## [1] 1 2
## 
## [[2]][[2]]
## [1] 1 3
## 
## 
## [[3]]
## [[3]][[1]]
## [1] 1 2
## 
## [[3]][[2]]
## [1] 1 3
```

```r
initpars <- c(pars[[1]][1],pars[[1]][2],pars[[2]][2])
initpars
```

```
## [1] 259.42034  44.45759  44.45759
```

```r
idpars <- c(0,1,1)
```

Here we just made up a data set where the two guilds have identical species abundances. We set the I parameters for the two guilds to be different, but they are equal across samples. We fix theta, but optimize these two I values. We now simply run the same script as before, but with option mult = 'both':


```r
result <- SADISA_ML(abund = abund2,initpars = initpars,idpars = idpars, labelpars = labelpars,model = model,mult = 'both')
```

```
## The loglikelihood for the initial parameter values is -2232.24.
## 
## Parameters after likelihood maximization:
## [[1]]
## [[1]][[1]]
## [1] 259.42034  44.45759
## 
## [[1]][[2]]
## [1] 259.42034  44.45759
## 
## 
## [[2]]
## [[2]][[1]]
## [1] 259.42034  44.45759
## 
## [[2]][[2]]
## [1] 259.42034  44.45759
## 
## 
## [[3]]
## [[3]][[1]]
## [1] 259.42034  44.45759
## 
## [[3]][[2]]
## [1] 259.42034  44.45759
## 
## 
## 
## Maximum loglikelihood:
##  -2232.24
```

Have fun with SADISA!

Haegeman, B. & Etienne, R.S. (2016). A general sampling formula for community abundance data. Methods in Ecology & Evolution. In review.
