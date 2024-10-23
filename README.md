# Confidence-Region-Optimiser


## Setup
In general settings, how does one construct a confidence region? We can make use of an MLE-based Wald's test, using the following theorem:

Let $X_1,X_2, ...$ be iid observations with pdf (or pmf) $f_{\theta}(x)$, where $\theta\in\Theta\in\mathrm{R}^d$. Let $\theta_0\in\mathrm{R}^d$ denote the true parameter. Under regularity conditions, the following holds: \
suppose $\hat{\theta}_n$ is a consistent sequence of MLEs. Then, 
$$\sqrt{n}(\hat{\theta}_n - \theta_0) \to^d MVN(0, I_f(\theta_0)^{-1})$$
Where $I_f(\theta_0)^{-1}$ is the Fisher Information of the sample.

In practise, we may assume this asymptotic relation to be an exact one, and we may thus construct our $1-\alpha$ confidence regions using the normal distribution. For $n$ large, this approximation is good. Can we improve it, and make it more effective for smaller $n$?

It turns out we can make use of the likelihood ratio test.

Let $Y_1, ..., Y_n$ be a random sample. Let $\mathbf{Y} = (Y_1, ..., Y_n)$. Under certain regularity conditions, $$2\log\frac{\sup_{\theta\in\Theta} L(\theta;\mathbf{y})}{\sup_{\theta\in\Theta_0} L(\theta;\mathbf{y})}\to^d\chi_r^2$$


In theory, the LRT, and by extension the likelihood function, is finer than WALD's theorem. We can take advantage of the asymmetric nature of the likelihood function to optimise an initial confidence region given by the WALD estimate.


## What the code does
Given number of covariates (dimension) and number of simulates, aims to optimise a "naive" confidence regions given by the Wald's test, and then optimise it using the likelihood ratio test.