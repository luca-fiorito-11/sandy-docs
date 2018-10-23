***********************
Uncertainty propagation
***********************

Below is reported the mathematical basis to perform a nuclear data uncertainty 
propagation using sets of perturbed files generated with SANDY.


Monte Carlo uncertainty propagation
===================================
It is common to refer with the term *Monte Carlo uncertainty propagation* to 
computational approaches that rely on repeated calculation of a quantity 
:math:`y`, each time varying the input data randomly according to a given 
multivariate PDF.

It is common to refer with the term *Monte Carlo* to computational approaches 
that rely on a repeated random sampling of chosen parameters to obtain 
statistical outcomes of selected responses.

Now, lets us take :math:`x=\left[x_1, x_2, \dots, x_m\right]^T` as the input 
parameters for a model :math:`y = f(x)`.
If :math:`p(x, \Sigma)` is the joint PDF of the input data and :math:`\Sigma` 
is the covariance matrix, then :math:`X` is the sample matrix containing on 
its columns :math:`n` independent sets of samples of the type 
:math:`x^{(k)}=\left[x_1^{(k)}, x_2^{(k)}, \dots, x_m^{(k)}\right]^T` with 
:math:`k \in 1,\dots,n`.

.. math::
   X=\begin{bmatrix}
    x_1^{(1)} & x_1^{(2)} & ... & x_1^{(n)} \\
    x_2^{(1)} & x_2^{(2)} & ... & x_2^{(n)} \\
    \vdots    & \vdots    & \ddots & \vdots \\
    x_m^{(1)} & x_m^{(2)} & ... & x_m^{(n)} \\
   \end{bmatrix}\,.
   :label: x-samples

The samples are taken *with replacement*, that is, each point in the model input 
domain is not excluded after being sampled.

One can prove that for any given parameter :math:`i`, when 
:math:`n\rightarrow\infty`, then 

.. math::
   \overline{x}_i \equiv & \frac{1}{n}\sum_{k=1}^n x_i^{(k)} \rightarrow E\left[x_i\right]\,, \\
   s^2_{x_i} \equiv & \frac{1}{n}\sum_{k=1}^n \left( x_i^{(k)} - E\left[x_i\right]\right)^2 \rightarrow V(x_i)\,.
   :label: x-mean-variance

Analogously, as :math:`n\rightarrow\infty`, the samples' distribution and 
covariance matrix converge to :math:`p(x, \Sigma)`.

If the user wants to quantify the uncertainty of :math:`y` produced by 
:math:`x` --- assuming that the model :math:`f` does not carry any error --- 
then, one could symply use

.. math::
   E\left[f(x)\right] &= \int f(x) p(x,\Sigma) dx\,, \\
   E\left[ (f(x) - E\left[f(x)\right])^2 \right] &= \int f(x)^2 p(x,\Sigma) dx - E\left[f(x)\right] \,.
   :label: y-sample-mean-variance

The Monte Carlo sampling approach for uncertainty quantification uses 
a straightforward technique to reproduce the integrals reported above.
For any given set :math:`k` of sampled input parameters :math:`x^{(k)}` 
a model *perturbed* response is quantified as :math:`y^{(k)} = f(x^{(k)})`.


Then, the mean and variance of the population of :math:`n` perturbed responses 
can be calculated as

.. math::
   \overline{y} \equiv & \frac{1}{n}\sum_{k=1}^n y^{(k)}\,, \\
   s^2_{y} \equiv & \frac{1}{n}\sum_{k=1}^n \left( y^{(k)} - \overline{y}\right)^2\,,
   :label: y-mean-variance

which converge to the real solutions for :math:`n\rightarrow\infty`.

This procedure can be further extended to multi-response models.

Conclusions
===========
In short, the Monte Carlo sampling for uncertainty propagation simply 
consists in running the same model multiple times, every time replacing the input 
parameters with a new set of samples.
This approach brings numerous advantages compared to other uncertainty propagation 
techniques, such as perturbation theory:

 * its application to any model is straightforward, as it does not interact with the 
   model itself but only with the model inputs;

 * there is no need for developing complicated algorithms to calculate adjoint 
   functions for the model and responses under study;

 * it does not apply only locally around the input best estimates, but it can cover 
   the whole input and output phase space;

 * it can represent any-order effect of the model, as well as interactions between 
   model inputs.

The major drawback of this method lies in the large amount of sets of samples 
that must be drawn in order to grant the convergence of the response statistics, 
that is, to reduce the statistical error :math:`\epsilon` on the response best estimate. 
The central limit theorem tells us that this error is proportional to 
:math:`1/\sqrt{n}` and :math:`lim_{n\rightarrow\infty} \epsilon = 0`

The huge improvements of computer performances in the last decades, combined with 
model simplifications and dimensionality reductions help reduce the computational 
time of the solvers, thus making Monte Carlo sampling a practical option for 
nuclear data uncertainty propagation.