********
Sampling
********

``SANDY`` can work as a nuclear data sampling tool applicable to sensitivity 
analysis (SA) and uncertainty quantification (UQ) problems.
The code exploits the basic theory of Monte Carlo sampling to 
propagate nuclear data covariances through the nuclear models under study.

Definition of best-estimates and covariances
============================================

Nuclear parameters :math:`\mathbf{x}=\left[x_1, x_2, \dots, x_m\right]^T` --- e.g. 
cross sections, resonances, decay constants, particle energy-angular 
distributions, ... --- 
are stored in the ``ENDF-6`` nuclear data libraries as best-estimate or 
evaluated data

.. math::
   x_{i,BE} = E \left[ x_i \right] = \int x_i p(\mathbf{x}) d\mathbf{x}\,,

where :math:`E` is the expectation operator which performs the integral of 
the selected *i-th* parameter :math:`x_i` over the joint probability 
distribution of :math:`x`, :math:`p(\mathbf{x})`.
Best-estimate parameters are often used as inputs of numerical models.

UQ studies aim at propagating the model input uncertainties and/or covariances 
to the model responses.
Along this working line, the ``ENDF-6`` format allowed storing the covariance 
matrices :math:`\mathbf{\Sigma}` of the best-estimate nuclear data provided.

A covariance matrix :math:`\mathbf{\Sigma}` of the best-estimate data :math:`\mathbf{x}_{BE}` 
contains the so-called variance terms on the diagonal

.. math::
   \Sigma_{i,i} \equiv cov(x_i, x_i) \equiv V(x_i) = E \left[ (x_i - x_{i,BE})^2 \right]\,.

From the variance term of a parameter :math:`x_i`, one can extract its 
standard deviation (or uncertainty)

.. math::
   \Delta x_i \equiv \sqrt{V(x_i)} = cov(x_i, x_i)^{1/2}\,.

The off-diagonal terms contain the covariances, a generalization of the 
variance integral to two different variables

.. math::
   \Sigma_{i,j} \equiv cov(x_i, x_j) = E \left[ (x_i - x_{i,BE})(x_j - x_{j,BE}) \right]\,.

The covariance is a measure of the joint variability of two variables to the 
same source of error.

A covariance :math:`\Sigma_{i,j}` can be expressed in terms of the standard 
deviations of its parameters introducing the Pearson correlation coefficient 
:math:`\rho_{i,j}`:

.. math::
   cov(x_i, x_j) = \rho_{i,j} \Delta x_i \Delta x_j\,.

The Pearson correlation coefficient is a parameter that varies between -1 and 
1 and that can be interpreted as it follows:

 * :math:`\rho_{i,j}=1` implies that :math:`x_i` linearly increases as 
   :math:`x_j` increases;
 * :math:`\rho_{i,j}=-1` implies that :math:`x_i` linearly increases as 
   :math:`x_j` decreases;
 * :math:`\rho_{i,j}=0` implies that there is no linear correlation between 
   the variables.

A last important concept, which is often used in the rest of the manual, is 
the relative covariance,

.. math::
   rcov(x_i, x_j) = \frac{cov(x_i, x_j)}{x_{i,BE}x_{j,BE}}\,.

Monte Carlo Sampling
====================

It is common to refer with the term *Monte Carlo sampling* to 
computational approaches that rely on a repeated random sampling of chosen 
parameters to obtain statistical outcomes of selected responses.

Monte Carlo sampling methods are based on drawing sets of samples from the model 
input joint probability density function (PDF), as

.. math::
   \mathbf{X}=\begin{bmatrix}
    x_1^{(1)} & x_1^{(2)} & ... & x_1^{(n)} \\
    x_2^{(1)} & x_2^{(2)} & ... & x_2^{(n)} \\
    \vdots    & \vdots    & \ddots & \vdots \\
    x_m^{(1)} & x_m^{(2)} & ... & x_m^{(n)} \\
   \end{bmatrix}\,.
   :label: x-samples

Given input parameters :math:`\mathbf{x}=\left[x_1, x_2, \dots, x_m\right]^T` 
with a joint PDF :math:`p(\mathbf{x})` and covariance matrix :math:`\mathbf{\Sigma}`, 
:math:`\mathbf{X}` is the sample matrix containing on its columns :math:`n` 
independent sets of samples of the type 
:math:`\mathbf{x}^{(k)}=\left[x_1^{(k)}, x_2^{(k)}, \dots, x_m^{(k)}\right]^T` with 
:math:`k \in 1,\dots,n`. 
The samples are taken *with replacement*, that is, each point in the model input 
domain is not excluded after being sampled.

One can prove that for any given parameter :math:`i`, when 
:math:`n\rightarrow\infty`, then 

.. math::
   \overline{x}_i \equiv & \frac{1}{n}\sum_{k=1}^n x_i^{(k)} \rightarrow E\left[x_i\right]\,, \\
   s^2_{x_i} \equiv & \frac{1}{n}\sum_{k=1}^n \left( x_i^{(k)} - E\left[x_i\right]\right)^2 \rightarrow V(x_i)\,.
   :label: x-mean-variance

Analogously, as :math:`n\rightarrow\infty`, the samples' distribution and 
covariance matrix converge to :math:`p(\mathbf{x})` and :math:`\mathbf{\Sigma}`, 
respectively.

Now, let us assume that the user is interested in a model of the type 
:math:`y =f(\mathbf{x})`, where :math:`y` is the model response value, or 
simply *response*.
More specifically, the user wants to quantify the uncertainty of :math:`y` 
produced by :math:`\mathbf{x}` and assuming that the model :math:`f` does 
not carry any error.

The average - or best estimate - and uncertainty of the response are simply

.. math::
   y_{BE} = & E\left[f(\mathbf{x})\right] = \int f(\mathbf{x}) p(\mathbf{x}) d\mathbf{x}\,, \\
   \Delta y = & E \left[ (f(\mathbf{x}) - y_{BE})^2 \right] ^{1/2}\,.
   :label: y-sample-mean-variance

The Monte Carlo sampling approach for uncertainty quantification uses 
a straightforward approach to reproduce the integrals in equation 
:eq:`y-sample-mean-variance`.
For any given set :math:`k` of sampled input parameters :math:`\mathbf{x}^{(k)}` 
the model response is quantified as :math:`y^{(k)} = f(\mathbf{x}^{(k)})`.

.. Note::
   In the rest of the manual we often refer to :math:`y^{(k)}` as to a *perturbed* response, 
   since it reflects the variation of :math:`y` generated by the perturbations on 
   :math:`\mathbf{x}` introduced by the Monte Carlo sampling.

Then, the :math:`n` perturbed responses can be grouped into a matrix

.. math::
   \mathbf{Y}=\begin{bmatrix}
    y^{(1)} & y^{(2)} & ... & y^{(n)} \\
   \end{bmatrix}\,,
   :label: y-samples-1response

Notice that this matrix contains only one row, since so far we only considered models 
with a single response value.

Following what was previously done in equation :eq:`x-mean-variance` for an input 
parameter :math:`x_i`, the best estimate and uncertainty of :math:`y` can be 
approximated with

.. math::
   \overline{y} \equiv & \frac{1}{n}\sum_{k=1}^n y^{(k)}\,, \\
   s^2_{y} \equiv & \frac{1}{n}\sum_{k=1}^n \left( y^{(k)} - E\left[y\right]\right)^2\,,
   :label: y-mean-variance

which converge to the real solutions for :math:`n\rightarrow\infty`.

Should the model response be a vector :math:`\mathbf{y}=\mathbf{F}(\mathbf{x})`, 
where :math:`\mathbf{y}=\left[ y_1, y_2, \dots, y_l \right]`, then the matrix of 
perturbed responses becomes

.. math::
   \mathbf{Y}=\begin{bmatrix}
    y_1^{(1)} & y_1^{(2)} & ... & y_1^{(n)} \\
    y_2^{(1)} & y_2^{(2)} & ... & y_2^{(n)} \\
    \vdots    & \vdots    & \ddots & \vdots \\
    y_l^{(1)} & y_l^{(2)} & ... & y_l^{(n)} \\
   \end{bmatrix}\,.
   :label: y-samples

Any two values :math:`y_i` and :math:`y_j` could be two completely different model 
responses, or the same response evaluated in two different points of the output 
phase-space, e.g. different time, energy or space.
Equations :eq:`y-mean-variance` still apply to each of the :math:`l` components 
of :math:`\mathbf{y}`.
In addition, the covariance coefficient of :math:`y_i` with :math:`y_j` can 
be evaluated as

.. math::
   lim_{n\rightarrow\infty} \frac{1}{n}\sum_{k=1}^n \left( y_i^{(k)} - E\left[y_i\right]\right)\left( y_j^{(k)} - E\left[y_j\right]\right) = cov(y_i,y_j)\,.

In short, the Monte Carlo sampling for uncertainty propagation simply 
consists in running the same model multiple times, every time replacing the input 
parameters with a set of samples.
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
The central limit theorem *[ref]* tells us that this error is proportional to 
:math:`1/\sqrt{n}` and :math:`lim_{n\rightarrow\infty} \epsilon = 0`

The huge improvements of computer performances in the last decades, combined with 
model simplifications and dimensionality reductions help reduce the computational 
time of the solvers, thus making Monte Carlo sampling a practical option for nuclear 
data uncertainty propagation.


Sampling ENDF-6 nuclear data using SANDY
========================================

As we reported in the previous chapter, the ``ENDF-6`` files contain best 
estimate and covariance data for many of the data important for nuclear 
technology.

``SANDY``'s sampling is consistent with the work of the nuclear data evaluators, 
as it draws random nuclear data samples from the covariance matrices that 
are located in the same ``ENDF-6`` files.
Generally, a covariance matrix is strictly related to its corresponding 
best-estimate data evaluation.
The sampling procedure should make sure to maintain this relation, therefore 
it should be avoided sampling from covariance matrices that come from different 
evaluations.

However, to proceed with a Monte Carlo sampling approach, one must also know 
the probability distribution of the input nuclear data, but there is 
currently no defined way on how to store such piece of information 
is ``ENDF-6`` files.
A thorough study of the evaluated uncertainties and of the evaluation 
procedures could give some prior estimation of the distribution.
However, insights on the uncertainties are rarely available in the data files.
Hence, the sensible way to proceed would be for the user to get directly in 
touch with the data evaluators, although this step is often cumbersome, if not 
impossible.
Most of the times the nuclear data user is left with the only option of 
taking assumptions on the input PDF.

A popular assumption in the nuclear community is that nuclear data behave 
according to multivariate Normal distributions \cite{Zwermann_2011}.
Such an assumption was often inherent to the ENDF-6 format and the 
covariance processing \cite{Rochman_2011} and it is generally justified by 
the following theorems:

 * `Central limit theorem` : the sample mean of a number of samples of 
   a PDF with a defined mean and a finite variance, approaches a Normal 
   distribution as the number of samples approaches infinity, regardless of the 
   underlying choice of the PDF \cite{Mood_1974};

 * `Principle of maximum entropy` : if only the mean and variance are 
   given, the Normal distribution is the one with maximum entropy \cite{Cover_2006}.

However, sometimes the Normal assumption should be 
reconsidered \cite{Rochman_2011}.
For examples, the Gaussian distribution would no longer be justified if most 
of the uncertainty came from systematic errors 
\cite{Shlyakhter_1993}.

At its current state, ``SANDY`` assumes that all the nuclear data are 
Normally distributed.
Some development work is currently ongoing to add the capability to 
sample from different distributions, such as LogNormal and Uniform.
However, this feature is not yet available in the latest release of 
the code.




Sampling methodology in details
===============================

.. Let us take a vector of input parameters 
 :math:`\mathbf{x}=\left[x_1, x_2, \dots, x_m\right]^T` 
 defined by a multivariate Normal distribution

.. math:
   p(\mathbf{x}) = N(\mathbf{x}_{BE}, \mathbf{\Sigma})
   :label: pdf

.. centered on the best estimates :math:`\mathbf{x}_{BE}` and with
 covariance matrix :math:`\mathbf{\Sigma}`.

This section describes how ``SANDY`` transfers the best-estimate 
(:math:`\mathbf{x}_{BE}`) and covariance (:math:`\mathbf{\Sigma}`) 
information stored in the ``ENDF-6`` file to a matrix of 
samples :math:`\mathbf{X}` for which the sample distribution tends to 
the multivariate *Normal* distribution 
:math:`p(\mathbf{x}) = N(\mathbf{x}_{BE}, \mathbf{\Sigma})` 
as the number :math:`n` of set of samples goes to infinity.

First, ``SANDY`` draws a sample matrix :math:`\mathbf{Z}` from a standard 
distribution :math:`N(\mathbf{0},\mathbf{I})`, with :math:`\mathbf{0}` 
being a null vector and :math:`\mathbf{I}` the identity matrix.

.. math::
   \mathbf{Z}=\begin{bmatrix}
    z_1^{(1)} & z_1^{(2)} & ... & z_1^{(n)} \\
    z_2^{(1)} & z_2^{(2)} & ... & z_2^{(n)} \\
    \vdots    & \vdots    & \ddots & \vdots \\
    z_m^{(1)} & z_m^{(2)} & ... & z_m^{(n)} \\
   \end{bmatrix}\,.
   :label: x-samples

The samples in :math:`\mathbf{Z}` are independent and uncorrelated. 

At this point, we define a linear operator :math:`\mathbf{L}` that, when 
applied to uncorrelated standard samples :math:`\mathbf{Z}`, converts them into 
:math:`N(\mathbf{0},\mathbf{\Sigma})`-distributed samples.

..
  .. math::
  \mathbf{X}_{\normc} = \mathbf{L}\overline{\mathbf{X}}_{\norm}
  \label{eq:corr_eq}

Operator :math:`\mathbf{L}` does not affect the mean of the distribution, since 

.. math::
   E \left[ \mathbf{L} \mathbf{Z} \right] = \mathbf{L} E \left[ \mathbf{Z} \right] = \mathbf{0}\,.
   :label: preserve-mean

Also, the linearity of the operator preserves the shape of the distribution.

Then, the covariance becomes

.. math::
   E \left[ \left(\mathbf{L}\mathbf{Z}\right) \left(\mathbf{L}\mathbf{Z}\right)^T \right] = \mathbf{L} E \left[ \mathbf{Z} \mathbf{Z}^T \right] \mathbf{L}^T 
   = \mathbf{L} \mathbf{L}^T\,,
   :label: preserve-cov

since by definition :math:`E \left[ \mathbf{Z} \mathbf{Z}^T \right] = \mathbf{I}`.

The problem of sampling Normally distributed variables with a prescribed 
covariance matrix reduces to finding :math:`\mathbf{L}` such that

.. math::
   \mathbf{L} \mathbf{L}^T = \mathbf{\Sigma}\,.
   :label: chol_condition

By default ``SANDY`` calculates :math:`\mathbf{L}` by performing the 
Cholesky decomposition of :math:`\mathbf{\Sigma}`.
The Cholesky decomposition applies only to hermitian, positive definite 
matrices and returns a lower triangular matrix :math:`\mathbf{L}` that is 
unique, with real and positive diagonal elements, 
so that :eq:`chol_condition` is satisfied.

By definition, the properties of a covariance matrix are:

 * it is symmetric;
 * its diagonal elements are real and non-negative;

Consequently, only covariance matrices that are 
positive-definite satisfy the requirements for a Cholesky factorization.

.. Important::
   It is not rare to find evaluations for which the covariance 
   matrices have some negative eigenvalues.
   Such situations are in contradiction with the definition of a covariance 
   matrix and must be handled with care.

   ``SANDY`` cannot investigate into the covariance evaluation process to 
   identify the sources of this issue.
   Therefore, it reconstructs an approximate covariance matrix 
   :math:`\widetilde{\mathbf{\Sigma}}` with the negative eigenvalues set to zero.
   The assumption is acceptable for small negative eigenvalues that are likely 
   to sprout from round-off or truncation processes.
   For such cases it follows that :math:`\widetilde{\mathbf{\Sigma}} \approx \mathbf{\Sigma}`.
   When the negative eigenvalues are large, :math:`\widetilde{\mathbf{\Sigma}}` 
   is not anymore a good approximation of :math:`\mathbf{\Sigma}`.
   Then, the covariance matrix evaluation should be reconsidered as it does not 
   satisfy the requirements for the Monte Carlo sampling.

Applying the operator :math:`\mathbf{L}` to the uncorrelated samples 
:math:`\mathbf{Z}` one obtains a samples matrix :math:`\mathbf{L} \mathbf{Z}` 
that is :math:`N(\mathbf{0}, \mathbf{\Sigma})`-distributed, 
centered in :math:`\mathbf{0}`, with covariance matrix :math:`\mathbf{\Sigma}`.
Then, by translating the samples to the desired mean vector :math:`\mathbf{x}_{BE}` 
--- we can do it without changing neither the covariance, nor the shape of the 
distribution --- one obtains :math:`N(\mathbf{x}_{BE}, \mathbf{\Sigma})`-distributed 
samples

.. math::
   \mathbf{X} = \mathbf{x}_{BE} + \mathbf{L} \mathbf{Z}\,.


These results are valid when :math:`\mathbf{\Sigma}` contains absolute 
covariances.
However, most of the ``ENDF-6`` sections provide relative covariances, which 
implies that the samples :math:`\mathbf{L} \mathbf{Z}` are in relative units.
Then, samples :math:`\mathbf{X}` can still be obtained as 

.. math::
   \mathbf{X} &= \mathbf{x}_{BE} + \mathbf{x}_{BE} \circ \mathbf{L} \mathbf{Z} \\
              &= \mathbf{x}_{BE} \circ \left(\mathbf{1} + \mathbf{L} \mathbf{Z} \right)\\
              &= \mathbf{x}_{BE} \circ \mathbf{P}\,,

where :math:`\mathbf{1}` is a unitary vector and :math:`\circ` is the element-wise 
product.
We generally refer to :math:`\mathbf{P} = \mathbf{1} + \mathbf{L} \mathbf{Z}` as 
to the perturbation matrix of the system, since any :math:`P_{i,j}` element 
introduces a perturbation to the corresponding best-estimate :math:`x_{i,BE}` 
via a simple multiplication.

.. math::
   X_{i,j} = x_{i,BE} P_{i,j}\,.






Cross sections (MF33)
=====================

Best-estimates
--------------

``ENDF-6`` files contain evaluated, best-estimates reaction cross sections 
:math:`\sigma(E)` in section ``MF3`` as a function of the incident neutron 
energy :math:`E`.
In addition to tabulated energy-cross section pairs 
:math:`\left[ E_i,\sigma(E_i) \right]`, interpolation schemes are 
given that specify the energy variation of the data for incident energies 
between two adjacent energy points.
Using that information, linearized cross sections can be reconstructed over an 
unionized energy grid.
The linearized cross section format is called pointwise-``ENDF`` (``PENDF``) 
and it is particularly suitable for codes working with continuous energies 
since a cross section value can be extracted at any energy 
point via a simple linear interpolation.

The ``MT`` numbers that define the most common reaction cross sections are:
 * :``MT1``: : total cross section;
 * :``MT2``: : elastic cross section;
 * :``MT4``: : inelastic cross section;
 * :``MT18``: : fission cross section;
 * :``MT51-91``: : discrete and continuous levels of the inelastic cross section;
 * :``MT102``: : radiative capture cross section.

Other ``MT`` numbers are reserved for more exotic reactions, for which the 
description is reported in the ``ENDF-6`` manual.

The information stored in the various ``MT`` sections is often *redundant*, as 
some cross sections can be derived from others via summation rules.

.. Note::
   A cross section is defined as *redundant* when it can be completely calculated 
   from the combination of other cross sections, according to the conservation 
   laws.

Example 1: the total cross section can be reconstructed as the sum of all the 
partial (non-redundant) cross sections::
  MF3/MT1 = MF3/MT2 + MF3/MT4 + MF3/MT18 + MF3/MT102 + ...

Example 2: the inelastic cross section can be reconstructed as the sum of 
all the discrete and continuous inelastic levels::
  MF3/MT4 = MF3/MT51 + MF3/MT52 + ... + MF3/MT90 + MF3/MT91

Covariances
-----------
..
 In the last decades the nuclear data evaluators moved their attention 
 towards the production of adequate cross section uncertainties and covariance 
 matrices.
 Such uncertainties are meant to represent the current degree of knowledge 
 that we can acquire from experimental measurements and the laws of physics.
 Since cross section evaluations are derived from nuclear models that are based 
 on the experimental measurements of some model parameters, cross section 
 covariances were produced from the propagation of the uncertainties 
 associated to those parameters.
 The covariance matrices were tuned to integral experiments to be suitable 
 for current reaction applications.

.. %%% DESCRIPTION OF MF33

Section ``MF33`` of the ``ENDF-6`` format files contains the covariances
of the neutron cross sections appearing in section ``MF3``.

Covariances are given for:

 * the cross section of a reaction :math:`r` evaluated at different energies :math:`E_i` and :math:`E_j`: 
   :math:`cov(\sigma_r(E_i),\sigma_r(E_j))`;

 * the cross sections of different reactions :math:`r` and :math:`p` evaluated at different energies :math:`E_i` and :math:`E_j`: 
   :math:`cov(\sigma_r(E_i),\sigma_p(E_j))`;

 * the cross sections of different materials.

The latter option is not yet implemented in ``SANDY``.
However, it has to be said that only few are the files that 
provide such covariances in the current general-purpose nuclear data 
libraries.

The covariance information is expressed either explicitly in the 
so-called *NI-Type* sections, or as a combination of explicit 
covariances in the *NC-Type* sections.


NI-Type section
~~~~~~~~~~~~~~~

The *NI-Type* sections provide group-averaged covariances that can 
be arranged into a multigroup covariance matrix :math:`\mathbf{\Sigma}`. 

.. figure:: figures/corr_bi209_ngamma.pdf
   :scale: 60 %

   Multigroup correlation matrix.

Although their multigroup structures, covariance matrices can also apply to 
continuous energy data.
Let us take a cross section :math:`\sigma(E)` that is a continuous 
function of the incident-neutron energy :math:`E`.
:math:`\mathbf{\Sigma}` is the multigroup covariance matrix of 
:math:`\sigma(E)`, with each term :math:`\Sigma_{g,g'}` being a constant 
value defined over the energy groups :math:`g` and :math:`g'`.
Then, taken two evaluations of the cross section at energies :math:`E_i` and 
:math:`E_j`,

.. math::
   cov(\sigma(E_i),\sigma(E_j))=\Sigma_{g,g'} \, \iff \, E_i \in g,\, E_j \in g'\,.

Each *NI-Type* section consists of tabulated functions and a 
parameter ``LB``.
The parameter ``LB`` defines how the tabulated functions must be used 
to directly reconstruct either the relative covariance matrix 
(``LB=1,2,3,4,5,6``) or the absolute covariance matrix (``LB=0,8,9``).

 * :``LB=0``: : only the absolute variances are given for one group structure assuming all energy groups to be uncorrelated;

 * :``LB=1``: : only the relative variances are given for one group structure assuming all energy groups to be uncorrelated;

 * :``LB=2``: : relative components are given for one group structure assuming all energy groups to be fully correlated;

 * :``LB=3``: : relative components are given for two group structures assuming all energy groups to be fully correlated;

 * :``LB=4``: : relative components are given for one group structure. A second group structure is provided to define ??? 

 * :``LB=5``: : all the components of the covariance matrix are explicitly given. The matrix is square and defined for a single group structure. If the matrix is symmetric, only the upper triangular part is provided.

 * :``LB=6``: : all the components of the covariance matrix are explicitly given. The matrix is not square square and two group structures are provided, one per matrix dimension. This option is often used to provide a covariance matrix for cross sections of two different reaction types or materials, which might be defined on different energy grids.

 * :``LB=8,9``: : absolute variance components are given. These formats are dedicated for multigroup processing codes, since the actual cross section variance depends on the user's choice of the energy intervals along which the cross sections are averaged. In continuous energy, the cross section are pointwise and the energy intervals tend to zero.

A covariance matrix can be given from a single *NI-Type* section or from a 
combination of several.
For extra information on the covariance reconstruction rules one 
could refer to the ``ENDF-6`` or to the ``NJOY`` manual.

.. Warning::
   At the current stage, ``SANDY`` does not process *NI-Type* sections 
   with ``LB=0,8,9``, that is, where the information to reconstruct 
   absolute covariance matrices is given.
   Hence, for the rest of the chapter we refer only to relative 
   covariance matrices.



NC-Type section
~~~~~~~~~~~~~~~

The *NC-Type* sections do not contain information to reconstruct the 
covariance matrices, but they describe how *derived* cross 
sections may be obtained in given energy ranges in terms of:

 * the summation of other cross sections for which the covariances
   are explicitly given;
 * the ratio to other cross sections fro which the covariances
   are explicitly given.

Correlations between cross section arise implicitly from the 
"derivation relations" described in the *NC-Type* cross sections.

``SANDY`` does not process any *NC-Type* section, but it rather applies 
pre-defined conservation laws to reconstruct the so-called 
*redundant* cross sections.

.. note:
   The neutron total cross section is defined as the sum of all the neutron 
   partial cross sections.
 
   The neutron fission cross section is defined as the sum of the 
   multi-chance fission cross sections


The process of *redundant* cross section reconstruction is implemented 
automatically in ``SANDY`` and applies to each set of unperturbed and 
perturbed cross sections, independently.

Lumped reaction covariances
---------------------------
Lumped reaction covariances are not covered by ``SANDY``.


Sampling
--------

``SANDY`` draws :math:`N` samples from the :math:`M`-group relative 
covariance matrix using a multivariate normal distribution centered 
in zero:

.. math::
    N(1, RCOV) \rightarrow 
    p = 
    \begin{bmatrix}
    p_1^{(1)}  & p_1^{(2)} & p_1^{(3)} & \dots  & p_1^{(N)} \\
    p_2^{(1)}  & p_2^{(2)} & p_2^{(3)} & \dots  & p_2^{(N)} \\
    \vdots     & \vdots    & \vdots    & \ddots & \vdots \\
    p_M^{(1)}  & p_M^{(2)} & p_M^{(3)} & \dots  & p_M^{(N)}
    \end{bmatrix}

Such samples are often referred as perturbation factors as they directly 
apply a perturbation to the cross section with a simple multiplication.

A perturbation factor :math:`p_i^{(k)}` refers to a given energy group 
of the covariance matrix :math:`E_i \leq E \leq E_{i+1}`.
Hence, it perturbs only the tabulated cross section points 
:math:`x(E_j)^{(k)}` for any energy :math:`E_j` that belongs to the 
covariance energy-group 

.. math::
   x(E_j)^{(k)} = x(E_j) \sum_i p_i^{(k)} \delta_{i,j}

where:

 * :math:`x(E_j)` is the tabulated cross section at energy :math:`E_j`;
 * :math:`x(E_j)^{(k)}` is the perturbed tabulated cross section at energy :math:`E_j`;
 * :math:`p_i^{(k)}` is the perturbation coefficient for the energy group :math:`E_i \leq E \leq E_{i+1}`;
 * :math:`\delta_{i,j}=1` when :math:`E_i \leq E_j \leq E_{i+1}`, otherwise :math:`\delta_{i,j}=0`.

As a consequence, any two tabulated cross section points :math:`x(E_i)` and 
:math:`x(E_j)`, for which :math:`E_i` and :math:`E_j` belong to the same 
covariance group, are assumed 100% correlated.


Subtraction rules
^^^^^^^^^^^^^^^^^
Often, the summation rules given in a *NC-Type* section are identical 
to the pre-defined laws used by ``SANDY``.
However, sometimes they are defined as the difference between 
other cross sections.
When this happens, there are chances that negative cross section values arise.
It is the case of :math:`^{16}O`, for which the inelastic cross section 
:math:`\sigma_{inel}` must be derived by subtracting the elastic 
cross section :math:`\sigma_{el}` and all the other partial cross 
sections (:math:`\sigma_r`) from the total cross sections :math:`\sigma_{tot}`, 
according to the *NC-Type* section in the 
``ENDF/B-VII.1`` file.

.. math::
    \sigma_{inel} = \sigma_{tot} - \sigma_{el} - \sum_r \sigma_r
    :label: sumO16
 
At thermal energy (0.0253 eV), the total and elastic cross sections are 
respectively

.. math::
   \sigma_{tot} = 3.9620963\,\,b \,\,\,\,\,\,\,\,\,\,\,\,\,\,
   \sigma_{el} = 3.961907\,\,b.

The two values differ by only 0.005%, since the elastic scattering is by far 
the dominant interaction at that energy.
According to the ``MF33`` section of ``ENDF/B-VII.1``, the two reactions 
are not correlated (no *NC-Type* section, nor an explicit covariance 
between the two reactions is given) and have both an uncertainty of 2%.

By drawing independent samples for the two reactions, there is a high probability 
that the perturbed elastic cross section be larger than the total cross section.
Then, according to the equation :eq:`sumO16`, the perturbed inelastic cross 
section :math:`\sigma_{inel}` would be negative.

It must be said that the ``ENDF-6`` rules advise against evaluating 
*NC-Type* sections for a given parameter as the difference of two other 
parameters with a similar larger order of magnitude.
However, since at the current stage there is no algorithm implemented in 
``SANDY`` to distinguish and exclude *NC-Type* sections such as 
the one evaluated for the :math:`^{16}O` inelastic cross section 
by ``ENDF/B-VII.1``, the code simply skips all the *NC-Type* sections.

Correlation between a *redunant* cross section and its components
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Sometimes a covariance matrix is given explicitly for a *redundant* 
cross section, but none is provided for its components, neither 
explicitly nor using a *NC-Type* section.

By definition, a redundant cross section :math:`\sigma_R` can be obtained 
from the sum of its components :math:`\sigma_c`, as 

.. math::
   \sigma_R = \sum_c \sigma_c\,,
   :label: sum

which implies correlations between the reactions.

By sampling the redundant cross sections and leaving its components 
unperturbed, the correlations between channels would be lost, as well 
as the physics expressed in equation :eq:`sum`
``SANDY`` reintroduces the correlations into the system by 
assigning the same perturbation coefficients used for 
:math:`\sigma_R` to all the :math:`\sigma_c` cross sections.
By doing so, every set of sample comply with the summation rule 
in equation :eq:`sum`.

.. Important::
   ``SANDY`` apportions the samples of a *redundant* cross section to 
   its components only when of its components has got its own samples.



Sampling options
~~~~~~~~~~~~~~~~

Use command MF : 33

cross sections are linearized with tolerance 0.1% using RECONR
option
linearization is important to have enough cross section points to cover the whole cross section energy domain, and not lose information from the multigroup 
covariance matrix.

Whole cross sections are perturbed, unless MT option is defined

All correlations explicitly given, except NC-Type are included, unless option 
is used

If the user whishes to erase all cross section correlations, including the 
correlations find in each *NI-Type* section for a single reaction... 
In the testing of ``SANDY``, we found this option particularly useful to bypass 
the problem of working with non positive-definite matrices (a covariance 
matrix without off-diagonal terms and with positive variances is always 
positive-definite).
Another important application for this option could be in linear regression 
analysis, to prevent problems related to multicollinearity.



Fission Neutron multiplicities (MF31)
=====================================



Resonance parameters (MF32)
===========================

Best-estimates
--------------
..
  Reaction cross section exhibit a strong dependence on the incoming neutron 
  energy.
  In particular, at certain energies the cross section suddenly rises of several 
  order of magnitudes creating peaks in the energy dependent function.
  Such energies correspond to excited states of the compound nucleus generated 
  by the interaction of the neutron with the target nucleus.
  The local increase of interaction probability is called resonance.
  At low neutron energies the resonances are well separated and can be observed 
  directly, or *resolved*, by experimental measurement.
  These resonances are found in the so called resolved resonance region (RRR).
  As the neutron energy increases the distance between resonances becomes 
  narrower, until they overlap and cannot be resolved anymore due to 
  the limits on the instrumental resolution.
  This region is called unresolved resonance region (URR).
  Only averages of resolved resonance parameters over energy are given in the 
  URR.
..
  Resonances are parts of the cross sections in all respects.
  Thus, resonance parameters impact on the same model calculations reported in 
  \autoref{sec:cross-sections}.

..
  The most evident effect of the resonances is visible on the neutron 
  spectra in reactor applications.
  Because of the sudden increase of probability for a neutron to be scattered 
  or captured, the neutron population in the reactor decreases in 
  correspondence of the resonance energies.
  Such effects result in localized dips in the neutron flux.

Resonance parameters are allowed by the ``ENDF-6`` format in section ``MF2``
Only one ``MT`` section --- ``MT151`` --- is given, which contains all the 
resonance data.
The resonances are divided into several ranges according to the 
incident neutron energy.
A given range corresponds either to the resolved resonance region (RRR) --- 
where resonance parameters for individual resonances are given --- or to 
the unresolved resonance region (URR) --- where the experimental resolution 
is not adequate to separate two non-oberlapping resonances.
No capability for the processing of the URR is currently implemented in 
``SANDY``, therefore this URR will not be described in this manual.

.. Important::
   The ``ENDF-6`` format allows storing resonance parameters for a maximum 
   of 10 isotopes in the same section ``MF2``, each isotope defined by a 
   given abundance.
   Notice that ``SANDY`` can only work with ``MF2`` sections were 
   resonance parameters are given only for one isotope.

Historically, six were the accepted representations for the resonance 
parameters in the RRR: 

 * single or multi-level Breit-Wigner (BW);
 * Reich-Moore (RM);
 * Adler-Adler;
 * General R-matrix;
 * Hybrid R-function;
 * R-matrix limited format.

However, the General R-matrix and the Hybrid R-function formats are no longer 
available in ``ENDF-6``.

Currently, ``SANDY`` can only sample resolved resonance parameters given in 
BW or RM formalism, which cover the majority of the current evaluations.
However, because of the increasing popularity of the R-matrix limited 
formalism, it is foreseen to add extra capabilities to the code to manage 
resonance parameters expressed in such a representation.

Resolved resonances are divided in further subsections according to 
the orbital angular momentum (or *l*-value).
For each subsection, the resolved resonances are listed in consecutive 
lines where the following parameters are explicitly given.

 * single- or multi-level Breit-Wigner formalism:

   * :math:`E_r` : the resonance energy in the laboratory system;
   * :math:`J` : the spin of the resonance;
   * :math:`\Gamma` : total resonance width evaluated at the resonance energy :math:`E_r`
   * :math:`\Gamma_n` : neutron width evaluated at the resonance energy :math:`E_r`
   * :math:`\Gamma_g` : radiation width evaluated at the resonance energy :math:`E_r`
   * :math:`\Gamma_f` : fission width evaluated at the resonance energy :math:`E_r`

   According to the conservation laws, the total width is equal to the sum of the 
   partial widths
   
   .. math::
      \Gamma = \Gamma_n + \Gamma_g + \Gamma_f + \Gamma_x\,,
      :label: sum-widths
   
   where :math:`\Gamma_x` is the competitive width (normally inelastic 
   scattering to the first excited state), which is not explicitly given.
   In case a competitive width exist, it must be recovered from equation 
   :eq:`sum-widths`.

 * Reich-Moore formalism:

   * :math:`E_r` : the resonance energy in the laboratory system;
   * :math:`J` : the spin of the resonance;
   * :math:`\Gamma` : total resonance width evaluated at the resonance energy :math:`E_r`
   * :math:`\Gamma_n` : neutron width evaluated at the resonance energy :math:`E_r`
   * :math:`\Gamma_{f_A}` : first partial fission width evaluated at the resonance energy :math:`E_r`
   * :math:`\Gamma_{f_B}` : second partial fission width evaluated at the resonance energy :math:`E_r`

   For the Reich-Moore formalisms, competitive reactions are not used.

Both formalism provide a scattering radius :math:`a_p` that can be constant 
or energy-dependent.
Also, in the Reich-Moore formalism, different radii are accepted for the 
several *l*-values.

The best-estimate parameters given in ``MF2`` is also repeated in 
``MF32`` for most of the resolved resonance.
As a general rule, the best-estimate parameters of a given resonance in 
section ``MF2`` are also reported in section ``MF32`` as long as they have 
a non-zero variance.
If this is not the case, the resonance is omitted from ``MF32``, thus 
preventing from storing zero-terms in the covariance section.

.. Notice::
   The list of resolved resonances in ``MF32`` can be contain 
   less resonances than the list in ``MF2``.


Covariances
-----------
Section ``MF32`` of the ``ENDF-6`` format contains covariances for most of 
(or all) the resonance parameters given in section ``MF2``.
Again, the only ``MT`` section available is ``MT151``.
The covariances contain either

 * *long-range* components, which follow the energy dependence of some 
   parameters, thus spanning over several resonances;
 * *short-range* components, which refer to parameters the neighborhood of 
   individual resonances.

In the process of testing ``SANDY``, we could not find any long-range covariance 
section in the most common general purpose nuclear data files in ``ENDF-6`` format.
.. No long-range in JEFF-3.2
   No long-range in ENDF/B-VII.1
As a consequence, no development effort was made to include in ``SANDY`` some 
capabilities for the processing of the long-range covariances.

FOOTNOTE
Many data allowed by the ``ENDF-6`` format are only rarely used by the 
evaluators.
If 

.. Note::
   ``SANDY`` raises an error message when long-range covariances are found.

In general, the format for the short-range covariances depends on the 
resonance formalism, however those used for BW and RM can be considered 
equivalent.
A compatibility flag ``lcomp`` defines the format in which the 
covariances are given:

 * :``lcomp=0``: : covariances can only be given between parameters for 
   the same resonance. This format is applicable only for the BW formalisms.
   Also, no long-range covariances can be defined.

 * :``lcomp=1``: : all the resonances for which covariances are to be 
   included are divided into blocks --- or short-range sections.
   Covariances between parameters can only be included for resonances 
   in the same block.
   ``SANDY`` returns an error if more than one block is given.
   Long-range covariances can be defined, but they are not covered 
   by ``SANDY``.

 * :``lcomp=2``: : standard deviations are given together with the 
   best estimates.
   Correlations coefficients (in compact form) follow.
   This format is suitable for sparse matrices, which would be impractical 
   to store in ``lcomp=1`` format.

.. Note::
   Contrary to most of the other covariance sections, absolute covariances 
   are used in ``MF32``.

In contrast with reaction cross section covariances in section ``MF33``, 
the covariances in section ``MF32`` have a bijective relation with their 
corresponding resonance parameters and do not introduce any multigroup 
assumption.


scattering radius



Sampling
--------
SANDY draws perturbation coefficients from the data and covariances 
in section ``MF32``.
The perturbation coefficients are applied to the best-estimate in ``MF32`` 

Recall that MF2 and MF32 do not have the same number of resonances.
To produce the perturbed output files, the parameter in ``MF2`` is replaced 
by samples when its analogous is found ``MF32`` 
Equivalence determined by equal spin and energy (tolerance 1e-5%)

.. Note::
   In contrast with reaction other energy-dependent data, the covariances in 
   section ``MF32`` have a bijective relation with their corresponding resonance 
   parameters.

Sometimes, the best-estimate resonance parameters differ between ``MF2`` 
and ``MF32``.
It is evident that such a situation is a mistake in the evaluation.
``SANDY`` detects such discrepancies and raises a warning message.
However, from the mere reading of the ``ENDF-6`` data the code cannot 
make a sensible choice based on physics.
Hence, during the sampling process, the perturbation factors are 
applied to the ``MF32`` best-estimates, to keep consistency 
with the covariances given in the same section.


``SANDY``'s sampling routine accounts for any correlation between different 
channels and energy levels.
On top of that, when given explicitly (i.e. Breit-Wigner format), 
``SANDY`` reconstruct the total width samples 
according to the physical constraint


Output files
------------

HENDF


Sampling options
----------------

MF

MT

stdmax

no_correlations





Angular distributions (MF34)
============================

In a neutron-nucleus reaction, particles can be emitted to specific 
directions according to the physics of the interaction.
The probability that a particle is emitted into the interval :math:`d\mu` 
about an angle whose cosine is :math:`\mu` following an interaction between 
a nucleus and a neutron with incident energy :math:`E` can be expressed by 
normalized distributions, such as:

.. math::
   \int_{-1}^{1} f(\mu, E) d\mu = 1

where :math:`f(\mu, E)` is the probability density function by unit angle and 
incident energy.

Since the angular distribution of secondary particles is generally assumed to 
have azimuthal symmetry, the dependence on the outgoing angle and 
incident-neutron energy can be decoupled by representing the distribution as 
a Legendre polynomial series:

.. math::
   f(\mu, E) = \sum_{l=0}^{NL}\frac{2l+1}{2}a_l{E}P_l{\mu}

where:
 * :math:`l` is the order of the Legendre polynomial;
 * :math:`NL` is the highest order of the Legendre polynomial that is given at each energy;
 * :math:`a_l(E)` is the *l*-th energy-dependent Legendre polynomial coefficient;
 * :math:`P_l(\mu)` is the *l*-th Legendre polynomial.

The angular distributions for secondary particles are provided in the `ENDF-6` 
format file in section `MF4`.
Such data can be stored in a file section with one of the following four 
formats:

 * purely isotropic angular distributions;
 * Legendre polynomial coefficients;
 * tabulated probability distributions;
 * angular distribution over two energy ranges.

The covariance matrices for the angular distributions of secondary particles 
are provided in the ``ENDF-6`` format file in section ``MF34``.
Because of the simplicity of representing the covariances of Legendre 
coefficients rather than normalized probability components, only the former is 
considered in the ``ENDF-6`` format, even for cases where ``MF4`` has tabulated 
probability distributions.


Sampling correlations
---------------------

Covariances for the Legendre polynomial coefficients are provided as average 
values over a multi-group structure of the incident-neutron energy.
``SANDY``'s sampling routine can process any correlation between energy 
groups and polynomial coefficients.
However, it does not handle correlations between different reactions or 
materials.


Adding extra points to the incident-neutron energy grid
-------------------------------------------------------
Given the multi-group energy structure of the covariance matrices
angular distributions are provided according 
to multigroup structures along the incident-neutron energy
``SANDY``'s sampling routine includes correlations between of Legendre 
polynomial coefficients at different correlations between Legendre 
polynomial coefficients aaccounts for any correlation between different 
channels and energy levels.



Energy distributions (MF35)
===========================


Best estimates
--------------

Continuous function dependent on the incident neutron and outgoing particle 
energy.
Normalized (not always to 1)
Tabulated values.

Covariances
-----------
Only for tabulated values
Few incident-energies because of strong correlations

No correlation between between separate covariances, which are sampled independently.

Add points along the energy grid to capture all the covariance terms

The ENDF format manual notes these matrices are probability distributions
that must remain normalized to unity and therefore the elements in these sym-
metric matrices are constrained such that the sum of the elements in any row
(or column) must be zero. This is sometimes referred to as the `zero sum` rule.

renormalize to comply with the `zero sum` rule.
Renormalize to 1, although original might not be by 1-2%.
The renormalization could change a lot the effect of the covariance matrix if the original covariance matrix did not comply with the zero-sum rule.

Sampling options
~~~~~~~~~~~~~~~~

MF : 35

One could chose the MT, bu we only found covariance for MT18
MT : 18

Decay data and fission yields (MF8)
===================================
Best estimates stored in MF8
Covariance are not allowed by the format, only uncertainty.

Covariances for branching ratios and fission yield in development.

This part is not yet included in the current release.
