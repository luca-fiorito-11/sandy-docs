******
Theory
******

Sampling methodology
====================
To describe the sampling methodology implemented in SANDY we introduce some 
assumptions with the purpose of simplifying the explanation.
Let us assume that a given ENDF-6 file contains only one covariance matrix 
:math:`\Sigma` that:

 * is defined over :math:`M` energy groups;
 * is in relative units;
 * is representative for a continuous energy data type, say, cross section.

No loss of generality follows by imposing this assumptions as it happens that 
they apply to almost all the covariances allowed by the ENDF-6 format.

SANDY parses the ENDF-6 file, it extracts the covariance matrix and constructs 
a multivariate Normal distribution :math:`N(0,\Sigma)`.
To sample from a multivariate Normal distribution, first we must 
draw a matrix :math:`X` of :math:`m` independent and identically :math:`N(0,1)` 
distributed variables :math:`x=[x_1, x_2, \dots, x_m]`.

.. math::
   X=\begin{bmatrix}
    x_1^{(1)} & x_1^{(2)} & ... & x_1^{(n)} \\
    x_2^{(1)} & x_2^{(2)} & ... & x_2^{(n)} \\
    \vdots    & \vdots    & \ddots & \vdots \\
    x_m^{(1)} & x_m^{(2)} & ... & x_m^{(n)} \\
   \end{bmatrix}\,,
   :label: x-samples

where :math:`n` is the number of samples.

Then, we define a linear operator :math:`L` that, when applied to 
uncorrelated standard samples :math:`X`, converts them into 
:math:`N(0,\Sigma)`-distributed samples :math:`K`.

.. math::
  K = L X
  \label{eq:corr_eq}

The operator :math:`L` does not affect the mean of the distribution, since 

.. math::
   E \left[ K \right] = E \left[ L X \right] = L E \left[ X \right] = 0 \,.
   :label: preserve-mean

Also, the linearity of the operator preserves the shape of the distribution.

Then, the covariance becomes

.. math::
   E \left[ K K^T \right] =
   E \left[ \left( L X \right) \left( L X \right)^T \right] =
   L E \left[ X X^T \right] L^T =
   L L^T \,,
   :label: preserve-cov

since by definition :math:`E \left[ X X ^T \right]` is the identity matrix.

The problem of sampling Normally distributed variables with covariance matrix 
:math:`\Sigma` reduces to finding :math:`L` such that

.. math::
   L L^T = \Sigma \,.
   :label: chol_condition

SANDY calculates :math:`L` by performing an eigendecomposition of 
:math:`\Sigma`.

Eventually, :math:`K` can be converted into a :math:`N(1,\Sigma)`-distributed  
matrix of perturbation coefficients :math:`P` by shifting the distribution mean 
to the unit vector :math:`1`   

.. math::
   P =
   \begin{bmatrix}
   p_1^{(1)} & p_1^{(2)} & ... & p_1^{(n)} \\
   p_2^{(1)} & p_2^{(2)} & ... & p_2^{(n)} \\
   \vdots    & \vdots    & \ddots & \vdots \\
   p_m^{(1)} & p_m^{(2)} & ... & p_m^{(n)} \\
   \end{bmatrix} = 
   K + 1 \,.
   :label: p-samples


Negative eigenvalues
~~~~~~~~~~~~~~~~~~~~
By definition, covariance matrices are positive-definite, that is, their 
eigenvalues are all positive.
These requirements must also be satisfied to apply the decomposition in Eq.
However, it is not rare to find covariance matrices for which it is not the case.

SANDY cannot investigate into the covariance evaluation process to 
identify the sources of this issue.
Therefore, it reconstructs an approximate covariance matrix 
:math:`\widetilde{\Sigma}` with all negative eigenvalues set to zero.
This assumption is acceptable for small negative eigenvalues that are likely 
to sprout from round-off or truncation processes.
For such cases it follows that :math:`\widetilde{\Sigma} \approx \Sigma`.

To check whether the covariance eigenvalues are well represented by the 
perturbations, option ``--eig`` was implemented in SANDY.

.. hint:: to compare the first 20 eigenvalues, type

	.. code:: bash

		sandy  <endf6_file>  --samples 100  --eig 20

	By default, SANDY displays the first 10 eigenvalues.

	
		
How to apply multigroup perturbations to continuous-energy data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SANDY perturbation coefficients :math:`P` reflect the multigroup energy structure 
of the covariance matrix used for sampling.
On the contrary, the tabulated data in the ENDF-6 files are not defined for a energy 
group structure, but they rather apply to a continuous-energy domain using a number of 
explicitely given energy-value pairs :math:`(e_k,v_k)` and interpolation laws.

For any given coefficient :math:`p_i^{(j)}` defined over an energy group 
:math:`[e_i,e_{i+1}]`, SANDY perturbs all energy-value pairs 
:math:`(e_k,v_k)` for which :math:`e_i \leq e_k \leq e_{i+1}` using the 
following formula

.. math::
   v_k^{(j)} = v_k p_i^{(j)} \,.

That is, if a cross section must be perturbed by 10% between 1 and 10 eV,
then all the values of the energy-value pairs in the corresponding ``MF3`` 
section are multiplied by :math:`1.1`.

.. important::
	This procedure implies that all cross section points in the energy interval 
	of interest are 100% correlated.

To make sure that the actual covariance structure is represented in the 
perturbed files, SANDY adds additional energy-value pairs to the tabulated data.
This is particularly important when the covariance energy structure is finer 
than the energy-value pair density, as it is often found in ``MF4`` and ``MF5`` 
for energies below 1 KeV. 


How to handle negative samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When sampling from a Normal distribution, and in particular when the standard 
deviations are large, it is likely to draw perturbation :math:`<=0` that, when 
applied to the evaluated data they will make them change sign.
Physically, many quantities such as cross sections or energy ditributions are 
intrinsically positive.
As a consequence, SANDY proposes two methods to handle negative perurbations:

* method 1:
	perturbation coefficients outside the range :math:`[0,2]` are set to 1;
* method 2: 
	perturbation coefficients outside the range :math:`[0,2]` are set either to 0 or 2 if they fall respectively below or above the defined range.
				
A comparison of the two methods to represent different level of uncertainty is 
reported in figure.

Given the strong nuclear data uncertainty reduction using the 1st method 
for standard deviations larger than 40%, it was decided to implement the 2nd 
method in the code.



How to perturb cross sections
=============================
To correctly perturb cross sections, SANDY needs the tabulated data in ``MF3`` to 
be reconstructed from the resonance parameters in ``MF2`` and to be linearized, 
so that for any energy point a cross section value can be retrieved via 
linear interpolation.
ENDF-6 files that satisfy these conditions are called PENDF, or pointwise-ENDF, 
files and are recognizible by a special flag in section ``MF1/MT451``.
To produce such files it is common to utilize processing codes such as 
NJOY_ (module RECONR) or PREPRO_ (modules RECENT and LINEAR).
Notice that SANDY does not perturb cross sections for files that are not in 
PENDF format.

.. _NJOY: http://www.njoy21.io/NJOY2016/

.. _PREPRO: https://www-nds.iaea.org/public/endf/prepro/


Covariances for resonance parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To perturb cross sections, one can choose to extract perturbations from 
the ``MF33`` section in the original ENDF-6 file.

.. code:: bash

	sandy  <pendf_file>  --cov <endf6_file>  --samples 100  --mf 33

By doing so, the user must be aware that the covariance information for the 
resonance region will not be included, if the dedicated section ``MF32`` is 
provided.

If that is the case, another choice to properly include both ``MF32`` and 
``MF33`` contributions is to:

 1. process the ENDF-6 file with the NJOY module ERRORR for a given 
    multigroup structure;
 2. extract perturbations from the processed covariance matrix written 
    in the ERRORR output file.

.. code:: bash

	sandy  <pendf_file>  --cov <errorr_file>  --samples 100  --mf 33

The ERRORR output contains a derived ``MF33`` sections that includes 
both the cross sections and resonance parameters covariances.

.. Note:: This is the preferred methodology used by the authors to produce perturbed files
	  with SANDY.
	
A template NJOY input file to produce a ERRORR output file is reported 
below.


Redundant cross sections
~~~~~~~~~~~~~~~~~~~~~~~~
The information stored in the various ``MT`` sections is often redundant [1]_, 
in the sense that some cross sections can be derived from others via summation rules.
For example, the total cross section can be reconstructed as the sum of all the 
partial (non-redundant) cross sections, or the total inelastic cross section can 
be reconstructed as the sum of all the discrete and continuous inelastic levels.

.. [1] A cross section is defined as redundant when it can be completely calculated 
   from the combination of other cross sections, according to the conservation 
   laws.

SANDY automatically reconstruct redundant cross sections for each perturbed file.

Also, in the case pertrubations exist only for a redundant cross section, say 
``MT4``, and not for its components, say from ``MT51`` to ``MT91``.
Then, SANDY applies the ``MT4`` perturbations also to ``MT51``, ``MT52``, and so 
on, to make sure that the perturbation is taken into account and that the sum of 
the components equals the redundant cross section.

.. important:: by imposing conservation laws such as summation rules, extra 
               levels of correlations are forced upon the random samples, 
               which might not be present in the original covariances.




How to perturb fission neutron multiplicities
=============================================
Covariances for average fission neutron multiplicities are given in ``MF31`` for:

 - ``MT452`` total fission neutrons;
 - ``MT455`` delayed fission neutrons;
 - ``MT456`` prompt fission neutrons.

To produce perturbed files where only the fission neutron mulitplicities are 
varied, type

.. code:: bash

	sandy  <endf6_file>  --samples 100  --mf 31

As for cross sections, the redundant total fission nubar is reconstructed from 
the prompt and delayed fission nubar.



How to perturb angular distributions
====================================
The ENDF-6 formats allows angular distributions to be stored in section 
``MF4`` in two forms:

 - by tabulating the normalized probability distribution as a function of incident energy;
 - by tabulating the Legendre polynomial expansion coefficientsas a function of incident neutron energy.

However, the corresponing covariances in ``MF34`` only accept the second form.

To run SANDY to perturb only the angular distributions, type

.. code:: bash

	sandy  <endf6_file>  --samples 100  --mf 34

If no option is specified, SANDY will perturb all Legendre polynomial 
coefficients up to any order [2]_, as long as covariance data are available.
However, to consider only covariances for Legendre polynomial coefficients 
up to a given order, say 2, one must type

.. [2] all correlations between different Legendre polynomial coefficients are 
	   also taken into account.

.. code:: bash

	sandy  <endf6_file>  --samples 100  --mf 34  --max-polynomial 2


	
How to perturb energy distributions
===================================

Covariances for outgoing energy distributions in ``MF35`` are mostly given for 
prompt fission neutron spectra (PFNS) and for few incident energy ranges, 
therefore assuming large correlations for distributions associated to incident 
energies in the same range.
In addition, the format does not allow correlations between spectra that do not 
belong to the same range.

To perturb only energy distributions, type

.. code:: bash

	sandy  <endf6_file>  --samples 100  --mf 35

SANDY draws samples from each covariance matrix independently.
Then, the perturbations are applied to each evaluated energy distribution only 
if the incident neutron energy belongs to the covariance range.

Normalization
~~~~~~~~~~~~~
Being probability distributions, all perturbed PFNS must be normalized to unity.
If the normalization was already included in the covariance matrix, the sums of 
the elements in any row (and in any column) would be equal to zero [3]_.

To make up for covariances that do not comply with this rule, SANDY normalizes 
all perturbed energy distributions.

.. [3] This constraint is also called the *zero-sum* rule.

.. important:: by imposing conservation laws such as data normalization, extra
	levels of correlations are forced upon the random samples, 
	which might not be present in the original covariances.
