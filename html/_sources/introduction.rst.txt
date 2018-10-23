************
Introduction
************

Preface
=======
Nuclear data are widely considered as one of the largest sources of uncertainty 
for several kinds of reactor analyses, including rector criticality and safety.
As a consequence, already several years ago the nuclear data evaluators have 
taken up the initiative to provide nuclear data uncertainty estimates that 
describe our data knowledge in terms of best fit to experimental measurements.
Recently, nuclear models uncertainties were also included to provide for missing 
experimental information.
All these quantites were included as covariance matrices into public 
evaluations along with their corresponding nuclear data.
Also, the standard format ENDF-6 used to store nuclear data into computer 
readable files was adjusted to accept the new covariance information.

As a second step, algorithms to assess linear sensitivity measures to nuclear 
data were implemented for a few neutron transport codes and mostly for 
:math:`\text{k}_{eff}`.
The resulting sensitivies could be used to propagate nuclear data covariances 
making use of perturbation theory or generalized perturbation theory (GPT).

Several shortcomings were identified for this technology:

 - sensitivities can only be calculated with few specific codes;
 - perturbation theory and GPT are often cumbersome to implement for models 
   other than criticality;
 - calculations become too constraining for a large number of selected 
   responses;
 - sometimes a linear surrogate does not correctly represent the physical model.


Consequently, several research institutes started investigating a different 
approach to propagate nuclear data covariances based on nuclear data sampling.
In this sense, some codes were developed to produce *random* nuclear data 
files with perturbed data.
However, again several limitations were found, since at least one of the 
following drawbacks was experienced:

 - only few nuclear data types were considered (mostly only cross sections);
 - the process was specific to a single transport code, either by the nuclear 
   data format or by the direct implementation of the sampling algorithms 
   in the code routines;
 - the sampling algorithm did not apply to continuous-energy data;
 - the sampling did not use the information stored in the covariance sections 
   of the evaluated file.

With these restrictions in mind, we set up the objective of developing a new 
nuclear data sampling code that would not be bounded by such constraints.
As a results, the SANDY code was created.


What is SANDY?
==============
SANDY is a python package that can read, write and perform a set of operations 
on nuclear data files in ENDF-6 format.
Its primay objective, as it was originally conceived, is to 
produce perturbed nuclear data files containing sampled parameters that 
represent the information stored in the evaluated nuclear data covariances.
Such files can be ultimately used to propagate uncertainties through any given 
compatible system using a brute force technique.

To parse ENDF-6 files, collect nuclear data, build covariance matrices, 
draw samples from probability distribution and perform other nuclear data-related 
tasks, a high-level ``python`` interface was created.

Currently, SANDY can draw samples for:

 * cross sections;
 * angular distrbutions of outgoing particles;
 * energy distrbutions of outgoing particles;
 * fission neutron multiplicities.

SANDY strongly relies on the work of the nuclear data evaluators as it draws 
random nuclear data samples from the covariance matrices that are located in 
the ENDF-6 files.
However, the ENDF-6 files do not provide information on the nuclear data 
probability distribution functions (PDF).
As an assumption, SANDY samples from a multivariate Normal PDF.

The outputs of SANDY are sets of *perturbed* ENDF-6 format files that can be 
used for brute force uncertainty propagation. This methodology to propagate 
uncertainties can be applied to any model, response or computer code, as long 
as they are compatible with ENDF-6 format files, either directly or via data 
conversion.
Also, the perturbed files are not constrained by any implicit effect due 
to processing, such as multigroup assumptions, self shielding or temperature 
effects.