*****************
The ENDF-6 format
*****************

To better understand how SANDY works it is important to have some knowledge 
of the standard format ENDF-6 used to store evaluated data into computer readable 
files.
The full ENDF-6 documentation can be found at 
https://www.oecd-nea.org/dbdata/data/manual-endf/endf102.pdf
.

The ENDF-6 format is a system developed for the storage and retrieval 
of evaluated nuclear data to be used for applications of nuclear technology.
Evaluations are processed from the combination of experimentally measured 
physical parameters and the predictions of nuclear model calculations 
in the attempt to extract the true values of such parameters.

Since they were constructed for nuclear data processing programs, the 
ENDF-6 files store collected evaluated data in computer readable 
format following constraining formatting rules, which render the data 
cumbersome to be processed without dedicated tools.

The ENDF-6 format provides representations for neutron-induced cross sections
and distributions, neutron, photon and charged-particle production data from 
neutron reactions, photo-atomic interaction data, thermal neutron scattering 
data, radionuclide production, decay data and fission products yields.
The rules to decrypt and process the ENDF-6 format are encoded in 
SANDY, which can read most of the sections of any ENDF-6 file for 
neutron-induced data.

An overview of the ENDF-6 structure is reported in the following sections.

The material number (``MAT``)
=============================

A ``ENDF-6`` file is divided into material sections, each defined by a unique material 
number ``MAT``.
A single file can contain nuclear data evaluations for one or many 
``MAT`` sections, which refer to different materials.
Generally, a material correspond to a single nuclide, a natural 
element containing several isotopes, or a mixture of several elements 
such as compounds, alloys or molecules.

.. Note::
    Most of the recent ENDF-6 files for neutron-induced data contain only one 
    ``MAT`` number.


The data type number (``MF``)
=============================

For each material, the evaluated nuclear data are provided in sections, 
specified by the section number ``MF``.
Amongst the several ``MF`` numbers, SANDY can process the following:

 * :``MF1``: : general information and fission multiplicities;
 * :``MF3``: : neutron cross sections;
 * :``MF4``: : angular distributions of secondary particles;
 * :``MF5``: : energy distributions of secondary particles;

The ``ENDF-6`` format also allows dedicated ``MF`` sections for the 
storage of evaluated nuclear data uncertainties and covariance matrices.
These evaluated data reflect the information coming from the measurements --- 
e.g. systematic errors, machine resolution --- and are fine-tuned from the 
inference of practical reactor applications.
The covariance ``MF`` sections processable by SANDY are

 * :``MF31``: : covariances for the average fission neutron multiplicities;
 * :``MF33``: : covariances for neutron cross sections;
 * :``MF34``: : covariances for angular distributions of secondary particles;
 * :``MF35``: : covariances for energy distributions of secondary particles.

.. Note::
   Each ``MF`` covariance section contains covariances for the data given in 
   section ``(MF - 30)``, e.g. ``MF33`` contains covariances for section ``MF3``.


The reaction number (``MT``)
=============================

Each ``MF`` section is divided in further subsections identified by the 
reaction number ``MT`` that uniquely defines a reaction type.
A list of the most common ``MT`` numbers for incident neutron reactions follows:

 * :``MT1``: : total cross section;
 * :``MT2``: : elastic scattering cross section;
 * :``MT4``: : inelastic scattering cross section;
 * :``MT18``: : fission cross section;
 * :``MT102``: : radiative capture crross section;
 * :``MT452``: : total average fission neutron multiplicity;
 * :``MT455``: : delayed average fission neutron multiplicity;
 * :``MT456``: : prompt average fission neutron multiplicity.