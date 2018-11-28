.. role:: bash(code)
   :language: bash

***************
Getting started
***************

These instructions will get you a copy of the project up and running on your 
local machine for development and testing purposes.

.. important::

    SANDY works with files in ENDF-6 format.
    If you are not familiar with the format, have a look at the ENDF-6_ documentation.

.. _ENDF-6: https://www.oecd-nea.org/dbdata/data/manual-endf/endf102.pdf



Prerequisites
=============

SANDY is developed in ``python3`` and does not support ``python2``.
In order to run SANDY, make sure that you have a version of ``python3`` installed.
In file ``requirements.txt`` you can find the python dependencies required to ensure 
the correct functioning of SANDY.

To read and write ENDF-6 formatted files, SANDY uses some fortran routines.
Make sure that you have a valid fortran compiler installed.

SANDY has been developed and tested on a Ubuntu machine.


Installation
============

To download SANDY, move to the folder where you want the source code and type

```git
git clone https://github.com/luca-fiorito-11/sandy.git
```

To install SANDY, run the following commands

.. code:: bash

    cd sandy
    python setup.py install

Make sure that :bash:`python` points to a correct ``python3`` version for which 
you have administration rights.



Running the tests
=================

Once the installation is completed, run ``pytest`` to automatically start SANDY's tests

.. code:: bash

    pytest

.. _pytest: https://docs.pytest.org/en/latest/

More pytest_ command line options can be added to customize the tests set.

More than 50 tests have been put in place in SANDY, and cover the whole range 
of data types that the code can process.
Their purpose is to check the correct functioning of each basic SANDY unit, as 
well as to ensure that the implemented physics is correct.

For a list of the existing SANDY test, type

.. code:: python

	pytest  --collect-only

Each name has a self-explanatory name for the function it performs.



Usage
=====

For an overview of SANDY's usage type

.. code:: bash

    sandy  --help

The code will print out the following message:

.. code::

    usage: python -m sandy.sampling [-h] [--cov COV] [--samples SAMPLES]
                                    [--outdir DIR] [--processes PROCESSES]
                                    [--max-polynomial MAX_POLYNOMIAL] [--eig N]
                                    [--mat {1,..,9999} [{1,..,9999} ...]]
                                    [--mf {31,33,34,35} [{31,33,34,35} ...]]
                                    [--mt {1,..,999} [{1,..,999} ...]]
                                    [--outname OUTNAME] [--debug]
                                    file

    Run sampling

    positional arguments:
      file                  ENDF-6 or PENDF format file

    optional arguments:
      -h, --help            show this help message and exit
      --cov COV, -C COV     file containing covariances
      --samples SAMPLES, -S SAMPLES
                            number of samples
                            (default = 200)
      --outdir DIR, -D DIR  target directory where outputs are stored
                            (default = current working directory)
                            if it does not exist it will be created
      --processes PROCESSES, -N PROCESSES
                            number of worker processes
                            (default = 1)
      --max-polynomial MAX_POLYNOMIAL, -P MAX_POLYNOMIAL
                            Maximum order of Legendre polynomial coefficients 
                            considered for sampling
                            (default = all)
      --eig N               print the first N eigenvalues of the evaluated 
                            covariance matrices
                            (default = do not print)
      --mat {1,..,9999} [{1,..,9999} ...]
                            draw samples only from the selected MAT sections 
                            (default = keep all)
      --mf {31,33,34,35} [{31,33,34,35} ...]
                            draw samples only from the selected MF sections 
                            (default = keep all)
      --mt {1,..,999} [{1,..,999} ...]
                            draw samples only from the selected MT sections 
                            (default = keep all)
      --outname OUTNAME, -O OUTNAME
                            basename for the output files 
                            (default is the the basename of <file>.)
      --debug               turn on debug mode



Examples
========


Data and covariances are in the same file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Produce 1000 perturbed copies of a ENDF-6 file :bash:`<tape>` that contains both 
evaluated data and covariances.

.. code:: bash

    sandy  <tape>  --samples 1000


Below are reported the ENDF-6 data sections that will be perturbed and the 
respective covariance sections.

+------------------------+--------------+--------------------+
| Data type              | Data section | Covariance section |
+========================+==============+====================+
| fission multiplicities | MF1          | MF31               |
+------------------------+--------------+--------------------+
| cross sections         | MF3          | MF33               |
+------------------------+--------------+--------------------+
| angular ditributions   | MF4          | MF34               |
+------------------------+--------------+--------------------+
| energy distributions   | MF5          | MF35               |
+------------------------+--------------+--------------------+

.. important:: cross sections will be perturbed **only** if they are linearized 
               and given in PENDF (pointwise-ENDF) format.
               To convert a ENDF-6 file into PENDF format, you can use nuclear data processing 
               codes such as NJOY_ or PREPRO_.

.. _NJOY: http://www.njoy21.io/NJOY2016/

.. _PREPRO: https://www-nds.iaea.org/public/endf/prepro/


Perturb only one or few data types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Add keyword option :bash:`--mf` to perturb only few data types.
For example, to produce 1000 perturbed copies of a file :bash:`<tape>` where 
only angular and energy distributions are perturbed, type

.. code:: bash

    sandy  <tape>  --samples 1000  --mf 34 35



Data and covariances are in different files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Produce 1000 perturbed copies of a file :bash:`<tape>` that contains evaluated 
data using covariances from file :bash:`<covtape>`.

.. code:: bash

    sandy  <tape>  --cov <covtape>  --samples 1000

.. important:: this command is often used for perturbing cross sections, where 
               the linearized data are in a PENDF file :bash:`<tape>` that might 
               not contain covariances and the covariance data are in the original 
               ENDF-6 file :bash:`<covtape>`.


Covariance data in ERRORR format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ERRORR is a NJOY_ module that processes the covariance information present in a 
ENDF-6 file into a given multigroup structure.
The resulting tabulated covariance is tabulated into an output file with a specific 
ERRORR format.
Not only does ERRORR process cross section covariances in ``MF33``, but it can 
also handle the resonance-resonance covariances in ENDF-6 covariance section 
``MF32``.

To produce 1000 perturbed copies of a PENDF file :bash:`<pendf_tape>` including the 
``MF32`` 
covariances for resonance parameters, type

.. code:: bash

    sandy  <pendf_tape>  --cov <errorr_tape>  --samples 1000

where :bash:`<errorr_tape>` is a ERRORR output file.
