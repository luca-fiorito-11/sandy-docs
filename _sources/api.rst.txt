*********************************
Ptyhon API for nuclear data files
*********************************

SANDY includes a Python API that allows parsing formatted files, analyzing and performin a series 
of operations on the different nuclear data types available.
Cross sections, neutron multiplicities, angular and energy ditributions, fission yields and the 
respective covariance data can be stored in user-friedly ``pandas.DataFrame``.
Each dataframe is customized to provide additional functions that are inherent to each 
specific data type.

The Python API documentation is available `here <https://luca-fiorito-11.github.io/sandy>`_.


Examples
========

Below is reported a list of examples where SANDY and its Python API are used to perform different 
operations:

- `extract and plot eigenvalues from cross sections/nubar covariance matrix <https://gist.github.com/ad352e80b01f81fe12599c079f7cb9d7.git>`_
