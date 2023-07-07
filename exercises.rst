*********
Exercises
*********

Exercise 1
==========
Generate 200 perturbed ACE files from the JEFF-3.3 evaluation for U-235 in file ``u235.jeff33``.
Use 20 CPUs and post-process each file at a temperature of 300 K. 
Perturb both cross sections and fission yields.


.. code::

	python -m sandy.sampling  u235.jeff33  --samples 200 --processes 20 --temperatures 300  --acer
	

Exercise 2
==========
Generate 5 perturbed copies of the ENDF-6 file ``pu239.endf80``.
Use 1 CPU and do not produce the ACE files, but only perturbed ENDF-6 and PENDF files.
Perturb both cross sections and fission yields.

.. code::

	python -m sandy.sampling  pu239.endfb80  --samples 5


More info
=========
	For more flexibility on the generation of random file, see this example_.

.. _example: https://github.com/luca-fiorito-11/sandy/blob/v1.0/notebooks/notebook_random_files.ipynb
