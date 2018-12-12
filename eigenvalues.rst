************************************************************
Extract the eigenvalues of a cross section covariance matrix
************************************************************

Generate 100 perturbed copies of the PENDF file ``sandy/data/U5/u235.pendf`` 
using the covariance data stored in the ERRORR file ``sandy/data/U5/u235.errorr``.
Perturb only the fission cross section.
Produce the files in output directory ``U235_outputs``.
Print the first 20 eigenvalues.

.. code:: python

	from sandy import XsCov, Endf6
	tape = Endf6.from_file("<FILE>").filter_by(mflist=[33])
        cov = XsCov.from_endf6(tape)
	eigs = cov.to_matrix().eig()[0]
	print(eigs.to_screen())

	import seaborn as sns
	import matpltolib.pyplot as plt
	
	sns.set_style("whitegrid")
	plt.plot(eigs.index, eigs.values, 'x')
	plt.gca().set_xlabel("#")
	plt.gca().set_ylabel("eigenvalues")
	plt.show()
