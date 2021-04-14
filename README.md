#Installation
This package requires [xraylib][xraylib]. On Windows, it can be installed through the Anaconda interface with <br> 
`conda install -c conda-forge xraylib=4.0.0` <br>
For further information on how to install xraylib on other operational systems, check xraylib [wiki][xlibwiki].
<br>
This module can be installed with:<br>
`pip install xrffitting`<br>

#Usage
This module provides the "Spectrum" class. It is possible to initialize this class with a numpy nd.array or loading an *.mca or *.spt file, by giving the path. <br>
##Example1
<br> 
```
import xfit
import numpy as np
import matplotlib.pyplot as plt
path = r"./test.mca"
pool_file = r"./pool.txt"
Spec = xfit.Spectrum(file_path=path)
Spec.calibrate("from_source")
Spec.fit_fano_and_noise()
Spec.create_pool(pool_file)
Spec.fit()
plt.plot(Spec.energyaxis, Spec.data)
for element in Spec.areas.keys():
	plt.plot(Spec.energyaxis, 
		Spec.plots[element],
		label=element, 
		color=ElementColors[element],
		linestyle="--")
plt.show() 
```
<br>
##Example2
<br>
```
import xfit
import numpy as np
ydata = np.arange(1024)
fit_pool = {}
fit_pool["elements"]["Cu"] = ["KA1","KA2","KB1","KB3"]
fit_pool["bg"] = 1 #Uses continuum estimation for the fit
Spec = xfit.Spectrum(array=path)
Spec.calibrate(x=channels, y=energies)
Spec.fit_fano_and_noise()
Spec.pool = fit_pool
Spec.fit()
plt.plot(Spec.energyaxis, Spec.data)
for element in Spec.areas.keys():
	plt.plot(Spec.energyaxis, 
		Spec.plots[element],
		label=element, 
		color=ElementColors[element],
		linestyle="--")
plt.show() 
```
<br>
[xraylib]: http://lvserver.ugent.be/xraylib
[xlibwiki]: https://github.com/tschoonj/xraylib/wiki/Installation-instructions
