# Installation

This package requires [xraylib][xraylib] and [compwizard][compwizard]. On Windows, xraylib can be installed through the Anaconda interface.
<br><br>
`conda install -c conda-forge xraylib=4.0.0`<br>
`pip instal compwizard`
<br><br>
For further information on how to install xraylib on other operational systems, check xraylib [wiki][xlibwiki].
<br>
This module can be installed with:
<br><br>
`pip install xfit`
<br><br>

# Usage
This module provides the "Spectrum" class. It is possible to initialize this class with a numpy nd.array or loading an *.mca or *.spt file, by giving the path.
For the continuum estimation algorithm, refer to the [Continuum.py][continuum] file.
<br>
Test.py file provides a common *.mca file to use with the example.py script provided.

## Example 1

```
import xfit
import numpy as np
import matplotlib.pyplot as plt
path = r"./test.mca"
pool_file = r"./pool.txt"
Spec = xfit.Spectrum(file_path=path)
Spec.calibrate() #if no arguments are passed, it gets the parameters from the mca or spt header
Spec.estimate_continuum(30, 11, 11, 3) #iterations, filter window, sav-gol window, sav-gol order
Spec.fit_fano_and_noise()
Spec.create_pool(pool_file)
Spec.fit()

#Plot ------
fig, ax = plt.subplots()
ax.plot(Spec.energyaxis, Spec.data, color="black", label="Data")
ax.plot(Spec.energyaxis, Spec.continuum, color="green", label="Continuum")
for element in Spec.areas.keys():
     ax.plot(Spec.energyaxis, 
            Spec.plots[element],
            label=element+" fit result", 
            color=ElementColors[element],
            linestyle="--")
ax.legend(loc=1, fancybox=1)
ax.set_yscale("log")
plt.show()
```
<br>

0.11400000005960464 80.00951851146041 (Fano and Noise values found, respectively)<br>
![Output](https://github.com/linssab/xfit/tree/master/xfit/Figure_1.png)

## Example 2

```
import xfit
import numpy as np
ydata = np.arange(1024)
fit_pool = {}
fit_pool["elements"] = {}
fit_pool["elements"]["Cu"] = ["KA1","KA2","KB1","KB3"]
fit_pool["bg"] = 1 #Forces the use of continuum estimation for the fit
Spec = xfit.Spectrum(array=ydata)
Spec.calibrate(x=channels, y=energies)
Spec.estimate_continuum(30, 11, 11, 3)
Spec.fit_fano_and_noise()
Spec.pool = fit_pool
Spec.fit()
#or simply: Spec.fit(pool=fit_pool)

#Plot ------
fig, ax = plt.subplots()
ax.plot(Spec.energyaxis, Spec.data, color="black", label="Data")
ax.plot(Spec.energyaxis, Spec.continuum, color="green", label="Continuum")
for element in Spec.areas.keys():
     ax.plot(Spec.energyaxis, 
            Spec.plots[element],
            label=element+" fit result", 
            color=ElementColors[element],
            linestyle="--")
ax.legend(loc=1, fancybox=1)
ax.set_yscale("log")
plt.show()
```
<br>

[xraylib]: http://lvserver.ugent.be/xraylib
[xlibwiki]: https://github.com/tschoonj/xraylib/wiki/Installation-instructions
[compwizard]: https://pypi.org/project/compwizard/#description
[continuum]: https://github.com/linssab/xfit/blob/master/xfit/Continuum.py
