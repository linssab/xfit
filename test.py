
if __name__.endswith("__main__"):
    import sys, os
    import xfit
    import matplotlib.pyplot as plt
    from Elements import *
    
    chan = [597,941,1875]
    kev = [8.04, 12.641, 25.3]
    data = np.arange(1024)

    path = sys.argv[1]
    Spec = xfit.Spectrum(file_path=path)
    #Spec = xfit.Spectrum(array=data)

    #Spec.calibrate(x=chan,y=kev)
    Spec.calibrate()

    Spec.estimate_continuum(30, 11, 11, 3)
    Spec.fit_fano_and_noise()

    fit_pool = {}
    fit_pool["elements"] = {}
    fit_pool["elements"]["Cu"] = ["KA1","KA2","KB1","KB3"]
    fit_pool["bg"] = 1 #Forces the use of continuum estimation for the fit
    Spec.fit(pool=fit_pool)

    #Spec.create_pool(sys.argv[2])
    #Spec.fit()

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
    sys.exit(0)
