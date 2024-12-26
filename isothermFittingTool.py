##############################################################################
# 
#  Imperial College London, United Kingdom
#  Multifunctional Nanomaterials Laboratory / Complex Porous Media Laboratory
# 
#  Project:  PhD
#  Year:     2021
#  MATLAB:   R2020a
#  Authors:  Hassan Azzan (HA)
# 
#  Purpose:
#  Fits input experimental adsorption isotherm data to dual-site
#  Langmuir model using MLE using differential evolution method
# 
#  Last modified:
#  - 2021-04-07, HA: Change minimizer method and add options and bounds   
#  - 2021-03-31, HA: Initial creation
# 
#  Input arguments:
# 
#  Output arguments:
# 
##############################################################################

def isothermFittingTool():
    # from scipy.optimize import basinhopping
    from scipy.optimize import differential_evolution
    import numpy as np
    optBounds = np.array([[0, 10],
                         [0, 10],
                         [0, 1],
                         [0, 1],
                         [0, 8e4],
                         [0, 8e4]])
    
    # minimizer_kwargs = {"method":"SLSQP","bounds":optBounds, "tol":1e-6}
    # x0 = [3,3,1e-5,1e-5,4e4,4e4]    
    # ret = basinhopping(generateObjectiveFunction, x0, 
    #                    minimizer_kwargs=minimizer_kwargs, 
    #                    disp = True,niter=10000, T = 10)
    ret = differential_evolution(generateObjectiveFunction, optBounds,
                                 disp=True)
    return ret

def generateObjectiveFunction(x):
    from scipy.io import loadmat
    import numpy as np
    matData = loadmat('zif8Data')
    rawData = matData["zif8Data"]
    pressure = np.array(rawData[:,0])
    qamount = np.array(rawData[:,1])
    temperature = np.array(rawData[:,2])
    err = 0
    Nvals = pressure.size;
    for kk in range(Nvals):
        qfun = x[0]*(x[2]*pressure[kk]*np.exp(x[4]/(8.314*temperature[kk])))/(1+(x[2]*pressure[kk]*np.exp(x[4]/(8.314*temperature[kk])))) + x[1]*(x[3]*pressure[kk]*np.exp(x[5]/(8.314*temperature[kk])))/(1+(x[3]*pressure[kk]*np.exp(x[5]/(8.314*temperature[kk]))))
        err += np.power((qamount[kk] - qfun),2)
    f =(Nvals/2 * np.log(err))
    return f