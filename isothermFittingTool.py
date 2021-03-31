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
#  Fits input experimental adsorption isotherm data to either dual-site
#  Langmuir model
# 
#  Last modified:
#  - 2021-03-31, HA: Initial creation
# 
#  Input arguments:
# 
#  Output arguments:
# 
##############################################################################

def isothermFittingTool():
    from scipy.io import loadmat
    import numpy as np
    from scipy.optimize import basinhopping
    minimizer_kwargs = {"method":"powell"}
    x0 = [3,3,1e-5,1e-5,4e4,4e4]
    mybounds = MyBounds()
    ret = basinhopping(generateObjectiveFunction, x0, minimizer_kwargs=minimizer_kwargs, niter = 10,T = 2, disp = True, accept_test=mybounds)
    return ret
    print("global minimum: x = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f], f(x0) = %.4f" % (ret.x[0],ret.x[1],ret.fun))
                   

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
        err = err + np.power((qamount[kk] - qfun),2)
    f =(Nvals/2 * np.log(err))
    return f


class MyBounds(object):
    def __init__(self, xmax=[10,10,1,1,8e4,8e4], xmin=[0,0,0,0,0,0] ):
        import numpy as np
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)
    def __call__(self, **kwargs):
        import numpy as np
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
    
        return tmax and tmin