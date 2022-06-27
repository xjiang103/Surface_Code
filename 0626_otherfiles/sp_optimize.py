from scipy.optimize import minimize
import subprocess
def print_callback(xs):
    print(xs)
def fun(delta):
    #For just delta_in
    #fid = float(subprocess.check_output(["./na_11_sp_1",'-delta_in',str(delta[0]),'-deltat',str(0.2165)]).split()[-1])
    #For deltat
    #fid = float(subprocess.check_output(["./na_11_sp",'-delta_in',str(-0.3),'-deltat',str(delta[0])]).split()[-1])
    #For both delta_in and deltat
    fid = float(subprocess.check_output(["./na_11_sp_1",'-delta_in',str(delta[0]),'-deltat',str(delta[1])]).split()[-1])
    return 1-fid
#For both delta_in and deltat
#res = minimize(fun,[-0.5,0.2],method="nelder-mead",callback=print_callback)
res = minimize(fun,[-0.8635,0.2165],method="cobyla",callback=print_callback)
#For just delta_in
#res = minimize(fun,[-0.8635],method="cobyla",callback=print_callback)
print("Final Fidelity: ",str(1-res.fun))
