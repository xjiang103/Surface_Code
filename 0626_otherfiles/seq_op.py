from scipy.optimize import minimize
import subprocess
def print_callback(xs):
    print(xs)
def fun(delta):
    fid = float(subprocess.check_output(["./na_seq_op_sp",'-delta_in',str(delta[0]),'-deltat',str(delta[1]),'-dmpos',str(2),'-dmstdpos',str(1)]).split()[-1])
    return 1-fid
res = minimize(fun,[-0.5,0.2],method="nelder-mead",callback=print_callback)
#For just delta_in
#res = minimize(fun,[-0.8635],method="cobyla",callback=print_callback)
print("Final Fidelity: ",str(1-res.fun))
