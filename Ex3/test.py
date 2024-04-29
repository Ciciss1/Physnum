import numpy as np
from scipy.optimize import root


d = 149598023e3
G = 6.674e-11
M_S = 1.98892e30
M_T = 5.9736e24
alpha = M_T / (M_S + M_T)
beta = M_S / (M_S + M_T)
X_S = -alpha * d
X_T = beta * d

a = (X_T - X_S)

#angular speed of the Earth around sun
Omega = np.sqrt(G * M_S / (X_T*d**2))


print(f"Omega nous = {Omega:.20e}")

def function(x) : 
    return (Omega**2)*(x**5) - 2*(Omega**2)*(x**4)*(X_S+X_T) + (Omega**2)*(x**3)*((X_S + X_T)**2 + 2*X_S*X_T) - (x**2)*(G*(M_S+ M_T) + 2*(Omega**2)*X_T*X_S*(X_S+X_T)) + x*(2*G*(M_S*X_T + M_T*X_S) + (Omega**2)*(X_T**2)*(X_S**2)) - G*(M_S*X_T**2 + M_T*X_S**2)


#L2 is the root of the function
#set the tolerance to 1e-15
L2 = root(function, 1.5e9, tol=1e-20)
L2x = L2.x[0]
print(f"L2 nous : {L2.x[0]:.20e}")


x=1.51099098533097e11

print("L2 prof : ",x)
def function(Omega) : 
    return (Omega**2)*(x**5) - 2*(Omega**2)*(x**4)*(X_S+X_T) + (Omega**2)*(x**3)*((X_S + X_T)**2 + 2*X_S*X_T) - (x**2)*(G*(M_S+ M_T) + 2*(Omega**2)*X_T*X_S*(X_S+X_T)) + x*(2*G*(M_S*X_T + M_T*X_S) + (Omega**2)*(X_T**2)*(X_S**2)) - G*(M_S*X_T**2 + M_T*X_S**2)

#L2 is the root of the function
#set the tolerance to 1e-15
Om = root(function, 1.5e9, tol=1e-20)
print(f"Omega prof = {Om.x[0]:.20e}")

Omega_prof = Om.x[0]



def norm(x, y) : 
    return np.sqrt(x**2 + y**2)

Y_S = 0
Y_T = 0
x = np.array([1.51099098533097e11, 0, 0, 0])
r_S = norm(x[0]-X_S, x[1]-Y_S)
r_T = norm(x[0]-X_T, x[1]-Y_T)

def function(Omega) : 
    return -G*M_S*(x[0]-X_S)/pow(r_S,3) - G*M_T*(x[0]-X_T)/pow(r_T,3) + pow(Omega, 2)*x[0] + 2*Omega*x[3]
    
Om = root(function, 1.5e9, tol=1e-20)
print(f"Omega opti (avec L2 prof) = {Om.x[0]:.20e}")



Omega = Omega_prof
def function(xx) : 
    return -G*M_S*(xx-X_S)/pow(r_S,3) - G*M_T*(xx-X_T)/pow(r_T,3) + pow(Omega, 2)*xx + 2*Omega*x[3] 
    
L2 = root(function, 1.5e9, tol=1e-20)
print(f"L2 opti (avec omega prof) = {L2.x[0]:.20e}")




