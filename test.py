import numpy as np
from scipy.optimize import root


d = 149.598023e9
G = 6.67430e-11
M_S = 1.989e30
M_T = 5.972e24
alpha = M_T / (M_S + M_T)
beta = M_S / (M_S + M_T)
X_S = -alpha * d
X_T = beta * d

a = (X_T - X_S)

#angular speed of the Earth around sun
omega = np.sqrt(G * M_S / a**3)

print(omega)

# omega = 1.99e-7

#angular speed around G center
Omega = omega 

def function(x) : 
    return (Omega**2)*(x**5) - 2*(Omega**2)*(x**4)*(X_S+X_T) + (Omega**2)*(x**3)*((X_S + X_T)**2 + 2*X_S*X_T) - (x**2)*(G*(M_S+ M_T) + 2*(Omega**2)*X_T*X_S*(X_S+X_T)) + x*(2*G*(M_S*X_T + M_T*X_S) + (Omega**2)*(X_T**2)*(X_S**2)) - G*(M_S*X_T**2 + M_T*X_S**2)

L2 = root(function, 1.5e9)
L2x = L2.x[0]
print(f"{L2.x[0]:.11e}")

#distance from the center of the Earth to L2
dist = L2x - X_T    
print(dist)