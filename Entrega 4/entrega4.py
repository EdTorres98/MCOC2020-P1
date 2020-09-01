from matplotlib.pylab import *
from scipy.integrate import odeint
import scipy as sp

#Notar que eulerint es inestable Nsubdivisiones=1....
##Notar que eulerint es estable Nsubdivisiones>1   5,10,100....

#a = -1.
m = 1. #[kg]
f = 1. #[Hz]
E = 0.2
phi = sp.pi
w = 2*phi*f
k = m*(w**2)
c = 2*E*w*m

def eulerint(zp, z0, t, Nsubdivisiones=1):
    Nt = len(t)
    Ndim = len(array(z0))

    z = sp.zeros((Nt, Ndim))
    z[0,:] = z0

    # z (i+1) = zp_i * dt + z_i
    for i in range(1, Nt):
        t_anterior = t[i-1]
        dt = (t[i] - t[i-1])/Nsubdivisiones
        z_temp = z[i-1, :].copy()
        for k in range(Nsubdivisiones):
            z_temp += dt * zp(z_temp, t_anterior + k*dt)
        z[i, :] = z_temp

    return z

# z'= a * z
def zp(z,t):
    #return a*z
    zp = sp.zeros(2)
    zp[0] = z[1]
    x = z[0]
    x_p = z[1]
    # ec diferencial armónico  ((d2y)/(dt2)) = - (cy/m) - k/m * (dy/dt)
    zp[1] = -((c*x_p) + (k*x))/m

    return zp

z0 = [1, 1]     #Condición inicial
t = sp.linspace(0, 4., 100) #enunciado

#--- odeint ---
z_odeint = odeint(zp, z0, t)
#z_odeint = sol[:, 0]

#--- euler para subdivisiones 1, 10 y 100 ---
z_euler1 = eulerint(zp, z0, t, Nsubdivisiones=1)
z_euler10 = eulerint(zp, z0, t, Nsubdivisiones=10)
z_euler100 = eulerint(zp, z0, t, Nsubdivisiones=100)
#z_euler = sol[:, 0]

#--- real calculada---      Solución analítica
z_real = ((2*(m**2)) / (sp.sqrt((c**2) - (4*k*m))*m - (c*m))) * sp.exp(((sp.sqrt((c**2)-(4*k*m)))/(2*m)-(c/(2*m)))*t) + \
         (sp.sqrt((c**2)-(4*k*m))*m-(c+2*m)*m)/(sp.sqrt(c**2-4*k*m)*m-c*m) * sp.exp((((-sp.sqrt((c**2) - \
         (4*k*m)))/(2*m)-(c/(2*m)))*t))

#z_real = (sp.exp(-E*w*t)*(m*cos(w*((1-E**2)**0.5)*t)) + ((1 + w*E*m)/(w*((1-E**2)**0.5)))*sin(w*((1-E**2)**0.5) *t))

#--- graficar ---
figure()
xlabel("t [s]")
ylabel("x(t)")
xlim(0, 4)
                #solicitados
plot(t, z_real, label="real (analítica)", color = "black", linewidth = 2)
plot(t, z_odeint[:, 0], label = "odeint", color= "blue")
plot(t, z_euler1[:, 0], linestyle = "--", label="eulerint 1", color ="green")
plot(t, z_euler10[:, 0], linestyle = "--", label="eulerint 10", color = "red")
plot(t, z_euler100[:, 0], linestyle = "--", label="eulerint 100", color= "orange")
#plot(t, z_odeint, label="z_odeint")
#plot(t, z_euler, label="z_eulerint")
#plot(t, z_real, "--", label="z_real")
legend()

savefig("entrega4.png")