import scipy as sp
from scipy.integrate import odeint
import matplotlib.pylab as plt

#Unidades base:
cm = 0.01  # m
inch = 2.54*cm
g = 9.81  # m/s/s

#Coeficiente de arrastre:
rho = 1.225  # kg / m**3
cd = 0.47
D = 8.5*inch
r = D/2
A = sp.pi * r **2
CD = 0.5*rho*cd*A

#Masa:
m = 15  # kg

#Función a integrar:
# z es el vector de estado

# z = [x,y,vx,vy]
# dz/dt = bala(z,t)

#         [  z2     ]
# dz/dt = [         ]       (modelo)
#         [FD/m  -g ]

#Vector de estado:
# z[0] -> x
# z[1] -> y
# z[2] -> vx
# z[3] -> vy

def bala(z, t):
    zp = sp.zeros(4)
    zp[0] = z[2]
    zp[1] = z[3]
    v = z[2:4]  #saca las últimas 2 componentes (veloc)
    v[0] = v[0] - V
    v2 = sp.dot(v,v)
    vnorm = sp.sqrt(v2)
    FD = -CD * v2 * (v / vnorm)
    zp[2] = FD[0] / m
    zp[3] = FD[1] / m - g

    return zp

#Vientos
V = [0, 10.0, 20.0]  # m/s

for i in V:
    V = i
    # vector de tiempo
    t = sp.linspace(0, 30, 1001)
    # parte en el origen y vx = vy = 100 km/h
    vi = 100*1000/3600
    z0 = sp.array([0, 0, vi, vi])
    sol = odeint(bala, z0, t)
    x = sol[:,0]
    y = sol[:,1]

    plt.plot(x, y, label=f"V = {i} m/s")

plt.figure(1)
plt.grid(True)
plt.ylim(0,50)
plt.xlim(0,150)
plt.ylabel("Y (m)")
plt.xlabel("X (m)")
plt.title("Trayectoria para distintos vientos")
plt.legend()
plt.tight_layout()

#Generar imagen png
plt.savefig("Gráfico - Trayectoria para distintos vientos.png")