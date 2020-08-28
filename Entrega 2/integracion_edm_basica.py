import scipy as sp
from scipy.integrate import odeint
import matplotlib.pylab as plt

#Función a integrar:
# z es el vector de estado

# z = [x,y,vx,vy]
# dz/dt = satelite(z,t)

#         [   z2   ]
# dz/dt = [        ]     (modelo)
#         [FD/m  -g]

#Vector de estado:
# z[0] -> x
# z[1] -> y
# z[2] -> vx
# z[3] -> vy

#Unidades base:
mt = 5.9722*(10**24)  #Masa tierra [kg]
G = 6.67*(10**(-11))*(10**(-9))*(3600**2)  #[kg-1 km3 h-2]
omega = 7.27*(10**(-5))/(1/3600) #[rad/h]
r = 6.371*(10**3) #radio tierra [km]

def satelite(z, t):

    zp = sp.zeros(6) #extensión del vector a uno de 6
    zp[0:3] = z[3:6] #tranf de coordenadas (zp[0:2] = z[2:4] video)
    z1 = z[0:3]
    z2 = z[3:6]
    d = sp.sqrt(z[0]**2 + z[1]**2 + z[2]**2) #distancia satelite -> tierra

    R_ = sp.array([[sp.cos(omega*t), -sp.sin(omega*t),     0],
                   [sp.sin(omega*t), sp.cos(omega*t),     0],
                   [0,                     0,             1]])

    R_p = (sp.array([[-sp.sin(omega*t), -sp.cos(omega*t),    0],
                    [sp.cos(omega*t),  -sp.sin(omega*t),    0],
                    [0,                 0,          0]]))*omega

    R_pp = (sp.array([[-sp.cos(omega*t), sp.sin(omega*t),    0],
                      [-sp.sin(omega*t), -sp.cos(omega*t),   0],
                      [0,                 0,         0]]))*(omega**2)

    zp[3:6] = ((-G*mt/(d**3)) * z1) - R_.T@(R_pp@z1) - 2*R_.T@(R_p@z2)  #EDM

    return zp

#Vector de tiempo
tiempo = 3.3  #[horas] Voy cambiando el tiempo segun el nº de orbitas**
t = sp.linspace(0, tiempo, 2000)

# Velocidad
vi = 24548 #[km/h]      #Voy cambiando la velocidad**

#Solución
inicio = 6371.0+700.0   #radio tierra + altura satelite
z0 = sp.array([inicio , 0, 0, 0, vi, 0])  #coordenadas iniciales
sol = odeint(satelite, z0, t)

#----------------- Graficar -------------------

x = sol[:,0]
y = sol[:,1]
z = sol[:,2]


#---- 1) Historias de tiempo de x(t), y(t), z(t) para dos órbitas completas -----
plt.figure(1)

#Para x(t)
plt.subplot(3,1,1)
plt.xlabel("t [Horas]", loc="center")
plt.ylabel("x(t) [km]", loc="center")
plt.xlim(0, tiempo)       #Para las 2 orbitas solicitadas
plt.grid(True)
plt.plot(t, x, "b--")

#Para y(t)
plt.subplot(3,1,2)
plt.xlabel("t [Horas]", loc="center")
plt.ylabel("y(t) [km]", loc="center")
plt.xlim(0, tiempo)       #Para las 2 orbitas solicitadas
plt.grid(True)
plt.plot(t, y, "b--")

#Para z(t)
plt.subplot(3,1,3)
plt.xlabel("t [Horas]", loc="center")
plt.ylabel("z(t) [km]", loc="center")
plt.xlim(0, tiempo)       #Para las 2 orbitas solicitadas
plt.grid(True)
plt.plot(t, z, "b--")

plt.tight_layout()
plt.savefig("Gráfico 1.png")


#---- 2) Gráfico Tierra, Atmósfera y 2 Órbitas completas -----
plt.figure(2)

#Circunferencia Tierra
alfa = sp.linspace(0,2*sp.pi,5000)
xT = sp.cos(alfa)*r
yT = sp.sin(alfa)*r
plt.plot(xT, yT, color="saddlebrown", label = "Tierra")

#Circunferencia Atmósfera
h_atm = r+80.0
xA = sp.cos(alfa)*h_atm
yA = sp.sin(alfa)*h_atm
plt.plot(xA, yA, color="deepskyblue", label = "Atmósfera")

#Órbitas del satélite
plt.plot(x, y, "b--", label = "Satélite")

plt.axis("Equal") #dim x=y para el circulo
plt.grid(True)
plt.legend(loc = "lower right")
plt.tight_layout()
plt.savefig("Gráfico 2.png")


#---- 3) Gráficos superficie de la Tierra, Atmósfera y Orbita del Satélite -----
plt.figure(3)

plt.xlabel("t [Horas]")
plt.ylabel("r(t) [km]")
plt.xlim(0,tiempo)
#plt.xlim(2,2.5)   #acercar

#Linea superficie tierra
plt.hlines(r, 0.0, tiempo, color="saddlebrown", label = "Tierra")

#Linea superficie atmósfera
h_atm = r+80.0    #asumir a 80 km de la tierra
plt.hlines(h_atm, 0.0, tiempo, color="deepskyblue", label = "Atmósfera")

#Linea órbita satélite
vect_sat = sp.sqrt((x**2)+(y**2)+(z**2))   #vector del satélite
plt.plot(t, vect_sat, "b--", label = "Satélite")

plt.grid(True)
plt.legend(loc = "lower right")
plt.tight_layout()
plt.savefig("Gráfico 3.png")