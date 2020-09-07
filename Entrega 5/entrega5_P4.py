from matplotlib.pylab import *
from scipy.integrate import odeint
import scipy as sp
from leer_eof import leer_eof
from time import perf_counter

ti = perf_counter()
#archivo E0F
fname = "S1B_OPER_AUX_POEORB_OPOD_20200813T110811_V20200723T225942_20200725T005942.EOF"
t, x, y, z, vx, vy, vz = leer_eof(fname)
z0 = [ x[0], y[0], z[0], vx[0], vy[0], vz[0] ]
vf_real = [ x[-1], y[-1], z[-1] ]
dT = t[-1]  #vector delta tiempo

#-----------------------------------------------------------------------------------------------
#-------                                                                                   -------

time = t/3600
xm = x/1000
ym = y/1000
zm = z/1000
m = 1000

#---------------------------------------------------------------------
#------- 2 MEJORADO) Comparación de la solución entre odeint y eulerint -------

#Unidades base:
mt = 5.9722*(10**24)  #Masa tierra [kg]
G = 6.67*(10**(-11))  #[kg-1 m3 s-2]
omega = 7.27*(10**(-5)) #[rad/s]
r = 6371000.0 #radio tierra [m]

J2 = (1.75553*(10**10))*(1000**5) #m5
J3 = -2.61913*(10**11)*(1000**6) #m6

def satelite(z, t):

    zp = sp.zeros(6)
    zp[0:3] = z[3:6]
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

    #Mejoras agregadas

    #J2
    zp[3] = zp[3]+ J2 *z[0]*((6*z[2]**2-3*(z[0]**2+z[1]**2)/2)/(r**7))
    zp[4] = zp[4]+ J2 *z[1]*((6*z[2]**2-3*(z[0]**2+z[1]**2)/2)/(r**7))
    zp[5] = zp[5]+ J2 *z[2]*((3*z[2]**2-9*(z[0]**2+z[1]**2)/2)/(r**7))
    #J3
    zp[3] = zp[3]+ J3 *z[0]*z[2]*((10*z[2]**2-15*(z[0]**2+z[1]**2)/2)/(r**9))
    zp[4] = zp[4]+ J3 *z[1]*z[2]*((10*z[2]**2-15*(z[0]**2+z[1]**2)/2)/(r**9))
    zp[5] = zp[5]+ J3 *(((4*z[2]**2*z[2]**2-3*(z[0]**2+z[1]**2))+1.5*(z[0]**2+z[1]**2)**2)/(r**9))

    return zp

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

# Deriva entre real y predicción
def deriva(t, z_ode, x, y, z):
    v_deriva = sp.zeros(len(t))
    for pos in range(len(t)):
        v_deriva[pos] = sp.sqrt(sp.dot((z_ode[pos, :3]- [x[pos],y[pos],z[pos]]),
                                       (z_ode[pos, :3]- [x[pos],y[pos],z[pos]])))

    v_d = v_deriva/1000

    dist_error = sp.sqrt(sp.dot([x[pos],y[pos],z[pos]],
                                [x[pos],y[pos],z[pos]]))

    mdistancia = sp.round_(sp.amax(v_d), 1)  # distancia máxima entre ambos
    err = sp.round_(v_deriva[-1] / dist_error, 2)  # % de error probando distintos N's
    #print(f"Error = {err * 100}%")

 #Gráfico P2
    figure(1)

    title(f"Distancia entre Eulerint y Odeint, d = {mdistancia} [km]")
    xlabel("Tiempo, t [Horas]")
    ylabel("Deriva, d [KM]")
    plot(time, v_d, color="dodgerblue")
    grid(False)

    tight_layout()
    savefig("P4.2.png")

#vector de tiempo
v_t = sp.linspace(0, dT, len(t))
v_tm = v_t/3600
#tiempos

t_1 = perf_counter()
z_ode = odeint(satelite, z0, v_t)
t_2 = perf_counter()
t_ode = t_2 - t_1
z_eul = eulerint(satelite, z0, v_t, Nsubdivisiones=1)
t_3 = perf_counter()
t_eul = t_3 - t_2

#print (f"Tiempo solución odeint = {t_ode} seg")
deriva(t, z_ode, x, y, z)  #Calcula y genera gráfico P2

#----------------------------------------------------------------------------------------------------------------------
#------- 1 MEJORADO) Grafico: Posición REAL Y PREDICHA (x,y,z) en el tiempo del vector de estado de Sentinel 1A/B -----

xode = z_ode[: ,0]
yode = z_ode[: ,1]
zode = z_ode[: ,2]

figure(2)
#X
subplot(3,1,1)
title("Posición Sentinel 1A/B")
ylabel("X(t) [KM]", loc="center")
plot(time,xm,color="dodgerblue")
plot(v_tm,xode/m,color="orange")
grid(False)
#Y
subplot(3,1,2)
ylabel("Y(t) [KM]", loc="center")
plot(time,ym,color="dodgerblue")
plot(v_tm,yode/m,color="orange")
grid(False)
#Z
subplot(3,1,3)
xlabel("Tiempo, t (horas)", loc="center")
ylabel("Z(t) [KM]", loc="center")
plot(time,zm,color="dodgerblue")
plot(v_tm,zode/m,color="orange")
grid(False)

tight_layout()
savefig("P4.1.png")

tf = perf_counter()
tiempo_total = tf - ti
print (f"Tiempo código {tiempo_total} seg")