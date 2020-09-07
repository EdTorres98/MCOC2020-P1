from matplotlib.pylab import *
from scipy.integrate import odeint
import scipy as sp
from time import perf_counter
from leer_eof import leer_eof

#archivo E0F
fname = "S1B_OPER_AUX_POEORB_OPOD_20200813T110811_V20200723T225942_20200725T005942.EOF"
t, x, y, z, vx, vy, vz = leer_eof(fname)
z0 = [ x[0], y[0], z[0], vx[0], vy[0], vz[0] ]
vf_real = [ x[-1], y[-1], z[-1] ]
dT = t[-1]  #vector delta tiempo

#-----------------------------------------------------------------------------------------------
#------- 1) Grafico: Posición (x,y,z) en el tiempo del vector de estado de Sentinel 1A/B -------

figure(1)
time = t/3600
xm = x/1000
ym = y/1000
zm = z/1000

#X
subplot(3,1,1)
title("Posición REAL Sentinel 1A/B")
ylabel("X(t) [KM]", loc="center")
plot(time,xm,color="dodgerblue")
grid(False)
#Y
subplot(3,1,2)
ylabel("Y(t) [KM]", loc="center")
plot(time,ym,color="dodgerblue")
grid(False)
#Z
subplot(3,1,3)
xlabel("Tiempo, t (horas)", loc="center")
ylabel("Z(t) [KM]", loc="center")
plot(time,zm,color="dodgerblue")
grid(False)

tight_layout()
savefig("P1.png")

#---------------------------------------------------------------------
#------- 2) Comparación de la solución entre odeint y eulerint -------

#Unidades base:
mt = 5.9722*(10**24)  #Masa tierra [kg]
G = 6.67*(10**(-11))  #[kg-1 m3 s-2]
omega = 7.27*(10**(-5)) #[rad/s]
r = 6371000.0 #radio tierra [m]

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

# Deriva entre odeint y eulerint
def deriva(t, z_ode, z_eul):
    v_deriva = sp.zeros(len(t))
    for pos in range(len(t)):
        v_deriva[pos] = sp.sqrt(sp.dot((z_ode[pos, :3]-z_eul[pos, :3]),
                                       (z_ode[pos, :3]-z_eul[pos, :3])))
    v_d = v_deriva/1000

    mdistancia = sp.round_(sp.amax(v_d), 1)  #distancia máxima entre ambos

 #Gráfico P2
    figure(2)

    title(f"Distancia entre Eulerint y Odeint, d = {mdistancia} [km]")
    xlabel("Tiempo, t [Horas]")
    ylabel("Deriva, d [KM]")
    plot(time, v_d, color="dodgerblue")
    grid(False)

    tight_layout()
    savefig("P2.png")

#vector de tiempo
v_t = sp.linspace(0, dT, len(t))
#tiempos

t_1 = perf_counter()
z_ode = odeint(satelite, z0, v_t)
t_2 = perf_counter()
t_ode = t_2 - t_1
z_eul = eulerint(satelite, z0, v_t, Nsubdivisiones=1)
t_3 = perf_counter()
t_eul = t_3 - t_2

print (f"Tiempo solución odeint = {t_ode} seg")
print (f"Tiempo solución eulerint = {t_eul} seg")
deriva(t, z_ode, z_eul)  #Calcula y genera gráfico P2

#---------------------------------------------------------------------
#------- 3) N para 1% de error en la preddicción con eulerint -------

def deriva2(t, z_ode, z_eul):
    v_deriva = sp.zeros(len(t))
    for pos in range(len(t)):
        v_deriva[pos] = sp.sqrt(sp.dot((z_ode[pos, :3]-z_eul[pos, :3]),
                                       (z_ode[pos, :3]-z_eul[pos, :3])))

    v_d = v_deriva/1000
    dist_error = sp.sqrt(sp.dot(z_ode[-1, :3],
                                z_ode[-1, :3]))

    mdistancia = sp.round_(sp.amax(v_d), 1)  #distancia máxima entre ambos
    err = sp.round_(v_deriva[-1] / dist_error, 1) #% de error probando distintos N's
    print(f"Error = {err*100}%")

 #Gráfico P3
    figure(3)

    title(f"Distancia entre Eulerint y Odeint, d = {mdistancia} [km]")
    xlabel("Tiempo, t [Horas]")
    ylabel("Deriva, d [KM]")
    plot(time, v_d, color="dodgerblue")
    grid(False)

    tight_layout()
    savefig("P3.png")

#vector de tiempo
v_t = sp.linspace(0, dT, len(t))
#tiempos

t_1 = perf_counter()
z_ode = odeint(satelite, z0, v_t)
t_2 = perf_counter()
t_ode = t_2 - t_1
z_eul = eulerint(satelite, z0, v_t, Nsubdivisiones=2)  #PARA REVISIÓN DE P3 CAMBIAR A N=1500 Para que resulte el gráfico utilizado
                                                          #Se puso este valor temporal para que P1 y P2 no se demoraran en obtener.
t_3 = perf_counter()
t_eul = t_3 - t_2

print (f"Tiempo solución odeint = {t_ode} seg")
print (f"Tiempo solución eulerint = {t_eul} seg")
deriva2(t, z_ode, z_eul)  #Calcula y genera gráfico P3

#N1 = 270% error
#N100 = 140% error
#N100 = 20% error
#N1500 = 10% error