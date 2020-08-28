import scipy as sp
from scipy.integrate import odeint
from datetime import datetime

'''
- Vector de estado inicial:

      <TAI>TAI=2020-07-23T23:00:19.000000</TAI>
      <UTC>UTC=2020-07-23T22:59:42.000000</UTC>
      <UT1>UT1=2020-07-23T22:59:41.786009</UT1>
      <Absolute_Orbit>+22605</Absolute_Orbit>
      <X unit="m">1584575.803741</X>
      <Y unit="m">-6709292.216083</Y>
      <Z unit="m">1590335.796951</Z>
      <VX unit="m/s">-1941.773309</VX>
      <VY unit="m/s">1267.169262</VY>
      <VZ unit="m/s">7235.465130</VZ>
      <Quality>NOMINAL</Quality>

- Vector de estado final:

      <TAI>TAI=2020-07-25T01:00:19.000000</TAI>
      <UTC>UTC=2020-07-25T00:59:42.000000</UTC>
      <UT1>UT1=2020-07-25T00:59:41.786602</UT1>
      <Absolute_Orbit>+22620</Absolute_Orbit>
      <X unit="m">-84713.447559</X>
      <Y unit="m">-3713900.787173</Y>
      <Z unit="m">-6029398.518733</Z>
      <VX unit="m/s">-2443.754947</VX>
      <VY unit="m/s">-6085.837364</VY>
      <VZ unit="m/s">3785.220040</VZ>
      <Quality>NOMINAL</Quality>
      
'''

#Unidades base:
mt = 5.9722*(10**24)  #Masa tierra [kg]
G = 6.67*(10**(-11))  #[kg-1 m3 s-2]
omega = 7.27*(10**(-5)) #[rad/s]
r = 6371000.0 #radio tierra [m]
hr = 3600 #s

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

#Vectores de tiempo

ti = "2020-07-23T22:59:42.000000"
ti = ti.split("T")
ti = "{} {}".format(ti[0], ti[1])
ti = datetime.strptime(ti, '%Y-%m-%d %H:%M:%S.%f')
tf = "2020-07-25T00:59:42.000000"
tf = tf.split("T")
tf = "{} {}".format(tf[0],tf[1])
tf = datetime.strptime(tf, '%Y-%m-%d %H:%M:%S.%f')

deltaT = (tf-ti).total_seconds()
#print("deltaT: ", deltaT)
#print("deltaT/hr: ",deltaT/hr)
#exit()

#Posición inicial
x_i = 1584575.803741 #m
y_i = -6709292.216083
z_i = 1590335.796951
vx_i = -1941.773309 #m/s
vy_i = 1267.169262
vz_i = 7235.465130

#Posición final
x_f = -84713.447559 #m
y_f = -3713900.787173
z_f = -6029398.518733
vx_f = -2443.754947 #m/s
vy_f = -6085.837364
vz_f = 3785.220040

#vector de tiempo
t = sp.linspace(0, deltaT, 9361)
#tiempo = 3.3  #[horas] Voy cambiando el tiempo segun el nº de orbitas**
#t = sp.linspace(0, tiempo, 2000)

# Velocidad
#vt = 24548 #[km/h]      #Voy cambiando la velocidad**

#Solución
z0 = sp.array([x_i, y_i, z_i, vx_i, vy_i, vz_i])
sol = odeint(satelite, z0, t)

#Diferencia entre vectores
def dif(vpredicho, vreal):
    return sp.sqrt(sp.sum((vpredicho - vreal)**2))

vf_predicho = sol[-1, 0:3] #primeras 3 coordenadas de sol predicha
vf_real = [x_f, y_f, z_f] #vector final real

print(dif(vf_predicho, vf_real))