# MCOC2020-P1

![Gráfico - Trayectoria para distintos vientos](https://user-images.githubusercontent.com/69275311/91117011-a271c980-e65b-11ea-8dc4-4e7900fc7d26.png)

# Primeras predicciones con la EDM básica del satélite

+ 1- Si se coloca el satélite a una altura de 700km (en la dirección x, ojo con el radio de la tierra) y se le da una velocidad inicial tangencial yp(0) = vt ¿Cuanto debe valer vt de modo que el satélite efectivamente orbite sin caer de vuelta dentro de la atmósfera? (asuma que esta comienza a una altura de 80km)
   + Para que el satélite realice 2 orbitas completas sin caer, vt debe ser aproximadamente mayor a 24548 km/h, ya que, cómo se puede apreciar en el gráfico 2, la órbita del satélite queda al límite de la atmósfera. Po lo tanto vt>24548 [km/h] para que no caiga.

+ 2- ¿Cómo encontró vt?
   + Principalmente el valor de vt se obtuvo mediante tanteo, se fueron probando valores para que la órbita del satélite en el gráfico 2 se fuera acercando o alejando (en caso de sobrepasar la atmósfera), por lo que se llegó a que para ese valor de vt (o mayores), la órbita no sobrepasaba la atmósfera, teniendo en cuenta las 2 órbitas mínimas que se pidieron.

+ Gráfico 1: Historias de tiempo x(t), y(t) y z(t) para 2 órbitas completas:

![Gráfico 1](https://user-images.githubusercontent.com/69275311/91509695-00dfb780-e8a9-11ea-8cf5-c79c46f96cf7.png)

   + Las horas fueron obtenidas de igual manera por tanteo, verificando en el gráfico 3 que el satélite realizara aproximadamente las dos órbitas completas, obteniendo un t = 3.3 Horas.

+ Gráfico 2: Órbita del satélite, superficie atmósfera y superficie de la tierra.

![Gráfico 3](https://user-images.githubusercontent.com/69275311/91509954-a561f980-e8a9-11ea-8b9d-4cbd991df125.png)

   + Gráfico realizado para vt = 24548 km/h, se puede apreciar que la órbita no logra sobrepasar la atmósfera. Posteriormente se presenta un gráfico mostrando el momento en el que la órbita queda al límite, pero no sobrepasa:

![Gráfico 3 cerca](https://user-images.githubusercontent.com/69275311/91510107-16a1ac80-e8aa-11ea-96ba-ac45bd89ecf5.png)

+ Gráfico 3: Órbitas del satélite junto a la superficie de la tierra y atmósfera completas:

![Gráfico 2](https://user-images.githubusercontent.com/69275311/91510189-55cffd80-e8aa-11ea-8fc9-d00526184a60.png)

   + Se logra apreciar de mejor manera que para las dos órbitas en un t = 3.3 Horas y una vt = 24548 km/h, estas no sobrepasan la atmósfera, quedando al límite en el extremo izquierdo, situación la cual fue descrita en el gráfico 1.

# Estudio de convergencia Método de Euler

+ Para el gráfico obtenido, se observa que al comparar las soluciones real, odeint, y eulerint para Nsubdivisiones(1,10,100), tenemos que las soluciones odeint y eulerint 100 son las que más se parecen y tienden a ser casi iguales a la solución real. Por lo que para la solución eulerint mientras menor sea la subdivisión, se alejará más de la solución real, siendo cada vez menos precisa, como se observa con eulerint 1.

   ![entrega4](https://user-images.githubusercontent.com/69275311/91844187-fe70bb00-ec24-11ea-8120-c3053bddd91a.png)


# Mejoras al modelo y estudio de convergencia

+ P1- El gráfico de la posición (x,y,z) real, en el tiempo del vector de estado de Sentinel 1A/B que me tocó, es el siguiente:

![P1](https://user-images.githubusercontent.com/69275311/92380952-5eada400-f0e0-11ea-90d1-3e7c2c630346.png)

+ P2- Usando la condición inicial (primer OSV) del archivo, se obtuvieron las soluciones de eulerint y odeint con N=1. A continuación se muestra el gráfico relacionado a la deriva con respecto al tiempo:

![P2](https://user-images.githubusercontent.com/69275311/92381236-d5e33800-f0e0-11ea-8af6-4864a66cf18e.png)

+ ¿Cuánto deriva eulerint de odeint en este caso al final del tiempo?             
+ Al final del tiempo (24 horas), la deriva es de 20868.8 km
+ ¿Cuanto se demora odeint y eulerint respectivamente en producir los resultados?          
+ Odeint = 1.01 seg, Eulerint = 2.63 seg

+ P3- ¿Cuantas subdivisiones hay que usar para que la predicción con eulerint al final del tiempo esté en menos de un 1% de error? 
+ Utilizando en este caso N=1500 se obtuvo un 10% de error, ya que para N mayores el tiempo de ejecución era mucho mayor y era inviable realizarlo. A continuación se muestra la deriva en el tiempo para eulerint:

![P3](https://user-images.githubusercontent.com/69275311/92382057-5b1b1c80-f0e2-11ea-97d7-2ed70e9285e0.png)

+ En este caso se obstuvo una deriva es de 750.2 km entre estos dos. Po otro lado el tiempo para producir los resultados de eulerint fueron 2947.57 seg o aproximadamente 49 minutos en ejecutar. Por lo que utilizar eulerint resulta ineficiente para estos casos.

+P4- Implementando las correcciones J2 y J3 y repitiendo los gráficos obtenidos se obtuvo lo siguiente:

![P4 1](https://user-images.githubusercontent.com/69275311/92382554-29ef1c00-f0e3-11ea-8fa2-40b233586736.png)
+ Se observa que la posición predicha para el satélite (color naranjo) con respecto a la posición obtenida anteriormente (color azul), es muy parecida, donde se nota una leve diferencia entre ambas curvas en su desplazamiento y llegando hacia el final de esta.

![P4 2](https://user-images.githubusercontent.com/69275311/92382596-370c0b00-f0e3-11ea-88d7-fa47b971799a.png)

+ ¿Cuánta deriva incurre al agregar las correcciones J2 y J3? 
+ Observando el gráfico obtenido, la deriva obtenida es de 660.9 km, muy inferior a los 20868.8 km que se obstuvieron para el gráfico sin las mejoras de J2 y J3, por lo que aplicarlas son de gran utilidad para obtener resultados más precisos.

+ ¿Cuanto se demora su código en correr?
+ El código con las mejoras implementadas demoró 5.91 segundos.
