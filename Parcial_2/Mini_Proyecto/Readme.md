El archivo Simulink Helicoptero_2DOF_EKF es un modelo de un helicóptero de 2 grados de libertado controlado por una red neuronal en la cual la matriz h de coeficientes 
utilicé cosenos así como en las ecuaciones de los estados lo cual encontramos que tuvo mejor eficiencia a largo plazo para la estabilización del sistema. 

Luego tenemos el modelo donde mediante el archivo HELI2DOF3.mat y el archivo MIN_CUADRADOS6.m obtenemos la función de transferencia mediante identificación por mínimos cuadrados. En ambas 
el seguimiento es bastante bueno pero en la red neuronal el error es menor, por lo cual concluimos que el control mediante la red neuronal es más eficiente.
