Simulado en Matlab 2024b

La simulación consiste en la planta del péndulo simple controlado por un control difuso con un controlador PD, el cual está programado mediante matlab function.
Entre sus entradas recibe la señal de entrada previamente afectada por el controlador PD, recibe el número de reglas que lleva nuestro controlador, que luego,
la raiz del número de reglas son las secciones en las que se divide nuestro error, derivada del error y salida. Recibe los rangos de cada uno de estos. 

Ahora bien, el objetivo en específico de esta simulación fue encontrar los valores de Kp y Kd ideales de nuestro sistema tales que se redujera al máximo el error 
a través del sistema y el tiempo de reacción. Para ello utilicé el método visto en clase llamado Optimización por enjambre de partículas. El cual mediante 2 archivos
.m uno proporcionado por el profesor que tiene la función de costo y el otro es para obtener y comparar las iteraciones y obtener las mejores locales y la mejor global.

Para ejecutar la simulación primero se pone un valor para Kp=5 y otro para Kd=1 en consola, luego corrí el simulink y finalmente ejecuté el siguiente código en la consola
que llama a mis 2 archivos .m hace 10 iteraciones en cada partícula (5 partículas), me las va escribiendo en consola y finalmente me escribe de cada partícula la 
mejor local y hasta abajo la mejor global misma que asigna directamente a mi valor de Kp y Kd de manera que si vuelvo a ejecutar mi simulink tengo en el scope la 
resultante del comportamiento con uso de el mejor global de Kp y Kd y el resultado fue óptimo debido a que se redujo la oscilación inicial, el tiempo de respuesta 
y en general el error respecto a la referencia.

% 5 partículas, 10 iteraciones
nP  = 5;
nIt = 10;

% (Opcional) parámetros: sin tope de velocidad => busca sin freno
params = struct('w',0.72,'c1',1.49,'c2',1.49,'vMax', [Inf Inf]);

[pbestX, pbestJ, gbestX, gbestJ, log] = pso_opt_kpkd(nP, nIt, init_lb, init_ub, params);

% Resultados:
disp('Mejores por partícula (Kp,Kd | costo):');
for i=1:nP
    fprintf('  p%02d: Kp=%8.5f  Kd=%8.5f   J=%g\n', i, pbestX(i,1), pbestX(i,2), pbestJ(i));
end
fprintf('\nMejor GLOBAL: Kp=%8.5f  Kd=%8.5f   J=%g\n', gbestX(1), gbestX(2), gbestJ);

% Aplicar los óptimos globales y volver a simular
assignin('base','Kp',gbestX(1));
assignin('base','Kd',gbestX(2));
sim('SimulacionPenduloSimpleFCN.slx');
