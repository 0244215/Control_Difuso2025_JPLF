% 5 partículas, 100 iteraciones
nP  = 5;
nIt = 10;

% Caja de INICIALIZACIÓN (no limita la búsqueda, solo el arranque)
init_lb = [0 0];      % arranque Kp >= 0, Kd >= 0
init_ub = [40 10];

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
