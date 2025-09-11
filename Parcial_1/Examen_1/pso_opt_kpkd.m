function [pbestX, pbestJ, gbestX, gbestJ, log] = pso_opt_kpkd(nParticles, nIter, init_lb, init_ub, params)

% - nParticles: # partículas 
% - nIter     : # iteraciones 
% - params    : struct opcional con campos w,c1,c2,vMax
% Devuelve:
% - pbestX (nParticles x 2): mejor [Kp,Kd] de CADA partícula
% - pbestJ (nParticles x 1): costo de cada partícula (su mejor)
% - gbestX (1x2) y gbestJ:   mejor global entre las 5
% - log: historial de gbest por iteración

    if nargin < 1 || isempty(nParticles), nParticles = 5; end
    if nargin < 2 || isempty(nIter),     nIter     = 10; end
    if nargin < 3 || isempty(init_lb),   init_lb   = [0 0]; end
    if nargin < 4 || isempty(init_ub),   init_ub   = [15 5]; end
    if nargin < 5, params = struct; end

    % Hiperparámetros
    if ~isfield(params,'w'),    params.w  = 0.72; end
    if ~isfield(params,'c1'),   params.c1 = 1.49; end
    if ~isfield(params,'c2'),   params.c2 = 1.49; end
    % vMax opcional: si NaN o vacío => sin tope de velocidad
    if ~isfield(params,'vMax') || isempty(params.vMax)
        params.vMax = 0.3*(init_ub - init_lb);   % razonable por defecto
    end

    rng('shuffle');
    dim = 2;

    % ---- Inicialización (SOLO usa init_lb/ub; luego no se recorta) ----
    X = init_lb + rand(nParticles,dim).*(init_ub - init_lb);  % posiciones
    V = zeros(nParticles,dim);                                 % velocidades

    % ---- Evaluación inicial ----
    pbestX   = X;
    pbestJ   = inf(nParticles,1);
    for i = 1:nParticles
        pbestJ(i) = sim_plant_pend_inv_pd(X(i,:)');   % tu función de costo (MSE)
    end
    [gbestJ, idx] = min(pbestJ);
    gbestX = pbestX(idx,:);

    % Historial
    log.gbestHist = zeros(nIter, dim);
    log.costHist  = zeros(nIter, 1);

    % ---- Bucle PSO ----
    for it = 1:nIter
        for i = 1:nParticles
            r1 = rand(1,dim);
            r2 = rand(1,dim);

            % Velocidad y posición (SIN recortar posición)
            V(i,:) = params.w*V(i,:) ...
                   + params.c1*r1.*(pbestX(i,:)-X(i,:)) ...
                   + params.c2*r2.*(gbestX   -X(i,:));

            % Tope opcional de velocidad
            if all(isfinite(params.vMax))
                V(i,:) = max(-params.vMax, min(params.vMax, V(i,:)));
            end

            X(i,:) = X(i,:) + V(i,:);    % <-- sin límites de posición

            % Evaluar
            J = sim_plant_pend_inv_pd(X(i,:)');

            % Mejor local de ESTA partícula
            if J < pbestJ(i)
                pbestJ(i) = J;
                pbestX(i,:) = X(i,:);
            end

            % Mejor global entre TODAS las partículas
            if J < gbestJ
                gbestJ = J;
                gbestX = X(i,:);
            end
        end

        % Log e impresión
        log.gbestHist(it,:) = gbestX;
        log.costHist(it)    = gbestJ;
        fprintf('Particle %d/%d | J=%.6g | Gbest: Kp=%.4g, Kd=%.4g\n', it, nParticles, gbestJ, gbestX(1), gbestX(2));
    end

    % Al final:
    %  pbestX/pbestJ -> "mejor por partícula" (5 resultados si nParticles=5)
    %  gbestX/gbestJ -> mejor global entre esas 5
end
