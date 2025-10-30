%% ============================================================
% 📘 Identificación por Mínimos Cuadrados (LS) con datos de Simulink (.mat)
% Método: Mínimos Cuadrados para Modelo ARX
% =============================================================
clear all; close all; clc;

%% --- 1. CARGA DE DATOS Y PARÁMETROS ---
load('HELI2DOF3.mat');   %HELI 3 COMPARA ENTRADAS CON VECTOR X, 
out_ok=out';
t = out_ok(:,1);             % tiempo
u_total = out_ok(:,2:3);     % dos señales de entrada (u1: pitch, u2: yaw)
y_total = out_ok(:,4:7);     % cuatro salidas (theta, theta_dot, psi, psi_dot)

% Mostrar dimensiones de los datos
fprintf('Datos cargados correctamente:\n');
fprintf('Tiempo: %d muestras\n', length(t));
fprintf('Entradas: %d señales (columnas 2–3)\n', size(u_total,2));
fprintf('Salidas: %d señales (columnas 4–7)\n\n', size(y_total,2));

% Seleccionar manualmente la entrada y salida para identificación
disp('--- Selección de señales ---');
disp('Entradas disponibles: 1=Pitch, 2=Yaw');
u_idx = input('Seleccione el número de la entrada a usar (1 o 2): ');
disp('Salidas disponibles: 1=theta, 2=theta_dot, 3=psi, 4=psi_dot');
y_idx = input('Seleccione el número de la salida a usar (1–4): ');

u = u_total(:,u_idx);
y = y_total(:,y_idx);

% Verificación de longitudes
N_muestras = length(u);
if N_muestras ~= length(y) || N_muestras ~= length(t)
    error('Los vectores de entrada, salida y tiempo deben tener la misma longitud.');
end

T_muestreo = t(2) - t(1);  % tiempo de muestreo

%% --- 2. PARÁMETROS DEL MODELO ---
disp(' ');
disp('--- Parámetros del Modelo a Identificar ---');
Na = input('Ingrese el grado del DENOMINADOR (n, e.g., 2): ');
Nb = input('Ingrese el grado del NUMERADOR (m, e.g., 2): ');

num_parametros = Na + Nb;
if N_muestras <= num_parametros
    error(['El número de muestras (', num2str(N_muestras), ...
        ') debe ser mayor que el número de parámetros a estimar (', ...
        num2str(num_parametros), ').']);
end

%% --- 3. CONSTRUCCIÓN DE LA MATRIZ DE REGRESORES (Phi) Y VECTOR Y ---
k_inicio = max(Na, Nb) + 1; 
Y = y(k_inicio:N_muestras); 
Phi = zeros(length(Y), num_parametros);

for k = k_inicio:N_muestras
    fila_actual = k - k_inicio + 1;

    % Parte autorregresiva (coeficientes negativos)
    for i = 1:Na
        Phi(fila_actual, i) = -y(k - i); 
    end

    % Parte de entrada (coeficientes positivos)
    for j = 1:Nb
        Phi(fila_actual, Na + j) = u(k - j);
    end
end

%% --- 4. SOLUCIÓN POR MÍNIMOS CUADRADOS (LS) ---
P = Phi' * Phi; 
R = Phi' * Y;   
theta_hat = P \ R; 

%% --- 5. RESULTADOS Y FUNCIÓN DE TRANSFERENCIA ---
disp(' ');
disp('==================================================');
disp('RESULTADOS DE LA IDENTIFICACIÓN POR MÍNIMOS CUADRADOS');
disp(['Entrada usada: u', num2str(u_idx)]);
disp(['Salida usada: y', num2str(y_idx)]);
disp(['Grado del Denominador (n_a): ', num2str(Na)]);
disp(['Grado del Numerador (n_b): ', num2str(Nb)]);
disp(['Número de Muestras (N): ', num2str(N_muestras)]);
disp(['Tiempo de Muestreo (T): ', num2str(T_muestreo, '%.4f'), ' segundos']);
disp('--------------------------------------------------');

a_hat = theta_hat(1:Na);
b_hat = theta_hat(Na+1:end);

disp('Vector de Coeficientes Estimados (theta_hat):');
disp(theta_hat);

% Crear función de transferencia discreta
A_tf = [1; a_hat];
B_tf = [0; b_hat];

% Ajustar longitud
len = max(length(A_tf), length(B_tf));
A_tf = [A_tf; zeros(len - length(A_tf), 1)]';
B_tf = [B_tf; zeros(len - length(B_tf), 1)]';

Gz = tf(B_tf, A_tf, T_muestreo);

disp(' ');
disp('FUNCIÓN DE TRANSFERENCIA IDENTIFICADA (G(z)):' );
disp(Gz);
disp('==================================================');

%% --- 6. VALIDACIÓN ---
y_predicha = Phi * theta_hat;
y_real = y(k_inicio:N_muestras);
t_plot = t(k_inicio:N_muestras); 
error_pred = y_real - y_predicha;

figure;
subplot(2,1,1);
plot(t_plot, y_real, 'b', 'LineWidth', 1.5); hold on;
plot(t_plot, y_predicha, 'r--', 'LineWidth', 1.5);
title(['Comparación: y', num2str(y_idx), ' real vs. predicha']);
legend('y real', 'y estimada');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

subplot(2,1,2);
plot(t_plot, error_pred, 'k', 'LineWidth', 1);
title('Error de Predicción e(k)');
xlabel('Tiempo (s)');
ylabel('Error');
grid on;

%% --- 7. MÉTRICA DE ERROR ---
ECM = mean(error_pred.^2);
fprintf('\nError Cuadrático Medio (ECM): %.6f\n', ECM);
