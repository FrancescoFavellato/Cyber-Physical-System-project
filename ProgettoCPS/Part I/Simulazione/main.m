% PROJECT - PART 1
clear all; close all; clc;

% Carico i dati
load('data9.mat');

n = size(A, 1);     % state dimension = 50
q = size(C, 1);     % number of sensors = 20   
h = sum(a ~= 0);    % number of attacked sensors = 4
T = 500;            % numero di campioni 

fprintf('  Numero di stati: %d\n',  n);
fprintf('  Numero di sensori: %d\n',  q);
fprintf('  Sensori attaccati (h): %d -> %s\n\n', h, mat2str(find(a ~= 0)'));

% Verifica se la matrice A è stocastica
isStochastic = all(abs(sum(A, 2) - 1) < 1e-5);

if isStochastic
    fprintf('La matrice A è stocastica');
    % Plot del grafo corrispondente 
    Adj = A;
    for i = 1:n
        Adj(i,i) = 0;
    end
    Adj = Adj ~= 0;
    leaders = find(abs(diag(A) - 1) < 1e-5);
    figure('Name', 'Topologia CPS', 'NumberTitle', 'off');
    p = plot(digraph(Adj'));    % traspongo perché indica archi uscenti
    p.NodeColor = 'r'; p.MarkerSize = 8; p.LineWidth = 1.5; p.ArrowSize = 12;
    highlight(p, leaders,'NodeColor', 'g', 'MarkerSize', 12);
    title('Visualizzazione del Grafo Orientato');
end

% ========== Generazioni misure y(k) =======
x_true = zeros(n, T);
y = zeros(q, T);
Apow_true = eye(n); 

for k = 1:T
    x_true(:, k) = Apow_true * x0;  % Per k=1 è x(0) = I * x0
    y(:, k) = C*x_true(:, k) + a; 
    Apow_true = A * Apow_true;      % A^k per il passo successivo
end

%% ========== ALGORITMO 1: O-PGD (Online Proximal Gradient Descent) =======
nu = 0.125;   
lambdaOPGD = 0.3;

fprintf('Parametri O-PGD:\n');
fprintf('  nu     = %.2e\n', nu);
fprintf('  lambda = %.4f\n\n', lambdaOPGD);

ris_opgd = OPGD(A, C, y, a, x_true, lambdaOPGD, nu, T);

fprintf('Tempo esecuzione O-PGD: %.4f\n\n', ris_opgd.time);

% --- Plot metriche ---
figure('Name', 'O-PGD - Metriche', 'NumberTitle', 'off');
 
subplot(2, 2, 1);
semilogy(1:T, ris_opgd.state_err, 'b', 'LineWidth', 1.5);
xlabel('k'); ylabel('Errore');
title('State Estimation Error');
grid on;
 
subplot(2, 2, 2);
plot(1:T, ris_opgd.support_err, 'r', 'LineWidth', 1.5);
xlabel('k'); ylabel('# sensori');
title('Support Attack Error');
ylim([-0.5, q + 0.5]);
grid on;
 
subplot(2, 2, 3);
plot(1:T, ris_opgd.false_pos, 'm', 'LineWidth', 1.5);
xlabel('k'); ylabel('# sensori');
title('False Positives');
ylim([-0.5, q + 0.5]);
grid on;
 
subplot(2, 2, 4);
plot(1:T, ris_opgd.false_neg, 'g', 'LineWidth', 1.5);
xlabel('k'); ylabel('# sensori');
title('False Negatives');
ylim([-0.5, q + 0.5]);
grid on;
 
sgtitle('O-PGD: metriche nel tempo');

%% ========== ALGORITMO 2: O-PAM (Online Proximal Alernating Minimization) =======
mu = 0.1;      
lambdaOPAM = 0.3;

fprintf('Parametri O-PAM:\n');
fprintf('  mu     = %.2e\n', mu);
fprintf('  lambda = %.4f\n\n', lambdaOPAM);

ris_opam = OPAM(A, C, y, a, x_true, lambdaOPAM, mu, T);

fprintf('Tempo esecuzione O-PAM: %.4f\n\n', ris_opam.time);

% --- Plot metriche ---
figure('Name', 'O-PAM - Metriche', 'NumberTitle', 'off');
 
subplot(2, 2, 1);
semilogy(1:T, ris_opam.state_err, 'b', 'LineWidth', 1.5);
xlabel('k'); ylabel('Errore');
title('State Estimation Error');
grid on;
 
subplot(2, 2, 2);
plot(1:T, ris_opam.support_err, 'r', 'LineWidth', 1.5);
xlabel('k'); ylabel('# sensori');
title('Support Attack Error');
ylim([-0.5, q + 0.5]);
grid on;
 
subplot(2, 2, 3);
plot(1:T, ris_opam.false_pos, 'm', 'LineWidth', 1.5);
xlabel('k'); ylabel('# sensori');
title('False Positives');
ylim([-0.5, q + 0.5]);
grid on;
 
subplot(2, 2, 4);
plot(1:T, ris_opam.false_neg, 'g', 'LineWidth', 1.5);
xlabel('k'); ylabel('# sensori');
title('False Negatives');
ylim([-0.5, q + 0.5]);
grid on;
 
sgtitle('O-PAM: metriche nel tempo');

%% ========== ALGORITMO 3: AO-PGD (Aggregated Online Proximal Gradient Descent) =======
nu = 0.125;   
lambdaAOPGD = 0.3;

fprintf('Parametri AO-PGD:\n');
fprintf('  nu     = %.2e\n', nu);
fprintf('  lambda = %.4f\n\n', lambdaAOPGD);

ris_aopgd = AOPGD(A, C, y, a, x_true, lambdaAOPGD, nu, T);

fprintf('Tempo esecuzione AO-PGD: %.4f\n\n', ris_aopgd.time);

% --- Plot metriche ---
figure('Name', 'AO-PGD - Metriche', 'NumberTitle', 'off');
 
subplot(2, 2, 1);
semilogy(1:T, ris_aopgd.state_err, 'b', 'LineWidth', 1.5);
xlabel('k'); ylabel('Errore');
title('State Estimation Error');
grid on;
 
subplot(2, 2, 2);
plot(1:T, ris_aopgd.support_err, 'r', 'LineWidth', 1.5);
xlabel('k'); ylabel('# sensori');
title('Support Attack Error');
ylim([-0.5, q + 0.5]);
grid on;
 
subplot(2, 2, 3);
plot(1:T, ris_aopgd.false_pos, 'm', 'LineWidth', 1.5);
xlabel('k'); ylabel('# sensori');
title('False Positives');
ylim([-0.5, q + 0.5]);
grid on;
 
subplot(2, 2, 4);
plot(1:T, ris_aopgd.false_neg, 'g', 'LineWidth', 1.5);
xlabel('k'); ylabel('# sensori');
title('False Negatives');
ylim([-0.5, q + 0.5]);
grid on;
 
sgtitle('AO-PGD: metriche nel tempo');

%% ========== ALGORITMO 4: AO-PAM (Aggregated Online Proximal Alternating Minimization) =======
mu = 0.1;   
lambdaAOPAM = 0.3;

fprintf('Parametri AO-PAM:\n');
fprintf('  mu     = %.2e\n', mu);
fprintf('  lambda = %.4f\n\n', lambdaAOPAM);

ris_aopam = AOPAM(A, C, y, a, x_true, lambdaAOPAM, mu, T);

fprintf('Tempo esecuzione AO-PAM: %.4f\n\n', ris_aopam.time);

% --- Plot metriche ---
figure('Name', 'AO-PAM - Metriche', 'NumberTitle', 'off');
 
subplot(2, 2, 1);
semilogy(1:T, ris_aopam.state_err, 'b', 'LineWidth', 1.5);
xlabel('k'); ylabel('Errore');
title('State Estimation Error');
grid on;
 
subplot(2, 2, 2);
plot(1:T, ris_aopam.support_err, 'r', 'LineWidth', 1.5);
xlabel('k'); ylabel('# sensori');
title('Support Attack Error');
ylim([-0.5, q + 0.5]);
grid on;
 
subplot(2, 2, 3);
plot(1:T, ris_aopam.false_pos, 'm', 'LineWidth', 1.5);
xlabel('k'); ylabel('# sensori');
title('False Positives');
ylim([-0.5, q + 0.5]);
grid on;
 
subplot(2, 2, 4);
plot(1:T, ris_aopam.false_neg, 'g', 'LineWidth', 1.5);
xlabel('k'); ylabel('# sensori');
title('False Negatives');
ylim([-0.5, q + 0.5]);
grid on;
 
sgtitle('AO-PAM: metriche nel tempo');

%% --- Confronto attacco vero vs stimato ---
figure('Name', 'Attacco: vero vs stimato', 'NumberTitle', 'off');
bar([a, ris_opgd.a_hat, ris_opam.a_hat, ris_aopgd.a_hat, ris_aopam.a_hat]);
legend('a vero', 'a stimato (O-PGD)', 'a stimato (O-PAM)', 'a stimato (AO-PGD)', 'a stimato (AO-PAM)');
xlabel('Indice sensore (1..20)');
ylabel('Valore');
title('Confronto vettore di attacco');
grid on;

% --- Confronto errore di stima tra i diversi algoritmi ---
figure('Name', 'Confronto Errore di Stato', 'NumberTitle', 'off');
semilogy(1:T, ris_opgd.state_err, 'r', 'LineWidth', 1.5); hold on;
semilogy(1:T, ris_opam.state_err, 'b', 'LineWidth', 1.5);
semilogy(1:T, ris_aopgd.state_err, 'm', 'LineWidth', 1.5);
semilogy(1:T, ris_aopam.state_err, 'c', 'LineWidth', 1.5);
xlabel('Passi temporali (k)'); ylabel('Errore quadratico medio');
legend('O-PGD', 'O-PAM', 'AO-PGD', 'AO-PAM');
title('Confronto State Estimation Error');
grid on;

% --- Confronto errore di stima del supporto dell'attacco ---
figure('Name', 'Confronto Errore di Stima', 'NumberTitle', 'off');
plot(1:T, ris_opgd.support_err, 'r', 'LineWidth', 1.5); hold on;
plot(1:T, ris_opam.support_err, 'b', 'LineWidth', 1.5);
plot(1:T, ris_aopgd.support_err, 'm', 'LineWidth', 1.5);
plot(1:T, ris_aopam.support_err, 'c', 'LineWidth', 1.5);
xlabel('Passi temporali (k)'); ylabel('Errore di Stima');
legend('O-PGD', 'O-PAM', 'AO-PGD', 'AO-PAM');
title('Confronto Support Attack Error');
grid on;