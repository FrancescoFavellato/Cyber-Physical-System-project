% PROJECT - PART 2
clear all; close all; clc

%% ========================= SELEZIONE PARAMETRI  ========================= 

% Scelta del coupling gain
% 1: c = cl -> Forza uguale
% 2: c != cl -> Tracking forte
gain_choice = 1;

% Scelta del riferimento:
% 1: Riferimento con velocità costante
% 2: Riferimento con traiettoria circolare
% 3: Riferimento con traiettoria sinusoidale
traj_choice = 2;

% Scelta dell rete di comunicazione:
% 1: (Con albero ricoprente) 1 -> 2 -> 3
% 2: (Con albero ricoprente) 3 <- 1 -> 2
connection_choice = 1;

% Scelta della saturazione
% 0: Senza Saturazione 
% 1: Con Saturazione
u_choice = 0;

% Scelta del modello da simulare
% 1: Sistema linearizzato
% 2: Sistema NON linearizzato
model_choice = 1;

% IMPORTANTE: La simulazione Dinamica e quella Robusta si escludono a vicenda. 
% Se dynamic_choice = 1, robustness_choice e le scelte precedenti sulla
% topologia  vengono ignorate. Se dynamic_choice = 0 e robustness_choice =
% 1, usa solo il modello con incertezza. Altrimenti tiene conto delle
% scelte precedenti

% Scelta del sistema dinamico: prima perde la connessione 3->2, dopo 2->3 => NO CONSENSUS !!
% Simulazione fatta solo con modello linearizzato e rete: 1 -> 2 <-> 3
% 0: N0
% 1: Si
dynamic_choice = 0;

% Scelta del sistema robusto: il modello dei tre agenti non è perfetto
% Simulazione fatta solo con modello linearizzato e rete: 1 -> 2 <-> 3
% 0: No
% 1: Si
robustness_choice = 0;

%% ========================================================================

% Definizione delle Posizioni relative 
h1 = [0; 0; 0; 0];                % Leader
h2 = [1; 0; 0; 0];                % pos. relativa del nodo 2 rispetto il nodo 1
h3 = [0.5; 0; sqrt(3)/2; 0];      % pos. relativa del nodo 3 rispetto il nodo 1

% Dinamica del sistema LINEARIZZATO: doppio integratore
A = [0 1 0 0;
     0 0 0 0;
     0 0 0 1;
     0 0 0 0];

B = [0 0;
     1 0;
     0 0;
     0 1];

% Feasibility Conditions
tol = 1e-10;
if all(abs(A*(h3-h1)) <= tol) && all(abs(A*(h2-h1)) <= tol)
    fprintf('Verifica superata: La formazione è fattibile.\n');
else
    error('Errore: Condizioni di fattibilità non rispettate.\n');
end

% Scelta dinamica delle matrici in base alla tipo di connessione scelta
if connection_choice == 1           % 1 -> 2 -> 3
    D_Graph = diag([0 3 2]);
    A_adj = [0 0 0;     
             3 0 0;
             0 2 0];
elseif connection_choice == 2       % 3 <- 1 -> 2
    D_Graph = diag([0 3 2]);
    A_adj = [0 0 0;     
             3 0 0;
             2 0 0];
end

% Scelta dinamica delle matrici in base alla scelta Robusta
if robustness_choice == 1
    D_Graph = diag([0 3 2]);
    A_adj = [0 0 0;     
             3 0 0;
             0 2 0];
end

% Scelta dinamica delle matrici in base alla scelta Dinamica
if dynamic_choice == 1
    D_Graph = diag([0 5 2]);
    A_adj = [0 0 0;     
             3 0 2;
             0 2 0];
end

L = D_Graph - A_adj;

% Scelta dinamica del valore della saturazione
if u_choice == 0            % No Saturazione
    u_max = realmax;
    u_min = -realmax;
elseif u_choice == 1        % Si Saturazione 
    if traj_choice ~= 1
        u_max = 2;
        u_min = -2;
    else %devo usare valori leggermente più alti a causa di problemi numerici con traiettoria rettilinea
        u_max  = 2.5;
        u_min = -2.5;
    end
end

% STATE-FEEDBACK PROTOCOL

% Parametri QR
weightQ = 10;
Q = weightQ*diag([1 1 1 1]);    
R = diag([1 1]);    

PI1 = diag([1 0 0]);

% Calcolo Guadagni (c, K)
P_SVFB = are(A, B*inv(R)*B', Q);

K = inv(R)*B'*P_SVFB;  

eig_L_PI1 = eig(L+PI1);
c_min = 1/(2*min(real(eig_L_PI1))); 
c = c_min*3;    

% Calcolo (cl, Kl)
if gain_choice == 1         % forza uguale
    cl = c;
elseif gain_choice == 2     % tracking forte
    cl = 1.5*c;
end
K_l = cl*K;

% CONDIZIONI INIZIALI

% Condizioni casuali (scommentare per utilizzarle e commentare le
% asseganzioni fatte dopo)

% min_pos = -5;
% max_pos = 5;
% 
% p1x_0 = min_pos + (max_pos - min_pos) * rand();
% p1y_0 = min_pos + (max_pos - min_pos) * rand();
% 
% p2x_0 = min_pos + (max_pos - min_pos) * rand();
% p2y_0 = min_pos + (max_pos - min_pos) * rand();
% 
% p3x_0 = min_pos + (max_pos - min_pos) * rand();
% p3y_0 = min_pos + (max_pos - min_pos) * rand();

% x1_0 = [p1x_0; 0; p1y_0; 0];
% x2_0 = [p2x_0; 0; p2y_0; 0];
% x3_0 = [p3x_0; 0; p3y_0; 0];

% Condizioni iniziali usate per i grafici del report
x1_0 = [3; 0; 1; 0];
x2_0 = [-4; 0; 7; 0];
x3_0 = [-2; 0; -4; 0];

%% SIMULAZIONE
Tsim = 100; 

fprintf('Avvio Simulazione ... ');
if dynamic_choice == 1                                      % Con interruziuone dei collegamenti
    fprintf('del caso Dinamico\n');
    % open('Simulink_DINAMICO.slx');
    out = sim('Simulink_DINAMICO.slx', 'StopTime', num2str(Tsim));

elseif robustness_choice == 1                               % Incertezza sul modello
    fprintf('del caso Robusto\n');
    % open('Simulink_ROBUSTO.slx');
    out = sim('Simulink_ROBUSTO.slx', 'StopTime', num2str(Tsim));

elseif model_choice == 1 && connection_choice == 1          % Lineare && 1 -> 2 -> 3
    fprintf('del modello lineare e connessione 1 -> 2 -> 3\n');
    % open('Simulink1.slx');
    out = sim('Simulink1.slx', 'StopTime', num2str(Tsim));

elseif model_choice == 1 && connection_choice == 2          % Lineare && 3 <- 1 -> 2
    fprintf('del modello lineare e connessione 3 <- 1 -> 2\n');
    % open('Simulink2.slx');
    out = sim('Simulink2.slx', 'StopTime', num2str(Tsim));

elseif model_choice == 2 && connection_choice == 1          % NON Lineare && 1 -> 2 -> 3
    fprintf('del modello non lineare e connessione 1 -> 2 -> 3\n');
    % open('Simulink_NL1.slx');
    out = sim('Simulink_NL1.slx', 'StopTime', num2str(Tsim));

elseif model_choice == 2 && connection_choice == 2          % NON Lineare && 3 <- 1 -> 2
    fprintf('del modello non lineare e connessione 3 <- 1 -> 2\n');
    % open('Simulink_NL2.slx');
    out = sim('Simulink_NL2.slx', 'StopTime', num2str(Tsim));
end
fprintf('Simulazione terminata\n');

t = out.tout; 

p1x = squeeze(out.out_p1x); p2x = squeeze(out.out_p2x); p3x = squeeze(out.out_p3x);
p1y = squeeze(out.out_p1y); p2y = squeeze(out.out_p2y); p3y = squeeze(out.out_p3y);

p1x = p1x(:); p2x = p2x(:); p3x = p3x(:);
p1y = p1y(:); p2y = p2y(:); p3y = p3y(:);

px = [p1x, p2x, p3x];
py = [p1y, p2y, p3y];

% Estrazione degli errori di tracking e degli ingressi fisici
err_data = squeeze(out.out_err_track);
u1_data  = squeeze(out.out_u1);
u2_data  = squeeze(out.out_u2);
u3_data  = squeeze(out.out_u3);

% Controllo finale sull'orientamento delle matrici rispetto al vettore tempo
if size(px, 1) ~= length(t), px = px'; end
if size(py, 1) ~= length(t), py = py'; end
if size(err_data, 1) ~= length(t), err_data = err_data'; end
if size(u1_data, 1) ~= length(t), u1_data = u1_data'; end
if size(u2_data, 1) ~= length(t), u2_data = u2_data'; end
if size(u3_data, 1) ~= length(t), u3_data = u3_data'; end

%% Plot e animazione

% Array dei colori: 1 = Rosso (Leader), 2 = Verde (Robot2), 3 = Blu (Robot3)
colors = ["r", "g", "b"]; 

% ----- GRAFICO 1: TRAIETTORIE NEL PIANO XY -----
figure('Name', 'Traiettorie nel Piano XY', 'NumberTitle', 'off');
hold on; grid on;
for i = 1:3
    plot(px(:, i), py(:, i), colors(i), 'LineWidth', 2, 'DisplayName', ['Robot ', num2str(i)]);
    plot(px(1, i), py(1, i), strcat(colors(i), 'o'), 'MarkerFaceColor', colors(i), 'MarkerSize', 8, 'HandleVisibility', 'off');
    plot(px(end, i), py(end, i), strcat(colors(i), 'x'), 'LineWidth', 3, 'MarkerSize', 12, 'HandleVisibility', 'off');
end

istanti_formazione = [1, length(t)]; 
for idx = istanti_formazione
    tx = [px(idx, 1), px(idx, 2), px(idx, 3), px(idx, 1)];
    ty = [py(idx, 1), py(idx, 2), py(idx, 3), py(idx, 1)];
    if idx == 1
        plot(tx, ty, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Formazione Iniziale');
    else
        plot(tx, ty, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Formazione Finale');
    end
end
xlabel('Posizione X [m]'); ylabel('Posizione Y [m]');
title('Traiettoria Multi-Agente nel Piano'); legend('Location', 'best'); axis equal;

% ---- GRAFICO 2: ERRORE DI TRACKING DEL LEADER (ROBOT 1 - ROSSO) -----
figure('Name', 'Errore di Tracking del Leader', 'NumberTitle', 'off');
subplot(2,1,1);
plot(t, err_data(:, 1), 'r', 'LineWidth', 1.5); hold on; 
plot(t, err_data(:, 3), 'm', 'LineWidth', 1.5); grid on;  
xlabel('Tempo [s]'); ylabel('Errore di Posizione [m]');
title('Errore di Inseguimento del Riferimento (Robot 1 - Leader)'); legend('Errore X', 'Errore Y');

subplot(2,1,2);
plot(t, err_data(:, 2), 'r--', 'LineWidth', 1.5); hold on; 
plot(t, err_data(:, 4), 'm--', 'LineWidth', 1.5); grid on;  
xlabel('Tempo [s]'); ylabel('Errore di Velocità [m/s]'); legend('Errore v_x', 'Errore v_y');

% ----- GRAFICO 3: INGRESSI FISICI REALI (PHYSICAL INPUTS u_i) -----
figure('Name', 'Ingressi Fisici Reali degli Attuatori', 'NumberTitle', 'off');
subplot(2,1,1);
plot(t, u1_data(:, 1), 'r', 'LineWidth', 1.5, 'DisplayName', 'U1_x (Leader)'); hold on;
plot(t, u2_data(:, 1), 'g', 'LineWidth', 1.5, 'DisplayName', 'U2_x');
plot(t, u3_data(:, 1), 'b', 'LineWidth', 1.5, 'DisplayName', 'U3_x');
grid on; xlabel('Tempo [s]'); ylabel('Ingresso Fisico u_x');
title('Ingressi Fisici Reali (Forze/Comandi Attuatori lungo X)'); legend('Location', 'best');

subplot(2,1,2);
plot(t, u1_data(:, 2), 'r--', 'LineWidth', 1.5, 'DisplayName', 'U1_y (Leader)'); hold on;
plot(t, u2_data(:, 2), 'g--', 'LineWidth', 1.5, 'DisplayName', 'U2_y');
plot(t, u3_data(:, 2), 'b--', 'LineWidth', 1.5, 'DisplayName', 'U3_y');
grid on; xlabel('Tempo [s]'); ylabel('Ingresso Fisico u_y');
title('Ingressi Fisici Reali (Forze/Comandi Attuatori lungo Y)'); legend('Location', 'best');

% ----- GRAFICO 4: ERROIRI DI POSIZIONE RELATIVA TRA I NODI -----

% Estrazione degli offset desiderati di posizione
h1_pos = [h1(1); h1(3)];
h2_pos = [h2(1); h2(3)];
h3_pos = [h3(1); h3(3)];

% Calcolo degli errori di posizione relativi reali rispetto al target

% Errore Nodo 1 - Nodo 2
err_p12_x = (px(:,1) - px(:,2)) - (h1_pos(1) - h2_pos(1));
err_p12_y = (py(:,1) - py(:,2)) - (h1_pos(2) - h2_pos(2));

% Errore Nodo 1 - Nodo 3
err_p13_x = (px(:,1) - px(:,3)) - (h1_pos(1) - h3_pos(1));
err_p13_y = (py(:,1) - py(:,3)) - (h1_pos(2) - h3_pos(2));

% Errore Nodo 2 - Nodo 3
err_p23_x = (px(:,2) - px(:,3)) - (h2_pos(1) - h3_pos(1));
err_p23_y = (py(:,2) - py(:,3)) - (h2_pos(2) - h3_pos(2));

figure('Name', 'Errori di Formazione Relativi tra i Nodi', 'NumberTitle', 'off');

% Subplot per gli errori lungo l'asse X
subplot(2,1,1);
plot(t, err_p12_x, 'g', 'LineWidth', 1.5, 'DisplayName', 'Errore 1-2 (X)'); hold on;
plot(t, err_p13_x, 'b', 'LineWidth', 1.5, 'DisplayName', 'Errore 1-3 (X)');
plot(t, err_p23_x, 'm', 'LineWidth', 1.5, 'DisplayName', 'Errore 2-3 (X)');
grid on; xlim([0 t(end)]);
xlabel('Tempo [s]'); ylabel('Errore [m]');
title('Errori di Formazione Relativi lungo l''asse X');
legend('Location', 'best');

% Subplot per gli errori lungo l'asse Y
subplot(2,1,2);
plot(t, err_p12_y, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Errore 1-2 (Y)'); hold on;
plot(t, err_p13_y, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Errore 1-3 (Y)');
plot(t, err_p23_y, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Errore 2-3 (Y)');
grid on; xlim([0 t(end)]);
xlabel('Tempo [s]'); ylabel('Errore [m]');
title('Errori di Formazione Relativi lungo l''asse Y');
legend('Location', 'best');

% =========================================================================
% ANIMAZIONE SPAZIO-TEMPORALE 3D DELLA FORMAZIONE
% =========================================================================
fig_anim = figure('Name', 'Animazione 3D Rallentata', 'NumberTitle', 'off');
hold on; grid on;
xlabel('Posizione X [m]'); ylabel('Posizione Y [m]'); zlabel('Tempo [s]');
title('Evoluzione della Formazione'); view(3); 
axis([min(px(:))-2, max(px(:))+2, min(py(:))-2, max(py(:))+2, 0, t(end)]);

h_robot1 = plot3(px(1,1), py(1,1), t(1), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
h_robot2 = plot3(px(1,2), py(1,2), t(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 6);
h_robot3 = plot3(px(1,3), py(1,3), t(1), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
h_triangolo = plot3([px(1,1) px(1,2) px(1,3) px(1,1)], ...
                    [py(1,1) py(1,2) py(1,3) py(1,1)], ...
                    [t(1) t(1) t(1) t(1)], 'k-', 'LineWidth', 1.5);

if dynamic_choice == 1
    % Flag per evitare di stampare il testo ripetutamente
    msg20_stampato = false;
    msg35_stampato = false;
    
    % Coordinata spaziale fissa in cui far apparire i testi
    x_text_pos = min(px(:));
    y_text_pos = max(py(:));
end

passo = 5; 
fprintf('\nAvvio animazione... \n');
for k = (passo + 1):passo:length(t)

    if dynamic_choice ~= 1
        view(0,90);
    end
    idx_inizio = k - passo;
    
    if ~ishandle(fig_anim) || ~isvalid(h_robot1)
        fprintf('Animazione chiusa dall''utente.\n');
        break;
    end
    
    plot3(px(idx_inizio:k, 1), py(idx_inizio:k, 1), t(idx_inizio:k), 'r-', 'LineWidth', 1.5, 'HandleVisibility','off');
    plot3(px(idx_inizio:k, 2), py(idx_inizio:k, 2), t(idx_inizio:k), 'g-', 'LineWidth', 1.5, 'HandleVisibility','off');
    plot3(px(idx_inizio:k, 3), py(idx_inizio:k, 3), t(idx_inizio:k), 'b-', 'LineWidth', 1.5, 'HandleVisibility','off');
    
    set(h_robot1, 'XData', px(k, 1), 'YData', py(k, 1), 'ZData', t(k));
    set(h_robot2, 'XData', px(k, 2), 'YData', py(k, 2), 'ZData', t(k));
    set(h_robot3, 'XData', px(k, 3), 'YData', py(k, 3), 'ZData', t(k));
    set(h_triangolo, 'XData', [px(k,1) px(k,2) px(k,3) px(k,1)], ...
                     'YData', [py(k,1) py(k,2) py(k,3) py(k,1)], ...
                     'ZData', [t(k) t(k) t(k) t(k)]);
    
    % --- LOGICA PER L'INSERIMENTO DEI TESTI AI TEMPI SEGNALATI ---
    if dynamic_choice == 1
        if t(k) >= 20 && ~msg20_stampato
            text(x_text_pos, y_text_pos, 20, ' \leftarrow connessione 3->2 disattivata', ...
                'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
            msg20_stampato = true;
        end
    
        if t(k) >= 35 && ~msg35_stampato
            text(x_text_pos, y_text_pos, 35, ' \leftarrow connessione 2->3 disattivata, no albero ricoprente', ...
                'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
            msg35_stampato = true;
        end
    end
    % -------------------------------------------------------------
    
    drawnow;
    pause(0.02); 
end

if ishandle(fig_anim) && k < length(t)
    plot3(px(k:end, 1), py(k:end, 1), t(k:end), 'r-', 'LineWidth', 1.5, 'HandleVisibility','off');
    plot3(px(k:end, 2), py(k:end, 2), t(k:end), 'g-', 'LineWidth', 1.5, 'HandleVisibility','off');
    plot3(px(k:end, 3), py(k:end, 3), t(k:end), 'b-', 'LineWidth', 1.5, 'HandleVisibility','off');
    drawnow;
    fprintf('Animazione completata.\n');
end