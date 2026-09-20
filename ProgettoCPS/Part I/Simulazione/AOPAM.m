function ris_aopam = AOPAM(A, C, y, a, x_true, lambda, mu, T)
    % OUTPUT:
    % ris_aopgd  - struct con campi:
    %     .x0_hat      stima finale di x(0)             [50 x 1]
    %     .a_hat       stima finale dell'attacco a      [20 x 1]
    %     .x_hat       stime stato corrente             [50 x T]
    %     .state_err   state estimation error           [1  x T]
    %     .support_err support attack error             [1  x T]
    %     .false_pos   false positives                  [1  x T]
    %     .false_neg   false negatives                  [1  x T]
    %     .time        tempo esecuzione                 [1  x 1]

    n = size(A, 1);
    q = size(C, 1);
    
    % Inizializzazione: x_hat(0) = 0, a_hat(0) = 0
    x0_hat = zeros(n, 1);
    a_hat = zeros(q, 1);
    x_hat = zeros(n, T);

    state_err = zeros(1, T);
    support_err = zeros(1, T);
    false_pos = zeros(1, T);
    false_neg = zeros(1, T);

    supp_true = (a ~= 0);           % supporto vero logical
    Apow = eye(n);
    tempo_complessivo = 0;

    % Inizializzazione degli accumulatori dei singoli termini della sommatoria
    sum_CAj_CAj = zeros(n, n);
    sum_CAjt    = zeros(n, q);
    sum_CAj_yj  = zeros(n, 1);
    sum_yj      = zeros(q, 1);

    for idx = 1:T               % idx = [1, T]   
        k = idx - 1;            % k = [0, T-1] 
        yj = y(:, idx);         % misura corrente al tempo j
        
        tic;                    % Inizio a calcolare il tempo
        
        CAj = C * Apow;         % Calcolo di C * A^j per il j corrente 
        
        % Aggiornamento delle sommatorie (da j = 0 a k)
        sum_CAj_CAj = sum_CAj_CAj + CAj' * CAj;
        sum_CAjt    = sum_CAjt + CAj';
        sum_CAj_yj  = sum_CAj_yj + CAj' * yj;
        sum_yj      = sum_yj + yj;

        inv_k1 = 1 / (k + 1);
        
        % Aggiornamento della stima dello stato iniziale
        LHS_x = inv_k1 * sum_CAj_CAj + mu * eye(n);
        RHS_x = mu * x0_hat + inv_k1 * (sum_CAj_yj - sum_CAjt * a_hat);
        x0_hat = LHS_x \ RHS_x;
        
        % Aggiornamento della stima dell'attacco a_temp
        a_temp = (1 / (1 + mu)) * (mu * a_hat + inv_k1 * (sum_yj - sum_CAjt' * x0_hat));
        a_hat  = soft_threshold(a_temp, lambda / (1 + mu));
    
        % Stima dello stato corrente al tempo k
        x_hat(:, idx) = Apow * x0_hat;
        
        % Aggiornamento della potenza di A per il prossimo passo
        Apow = A * Apow;
        
        tempo_calcoli = toc;
        tempo_complessivo = tempo_complessivo + tempo_calcoli;
    
        % === Calcolo Metriche ===
        x_true_k = x_true(:, idx);    
        state_err(idx) = norm(x_hat(:, idx) - x_true_k)^2 / n;
     
        supp_hat  = (a_hat ~= 0);       
        diff_supp = double(supp_hat) - double(supp_true);   
     
        support_err(idx) = sum(abs(diff_supp));   
        false_pos(idx)   = sum(diff_supp == 1);
        false_neg(idx)   = sum(diff_supp == -1); 
    end
    
    % Output
    ris_aopam.x0_hat = x0_hat;
    ris_aopam.a_hat = a_hat;
    ris_aopam.x_hat = x_hat;
    ris_aopam.state_err = state_err;
    ris_aopam.support_err = support_err;
    ris_aopam.false_pos = false_pos;
    ris_aopam.false_neg = false_neg;
    ris_aopam.time = tempo_complessivo;
end