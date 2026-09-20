function ris_aopgd = AOPGD(A, C, y, a, x_true, lambda, nu, T)
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

    supp_true = (a ~= 0);       % supporto vero logical
    Apow = eye(n);
    tempo_complessivo = 0;

    % Inizializzazione degli accumulatori dei singoli termini della sommatoria
    sum_CAj_CAj = zeros(n, n);
    sum_CAj     = zeros(n, q);
    sum_CAj_yj  = zeros(n, 1);
    sum_yj      = zeros(q, 1);

    for idx = 1:T               % idx = [1, T]   
        k = idx - 1;            % k = [0, T-1] 
        yj = y(:, idx);         % misura corrente al tempo j
        
        tic;                    % Inizio a calcolare il tempo
        
        CAj = C * Apow;         % Calcolo di C * A^j per il j corrente 
        
        % Aggiornamento delle sommatorie (da j = 0 a k)
        sum_CAj_CAj = sum_CAj_CAj + CAj' * CAj;
        sum_CAj     = sum_CAj + CAj';
        sum_CAj_yj  = sum_CAj_yj + CAj' * yj;
        sum_yj      = sum_yj + yj;
        
        % Calcolo dei gradienti aggregati valutati in x0_hat e a_hat
        grad_x = (sum_CAj_CAj * x0_hat + sum_CAj * a_hat - sum_CAj_yj) / (k + 1);
        grad_a = (sum_CAj' * x0_hat + (k + 1) * a_hat - sum_yj) / (k + 1);
        
        % Aggiornamento delle stime
        x0_hat = x0_hat - nu * grad_x;

        a_temp = a_hat - nu * grad_a;
        a_hat  = soft_threshold(a_temp, nu * lambda);
    
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
    ris_aopgd.x0_hat = x0_hat;
    ris_aopgd.a_hat = a_hat;
    ris_aopgd.x_hat = x_hat;
    ris_aopgd.state_err = state_err;
    ris_aopgd.support_err = support_err;
    ris_aopgd.false_pos = false_pos;
    ris_aopgd.false_neg = false_neg;
    ris_aopgd.time = tempo_complessivo;
end