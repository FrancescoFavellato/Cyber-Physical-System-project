function ris_opam = OPAM(A, C, y, a, x_true, lambda, mu, T)
    % OUTPUT:
    % ris_opam  - struct con campi:
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
    
    for k = 1:T
        
        yk = y(:, k);               % misura y(k) al tempo corrente
        tic
        CAk = C * Apow; % C * A^k

        % Passo di aggiornamento su x0_hat
        x0_hat = (CAk'*CAk + mu*eye(n)) \ (mu*x0_hat + CAk'*(yk - a_hat));
        
        % Passo di aggiornamento su a_hat
        z = yk - CAk*x0_hat + mu*a_hat;
        a_hat = (1/(1+mu))*soft_threshold(z, lambda);

        % Stima dello stato corrente al tempo k
        x_hat(:, k) = Apow * x0_hat;
        
        % Aggiornamento di Apow per il prossimo passo
        Apow = A * Apow;

        tempo_calcoli = toc;
        tempo_complessivo = tempo_complessivo + tempo_calcoli;
        
        % === Calcolo delle metriche ===
        x_true_k = x_true(:, k);
        state_err(k) = norm(x_hat(:, k) - x_true_k)^2 / n;
        
        % Supporto stimato
        supp_hat = (a_hat ~= 0);
        diff_supp = double(supp_hat) - double(supp_true);
        
        % Support attack error e tassi di errore
        support_err(k) = sum(abs(diff_supp));
        false_pos(k) = sum(max(diff_supp, 0));
        false_neg(k) = sum(max(-diff_supp, 0));
        
    end
    
    % Output
    ris_opam.x0_hat = x0_hat;
    ris_opam.a_hat = a_hat;
    ris_opam.x_hat = x_hat;
    ris_opam.state_err = state_err;
    ris_opam.support_err = support_err;
    ris_opam.false_pos = false_pos;
    ris_opam.false_neg = false_neg;
    ris_opam.time = tempo_complessivo;
end
