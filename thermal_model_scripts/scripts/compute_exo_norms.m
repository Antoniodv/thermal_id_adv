function exo_norms = compute_exo_norms(global_matrix_avg)
    % Inizializza la matrice delle norme
    exo_norms = zeros(size(global_matrix_avg, 1), size(global_matrix_avg, 2)/2);
    
    % Calcola la norma per ogni riga (core)
    for line_idx = 1:size(global_matrix_avg, 1)
        % Prendi la riga corrente
        current_row = global_matrix_avg(line_idx, :);
        
        % Dividi in coppie
        b1 = current_row(1:2:end);   % primi di ogni coppia (dispari)
        b2 = current_row(2:2:end);   % secondi di ogni coppia (pari)
        
        % Calcola la norma di ogni coppia
        exo_norms(line_idx, :) = sqrt(b1.^2 + b2.^2);
    end
end