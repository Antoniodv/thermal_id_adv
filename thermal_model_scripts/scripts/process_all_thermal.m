function [global_matrix, active_cores_indexes_matrix] = process_all_thermal(root_path, path_to_active_cores, mode)    
    if strcmp(mode, "all") || strcmp(mode, "patch")
        disp("mode: ")
        disp(mode)
    else
        error(['SELECTED MODE NOT ALLOWED:  ', mode]);
    end

    addpath("adv_functions\");
    
    if ~isfolder(root_path)
        error('The specified root path does not exist.');
    end
    
    % Trova tutte le cartelle 'power_model'
    power_model_folders = dir(fullfile(root_path, '**', 'power_model'));

    % Raccogli i percorsi univoci
    unique_paths = unique({power_model_folders.folder});
    
    disp("Cartelle trovate in MATLAB:");
    for k = 1:length(unique_paths)
        disp(unique_paths{k});
    end    

    for k = 1:length(unique_paths)
        disp(k)
        try
            power_model_path = unique_paths{k}; % Percorso corretto della cartella 'power_model'
            parent_path = fileparts(power_model_path); % Cartella padre di 'power_model'
            
            % Definisce la cartella di output 'thermal_model'
            thermal_model_path = fullfile(parent_path, 'thermal_model');
            if ~exist(thermal_model_path, 'dir')
                mkdir(thermal_model_path);
            end
            
            % Trova il file .pkl dentro 'power_model'
            pkl_file = fullfile(power_model_path, 'gb_core_uncore_tot_temp.pkl');
            
            % Controlla se il file esiste prima di procedere
            if exist(pkl_file, 'file')
                if strcmp(mode, "all")
                    [tmp_coeff_matrix, active_cores_indexes] = thermal_model_est(pkl_file, "./", thermal_model_path, mode);
                    global_matrix(:, :, k) = tmp_coeff_matrix;
                    active_cores_indexes_matrix(k,:) = active_cores_indexes;

                elseif strcmp(mode, "patch")    
                    disp("patch")
                    [tmp_coeff_matrix, active_cores_indexes] = thermal_model_est(pkl_file, path_to_active_cores, thermal_model_path, mode);
                    global_matrix(:, :, k) = tmp_coeff_matrix;
                    active_cores_indexes_matrix(k,:) = active_cores_indexes;
                    
                end
            else
                warning(['File non trovato: ', pkl_file]);
            end
        catch ME
            disp(['Errore elaborando la cartella: ', unique_paths{k}]);
            disp(['Messaggio errore: ', ME.message]);

            % Stampa dettagli della traccia dello stack
            for i = 1:length(ME.stack)
                disp(['Errore in file: ', ME.stack(i).file]);
                disp(['Funzione: ', ME.stack(i).name]);
                disp(['Riga: ', num2str(ME.stack(i).line)]);
            end
        end
    end

% % HEATMAP GLOBALE MEDIATA CON ELABORAZIONE COEFFICIENTI
% global_matrix_sum = sum(global_matrix, 3);
% global_matrix_avg = global_matrix_sum ./ 6; 
% 
% % exogen 
% exog_matrix = global_matrix_avg(:,3:end);
% 
% % media esogeni a coppie
% b1 = exog_matrix(:, 1:2:end);  % colonne dispari
% b2 = exog_matrix(:, 2:2:end);  % colonne pari
% b_mat = (b1 + b2) / 2;
% 
% figure('Position', [100, 100, 800, 600]);
% imagesc(b_mat); 
% hold on; 
% title("Global coefficients matrix " + mode);
% yticks(1:size(b_mat, 1)); 
% [m, n] = size(b_mat);
% plot(1:min(m,n), 1:min(m,n), 'w-', 'LineWidth', 1); 
% colorbar;
% colormap(flipud(autumn));
end


