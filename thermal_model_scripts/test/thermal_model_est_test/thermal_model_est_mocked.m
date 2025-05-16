function [theta_matrix, active_cores_indexes , U_active_cores, round_index_matlab] = thermal_model_est_mocked(path_to_metrics, path_to_active_cores, output_folder, mode)

DBG = 0;

disp("    ")
disp("start processing")
disp(path_to_metrics)
disp("deliv folder")
disp(output_folder)
disp("path_to_active_cores")
disp(path_to_active_cores)
disp("                  ")


py.importlib.import_module('pandas');
py.importlib.import_module('numpy'); 


N_CORE = 56;
thrm = THID(2,N_CORE);
thrm.Nmc = 4;
thrm.Nmp = 2;
thrm.test_sampling_time = 0.5;

TEST_TYPE = 1; % Gathering the PRBS data (only one available)
dfmat_idx = 1; %pkl
dtype_idx = 2; %list of dictionaries 1, dict 2

path_to_pm_file = string(path_to_metrics);

thrm.model_order = 2;
thrm.num_add_eq = 10;

if exist(fullfile(path_to_metrics, 'gb_core_uncore_tot_temp.pkl'), 'file') ~= 2
    error('The path does not bring to the right pkl file');
end


% add algorithms functions from Diversi
addpath("Diversi/");

if strcmp(mode, "patch")
    round_str = regexp(path_to_pm_file, 'round(\d+)', 'tokens');
    if isempty(round_str)
        error('No "round" word in this path');
    end
    
    % Estrae il numero come stringa e poi lo converte in numero
    round_index = (str2double(round_str{1}{1}));    
    round_index_matlab = round_index + 1;
    disp("round_index_matlab")
    disp(round_index_matlab)
    active_cores_indexes = sort(get_row_from_file(round_index_matlab, path_to_active_cores));
    disp("cores read from file")
    disp(active_cores_indexes)

    % from python/C-like indexes to MATLAB indexes
    active_cores_indexes = active_cores_indexes+1;
elseif strcmp(mode, "all")
    active_cores_indexes = (1:N_CORE);
end

n_active_cores = length(active_cores_indexes);

if DBG == 1
    % HEATMAP MATRICE POTENZE COMPLESSIVA DBG----------------------------------
    % Creazione della hean_active_corestmap
    figure('Position', [100, 100, 800, 600]); % Imposta una finestra più grande
    
    imagesc(U'); % Plotta la heatmap
    title(["round" + round_index]);
    
    
    % Imposta i tick su tutti i valori possibili dell'asse Y
    yticks(1:size(U, 1)); 
    ylabel('Indice Y');
    
    % Aggiunge una barra dei colori per maggiore leggibilità
    colorbar;
    colormap(jet); % Imposta una mappa colori (opzionale)
    % -------------------------------------------------------------------------
end

M = thrm.parser_script(0, 0,[], [], dfmat_idx, dtype_idx, path_to_pm_file);

% extraction of Power and Temperature
U = M(:,1:thrm.num_cores); % Cores Power
tt = (thrm.num_cores)+(1:2);
PP = [M(:,tt)]; %CPU Power
tt = (thrm.num_cores+3)+(1:thrm.num_cores);
T = M(:,tt); %Cores Temperature

% allign power and temperature according to model
U = U(2:end,:);
T = T(1:end-1,:);

% name change
Y = T;
cpu=1;

% reduce power matrix to the matrix of active cores
U_active_cores = U(:, active_cores_indexes);
dict_list = py.list();
theta_matrix = zeros(N_CORE, 2+2*N_CORE);
disp('active_cores')
disp(active_cores_indexes)
for core = 1:thrm.num_cores
     if strcmp(mode, "all")
        % Find thermal model
        [THETA_TH_vector, ~, ~, ~, ~, ~] = ...
            misoarxbls_mocked(U_active_cores,Y(:,core),thrm.model_order,thrm.num_add_eq);
        
        THETA_TH_vector = THETA_TH_vector';
        theta_matrix(core, :) = THETA_TH_vector;
        THETA_TH_vector_deliv = THETA_TH_vector;

        figure;
        data = THETA_TH_vector;
        bar(data)
        set(gcf, 'Position', [100, 100, 800, 600])
        title('theta vector debug')
        xlabel('Cores id')
        xticks(1:N_CORE)                 
        xticklabels(1:N_CORE)

        
     elseif strcmp(mode, "patch") && ismember(core, active_cores_indexes)
        
        % Find thermal model
        [THETA_TH_vector, ~, ~,~,~,~] = ...
            misoarxbls_mocked(U_active_cores,Y(:,core),thrm.model_order,thrm.num_add_eq);
    
        disp("theta_vector size before transponse: ")
        disp(size(THETA_TH_vector))
        THETA_TH_vector = THETA_TH_vector';
        
        theta_vector_remapped = remap_theta_vector(THETA_TH_vector, thrm.num_cores, active_cores_indexes);
        theta_matrix(core, :) = theta_vector_remapped;

    end
        
end

end


