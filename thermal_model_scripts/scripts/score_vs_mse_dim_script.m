close all;
clear;

mode = "patch";
plot_dir_path = './score_vs_mse_plot';

% APR5
path_to_pw_temp = "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification\APR5\prbs_random\112\";

% % FEB5
% path_to_pw_temp = "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification\prbs_random_FEB5\112\round0\power_model\gb_core_uncore_tot_temp.pkl";

path_to_active_cores = "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification\active_cores_APR5.txt";
load('global_data.mat', 'b_mat', 'global_matrix', 'active_cores_indexes_matrix');
N_CORE = 56;

thrm = THID(2,N_CORE);
thrm.Nmc = 4;
thrm.Nmp = 2;
thrm.test_sampling_time = 0.5;
thrm.model_order = 2;
thrm.num_add_eq = 10;
TEST_TYPE = 1; % Gathering the PRBS data (only one available)
dfmat_idx = 1; %pkl
dtype_idx = 2; %list of dictionaries 1, dict 2
path_to_pm_file = path_to_pw_temp;



%% Evalutate coefficient distribution for each evaluated core (High to Low)
arg_sorted_b = zeros(size(b_mat));
sorted_b = zeros(size(b_mat));
for line_idx = 1:size(b_mat, 1)
    [sorted_b(line_idx, :), arg_sorted_b(line_idx, :)] = sort(b_mat(line_idx, :), 'descend');
end
sorted_influencers = arg_sorted_b(:, 1:end);



%% Best round for core
chosen_cores_array = (0:55);

for core_to_evaluate_idx = chosen_cores_array
    
    core_to_evaluate_index = core_to_evaluate_idx +1; % adapt to MATLAB index
    
    % find the round with most of the top influencers cores for chosen core 
    % from now on sorted_influencers only has cols of active cores in this
    % round !!!!
    [best_round, sorted_influencers_active] = find_best_round_idx(core_to_evaluate_index, active_cores_indexes_matrix, sorted_influencers);

    % EXTRACT NORMS OF ACTIVE CORES IN THE SINGLE ROUND
    theta_from_glb = global_matrix(core_to_evaluate_index, :, best_round);
    [~, exo_norms] = separate_coefficients_norm(theta_from_glb);

    % CORE CHOSING CRITERIA 
    
    % 1 sort based on global matrix score
    ordered_cores = sorted_influencers_active;
    % --------------------------------------------

    % 2 sort based on score in the chosen round
    % [ordered_norms, ordered_cores] = sort(exo_norms, 'descend'); 
    %  ordered_cores = ordered_cores(1:14); % altrimenti fai car. aggiungendo anche core spenti dopo i primi 14
    % -------------------------------------------- 

    % 3 random, both in cores choice and in their order
    % ordered_cores = randsample(N_CORE, 14);
    % -------------------------------------------- 
    
    % 4 neig under devel
    % ordered_cores = get_neighbor_ids(core_to_evaluate_index, 4, 4);
    % -------------------------------------------- 

%% Thermal Model Estimation and Inference for this core 

    % extraction of Power and Temperature
    i=1;
    input_configurations = {};
    
    for n_cores_for_evaluation = 1:length(ordered_cores)
        active_cores_nosort = ordered_cores(1:n_cores_for_evaluation);
        active_cores = sort(active_cores_nosort);
        
        % EXTRACT TEMPERATURE AND POWER
        % [U, T] = extract_power_temperature_from_pkl(path_to_pw_temp, num_cores, best_round);
        % U_active_cores = U(:, active_cores);
        % T_measured = T(:,core_to_evaluate_index);

        % PER APR5
        path_to_pm_file = strcat(path_to_pw_temp, 'round', int2str(best_round), '\', 'power_model', '\', 'gb_core_uncore_tot_temp.pkl');
        
        % % PER FEB5
        % path_to_pm_file = strcat(path_to_pw_temp); %strcat(path_to_pw_temp, 'round', int2str(best_round), '\', 'power_model', '\', 'gb_core_uncore_tot_temp.pkl');
        % 
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

        U_active_cores = U(:, active_cores);
        T_measured = T(:,core_to_evaluate_index);

        % COMPUTE THERMAL MODEL OF CHOSEN CORES
        [THETA_TH_vector, sigv,sigw, res, condn, svd_ident] = ...
            misoarxbls(U_active_cores, T_measured, thrm.model_order,  thrm.num_add_eq);

        % this vector has 2+2*n_active_cores elements only, not remapped N_CORE
        THETA_TH_vector = THETA_TH_vector';

        [INN, T_predicted, T_filtered] = infer_Kalman(THETA_TH_vector, T_measured, U_active_cores, sigv, sigw, thrm.model_order);

        % saving data in array for barplot
        mse(i) = var(INN(7:end));
        
        score_b_mat(i) = sum(b_mat(core_to_evaluate_index, active_cores));
        score_test(i) = sum(exo_norms(active_cores));
        input_configurations{end+1} = active_cores;      

        T_predicted = T_predicted(8:end);
        T_measured =  T_measured(7:end-1);


        % base dir for plots
        core_plot_dir = fullfile('./score_vs_mse_plot', num2str(core_to_evaluate_index));
        
        % create subfolders
        inference_dir = fullfile(core_plot_dir, 'inferences');
        mkdir_if_not_exist(inference_dir);  
        acf_dir = fullfile(core_plot_dir, 'acf');
        mkdir_if_not_exist(acf_dir);  
        res_dir = fullfile(core_plot_dir, 'res');
        mkdir_if_not_exist(res_dir);  
    
        % ===== 1. PLOT TEMPERATURE (inferences folder) =====
        fig_temp = figure('Visible', 'off', 'Position', [100, 100, 1200, 600]);
        
        % Plot
        h_meas = plot(T_measured, 'b', 'LineWidth', 1.8, 'DisplayName', 'Misurata');
        hold on;
        h_pred = plot(T_predicted, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Predetta');
        
        title_str = sprintf('Core %d - Config %d: %s\n', ...
               core_to_evaluate_index, i, arrayToString(active_cores));
        
        title(title_str, 'FontSize', 10, 'FontWeight', 'normal');
        xlabel('Campioni temporali', 'FontSize', 9);
        ylabel('Temperatura (°C)', 'FontSize', 9);
        legend([h_meas, h_pred], 'Location', 'best');
        grid on;
        set(gca, 'FontSize', 8);
        
        % Save PNG and FIG
        temp_filename_base = sprintf('temp_core%d_config%02d', core_to_evaluate_index, i);
        saveas(fig_temp, fullfile(inference_dir, [temp_filename_base '.png']));
        savefig(fig_temp, fullfile(inference_dir, [temp_filename_base '.fig']));
        close(fig_temp);
        
        
        % ===== 2. PLOT ACF =====     
        max_lag = 40;
        [acf, lags] = autocorr(res, 'NumLags', max_lag);
        
        fig_acf = figure('Visible', 'off', 'Position', [100, 100, 800, 400]);
        stem(lags, acf, 'filled', 'b');
        title_str = sprintf('ACF Core %d - Config %d: %s\n', ...
            core_to_evaluate_index, i, arrayToString(active_cores));
        
        title(title_str);
        xlabel('Lag', 'FontSize', 10);
        ylabel('ACF', 'FontSize', 10);
        grid on;
        set(gca, 'FontSize', 9);
        
        % Save PNG and FIG
        acf_filename_base = sprintf('acf_residuals_cycle%02d', i);
        saveas(fig_acf, fullfile(acf_dir, [acf_filename_base '.png']));
        savefig(fig_acf, fullfile(acf_dir, [acf_filename_base '.fig']));
        close(fig_acf);
        
        
        % ===== 3. PLOT RESIDUALS =====     
        fig_res = figure('Visible', 'off', 'Position', [100, 100, 1200, 400]);
        scatter(1:length(res), res, 25, 'k', 'filled');
        title_str = sprintf('Core %d - Config %d: %s\nResidui', ...
                core_to_evaluate_index, i, arrayToString(active_cores));
        title(title_str, 'FontSize', 10);
        xlabel('Samples', 'FontSize', 9);
        ylabel('Residuals', 'FontSize', 9);
        grid on;
        set(gca, 'FontSize', 8);
        
        % Save PNG and FIG
        res_filename_base = sprintf('residuals_core%d_config%02d', core_to_evaluate_index, i);
        saveas(fig_res, fullfile(res_dir, [res_filename_base '.png']));
        savefig(fig_res, fullfile(res_dir, [res_filename_base '.fig']));
        close(fig_res);



        i = i+1;  
            
    end
   
labels = cell(size(input_configurations));
for j = 1:length(input_configurations)
    labels{j} = mat2str(input_configurations{j});
end

% Percorso base per i plot
plot_dir_path = fullfile('./score_vs_mse_plot');
mkdir_if_not_exist(plot_dir_path);

% Prepara la figura
fig = figure('Visible', 'off', 'Position', [100, 100, 2000, 1400]);

% Passa all'asse Y sinistro prima di plottare
yyaxis left;
b = bar(1:length(mse), mse, 'FaceColor', 'b');
ylabel('MSE', 'FontSize', 10);

% Passa all'asse Y destro per plottare lo Score
yyaxis right;
hold on;
plot(1:length(score_b_mat), score_b_mat, '-or', 'LineWidth', 2);      
plot(1:length(score_test), score_test, '-sg', 'LineWidth', 2);              
ylabel('Score', 'FontSize', 10);

% Configura asse X
set(gca, 'XTick', 1:length(input_configurations));
set(gca, 'XTickLabelRotation', 45);
xlabel('Model Configurations', 'FontSize', 10);
xticklabels(string(1:length(input_configurations))); 

% Titolo e legenda
title_str = sprintf('MSE vs Score, eval. core: %d, %s', ...
                           core_to_evaluate_index, mat2str(active_cores_nosort));
title(title_str, 'FontSize', 12);
legend("MSE", "Score global", "Score test", 'Location', 'northeastoutside');

% Salva la figura
filename = sprintf('score_vs_mse_core_%d.png', core_to_evaluate_index);
fullpath = fullfile(plot_dir_path, filename);
disp("Saving combined MSE-Score plot to: " + fullpath);
saveas(fig, fullpath);
close(fig);

end














% FUNCTIONS 

function mkdir_if_not_exist(path)
    if ~exist(path, 'dir')
        mkdir(path);
    end
end

function reordered = reorder_by_reference(original, reference)
    [~, idx] = ismember(reference, original);
    valid_idx = idx > 0;
    reordered = original(idx(valid_idx));
end

function [best_round, higher_common_cores] = find_best_round_idx(core_to_evaluate_index, active_cores_indexes_matrix, sorted_influencers)
    best_round = -1;
    higher_common_cores = [];
    for line_idx = 1: size(active_cores_indexes_matrix, 1)
        %se il core in eval c'é allora si valuta quanti attivi ha questo round
        if ismember(core_to_evaluate_index, active_cores_indexes_matrix(line_idx, :)) 
            common_cores = intersect(active_cores_indexes_matrix(line_idx, :), sorted_influencers(core_to_evaluate_index, :));
            if best_round == -1
                if size(common_cores, 1) > 0
                    higher_common_cores = common_cores;
                    best_round = line_idx;
                end
            else
                if size(common_cores, 1) > size(higher_common_cores, 1)
                    higher_common_cores = common_cores;
                    best_round = line_idx;
                end
            end
        end
    end

    higher_common_cores = reorder_by_reference(higher_common_cores, sorted_influencers(core_to_evaluate_index, :));
    
    if best_round == -1
        disp("no best round found for core")
        disp(core_to_evaluate_index)
    end
end

function str = arrayToString(arr)
    % Crea una stringa con numeri separati da virgole
    str = strjoin(string(arr), ',');
end


function [Power, Temperature] = extract_power_temperature_from_pkl(path_to_pw_temp, num_cores, round)
    pickle = py.importlib.import_module('pickle');
    best_round = round;
    path_to_pw_temp_tmp = strcat(path_to_pw_temp, 'round', int2str(best_round), '\', 'power_model', '\', 'gb_core_uncore_tot_temp.pkl');
    gb_mat_handle = py.open(path_to_pw_temp_tmp, 'rb');
    gb_mat = pickle.load(gb_mat_handle);   
    gb_mat_handle.close();
    gb_mat = double(gb_mat.values);
    
    M = gb_mat;
    U = M(:,1:num_cores); % Cores Power
    tt = (num_cores)+(1:2);
    PP = [M(:,tt)]; %CPU Power
    tt = (num_cores+3)+(1:num_cores);
    T = M(:,tt); %Cores Temperature
    
    % allign power and temperature according to model
    Power = U(2:end,:);
    Temperature = T(1:end-1,:);
    
    end 


function [ar, exo] = separate_coefficients(input_coefficients)
    %ar
    ar = input_coefficients(1:2);
    
    % exogen 
    exo = input_coefficients(3:end);
end

function theta_subset = extract_theta_subset(theta_full, active_cores, N_CORE)
    % Primi 2 coefficienti AR
    theta_subset = theta_full(1:2);

    % Ogni core ha 2 coefficienti (es. delay 1 e 2)
    % Quindi la parte esogena inizia da indice 3
    for c = active_cores
        idx_start = 2 + (c - 1) * 2 + 1;
        idx_end = idx_start + 1;
        theta_subset = [theta_subset, theta_full(idx_start:idx_end)];
    end
end

% functions to move in the grid
function neighbors = get_neighbor_ids(core_id, n_rows, n_cols)
 
    [x_core, y_core] = id_to_coordinates(core_id, n_rows, n_cols);

    [x_north, y_north] = find_north_coordinate(x_core, y_core, n_rows, n_cols);
    id_north = coordinates_to_id(x_north, y_north, n_rows, n_cols);

    [x_south, y_south] = find_south_coordinate(x_core, y_core, n_rows, n_cols);
    id_south = coordinates_to_id(x_south, y_south, n_rows, n_cols);

    [x_east, y_east] = find_east_coordinate(x_core, y_core, n_rows, n_cols);
    id_east = coordinates_to_id(x_east, y_east, n_rows, n_cols);

    [x_west, y_west] = find_west_coordinate(x_core, y_core, n_rows, n_cols);
    id_west = coordinates_to_id(x_west, y_west, n_rows, n_cols);

    neighbors = [ id_north, id_south, id_east, id_west ]; 
    neighbors = neighbors(neighbors ~= -1); 
end

function [x, y] = id_to_coordinates(core_id, n_rows, n_cols)
    % this function returns x, y coordinates on a grid given core id
    % x is the horizontal coordinate
    % y is the vertical coordinate
    % numeration of slots (core_id) has the lower number at the beginning of
    % each row on the left and increases from top to bottom

    x = mod(core_id - 1, n_cols) + 1;
    y = floor((core_id - 1) / n_cols) + 1;
end

function core_id = coordinates_to_id(x_core, y_core, n_rows, n_cols)
    % this function returns core_id on a grid given its coordinates
    % x is the horizontal coordinate
    % y is the vertical coordinate

    if (x_core ~= -1) && (y_core~= -1)
        y_contrib = (y_core - 1) * n_cols;
        x_contrib = x_core;
    
        core_id = x_contrib + y_contrib;    
    else
        core_id = -1;
    end
end

function [x, y] = find_north_coordinate(x_core, y_core, ~, ~)
    if y_core > 1
        x = x_core;
        y = y_core - 1;
    else
        x = -1;
        y = -1;
    end
end

function [x, y] = find_south_coordinate(x_core, y_core, n_rows, ~)
    if y_core < n_rows
        x = x_core;
        y = y_core + 1;
    else
        x = -1;
        y = -1;
    end
end

function [x, y] = find_east_coordinate(x_core, y_core, ~, n_cols)
    if x_core < n_cols
        x = x_core + 1;
        y = y_core;
    else
        x = -1;
        y = -1;
    end
end

function [x, y] = find_west_coordinate(x_core, y_core, ~, ~)
    if x_core > 1
        x = x_core - 1;
        y = y_core;
    else
        x = -1;
        y = -1;
    end
end


