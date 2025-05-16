close all
clear all
clc
% path_to_pm_file =  "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification\MAR17\prbs_random\112\round0\power_model\gb_core_uncore_tot_temp.pkl";
% path_to_active_cores = "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification\MAR17\active_cores.txt";
path_to_pm_file = "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification\APR5\prbs_random\112\round0\power_model\gb_core_uncore_tot_temp.pkl";
path_to_active_cores = "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification\active_cores_APR5.txt";
mode = "patch";

abs_or_signed = "abs";

disp("    ")
disp("start processing")
disp(path_to_pm_file)
disp("path_to_active_cores")
disp(path_to_active_cores)
disp("                  ")

py.importlib.import_module('pandas');
py.importlib.import_module('numpy'); 

% active cores from file
round_str = regexp(path_to_pm_file, 'round(\d+)', 'tokens');
if isempty(round_str)
    error('No "round" word in this path');
end

% Estrae il numero come stringa e poi lo converte in numero
round_index = (str2double(round_str{1}{1}));    
round_index_matlab = round_index + 1;
active_cores_indexes = sort(get_row_from_file(round_index_matlab, path_to_active_cores));

% from python/C-like indexes to MATLAB indexes
active_cores_indexes = active_cores_indexes+1;

N_CORE = 56;
thrm = THID(2,N_CORE);
thrm.Nmc = 4;
thrm.Nmp = 2;
thrm.test_sampling_time = 0.5;

TEST_TYPE = 1; % Gathering the PRBS data (only one available)
dfmat_idx = 1; %pkl
dtype_idx = 2; %list of dictionaries 1, dict 2

thrm.model_order = 2;
thrm.num_add_eq = 10;

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

% normalize 
Y = normalize(Y, 'range');
U = normalize(U, 'range');

% reduce power matrix to the matrix of active cores
U_active_cores = U(:, active_cores_indexes);
% noise_levels = [0.01, 0.1, 1];
noise_levels = [0.001, 0.01, 0.1, 1];
n_added_dummy_inputs = 42;



% uniform distr
U_rnd = 2 * rand(size(U,1), n_added_dummy_inputs) - 1;


mse_mat = zeros(length(noise_levels), n_added_dummy_inputs+1);
condn_mat =  zeros(length(noise_levels), n_added_dummy_inputs+1);

for core = active_cores_indexes
    disp(core)
    % Testing adding power traces of noise
    for i=1:length(noise_levels)
        
        U_rnd_scaled = noise_levels(i)*U_rnd;
        mse = zeros(1, n_added_dummy_inputs + 1);

        for j=1:n_added_dummy_inputs+1
            if j == 1
                U_dummy = U_active_cores;
            else
                U_dummy = [U_dummy, U_rnd_scaled(:,j-1)];
            end
        
            [THETA_TH_vector, sigv,sigw, res, condn, svd_ident] = ...
            misoarxbls(U_dummy,Y(:,core),thrm.model_order,thrm.num_add_eq);

            
            % inference 
            [INN, T_predicted, T_filtered] = infer_Kalman(THETA_TH_vector, Y(:,core), U_dummy, sigv, sigw, thrm.model_order);
            
            % align and cut initial overshoot
            % T_predicted = T_predicted();
            % T_filtered  = T_filtered();

            % compute mse
            mse(j) = var(INN);

            condn_mat(j) = condn;

        end
        mse_mat(i, :) = mse;
    end

    
    figure('visible', 'off', 'Position', [100, 100, 2000, 1400]);
    for i = 1:size(mse_mat, 1)
        subplot(ceil(length(noise_levels)/2), 2, i);
        
        % Create the left y-axis for MSE bars
        yyaxis left;
        bar(mse_mat(i, :), 'FaceColor', [0.3010, 0.7450, 0.9330]);
        ylabel('MSE');
        % Set appropriate y-limits based on your data
        % ylim([0, max(mse_mat(i, :))*1.1]);
        
        % Create the right y-axis for condition number line
        yyaxis right;
        scaled_condn = condn_mat(i, :);
        plot(1:length(scaled_condn), scaled_condn, '-o', 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980], 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
        ylabel('Condition Number');
        
        % Set title and other properties
        title(sprintf('Noise = %.3f', noise_levels(i)));
        xlabel('Number of added inputs');
        xticks(1:(n_added_dummy_inputs + 1));
        grid on;
        
        % Add legend
        legend('MSE', 'Condition Number');
    end
    
    sgtitle(sprintf('Core %d: MSE and Condition Number vs Added Inputs', core-1)); 

    output_dir = 'stress_misoarxbls_plot';
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
        disp(['Created output directory: ' output_dir]);
    end

    filename = sprintf('%s/core%d_plot.png', output_dir, core-1);
    saveas(gcf, filename);
    disp(['Saved combined plot: ' filename]);

end
