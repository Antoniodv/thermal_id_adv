clear;

% 
% % % MAR17
% path_to_pw_temp = "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification\MAR17\prbs_random\112\";
% path_to_active_cores = "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification\MAR17\active_cores.txt";

% APR5
path_to_pw_temp = "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification\APR5\prbs_random\112\";
path_to_active_cores = "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification\active_cores_APR5.txt";


mode = "patch";

N_CORE = 56;
thrm = THID(2,N_CORE);
thrm.Nmc = 4;
thrm.Nmp = 2;
thrm.test_sampling_time = 0.5;

% EXTRACT GLOBAL MATRIX
[global_matrix, active_cores_indexes_matrix] = process_all_thermal(path_to_pw_temp, path_to_active_cores, mode);

%% 

% EXTRACT EXO
global_matrix_exo = global_matrix(:,3:end,:);

% NORMALIZE EXO MATRIX
global_mat_exo_normalized = normalize_3d_by_row_max(global_matrix_exo);

% SUM EXO MATRIX
global_mat_exo_sum = sum(global_mat_exo_normalized, 3); 

% AVERAGE 
occurence_matrix = zeros(N_CORE, N_CORE); 
for i = 1:size(active_cores_indexes_matrix, 1)
    current_round_cores = active_cores_indexes_matrix(i, :);
    for j=1:length(current_round_cores)
        for k=1:length(current_round_cores)
                row_core_id = current_round_cores(j);
                col_core_id = current_round_cores(k);
                occurence_matrix(row_core_id, col_core_id) = occurence_matrix(row_core_id, col_core_id) + 1;
        end
    end
end
occurence_matrix(occurence_matrix == 0) = 1; 
occurrence_double = zeros(size(occurence_matrix, 1), 2*size(occurence_matrix, 2));
occurrence_double(1:2:end-1) = occurence_matrix;
occurrence_double(2:2:end) = occurence_matrix;

global_matrix_avg = global_mat_exo_sum ./ occurrence_double;

% VECTOR NORM
b_mat = compute_exo_norms(global_matrix_avg);

% SAVE DATA
save('global_data.mat', 'b_mat', 'global_matrix', 'active_cores_indexes_matrix');

% PLOT HEATMAP
fig = figure('Position', [100, 100, 800, 600]);
imagesc(b_mat);
hold on;
title("Global coefficients matrix " + mode);
yticks(1:size(b_mat, 1));
[m, n] = size(b_mat);

% A diagonal line
plot(1:min(m, n), 1:min(m, n), 'w-', 'LineWidth', 1);
colorbar;
colormap(flipud(autumn));

% Log color scale
set(gca, 'ColorScale', 'log');

% Find max elem for each row
[~, max_indices] = max(b_mat, [], 2);

% Add X for each maximum element in row
for i = 1:m
    text(max_indices(i), i, 'x', 'Color', 'w', 'FontSize', 8, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

%%
% SAVE HEATMAP
plot_dir_path = './score_vs_mse_plot'; 
if ~exist(plot_dir_path, 'dir')
    mkdir(plot_dir_path);
end


% filename = 'global_coefficients_matrix_log.png';
% fullpath = fullfile(plot_dir_path, filename);
% disp("Saving figure to: " + fullpath);
% saveas(fig, fullpath);
% close(fig);


function normalized_mat = normalize_3d_by_row_max(input_mat)
    normalized_mat = zeros(size(input_mat));  
    for round_idx = 1:size(input_mat, 3)
        for row_idx = 1:size(input_mat, 1)
            curr_line = input_mat(row_idx, :, round_idx);
            curr_line_max_value = max(curr_line);
            if curr_line_max_value > 0
                normalized_mat(row_idx, :, round_idx) = curr_line ./ curr_line_max_value;
            end
        end
    end
end
