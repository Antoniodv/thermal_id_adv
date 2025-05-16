clear;

% APR5
path_to_pw_temp = "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification\APR5\prbs_random\112\";
path_to_active_cores = "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification\active_cores_APR5.txt";

% % MAR17
% path_to_pw_temp = "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification\MAR17\prbs_random\112\";
% path_to_active_cores = "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification\MAR17\active_cores.txt";


mode = "patch";
N_CORE = 56;


% EXTRACT GLOBAL MATRIX
[global_matrix, active_cores_indexes_matrix] = process_all_thermal(path_to_pw_temp, path_to_active_cores, mode);

%%
% SUM GLOBAL MATRIX
global_matrix_sum = sum(global_matrix, 3);

% create avg matrix of exo coefficients
b_mat = zeros(N_CORE, N_CORE);
for i = 1:size(global_matrix_sum, 1)
    [~, b_mat_arr] = separate_coefficients_norm(global_matrix_sum(i,:));
    b_mat(i,:) = b_mat_arr;
end


% AVERAGE 
occurence_matrix = zeros(N_CORE, N_CORE); 
for i = 1:size(active_cores_indexes_matrix, 1)
    current_round_cores = active_cores_indexes_matrix(i, :);
    for j=1:length(current_round_cores)
        for k=1:length(current_round_cores)
            if j ~= k
                row_core_id = current_round_cores(j);
                col_core_id = current_round_cores(k);
                occurence_matrix(row_core_id, col_core_id) = occurence_matrix(row_core_id, col_core_id) + 1;
            end
        end
    end
end
occurence_matrix(occurence_matrix == 0) = 1; 
b_mat = b_mat ./ occurence_matrix;


% % NORMALIZATION
% b_mat_norm = b_mat; 
% for i = 1:size(b_mat, 1)
%     curr_line = b_mat(i, :);
%     max_of_curr_line = max(curr_line);
%     if max_of_curr_line > 0
%         b_mat_norm(i, :) = curr_line ./ max_of_curr_line;
%     end
% end



% Save computed data
save('global_data.mat', 'b_mat', 'global_matrix', 'active_cores_indexes_matrix');

% Plot theta matrix
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

% Imposta la scala dei colori su logaritmica
set(gca, 'ColorScale', 'log');

% Find max elem for each row
[~, max_indices] = max(b_mat, [], 2);

% Add X for each maximum element in row
for i = 1:m
    text(max_indices(i), i, 'x', 'Color', 'w', 'FontSize', 8, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end


plot_dir_path = './score_vs_mse_plot'; 
if ~exist(plot_dir_path, 'dir')
    mkdir(plot_dir_path);
end


filename = 'global_coefficients_matrix_log.png';
fullpath = fullfile(plot_dir_path, filename);
disp("Saving figure to: " + fullpath);
saveas(fig, fullpath);
% close(fig);
