function thermal_model_est_test_mngr

ROOT_PATH = "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification\thermal_model_scripts";

path_to_original_script = fullfile(ROOT_PATH, "\scripts\thermal_model_est.m");

% Open original script
fid = fopen(path_to_original_script, "r");
original_script_code = fileread(path_to_original_script);
fclose(fid);

% Sobstitute real function with mocked function in test script
mocked_script_code = strrep(original_script_code, ...
    "misoarxbls",...
    "misoarxbls_mocked");

mocked_script_code = strrep(mocked_script_code, ...
    "function [theta_matrix, active_cores_indexes] = thermal_model_est(path_to_metrics, path_to_active_cores, output_folder, mode)",...
    "function [theta_matrix, active_cores_indexes , U_active_cores, round_index_matlab] = thermal_model_est_mocked(path_to_metrics, path_to_active_cores, output_folder, mode)");


% Save new script 
mocked_script_name = fullfile(ROOT_PATH, "\test\thermal_model_est_test\thermal_model_est_mocked.m");
fid = fopen(mocked_script_name, "w");
fwrite(fid, mocked_script_code);
fclose(fid);

% Create mocked function
mocked_function_code = [...
    'function [ theta_th_array, a, b, c, d, e ] = misoarxbls_mocked(U, T_single_core, thermal_model_order, num_add_eq)', newline, ...
    '    % Mock implementation for testing', newline, ...
    '    % Return test values', newline, ...
    '    [n_rows, n_cols] = size(U);', newline, ...
    '    theta_th_array = ones(2+2*n_cols, 1);', newline, ... % mock output
    '    a = 1;', newline, ... % mock output
    '    b = 1;', newline, ... % mock output
    '    c = 1;', newline, ... % mock output
    '    d = 1;', newline, ... % mock output
    '    e = 1;', newline, ... % mock output
    'end' ...
];

mocked_function_name =fullfile(ROOT_PATH, "\test\thermal_model_est_test\misoarxbls_mocked.m");
fid = fopen(mocked_function_name, "w");
fwrite(fid, mocked_function_code);
fclose(fid);

% Prepare different input sets
pkl_file_inputs = {
    fullfile(ROOT_PATH,'test\mock_root_good/round0/power_model/gb_core_uncore_tot_temp.pkl'),
    fullfile(ROOT_PATH,'test\mock_root_good/round1/power_model/gb_core_uncore_tot_temp.pkl'),
    fullfile(ROOT_PATH,'test\mock_root_good/round2/power_model/gb_core_uncore_tot_temp.pkl'),
};

ref_pickle_file = fullfile(ROOT_PATH, "test\thermal_model_est_test\gb_core_uncore_tot_temp.pkl");


% Create pkl files trees and copy reference test pickle
for i = 1:length(pkl_file_inputs)
    target_path = pkl_file_inputs{i};
    
    % Ensure directory exists
    target_dir = fileparts(target_path);
    if ~isfolder(target_dir)
        mkdir(target_dir);
    end

    % Copy reference pickle to target path
    copyfile(ref_pickle_file, target_path);
end


% Create test cases
input_test_cases = {
        % Recognizing right active cores, no errors
    	pkl_file_inputs{1}, ... %round0
        fullfile(ROOT_PATH, "test/thermal_model_est_test/active_cores_orig.txt"), ...
        fullfile(ROOT_PATH, "test/thermal_model_est_test/output_folder"), ...
        "patch";

        % % Recognizing right active cores, first numeric line in the file has both numbers and letters 
        pkl_file_inputs{1}, ... %round0
        fullfile(ROOT_PATH, "test/thermal_model_est_test/active_cores_mix.txt"), ...
        fullfile(ROOT_PATH, "test/thermal_model_est_test/output_folder"), ...
        "patch";

        % Recognizing right active cores, active cores file is empty
        pkl_file_inputs{1}, ... %round0
        fullfile(ROOT_PATH, "test/thermal_model_est_test/active_cores_empty.txt"), ...
        fullfile(ROOT_PATH, "test/thermal_model_est_test/output_folder"), ...
        "patch";

        % Recognizing right active cores, active cores file does not exist
        pkl_file_inputs{1}, ... %round0
        fullfile(ROOT_PATH, "test/thermal_model_est_test/active_cores_no_exist.txt"), ...
        fullfile(ROOT_PATH, "test/thermal_model_est_test/output_folder"), ...
        "patch";

        % Round recognition, correct
        pkl_file_inputs{2}, ... %round1
        fullfile(ROOT_PATH, "test/thermal_model_est_test/active_cores_orig.txt"), ...
        fullfile(ROOT_PATH, "test/thermal_model_est_test/output_folder"), ...
        "patch"; 

        % Round recognition, correct
        pkl_file_inputs{3}, ... %round2
        fullfile(ROOT_PATH, "test/thermal_model_est_test/active_cores_orig.txt"), ...
        fullfile(ROOT_PATH, "test/thermal_model_est_test/output_folder"), ...
        "patch";

        % Round recognition, no round word to parse in path
        strrep(pkl_file_inputs{1}, "round", "") ... %round0 but removing "round" word
        fullfile(ROOT_PATH, "test/thermal_model_est_test/active_cores_orig.txt"), ...
        fullfile(ROOT_PATH, "test/thermal_model_est_test/output_folder"), ...
        "patch";

        % % Wrong path to pkl even if process_all_thermal is doing this
        % check before calling thermal_model_est
        strrep(pkl_file_inputs{1}, "mock_", ""), ... %round0
        fullfile(ROOT_PATH, "test/thermal_model_est_test/active_cores_orig.txt"), ...
        fullfile(ROOT_PATH, "test/thermal_model_est_test/output_folder"), ...
        "patch";
        
        % % Test U active cores column selection
        % pkl_file_inputs{1}, ... %round0
        % fullfile(ROOT_PATH, "test/thermal_model_est_test/active_cores_orig.txt"), ...
        % fullfile(ROOT_PATH, "test/thermal_model_est_test/output_folder"), ...
        % "patch";
        % 
        % % Test U active cores column selection
        % pkl_file_inputs{2}, ... %round1
        % fullfile(ROOT_PATH, "test/thermal_model_est_test/active_cores_orig.txt"), ...
        % fullfile(ROOT_PATH, "test/thermal_model_est_test/output_folder"), ...
        % "patch";

% Check on Power columns selection in patch mode
% Check on Power columns selection in all mode
};

expected_errors = {
    '',    
    'active_cores_indexes has incorrect content',
    'Index out of range. The file contains 0 numeric rows',
    'File opening error',
    '',    
    '', 
    'No "round" word in this path'
    'The path does not bring to the right pkl file'
};

expected_outputs = {
    % round 0
    [5 8 12 15 17 18 27 31 34 39 43 49 52 54], 1;
    [5 8 12 15 17 18 27 31 34 39 43 49 52 54], 1;
    [5 8 12 15 17 18 27 31 34 39 43 49 52 54], 1;
    [5 8 12 15 17 18 27 31 34 39 43 49 52 54], 1;
    % round 1
    [6 7 9 10 16 22 23 28 33 37 40 45 53 56], 2;
    % round 2
    [2 3 14 19 20 26 29 32 41 44 46 47 51 55], 3;
    % round wrong
    [], 1;
    % pkl path wrong
    [], 1;
    
};

% Testing
num_cases = size(input_test_cases, 1);

for i = 1:num_cases
    path_to_metrics = input_test_cases{i,1};
    path_to_active_cores = input_test_cases{i,2};
    output_folder = input_test_cases{i,3};
    mode = input_test_cases{i,4};

    try
        [theta_matrix, active_cores_indexes , U_active_cores, round_index] = thermal_model_est_mocked(path_to_metrics, path_to_active_cores, output_folder, mode);
 
        fprintf('\nTest #%d function executed without errors\n', i);
        disp(['In Mode: ', mode]);
        disp(['Out theta_matrix size: ', mat2str(size(theta_matrix))]);
        disp(['Out Active cores shape: ', mat2str(size(active_cores_indexes))]);
        
        % expected vs computed
        if isequal(size(theta_matrix), [56, 114]) 
            disp('theta_matrix has correct size [56, 114]');
        else
            error('theta_matrix size or values are incorrect.');
        end

        if isequal(size(active_cores_indexes), [1, 14]) 
            disp('active_cores_indexes has correct size (1x14)');
        else
            error('active_cores_indexes size or values are incorrect.');
        end
        
        active_cores_indexes_expected = expected_outputs{i, 1};
        round_index_expected = expected_outputs{i, 2};
        
        if isequal(active_cores_indexes, active_cores_indexes_expected)
            disp('active_cores_indexes has correct content');
        else
            error('active_cores_indexes has incorrect content');
        end
        disp("PASSED");
        
        if isequal(round_index, round_index_expected)
            disp('round_index has correct content');
        else
            error('round_index has incorrect content');
        end

        disp("PASSED");
 
    catch ME

        fprintf('\nTest #%d Executed with errors\n', i);
        disp(['Error: ', ME.message]);
        disp('Error Stack:');
        for stackIdx = 1:length(ME.stack)
            fprintf('   File: %s > Function: %s (Line %d)\n', ...
                ME.stack(stackIdx).file, ...
                ME.stack(stackIdx).name, ...
                ME.stack(stackIdx).line);
        end


        expected_error = expected_errors{i};

        % if expected error is empy (but we are in catch case so there is
        % an error)
        if isempty(expected_error)
            fprintf('Test #%d failed unexpectedly.\n', i);
            disp(['Expected no error, but got: ', ME.message]);
        else
            % generated error vs expected
            if contains(ME.message, expected_error)
                fprintf('Test #%d correctly failed with expected error: %s\n', i, expected_error);
                fprintf('PASSED\n');
            else
                fprintf('Test #%d failed with an unexpected error.\n', i);
                disp(['Expected: ', expected_error]);
                disp(['Got: ', ME.message]);
                fprintf('FAILED\n');
            end
        end
    end
end



end %end thermal_model_est_test_mngr



