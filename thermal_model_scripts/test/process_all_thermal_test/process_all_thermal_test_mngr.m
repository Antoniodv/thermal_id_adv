function process_all_thermal_test_mngr

ROOT_PATH = "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification\thermal_model_scripts";

path_to_original_script = fullfile(ROOT_PATH, "\scripts\process_all_thermal.m");

% Open original script
fid = fopen(path_to_original_script, "r");
original_script_code = fileread(path_to_original_script);
fclose(fid);

% Sobstitute real function with mocked function in test script
mocked_script_code = strrep(original_script_code, ...
    "thermal_model_est",...
    "thermal_model_est_mocked");

% Save new script 
mocked_script_name = fullfile(ROOT_PATH, "\test\process_all_thermal_test\process_all_thermal_mocked.m");
fid = fopen(mocked_script_name, "w");
fwrite(fid, mocked_script_code);
fclose(fid);

% Create mocked function
mocked_function_code = [...
    'function [coeff_matrix, active_cores_indexes] = thermal_model_est_mocked(pkl_file, path_to_active_cores, thermal_model_path, mode)', newline, ...
    '    % Mock implementation for testing', newline, ...
    '    disp([''Mock thermal_model_est called with: '' pkl_file]);', newline, ...
    '    ', newline, ...
    '    % Return test values', newline, ...
    '    coeff_matrix = ones(56, 56);', newline, ... % mock output
    '    active_cores_indexes = ones(14);', newline, ... % mock output
    'end' ...
];

mocked_function_name = fullfile(ROOT_PATH, "\test\process_all_thermal_test\thermal_model_est_mocked.m");
fid = fopen(mocked_function_name, "w");
fwrite(fid, mocked_function_code);
fclose(fid);

% Prepare different input sets
pkl_file_inputs = {
    fullfile(ROOT_PATH, 'test', 'mock_root_good/round0/power_model/gb_core_uncore_tot_temp.pkl'),
    fullfile(ROOT_PATH, 'test', 'mock_root_good/round1/power_model/gb_core_uncore_tot_temp.pkl'),
    fullfile(ROOT_PATH, 'test', 'mock_root_good/round2/power_model/gb_core_uncore_tot_temp.pkl'),
};


% Create pkl files trees
for i = 1:length(pkl_file_inputs)
    create_dummy_file(pkl_file_inputs{i});
end

% Create test cases
%  all correct, mode all
%  all correct, mode patch
%  all correct, mode not recog
%  all correct, mode not recog
%  all correct, mode not recog
%  all correct but root path does not exist

input_test_cases = {
    'mock_root_good','active_cores_orig','all';
    'mock_root_good','active_cores_orig','patch';
    'mock_root_good','active_cores_orig','nope';
    'mock_root_good','active_cores_orig','patchall';
    'mock_root_good','active_cores_orig','allpatch';
    '_','active_cores_orig','all';
};

expected_errors = {
    '',    
    '',     
    'SELECTED MODE NOT ALLOWED',  
    'SELECTED MODE NOT ALLOWED',   
    'SELECTED MODE NOT ALLOWED',   
    'The specified root path does not exist.'   
};


% Testing
num_cases = size(input_test_cases, 1);

for i = 1:num_cases
    root_path = input_test_cases{i,1};
    path_to_active_cores = input_test_cases{i,2};
    mode = input_test_cases{i,3};
    
    try
        [thermal_coeff_mat, active_cores_mat] = process_all_thermal_mocked(root_path, path_to_active_cores, mode);

        fprintf('\nTest #%d function executed without errors\n', i);
        disp(['In Mode: ', mode]);
        disp(['Out Global matrix size: ', mat2str(size(thermal_coeff_mat))]);
        disp(['Out Active cores shape: ', mat2str(size(active_cores_mat))]);

        % expected vs computed
        if isequal(size(thermal_coeff_mat), [56, 56, 3]) && all(thermal_coeff_mat(:) == 1)
            disp('Global matrix has correct size (56x56x3) and all ones.');
        else
            error('Global matrix size or values are incorrect.');
        end

        if isequal(size(active_cores_mat), [3, 14]) && all(active_cores_mat(:) == 1)
            disp('Active cores matrix has correct size (3x14) and all ones.');
        else
            error('Active cores matrix size or values are incorrect.');
        end

        % Verifica se il test è stato superato
        disp("PASSED");
        
    catch ME

        fprintf('\nTest #%d Executed with errors\n', i);
        disp(['Error: ', ME.message]);

  
        expected_error = expected_errors{i};
        
        % Se non c'è un errore atteso (stringa vuota), il test è fallito in modo imprevisto
        if isempty(expected_error)
            fprintf('Test #%d failed unexpectedly.\n', i);
            disp(['Expected no error, but got: ', ME.message]);
        else
            % Altrimenti, verifica che l'errore generato sia uguale a quello atteso
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



end %end process_all_thermal_test_mngr



