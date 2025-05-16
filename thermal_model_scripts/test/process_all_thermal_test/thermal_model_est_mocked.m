function [coeff_matrix, active_cores_indexes] = thermal_model_est_mocked(pkl_file, path_to_active_cores, thermal_model_path, mode)
    % Mock implementation for testing
    disp(['Mock thermal_model_est called with: ' pkl_file]);
    
    % Return test values
    coeff_matrix = ones(56, 56);
    active_cores_indexes = ones(14);
end