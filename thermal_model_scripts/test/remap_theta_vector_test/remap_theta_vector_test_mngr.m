function remap_theta_vector_test_mngr

% Remapping works as described below:
% the first 2 elements are arx and are kept at the beginning of the array
% for each group of elements (2 elements for order 2 of the model) after the arx,
% they are moved the the indexes of the first available active core


input_test_cases = {
    % no errors
    [1 2 3 4 5 6 7 8], ... % coefficients
    6, ... % n cores
    [2 4 5]; ... % active cores indexes
    
    % no errors
    [1 2 3 4 5 6 7 8], ... % coefficients
    8, ... % n cores
    [2 4 5]; ... % active cores indexes

    % no errors
    [1 2 3 4 5 6 8 7], ... % coefficients
    8, ... % n cores
    [2 4 5]; ... % active cores indexes

    % no errors
    [1 2 3 4 5 6 8 7], ... % coefficients
    8, ... % n cores
    [2 4 8]; ... % active cores indexes

     % error, one active core has 0 as index
    [1 2 3 4 5 6 8 7], ... % coefficients
    8, ... % n cores
    [0 4 8]; ... % active cores indexes

    
    % active cores are more than the number of cores that can be covered
    % with the coefficients
    [1 2 3 4 ], ... % coefficients
    56, ... % n cores
    [2 7 8 9 10]; ...  % active cores indexes

    % there is more coefficients than needed
    [1 2 3 4 5 6 7 8], ... % coefficients
    56, ... % n cores
    [2 7]; ...  % active cores indexes

    % the highest core index is higher than the maximum allowd by number of
    [1 2 3 4 5 6 7 8], ... % coefficients
    6, ... % n cores
    [2 4 800]; ... % active cores indexes
    
    };

    
expected_outputs = {
    [1 2 0 0 3 4 0 0 5 6 7 8 0 0];
    [1 2 0 0 3 4 0 0 5 6 7 8 0 0 0 0 0 0];
    [1 2 0 0 3 4 0 0 5 6 8 7 0 0 0 0 0 0];
    [1 2 0 0 3 4 0 0 5 6 0 0 0 0 0 0 8 7];
    [1 2 0 0 3 4 0 0 5 6 0 0 0 0 0 0 8 7];
    [1 2 0 0 3 4 0 0 5 6 0 0 0 0 0 0 8 7];
    [1 2 0 0 3 4 0 0 5 6 0 0 0 0 0 0 8 7];
    [1 2 0 0 3 4 0 0 5 6 0 0 0 0 0 0 8 7];
};
    expected_errors = {
    '',
    '',
    '',    
    '',    
    'at least one core index is lower than 1',    
    'not enough exo coefficients to cover all active cores';
    'too many exo coefficients to cover all active cores';
    'at least one core index is over the maximum allowed from num_cores',
    };


num_cases = size(input_test_cases, 1);

for i = 1:num_cases
    theta_vector = input_test_cases{i, 1};
    num_cores = input_test_cases{i, 2};
    active_cores_indexes = input_test_cases{i, 3};
    expected_output = expected_outputs{i};
    expected_error = expected_errors{i};
    try 
        theta_vector_remapped = remap_theta_vector(theta_vector, num_cores, active_cores_indexes);
        if theta_vector_remapped == expected_output 
            fprintf('Test #%d PASSED\n', i);
        else
            error("theta_vector_remapped has incorrect content")
        end
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

end % end remap_theta_vector_test_mngr

