function [ theta_th_array, a, b, c, d, e ] = misoarxbls_mocked(U, T_single_core, thermal_model_order, num_add_eq)
    % Mock implementation for testing
    % Return test values
    [n_rows, n_cols] = size(U);
    theta_th_array = ones(2+2*n_cols, 1);
    a = 1;
    b = 1;
    c = 1;
    d = 1;
    e = 1;
end