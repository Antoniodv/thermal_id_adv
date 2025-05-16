function [ar, exo] = separate_coefficients_avg(input_coefficients)
    %ar
    ar = input_coefficients(1:2);
    
    % exogen 
    exog = input_coefficients(3:end);
    % exo avg by couple of coefficients
    b1 = exog(1:2:end);  % odd
    b2 = exog(2:2:end);  % even
    b_arr = (b1 + b2) / 2;
    
    exo = b_arr;
end