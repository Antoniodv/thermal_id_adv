function [ar_norm, exo_norm] = separate_coefficients_norm(input_coefficients)
    
    %ar
    ar = input_coefficients(1:2);
    ar1 = ar(1);
    ar2 = ar(2);
    ar1_sqr = ar1^2;
    ar2_sqr = ar2^2;
    ar_norm = sqrt(ar1_sqr + ar2_sqr);

    % exogen 
    exog = input_coefficients(3:end);

    % exo avg by couple of coefficients
    b1 = exog(1:2:end);  % odd
    b2 = exog(2:2:end);  % even

    b1_sqr = b1.^2;
    b2_sqr = b2.^2;

    b_arr = b1_sqr + b2_sqr;
    b_arr = sqrt(b_arr);
    
    exo_norm = b_arr;
end