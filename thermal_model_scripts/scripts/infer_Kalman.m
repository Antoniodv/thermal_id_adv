function [INN, T_predicted, T_filtered] = infer_Kalman(THETA_TH_vector, T_measured, Power, sigv, sigw, model_order)
    [A,B,C] = getSSmat(squeeze(THETA_TH_vector), model_order);
    G = C';
    Q = sigw^2; 
    R = sigv^2;

    [INN,K,Pp,Pf,T_predicted,T_filtered,Xp,Xf] = Kalm(A,B,C,G,Power,T_measured,Q,R);
end