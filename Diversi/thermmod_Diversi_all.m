num_cores = 8;
num_sock = 2;
num_cpufeat = 3;
num_corefeat = 5;

model_order = 2;
num_add_eq = 10;

lab_std = 80;
temp_std = 50;

num_win = 14;
num_smallwin = 20;  % per window


%% Dati Galileo di produzione

file_path = "../outputs/";
out_path = file_path + "KalmanIdentRes/";
listing = dir(char(file_path));
file_names = "";
for indf = 1:length(listing)
    file_names = [file_names, string(listing(indf).name)];
end
%files_data = file_path + file_names(contains(file_names,"PowPartTherm_node"))';
files_data = file_path + file_names(contains(file_names,"PowPartTherm_Antarex"))';
nodespkg = extractBetween(files_data,"_",".mat");
num_nodpkg = length(files_data);
%files_csv = file_path + file_names(contains(file_names,"PredLinRegThermMod_node"))';
files_csv = file_path + file_names(contains(file_names,"PredLinRegThermMod_Antarex"))';
opts = detectImportOptions(char(files_csv(1)),'NumHeaderLines',1,'FileType','text');
pow_pkg_all = cell(length(files_csv)*num_sock,1);
for indf = 1:length(files_csv)
    input_data = readtable(char(files_csv(indf)), opts);
    pow_curr = table2array(input_data(:,2+(1:num_sock))) * lab_std;
    for inds = 1:num_sock
        pow_pkg_all{(indf-1)*num_sock + inds} = pow_curr(:,inds);
    end
end


% Do the identification and tracking per each node and package
for ind_node = 1:num_nodpkg
    data = load(files_data(ind_node));
    disp(files_data(ind_node))

    pow_pkg = pow_pkg_all{ind_node};
    pow_part = [sum(data.pow_part(:,1:2),2), data.pow_part(:,3:end)];
    temperatures = data.temperatures - data.temp_amb';
    pow_part = pow_part(1:end-1,:) * lab_std;
    temperatures = temperatures(2:end,:) * temp_std;
    win_length = floor(size(temperatures,1) / num_win);
    win_length_small = floor(win_length / num_smallwin);

    % Normalisations
    pow_mean = mean(pow_part,1);
    temp_mean = mean(temperatures, 1);
    % pow_mean = [0.16, 0.06 * ones(1,num_cores)];
    % temp_mean = [0.55, 0.36 * ones(1,num_cores)];
    pow_part = pow_part - pow_mean;
    temperatures = temperatures - temp_mean;


    % Calculate parameters and noise matrices in all windows
    num_autocorr = 100;
    pow_touse = pow_part;
    theta_init = []; 
    theta = NaN(10*model_order, num_cores, num_win);
    sigv = NaN(num_cores, num_win);
    sigw = NaN(num_cores, num_win);
    eigenvalues = NaN(num_cores, model_order, num_win);
    resid_autocorr = NaN(num_autocorr, num_cores, num_win);
    cond_ident = NaN(num_cores, num_win);
    svd_ident = NaN(20, num_cores, num_win);
    for core = 1:num_cores
        for ind_win = 1:num_win
            ind_curr = (ind_win-1)*win_length + (1:win_length);
            [theta(:,core,ind_win), sigv(core,ind_win), sigw(core,ind_win), res, ...
                cond_ident(core,ind_win), svd_ident(:,core,ind_win)] = ...
                misoarxbls(pow_touse(ind_curr,:), temperatures(ind_curr, 1+core), model_order, theta_init);
            %resid_autocorr(:,core,ind_win) = autocorr(res, num_autocorr-1);
            for indc = 1:num_autocorr
                resid_autocorr(indc,core,ind_win) = corr(res(1:end-indc+1), res(indc:end));
            end

            A = [-theta(1:model_order,core,ind_win), [eye(model_order-1); zeros(1,model_order-1)]];
            eigenvalues(core,:,ind_win) = eig(A);
        end
    end



    % Use the Kalman filter (evaluate the filter with coefficients derived from the chosen time window)
    Ts = 2;
    T_pred = NaN(size(temperatures,1),num_cores,num_win);
    T_filt = T_pred;
    residuals = T_pred;
    for core = 1:num_cores
        for ind_win = 1:num_win
            A = [-theta(1:model_order,core,ind_win), [eye(model_order-1); zeros(1,model_order-1)]];
            B = reshape(theta(model_order+1:end,core,ind_win), model_order, num_cores+1);
            C = [1, zeros(1,model_order-1)];
            Q = sigw(core,ind_win) * (C'*C);
            R = sigv(core,ind_win);
            Xi = zeros(model_order,1);
            Pi = eye(model_order); % * var(temperatures(:, 1+core));

        %     sys = ss(A, B, C, 0, Ts);
        % %     [kalmf,L,P] = kalman(sys, sigw(core), sigv(core), 0);
        %     y = lsim(sys, pow_part_re, 0:Ts:Ts*(size(temperatures,1)-1), Xi);

            for ind_win_Kalman = 1:num_win
                ind_curr = (ind_win_Kalman-1)*win_length + (1:win_length);  %1:size(temperatures,1);  %
                ind_time = ind_curr(1);
                [Xi, Pi, zf, yk] = KalmanFilter(A, eye(model_order), C, -A*Pi*A.' + Pi, R, ...
                    Xi, Pi, temperatures(ind_curr(1:model_order),1+core), temperatures(ind_time,1+core));
                T_pred(ind_time,core,ind_win) = yk(1) + temperatures(ind_time,1+core);
                T_filt(ind_time,core,ind_win) = zf;
                residuals(ind_time,core) = yk(1);
                for ind_time = ind_curr(2):ind_curr(end)
                    [Xi, Pi, zf, yk] = KalmanFilter(A, B, C, Q, R, Xi, Pi, pow_touse(ind_time-1,:)', temperatures(ind_time,1+core));
                    T_pred(ind_time,core,ind_win) = yk(1) + temperatures(ind_time,1+core);
                    T_filt(ind_time,core,ind_win) = zf;
                    residuals(ind_time,core) = yk(1);
                end
            end
        end
    end
    temperatures = (temperatures + temp_mean) + data.temp_amb(2:end)' * temp_std;
    temperatures_all = temperatures;
    T_pred = (T_pred + temp_mean(2:end)) + data.temp_amb(2:end)' * temp_std;
    T_filt = (T_filt + temp_mean(2:end)) + data.temp_amb(2:end)' * temp_std;
    err_pred = temperatures(:,2:9) - T_pred;
    err_filt = temperatures(:,2:9) - T_filt;
    
    % Save results
    save(out_path + "KalmanIdentAll_WinSep_" + nodespkg(ind_node) + ".mat", ...
        'T_pred', 'eigenvalues', 'pow_pkg', 'nodespkg', 'temperatures_all', 'resid_autocorr', ...
        'T_filt', 'sigv', 'sigw', 'theta', 'cond_ident', 'svd_ident') %, '-v7.3')
    
end

