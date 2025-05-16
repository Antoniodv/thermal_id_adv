%% This script builds a Diversi's model combining data from all windows in the training set and evaluates it on the test set
% The input to the algorithm is the same of the LSTM: core estimated power,
% measured package power and core previous temperature

num_cores = 8;
num_sock = 2;
num_cpufeat = 3;
num_corefeat = 5;

model_order = 2;
num_add_eq = 10;

lab_std = 80;
temp_std = 50;

num_win = 25;
num_smallwin = 20;  % per window


%% Dati Galileo di produzione

file_path = "../outputs/";
out_path = file_path + "KalmanIdentRes/";
data = importdata(out_path + 'dataset_core_summary.mat');
nodwin_train = data.nodwin_train + 1;
nodwin_test = data.nodwin_test + 1;
pow_pkg = data.pow_pkg;
win_length_min = Inf;
for ind_node = 1:length(pow_pkg)
    win_length_min = min(win_length_min, floor(max(size(pow_pkg{ind_node})) / num_win));
end
clear pow_pkg
% data = importdata(out_path + 'dataset_core_summary_powT.mat');
% pow_part = data.pow_part;
% temperatures = data.temperatures_all;

listing = dir(char(file_path));
file_names = "";
for indf = 1:length(listing)
    file_names = [file_names, string(listing(indf).name)];
end
files_data = file_path + file_names(contains(file_names,"PowPartTherm"))';
nodespkg = extractBetween(files_data,"_",".mat");
num_nodpkg = length(files_data);
files_csv = file_path + file_names(contains(file_names,"PredLinRegThermMod"))';
opts = detectImportOptions(char(files_csv(1)),'NumHeaderLines',1,'FileType','text');
pow_pkg_all = cell(length(files_csv)*num_sock,1);
for indf = 1:length(files_csv)
    input_data = readtable(char(files_csv(indf)), opts);
    pow_curr = table2array(input_data(:,2+(1:num_sock))) * lab_std;
    for inds = 1:num_sock
        pow_pkg_all{(indf-1)*num_sock + inds} = pow_curr(:,inds);
    end
end


% Perform the identification with all data in the training set (using subsets)

num_nodwin = num_nodpkg * num_cores * num_win;
win_length = win_length_min;

num_sub = 10;
theta_all = NaN(2*model_order, num_sub);
sigv_all = NaN(num_sub, 1);
sigw_all = NaN(num_sub, 1);
for ind_sub = 1:num_sub
    nodwin_train_sub = downsample(nodwin_train, num_sub, ind_sub-1);
    temperatures_all = NaN(length(nodwin_train_sub) * win_length, 1);
    pow_touse = NaN(length(nodwin_train_sub) * win_length, 1);
    ind_wintrain = 0;
    for ind_node = 1:num_nodpkg
        data = load(files_data(ind_node));

        pow_pkg = pow_pkg_all{ind_node};
        pow_part = [sum(data.pow_part(:,1:2),2), data.pow_part(:,3:end)];
        temperatures = data.temperatures - data.temp_amb';
        pow_pkg = pow_pkg(1:end-1);
        pow_part = pow_part(1:end-1,:) * lab_std;
        temperatures = temperatures(2:end,:) * temp_std;
        %win_length = floor(size(temperatures,1) / num_win);
        win_length_small = floor(win_length / num_smallwin);

        % Normalisations
        pow_pkg_mean = mean(pow_pkg,1);
        pow_mean = mean(pow_part,1);
        temp_mean = mean(temperatures, 1);
        % pow_mean = [0.16, 0.06 * ones(1,num_cores)];
        % temp_mean = [0.55, 0.36 * ones(1,num_cores)];
        pow_pkg = pow_pkg - pow_pkg_mean;
        pow_part = pow_part - pow_mean;
        temperatures = temperatures - temp_mean;


        % Calculate parameters and noise matrices in all windows
        num_autocorr = 100;
        for core = 1:num_cores
            for ind_win = 1:num_win
                ind_curr = (ind_win-1)*win_length + (1:win_length);
                ind_nodwin = ((ind_node-1) * num_cores + core-1) * num_win + ind_win;

                if ismember(ind_nodwin, nodwin_train_sub)
                    temperatures_all(ind_wintrain + (1:win_length)) = temperatures(ind_curr, 1+core);
                    pow_touse(ind_wintrain + (1:win_length), :) = [pow_pkg(ind_curr)];%, pow_part(ind_curr, 1+core)];
                    ind_wintrain = ind_wintrain + win_length;
                end
            end
        end

    end

    % Identification
    theta_init = []; 
    [theta_all(:,ind_sub), sigv_all(ind_sub), sigw_all(ind_sub), res] = ...
        misoarxbls(pow_touse, temperatures_all, model_order, theta_init);
    clear pow_touse temperatures_all
    
end

% Define an average model
theta = mean(theta_all,2);
sigv = mean(sigv_all);
sigw = mean(sigw_all);


% Do the tracking on all nodes and packages
for ind_node = 1:num_nodpkg
    data = load(files_data(ind_node));
    disp(files_data(ind_node))

    pow_pkg = pow_pkg_all{ind_node};
    pow_part = [sum(data.pow_part(:,1:2),2), data.pow_part(:,3:end)];
    temperatures = data.temperatures - data.temp_amb';
    pow_pkg = pow_pkg(1:end-1);
    pow_part = pow_part(1:end-1,:) * lab_std;
    temperatures = temperatures(2:end,:) * temp_std;
    %win_length = floor(size(temperatures,1) / num_win);
    win_length_small = floor(win_length / num_smallwin);

    % Normalisations
    pow_pkg_mean = mean(pow_pkg,1);
    pow_mean = mean(pow_part,1);
    temp_mean = mean(temperatures, 1);
    % pow_mean = [0.16, 0.06 * ones(1,num_cores)];
    % temp_mean = [0.55, 0.36 * ones(1,num_cores)];
    pow_pkg = pow_pkg - pow_pkg_mean;
    pow_part = pow_part - pow_mean;
    temperatures = temperatures - temp_mean;


    % Calculate parameters and noise matrices in all windows
    num_autocorr = 100;



    % Use the Kalman filter (evaluate the filter with coefficients derived from the chosen time window)
    Ts = 2;
    T_pred = NaN(size(temperatures,1),num_cores);
    T_filt = T_pred;
    residuals = T_pred;
    for core = 1:num_cores
        A = [-theta(1:model_order), [eye(model_order-1); zeros(1,model_order-1)]];
        B = reshape(theta(model_order+1:end), model_order, 1);
        C = [1, zeros(1,model_order-1)];
        Q = sigw * (C'*C);
        R = sigv;
        Xi = zeros(model_order,1);
        Pi = eye(model_order); % * var(temperatures(:, 1+core));
        pow_touse = [pow_pkg];%, pow_part(:,1+core)];

    %     sys = ss(A, B, C, 0, Ts);
    % %     [kalmf,L,P] = kalman(sys, sigw(core), sigv(core), 0);
    %     y = lsim(sys, pow_part_re, 0:Ts:Ts*(size(temperatures,1)-1), Xi);

        %for ind_win_Kalman = 1:num_win
        ind_curr = 1:size(temperatures,1);  %(ind_win_Kalman-1)*win_length + (1:win_length);  %
        ind_time = ind_curr(1);
        [Xi, Pi, zf, yk] = KalmanFilter(A, eye(model_order), C, -A*Pi*A.' + Pi, R, ...
            Xi, Pi, temperatures(ind_curr(1:model_order),1+core), temperatures(ind_time,1+core));
        T_pred(ind_time,core) = yk(1) + temperatures(ind_time,1+core);
        T_filt(ind_time,core) = zf;
        residuals(ind_time,core) = yk(1);
        for ind_time = ind_curr(2):ind_curr(end)
            [Xi, Pi, zf, yk] = KalmanFilter(A, B, C, Q, R, Xi, Pi, pow_touse(ind_time-1,:)', temperatures(ind_time,1+core));
            T_pred(ind_time,core) = yk(1) + temperatures(ind_time,1+core);
            T_filt(ind_time,core) = zf;
            residuals(ind_time,core) = yk(1);
        %end
        end
    end
    temperatures = (temperatures + temp_mean) + data.temp_amb(2:end)' * temp_std;
    temperatures_all = temperatures;
    T_pred = (T_pred + temp_mean(2:end)) + data.temp_amb(2:end)' * temp_std;
    T_filt = (T_filt + temp_mean(2:end)) + data.temp_amb(2:end)' * temp_std;
    err_pred = temperatures(:,2:9) - T_pred;
    err_filt = temperatures(:,2:9) - T_filt;
    
    % Save results
    save(out_path + "KalmanIdentAllCombined_" + nodespkg(ind_node) + ".mat", ...
        'T_pred', 'pow_pkg', 'nodespkg', 'temperatures_all', ...
        'T_filt', 'sigv_all', 'sigw_all', 'theta_all') %, '-v7.3')
    
end

