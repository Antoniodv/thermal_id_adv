num_cores = 8;
num_sock = 2;
num_cpufeat = 3;
num_corefeat = 5;

model_order = 2;
num_add_eq = 10;

lab_std = 80;
temp_std = 50;


%% Dati Galileo di Diversi

% load('/Users/FedeAdm/Documents/Postdoc Bologna/HPC Applications/power-models/outputs/datiGalileoDiversi.mat')
% 
% Ts = 2;
% theta = NaN(20, num_cores);
% sigv = NaN(num_cores, 1);
% sigw = NaN(num_cores, 1);
% T_pred = zeros(size(Y,1),num_cores);
% T_pred2 = zeros(size(Y,1),num_cores);
% residuals = T_pred;
% eigenvalues = NaN(num_cores, model_order);
% for core = 1:num_cores
%     [theta(:,core), sigv(core), sigw(core)] = misoarxbls(U, Y(:,core), model_order, []);
%  
%     A = [-theta(1:model_order,core), [eye(model_order-1); zeros(1,model_order-1)]];
%     B = reshape(theta(model_order+1:end,core), model_order, num_cores+1);
%     C = [1, zeros(1,model_order-1)];
%     Q = sigw(core) * (C'*C);
%     R = sigv(core);
%     Xi = zeros(model_order,1);
%     Pi = eye(model_order);% * var(Y(:,core));
%     eigenvalues(core,:) = eig(A);
%     
% %     sys = ss(A, B, C, 0, Ts);
% % %     [kalmf,L,P] = kalman(sys, sigw(core), sigv(core), 0);
% %     y = lsim(sys, pow_part_re, 0:Ts:Ts*(size(temperatures,1)-1), Xi);
%     
%     ind_time = 1;
%     [Xi, Pi, yk, res] = KalmanFilter(A, eye(model_order), C, -A*Pi*A.' + Pi, R, Xi, Pi, Y(1:model_order,core), Y(ind_time,core));
%     T_pred(ind_time,core) = yk(1);
%     residuals(ind_time,core) = res(1);
%     for ind_time = 2:size(T_pred,1)
%         [Xi, Pi, yk, res] = KalmanFilter(A, B, C, Q, R, Xi, Pi, U(ind_time-1,:)', Y(ind_time,core));
%         T_pred(ind_time,core) = yk(1);
%         residuals(ind_time,core) = res(1);
%     end
%     
%     %[INN,K,Pp,Pf,Yp,Yf,Xp,Xf] = Kalm(A,B,C,C',U,Y(:,core),Q(1),R);
% end

%% Dati Galileo di produzione

%load('/Users/FedeAdm/Documents/Postdoc Bologna/HPC Applications/power-models/outputs/PowPartTherm_node062-0.mat')
%load('/Users/FedeAdm/Documents/Postdoc Bologna/HPC Applications/power-models/outputs/PowPartTherm_node063-0.mat')
load('/Users/FedeAdm/Documents/Postdoc Bologna/HPC Applications/power-models/outputs/PowPartTherm_node065-0.mat')
%load('/Users/FedeAdm/Documents/Postdoc Bologna/HPC Applications/power-models/outputs/PowPartTherm_node067-0.mat')
%load('/Users/FedeAdm/Documents/Postdoc Bologna/HPC Applications/power-models/outputs/PowPartTherm_node326-0.mat')
num_win = 25;
num_smallwin = 20;  % per window

ind_sock = 1;
%file_name = "../Results ThermMod backup Galileo reg1e-3 Tr3days/PredLinRegThermMod_node062.csv";
%file_name = "../Results ThermMod backup Galileo reg1e-3 Tr3days/PredLinRegThermMod_node063.csv";
file_name = "../Results ThermMod backup Galileo reg1e-3 Tr3days/PredLinRegThermMod_node065.csv";
%file_name = "../Results ThermMod backup Galileo reg1e-3 Tr3days/PredLinRegThermMod_node067.csv";
%file_name = "../Results ThermMod backup Galileo reg1e-3 Tr3days/PredLinRegThermMod_node326.csv";
opts = detectImportOptions(file_name,'NumHeaderLines',1,'FileType','text');
input_data = readtable(file_name, opts);
labels = table2array(input_data(:,2+(1:num_sock)));
predictions = table2array(input_data(:,2+num_sock+(1:num_sock)));

pow_part = [sum(pow_part(:,1:2),2), pow_part(:,3:end)];
pow_part_re = repmat(labels(:,ind_sock) ./ predictions(:,ind_sock), 1,9) .* pow_part;
temperatures = temperatures - temp_amb';
pow_part = pow_part(1:end-1,:) * lab_std;
pow_part_re = pow_part_re(1:end-1,:) * lab_std;
temperatures = temperatures(2:end,:) * temp_std;
win_length = floor(size(temperatures,1) / num_win);
win_length_small = floor(win_length / num_smallwin);

% Normalisations
pow_mean = mean(pow_part,1);
pow_re_mean = mean(pow_part_re,1); 
temp_mean = mean(temperatures, 1);
% pow_mean = [0.16, 0.06 * ones(1,num_cores)];
% pow_re_mean = [0.16, 0.06 * ones(1,num_cores)]; 
% temp_mean = [0.55, 0.36 * ones(1,num_cores)];
pow_part = pow_part - pow_mean;
pow_part_re = pow_part_re - pow_re_mean;
temperatures = temperatures - temp_mean;


% Calculate parameters and noise matrices in all windows
pow_touse = pow_part;
theta_init = []; %[-0.561393504094913;-0.375334504166025;-0.00732070967008602;0.0508844797574530;0.0402305015826753;-0.00498611610057687;-0.0381380143686845;0.0662688385933748;-0.0344786836329175;0.0714024208606131;-0.0910585379529112;0.141485847318916;-0.0272539224918401;0.0681830069080370;-0.112129484836070;0.146741377537944;-0.297567057824103;0.329840644135576;-0.227737248931772;0.255391121002097];
theta = NaN(10*model_order, num_cores, num_win);
sigv = NaN(num_cores, num_win);
sigw = NaN(num_cores, num_win);
eigenvalues = NaN(num_cores, model_order, num_win);
res_ident = NaN(win_length-model_order,num_cores, num_win);
cond_ident = NaN(num_cores, num_win);
svd_ident = NaN(20, num_cores, num_win);
corr_res = NaN(21, num_cores, num_win);
for core = 1:num_cores
    for ind_win = 1:num_win
        ind_curr = (ind_win-1)*win_length + (1:win_length);
        [theta(:,core,ind_win), sigv(core,ind_win), sigw(core,ind_win), res_ident(:,core,ind_win), ...
            cond_ident(core,ind_win), svd_ident(:,core,ind_win)] = ...
            misoarxbls(pow_touse(ind_curr,:), temperatures(ind_curr, 1+core), model_order, theta_init);
        
        A = [-theta(1:model_order,core,ind_win), [eye(model_order-1); zeros(1,model_order-1)]];
        eigenvalues(core,:,ind_win) = eig(A);
        corr_res(:,core,ind_win) = autocorr(res_ident(:,core,ind_win));
    end
end

% % Calculate temperature variance in small windows and choose for the
% % identification the window with highest third momentum and stable and real
% % eigenvalues
% varT = NaN(num_cores,num_win,num_smallwin);
% meanT = NaN(num_cores,num_win,num_smallwin);
% varT_win = NaN(num_cores,num_win);
% meanT_win = NaN(num_cores,num_win);
% moment3 = NaN(num_cores,num_win);
% for ind_win = 1:num_win
%     ind_curr = (ind_win-1)*win_length + (1:win_length);
%     moment3(:,ind_win) = moment(temperatures(ind_curr, 2:end), 3, 1);
%     varT_win(:,ind_win) = var(temperatures(ind_curr, 2:end), 0, 1);
%     meanT_win(:,ind_win) = mean(temperatures(ind_curr, 2:end), 1);
%     
%     ind_curr = reshape(ind_curr(1:win_length_small*num_smallwin), win_length_small, num_smallwin);
%     for ind_s = 1:num_smallwin
%         meanT(:,ind_win,ind_s) = mean(temperatures(ind_curr(:,ind_s), 2:end), 1);
%         varT(:,ind_win,ind_s) = var(temperatures(ind_curr(:,ind_s), 2:end), 0, 1);
%     end
% end
% sigmaT_avgwin = mean(sqrt(varT),3);
% %sigmaT_avgwin = sqrt(var(meanT, 0, 3));
% %[~,wins_order] = sort(mean(sigmaT_avgwin,1)', 'descend');
% [~,wins_order] = sort(mean(abs(moment3),1)', 'ascend');
% win_found = false;
% win_ind = 1;
% while ~win_found
%     win_ident_curr = wins_order(win_ind);
%     eig_curr = eigenvalues(:,:,win_ident_curr);
%     if nnz(abs(imag(eig_curr))) == 0 && nnz(abs(eig_curr) >= 1) == 0 && mean(abs(eig_curr(:,1))) > 0.8
%         win_ident = win_ident_curr;
%         win_found = true;
%     else
%         if win_ind == length(wins_order)
%             error('No window for identification found')
%         end
%         win_ind = win_ind + 1;
%     end
% end
% %[~,win_ident] = max(mean(sigmaT_win,1));


% Calculate package power average in small windows and autocorrelation, and
% then choose the best window for identification
peak_th = 20;
pow_pkg = labels(:,ind_sock) * lab_std;
meanP = NaN(num_win,num_smallwin);
corr_all = NaN(num_win,floor(win_length/2));
corr_allsw = NaN(num_win,floor(num_smallwin/2));
for ind_win = 1:num_win
    ind_curr = (ind_win-1)*win_length + (1:win_length);    
    
    % Power autocorrelation (entire trace)
    for indc = 1:floor(win_length/2)
        corr_all(ind_win,indc) = corr(pow_pkg(ind_curr(1:end-indc+1)), pow_pkg(ind_curr(indc:end)));
    end
    
    % Average power in small windows
    ind_curr = reshape(ind_curr(1:win_length_small*num_smallwin), win_length_small, num_smallwin);
    for ind_s = 1:num_smallwin
        meanP(ind_win,ind_s) = mean(pow_pkg(ind_curr(:,ind_s)), 1);
    end
    
    % Power autocorrelation (based on averages in small windows)
    for indc = 1:floor(num_smallwin/2)
        corr_allsw(ind_win,indc) = corr(meanP(ind_win,1:end-indc+1)', meanP(ind_win,indc:end)');
    end
end
wins_good = find(max(meanP,[],2) - min(meanP,[],2) > peak_th);
%[~,wins_order] = sort(abs(sum(corr_all, 2)/floor(win_length/2)), 'ascend');
%[~,wins_order] = sort(abs(sum(corr_allsw, 2)/floor(num_smallwin/2)), 'ascend');
[~,corr_ind] = min(abs(corr_all - 0.5), [], 2);
[~,wins_order] = sort(corr_ind, 'ascend');
wins_order = wins_order(ismember(wins_order, wins_good));
win_found = false;
win_ind = 1;
while ~win_found
    win_ident_curr = wins_order(win_ind);
    eig_curr = eigenvalues(:,:,win_ident_curr);
    if nnz(abs(imag(eig_curr))) == 0 && nnz(abs(eig_curr) >= 1) == 0 %&& mean(abs(eig_curr(:,1))) > 0.8
        win_ident = win_ident_curr;
        win_found = true;
    else
        if win_ind == length(wins_order)
            error('No window for identification found')
        end
        win_ind = win_ind + 1;
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
temperatures = (temperatures + temp_mean) + temp_amb(2:end)' * temp_std;
T_pred = (T_pred + temp_mean(2:end)) + temp_amb(2:end)' * temp_std;
T_filt = (T_filt + temp_mean(2:end)) + temp_amb(2:end)' * temp_std;
err_pred = temperatures(:,2:9) - T_pred;
err_filt = temperatures(:,2:9) - T_filt;


%% Plots

line_width = 2;
plot_color = 'w';
font_size_legend = 14;
font_size = 14;
colours_default = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

% Cumulative errors per core using one window
figure
temp = cell(0);
ind_win = 14; %win_ident;
for core = 1:num_cores
    if ~isnan(err_pred(1,core))
    [cumerrordistrib_prob, cumerrordistrib] = ecdf(abs(err_pred(:,core,ind_win)));
    semilogx(cumerrordistrib, cumerrordistrib_prob, '-', 'LineWidth',line_width, 'Color', colours_default(mod(core-1,7)+1,:))
    hold on
    temp{end+1} = num2str(core);
    end
end
lgnd=legend(temp,'Location','NorthWest');
set(lgnd,'Fontsize',font_size_legend,'Color','none','TextColor',plot_color,'EdgeColor',plot_color)
xlabel('Absolute core temperature error [°C]')
ylabel('Cumulative distribution')
grid on
set(gca,'FontSize',font_size, 'LineWidth',line_width)
set(gca,'Color','none','XColor',plot_color,'YColor',plot_color)
xlim([1e-2,1e1])
set(gcf,'Position',[10,1000,500,250])


% Temperature time traces in one core and window
figure
core = 1;
plot((1:size(T_pred,1))/1800, T_pred(:,core,ind_win), '-')
hold on
plot((1:size(T_pred,1))/1800, temperatures(:,1+core), '-')
lgnd=legend({'predicted','measured'},'Location','NorthWest');
set(lgnd,'Fontsize',font_size_legend,'Color','none','TextColor',plot_color,'EdgeColor',plot_color)
xlabel('Time [hours]')
ylabel('Temperature [°C]')
set(gca,'FontSize',font_size, 'LineWidth',line_width)
set(gca,'Color','none','XColor',plot_color,'YColor',plot_color)
set(gcf,'Position',[10,1000,500,250])


% Mean errors per window
figure
semilogy(1:num_win, squeeze(mean(mean(abs(err_pred),2),1)), '-', 'LineWidth',line_width, 'Color', colours_default(1,:))
hold on
semilogy(1:num_win, squeeze(mean(sqrt(var(abs(err_pred),1,1)),2)), '-.', 'LineWidth',line_width, 'Color', colours_default(1,:))
% semilogy(1:num_win, squeeze(mean(mean(abs(err_filt),2),1)), '-', 'LineWidth',line_width, 'Color', colours_default(2,:))
% semilogy(1:num_win, squeeze(mean(sqrt(var(abs(err_filt),1,1)),2)), '-.', 'LineWidth',line_width, 'Color', colours_default(2,:))
semilogy([1,num_win], [1.2, 1.2], '--', 'LineWidth',line_width, 'Color', colours_default(3,:))
xlabel('Window used for identification')
ylabel('Average temperature error [°C]')
set(gca,'FontSize',font_size, 'LineWidth',line_width)
set(gca,'Color','none','XColor',plot_color,'YColor',plot_color)
ylim([1e-1,1e3])
set(gcf,'Position',[10,1000,500,250])


% Eigenvalues characteristics per window
figure
plot(1:num_win, squeeze(max(max(abs(eigenvalues),[],2),[],1)), '-', 'LineWidth',line_width)
hold on
plot(1:num_win, squeeze(min(max(abs(eigenvalues),[],2),[],1)), '-', 'LineWidth',line_width)
%plot(1:num_win, squeeze(max(max(abs(imag(eigenvalues)),[],2),[],1)), '-', 'LineWidth',line_width)
semilogy([1,num_win], [1, 1], '--', 'LineWidth',line_width)
lgnd=legend({'max module','min module'},'Location','NorthWest');
set(lgnd,'Fontsize',font_size_legend,'Color','none','TextColor',plot_color,'EdgeColor',plot_color)
xlabel('Window used for identification')
ylabel('Eigenvalues min or max')
set(gca,'FontSize',font_size, 'LineWidth',line_width)
set(gca,'Color','none','XColor',plot_color,'YColor',plot_color)
set(gcf,'Position',[10,1000,500,250])


% Condition number
figure
semilogy(1:num_win, squeeze(mean(cond_ident,1)), '-', 'LineWidth',line_width)
hold on
semilogy(1:num_win, squeeze(sqrt(mean(cond_ident,1))), '-', 'LineWidth',line_width)
lgnd=legend({'average','sqrt'},'Location','NorthWest');
set(lgnd,'Fontsize',font_size_legend,'Color','none','TextColor',plot_color,'EdgeColor',plot_color)
xlabel('Window used for identification')
ylabel('Identification condition number')
set(gca,'FontSize',font_size, 'LineWidth',line_width)
set(gca,'Color','none','XColor',plot_color,'YColor',plot_color)
set(gcf,'Position',[10,1000,500,250])


% Residual autocorrelation
figure
plot(1:num_win, squeeze(mean(sum(corr_res(2:end,:,:),1),2)), '-', 'LineWidth',line_width)
xlabel('Window used for identification')
ylabel('Identification condition number')
set(gca,'FontSize',font_size, 'LineWidth',line_width)
set(gca,'Color','none','XColor',plot_color,'YColor',plot_color)
set(gcf,'Position',[10,1000,500,250])
