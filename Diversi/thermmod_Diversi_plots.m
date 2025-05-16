%% Plots of results processed with Spark on Galileo's backup for Diversi's power model

line_width = 2;
plot_color = 'k';
font_size_legend = 14;
font_size = 14;
colours_default = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

num_cores = 8;
num_sock = 2;
num_cpufeat = 3;
num_corefeat = 5;

freq_std = 2600;
lab_std = 80;
temp_std = 50;
ips_scale = 3e10;
avx_scale = 12;

model_order = 2;
num_add_eq = 10;


%% Import data from Python elaborations

file_path = "../outputs/";
listing = dir(char(file_path));
file_names = "";
for indf = 1:length(listing)
    file_names = [file_names, string(listing(indf).name)];
end
load('Kalman_Tpred.mat')
files_data = file_names(contains(file_names,"PowPartTherm"))';
files_data = files_data(1:length(T_pred));
nodespkg = extractBetween(files_data,"_",".mat");
num_nodpkg = length(T_pred);  %length(files_data);
temperatures = cell(num_nodpkg,1);
errors_abs = cell(num_nodpkg,1);
for indf = 1:num_nodpkg
    input_data = load(file_path + files_data(indf));
    temperatures{indf} = input_data.temperatures(2:end,2:end) * temp_std;
    T_pred{indf} = T_pred{indf};
    errors_abs{indf} = temperatures{indf} - T_pred{indf};
end


% % Exclude cores with unstable models
% for ind_nod = 1:num_nodpkg
%     unst_cores = abs(max(T_pred{ind_nod},[],1)) > 100;
%     T_pred{ind_nod}(:,unst_cores) = NaN;
%     errors_abs{ind_nod}(:,unst_cores) = NaN;
% end


% Count points below thresholds

term_th = [0.5, 1];  %[°C]
num_th = length(term_th);
percent_th = cell(num_nodpkg,1);

for ind_nod = 1:num_nodpkg
    percent_th{ind_nod} = zeros(num_cores, num_th);
    
    % Calculate errors per node and core
    for core = 1:num_cores
        percent_th{ind_nod}(core,:) = sum(abs(repmat(errors_abs{ind_nod}(:,core), 1,num_th)) ...
            <= term_th, 1) ./ sum(repmat(~isnan(errors_abs{ind_nod}(:,core)), 1,num_th), 1);
    end
end

% Calculate statistics on all nodes and packages
nod_part = [5, 10, 20, 30, 40]; %num_nodpkg];
num_part = length(nod_part);
pth_var = zeros(3, num_part, num_th);
ind_nod_part = cell(num_part,1);
for ind_part = 1:num_part
    ind_nod_part{ind_part} = randsample(1:num_nodpkg, nod_part(ind_part));
    pth_curr = [];
    for ind = 1:nod_part(ind_part)
        ind_nod = ind_nod_part{ind_part}(ind);
        cores_stable = ~isnan(percent_th{ind_nod}(:,1));
        pth_curr = [pth_curr; reshape(percent_th{ind_nod}(cores_stable,:), sum(cores_stable), num_th)];
        pth_var(1,ind_part,:) = prctile(pth_curr, 25, 1);
        pth_var(2,ind_part,:) = median(pth_curr, 1);
        pth_var(3,ind_part,:) = prctile(pth_curr, 75, 1);
    end
end


%% Plots

% Plot cumulated temperature errors

figure
temp = cell(0);
ind_nod = 64;
for core = 1:num_cores
    if ~isnan(errors_abs{ind_nod}(1,core))
        [cumerrordistrib_prob, cumerrordistrib] = ecdf(abs(errors_abs{ind_nod}(:,core)));
        semilogx(cumerrordistrib, cumerrordistrib_prob, '-', 'LineWidth',line_width, 'Color', colours_default(mod(core-1,7)+1,:))
        hold on
        temp{end+1} = num2str(core);
    end
end
% yl = ylim;
% semilogx(repmat(pow_th',1,2)', repmat(yl,2,1)', '-.', 'LineWidth',line_width, 'Color', plot_color)
lgnd=legend(temp,'Location','NorthWest');
set(lgnd,'Fontsize',font_size_legend,'Color','none','TextColor',plot_color,'EdgeColor',plot_color)
xlabel('Absolute core temperature error [°C]')
ylabel('Cumulative distribution')
grid on
set(gca,'FontSize',font_size, 'LineWidth',line_width)
set(gca,'Color','none','XColor',plot_color,'YColor',plot_color)
xlim([1e-2,1e1])


% Plot percentages of points below the thresholds varying number of nodes

figure
ind_th = 1;
pth_var_curr = reshape(pth_var(:,:,ind_th), 3, num_part) * 100;
h = errorbar(nod_part, pth_var_curr(2,:), pth_var_curr(2,:)-pth_var_curr(1,:), ...
    pth_var_curr(3,:)-pth_var_curr(2,:), '-o', 'LineWidth',line_width);
h.CapSize = 12;
hold on
ind_th = 2;
pth_var_curr = reshape(pth_var(:,:,ind_th), 3, num_part) * 100;
h = errorbar(nod_part, pth_var_curr(2,:), pth_var_curr(2,:)-pth_var_curr(1,:), ...
    pth_var_curr(3,:)-pth_var_curr(2,:), '--o', 'LineWidth',line_width);
h.CapSize = 12;
lgnd=legend({[num2str(term_th(1)),'°C'], [num2str(term_th(2)),'°C']},'Location','SouthEast');
set(lgnd,'Fontsize',font_size_legend,'Color','none','TextColor',plot_color,'EdgeColor',plot_color)
xlabel('Number of considered nodes and packages')
ylabel('Points below threshold [%]')
%grid on
ylim([0.3,1] * 100)
set(gca,'FontSize',font_size, 'LineWidth',line_width)
set(gca,'Color','none','XColor',plot_color,'YColor',plot_color)
set(gcf,'Position',[10,1000,500,250])
