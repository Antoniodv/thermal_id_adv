%% Plots of investigations on Diversi's thermal model for TCAD paper

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

file_path = "../outputs/KalmanIdentRes/";
listing = dir(char(file_path));
file_names = "";
for indf = 1:length(listing)
    file_names = [file_names, string(listing(indf).name)];
end
files_data = file_names(contains(file_names,"KalmanIdentAll"))';
nodespkg = extractBetween(files_data,"node",".mat");
files_algs = extractBetween(files_data,"KalmanIdentAll_","_node");
num_calc = length(files_data);  %length(files_data);
T_pred = cell(num_calc,1);
errors_abs = cell(num_calc,1);
theta_all = cell(num_calc,1);
sigw_all = cell(num_calc,1);
sigv_all = cell(num_calc,1);
for indf = 1:num_calc
    input_data = load(file_path + files_data(indf));
    errors_abs{indf} = input_data.temperatures_all(:,2:end) - input_data.T_pred;
    temperatures = input_data.temperatures_all(:,2:end);
    T_pred{indf} = input_data.T_pred;
    
    if contains(files_data(indf), "Combined")
        theta_all{indf} = input_data.theta_all;
        sigw_all{indf} = input_data.sigw_all;
        sigv_all{indf} = input_data.sigv_all;
    else
        theta_all{indf} = input_data.theta;
        sigw_all{indf} = input_data.sigw;
        sigv_all{indf} = input_data.sigv;
    end
end
num_win = size(errors_abs{indf},3);


% Calculate average errors per window
sepwins = (files_algs ~= "Combined");
Terr_mean_sepwin = zeros(num_calc, num_cores, num_win);
Terr_mean_comb = zeros(num_calc, num_cores);
for indf = 1:num_calc
    if sepwins(indf)
        Terr_mean_sepwin(indf,:,:) = mean(abs(errors_abs{indf}), 1, 'omitnan');
    else
        Terr_mean_comb(indf,:) = mean(abs(errors_abs{indf}), 1, 'omitnan');
    end
end


%% Import Diversi's IECON data

input_data = load(file_path + 'Diversi_IECON2016_prediction_error.mat');
sampl_time = 2;
Terr_IECON = input_data.PRED_ERR;
poles_IECON = [0.846, 0.024;
    0.841,0.007;
    0.861,0.005;
    0.880,0.020;
    0.863,0.021;
    0.904,0.031;
    0.882,0.006;
    0.879,0.011];
tau_IECON = -sampl_time ./ log(poles_IECON);


%% Plots for us

core = 1;

% Plot temperature errors

figure
temp = cell(0);
semilogy(squeeze(Terr_mean_sepwin(sepwins,core,:))', 'LineWidth',line_width)
hold on
semilogy([1,num_win], repmat(squeeze(Terr_mean_comb(~sepwins,core))',1,2), '--', 'LineWidth',line_width)
% yl = ylim;
% semilogx(repmat(pow_th',1,2)', repmat(yl,2,1)', '-.', 'LineWidth',line_width, 'Color', plot_color)
lgnd=legend([files_algs(sepwins); files_algs(~sepwins)],'Location','NorthWest', 'Interpreter', 'none');
set(lgnd,'Fontsize',font_size_legend,'Color','none','TextColor',plot_color,'EdgeColor',plot_color)
xlabel('Window number')
ylabel('Absolute core temperature error [°C]')
%grid on
set(gca,'FontSize',font_size, 'LineWidth',line_width)
set(gca,'Color','none','XColor',plot_color,'YColor',plot_color)
%xlim([1e-2,1e1])


% Error distributions for one window

win = 22;

figure
temp = cell(0);
for indf = 2:num_calc
    [cumerrordistrib_prob, cumerrordistrib] = ecdf(abs(errors_abs{indf}(:,core,win)));
    semilogx(cumerrordistrib, cumerrordistrib_prob, '-', 'LineWidth',line_width)
    hold on
end
indf = 1;
[cumerrordistrib_prob, cumerrordistrib] = ecdf(abs(errors_abs{indf}(:,core)));
semilogx(cumerrordistrib, cumerrordistrib_prob, '--', 'LineWidth',line_width)

lgnd=legend([files_algs(sepwins); files_algs(~sepwins)],'Location','NorthWest', 'Interpreter', 'none');
set(lgnd,'Fontsize',font_size_legend,'Color','none','TextColor',plot_color,'EdgeColor',plot_color)
xlabel('Absolute core temperature error [°C]')
ylabel('Cumulative distribution')
grid on
set(gca,'FontSize',font_size, 'LineWidth',line_width)
set(gca,'Color','none','XColor',plot_color,'YColor',plot_color)
xlim([1e-2,1e1])


% Temperature traces for the same window


figure
temp = cell(0);
for indf = 2:num_calc
    plot(T_pred{indf}(:,core,win), '-')
    hold on
end
indf = 1;
plot(T_pred{indf}(:,core), '--')
plot(temperatures(:,core), '-.')

lgnd=legend([files_algs(sepwins); files_algs(~sepwins)],'Location','NorthWest', 'Interpreter', 'none');
set(lgnd,'Fontsize',font_size_legend,'Color','none','TextColor',plot_color,'EdgeColor',plot_color)
xlabel('Time index')
ylabel('Temperature [°C]')
set(gca,'FontSize',font_size, 'LineWidth',line_width)
set(gca,'Color','none','XColor',plot_color,'YColor',plot_color)


% Error distributions for one window (comparison with IECON)

win = 22;

figure
temp = cell(0);
semilogx(1e-2,0, '-k', 'LineWidth',line_width)
hold on
semilogx(1e-2,0, '--k', 'LineWidth',line_width)
indf = 4;
for core = 1:num_cores
    [cumerrordistrib_prob, cumerrordistrib] = ecdf(abs(errors_abs{indf}(:,core,win)));
    semilogx(cumerrordistrib, cumerrordistrib_prob, '-', 'LineWidth',line_width)
    hold on
end
for core = 1:num_cores
    [cumerrordistrib_prob, cumerrordistrib] = ecdf(abs(Terr_IECON(:,core)));
    semilogx(cumerrordistrib, cumerrordistrib_prob, '--', 'LineWidth',line_width)
    hold on
end

lgnd=legend(["Current work", "Ref. [25]"],'Location','NorthWest', 'Interpreter', 'none');
set(lgnd,'Fontsize',font_size_legend,'Color','none','TextColor',plot_color,'EdgeColor',plot_color)
xlabel('Absolute core temperature error [°C]')
ylabel('Cumulative distribution')
grid on
set(gca,'FontSize',font_size, 'LineWidth',line_width)
set(gca,'Color','none','XColor',plot_color,'YColor',plot_color)
xlim([1e-2,1e1])
set(gcf,'Position',[10,1000,500,250])


%% Plots for TCAD review

% Plot temperature errors for all windows

figure
indf = 4;
semilogy(squeeze(mean(Terr_mean_sepwin(indf,:,:),2)), 'LineWidth',line_width)
hold on
indf = 2;
semilogy(squeeze(mean(Terr_mean_sepwin(indf,:,:),2)), 'LineWidth',line_width)
semilogy([1,num_win], [1.2, 1.2], '--', 'LineWidth',line_width, 'Color', colours_default(3,:))

lgnd=legend(["Powers", "Temperature and powers", "1.2°C threshold"], 'Interpreter', 'none');
set(lgnd,'Fontsize',font_size_legend,'Color','none','TextColor',plot_color,'EdgeColor',plot_color)
xlabel('Window number')
ylabel('Absolute core temperature error [°C]')
%grid on
set(gca,'FontSize',font_size, 'LineWidth',line_width)
set(gca,'Color','none','XColor',plot_color,'YColor',plot_color)
ylim([1e-1,1e3])


% Error distributions for one window (different models)

win = 22;

figure
temp = cell(0);
semilogx(1e-2,0, '-k', 'LineWidth',line_width)
hold on
semilogx(1e-2,0, '--k', 'LineWidth',line_width)
indf = 4;
for core = 1:num_cores
    [cumerrordistrib_prob, cumerrordistrib] = ecdf(abs(errors_abs{indf}(:,core,win)));
    semilogx(cumerrordistrib, cumerrordistrib_prob, '-', 'LineWidth',line_width)
    hold on
end
indf = 2;
for core = 1:num_cores
    [cumerrordistrib_prob, cumerrordistrib] = ecdf(abs(errors_abs{indf}(:,core,win)));
    semilogx(cumerrordistrib, cumerrordistrib_prob, '--', 'LineWidth',line_width)
    hold on
end

lgnd=legend(["Powers", "Temperature and powers"],'Location','NorthWest', 'Interpreter', 'none');
set(lgnd,'Fontsize',font_size_legend,'Color','none','TextColor',plot_color,'EdgeColor',plot_color)
xlabel('Absolute core temperature error [°C]')
ylabel('Cumulative distribution')
grid on
set(gca,'FontSize',font_size, 'LineWidth',line_width)
set(gca,'Color','none','XColor',plot_color,'YColor',plot_color)
xlim([1e-2,1e1])
set(gcf,'Position',[10,1000,500,250])


% Error distributions for two windows (comparison with IECON, all cores together)

figure
indf = 4;
win = 22;
num_data = size(errors_abs{indf}(:,:,win));
[cumerrordistrib_prob, cumerrordistrib] = ecdf(reshape(abs(errors_abs{indf}(:,:,win)), num_data(1)*num_data(2),1));
semilogx(cumerrordistrib, cumerrordistrib_prob, '-', 'LineWidth',line_width)
hold on
win = 16;
num_data = size(errors_abs{indf}(:,:,win));
[cumerrordistrib_prob, cumerrordistrib] = ecdf(reshape(abs(errors_abs{indf}(:,:,win)), num_data(1)*num_data(2),1));
semilogx(cumerrordistrib, cumerrordistrib_prob, '-', 'LineWidth',line_width)
num_data = size(Terr_IECON);
[cumerrordistrib_prob, cumerrordistrib] = ecdf(reshape(abs(Terr_IECON), num_data(1)*num_data(2),1));
semilogx(cumerrordistrib, cumerrordistrib_prob, '--', 'LineWidth',line_width)

lgnd=legend(["Current work, good win", "Current work, bad win", "Ref. [25]"],'Location','NorthWest', 'Interpreter', 'none');
set(lgnd,'Fontsize',font_size_legend,'Color','none','TextColor',plot_color,'EdgeColor',plot_color)
xlabel('Absolute core temperature error [°C]')
ylabel('Cumulative distribution')
grid on
set(gca,'FontSize',font_size, 'LineWidth',line_width)
set(gca,'Color','none','XColor',plot_color,'YColor',plot_color)
xlim([1e-2,1e1])
set(gcf,'Position',[10,1000,500,250])
