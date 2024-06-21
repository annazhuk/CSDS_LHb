%% colors
green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

%% Fig S1A: HC correlations control
for xc = 1
clear 
load('ALL_fig_data.mat')
figure('DefaultAxesFontSize',20)
pos5 = [0.87 0.6 0.05 0.1];
subplot('Position',pos5)
sz = 30;
% fit line to data
A = HC_CTRL_SI;
B = HC_CTRL_SW;
mdl = fitlm(A,B,'linear');
h = plot(mdl);
hold on
% Get handles to plot components
dataHandle = findobj(h,'DisplayName','Data');
fitHandle = findobj(h,'DisplayName','Fit');
cbHandles = findobj(h,'DisplayName','Confidence bounds');
cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
dataHandle.Color = 'w'; 
fitHandle.Color = 'k';
fitHandle.LineWidth = 1;
fitHandle.LineStyle = '--';
set(cbHandles, 'Color', 'k', 'LineWidth', 1,'LineStyle','-')
scatter(HC_CTRL_SI,HC_CTRL_SW,sz,'filled','MarkerEdgeColor','k','MarkerEdgeAlpha',1,'MarkerFaceColor','k','MarkerFaceAlpha',.2)
xlabel('S       SI Time (%)     R')
ylabel('% time investigating')
title('HC SW') %,'FontSize',20)
xlim([20 85])
ylim([-10 100])
box off
end

%% Fig S1B: EPM control correlation
for xc = 1
clear
load('EPM_susVScontrol_stress.mat')
figure('DefaultAxesFontSize',12)
pos3 = [0.35 0.88 0.05 0.1];
subplot('Position',pos3)

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

% fit line to data
sz = 30;
A = opencenter_control_SI;
B = opensum_control;
mdl = fitlm(A,B,'linear');
h = plot(mdl);
hold on
% Get handles to plot components
dataHandle = findobj(h,'DisplayName','Data');
fitHandle = findobj(h,'DisplayName','Fit');
cbHandles = findobj(h,'DisplayName','Confidence bounds');
cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
dataHandle.Color = 'w'; 
fitHandle.Color = 'k';
fitHandle.LineWidth = 1;
fitHandle.LineStyle = '--';
set(cbHandles, 'Color', 'k', 'LineWidth', 1,'LineStyle','-')
scatter(opencenter_control_SI,opensum_control,sz,'filled','MarkerEdgeColor','k','MarkerEdgeAlpha',1,'MarkerFaceColor','k','MarkerFaceAlpha',.2)
xlabel('SI Time (%)')
ylabel({'% time spent';'in open arms'})
title('EPM') %,'FontSize',20)
xlim([35 75])
ylim([-5 40])
box off
end

%% Fig S1C: NSF control correlation
for xc = 1
clear
load('NSF_susVcontrol_stress.mat')
figure('DefaultAxesFontSize',12)
pos3 = [0.35 0.88 0.05 0.1];
subplot('Position',pos3)

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

% fit line to data
sz = 30;
A = NSF_control_SI;
B = NSF_control;
mdl = fitlm(A,B,'linear');
h = plot(mdl);
hold on
% Get handles to plot components
dataHandle = findobj(h,'DisplayName','Data');
fitHandle = findobj(h,'DisplayName','Fit');
cbHandles = findobj(h,'DisplayName','Confidence bounds');
cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
dataHandle.Color = 'w'; 
fitHandle.Color = 'k';
fitHandle.LineWidth = 1;
fitHandle.LineStyle = '--';
set(cbHandles, 'Color', 'k', 'LineWidth', 1,'LineStyle','-')
scatter(NSF_control_SI,NSF_control,sz,'filled','MarkerEdgeColor','k','MarkerEdgeAlpha',1,'MarkerFaceColor','k','MarkerFaceAlpha',.2)
xlabel('SI Time (%)')
ylabel('Latency to feed (s)')
title('NSF') %,'FontSize',20)
xlim([35 75])
ylim([-70 650])
box off
end

%% Fig S1D: Immobility control
for xc = 1
clear 
load('ALL_fig_data.mat')
figure('DefaultAxesFontSize',15)
pos5 = [0.55 0.6 0.05 0.1];
subplot('Position',pos5)

% fit line to data
A = cropped_CTRL_SI_all_wfos ;
B = immobility_ctrl_all.velocity;
mdl = fitlm(A,B,'linear');
h = plot(mdl);
hold on
% Get handles to plot components
dataHandle = findobj(h,'DisplayName','Data');
fitHandle = findobj(h,'DisplayName','Fit');
cbHandles = findobj(h,'DisplayName','Confidence bounds');
cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
dataHandle.Color = 'w'; 
fitHandle.Color = 'k';
fitHandle.LineWidth = 1;
fitHandle.LineStyle = '--';
sz=30;
set(cbHandles, 'Color', 'k', 'LineWidth', 1,'LineStyle','-')
scatter(cropped_CTRL_SI_all_wfos,immobility_ctrl_all.velocity,sz,'filled','MarkerEdgeColor','k','MarkerEdgeAlpha',1,'MarkerFaceColor','k','MarkerFaceAlpha',.2)
xlabel('S       SI Time (%)     R')
ylabel({'% time spent';'immobile'})
title('Chamber Exploration') %,'FontSize',20)
xlim([20 85])
ylim([-5 20])
box off
end

%% Fig S1E: Immobility PRE correlation
for xc = 1
clear
load('immobility_PRE.mat')
figure('DefaultAxesFontSize',15)
pos5 = [0.55 0.6 0.05 0.1];
subplot('Position',pos5)

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

SI_immobility_sus = cropped_SI_all(cropped_SI_all<med); 
SI_immobility_res = cropped_SI_all(cropped_SI_all>=med); 

% fit line to data
A = [SI_immobility_res;SI_immobility_sus];
B = [immobility_res_all.velocity;immobility_sus_all.velocity];
mdl = fitlm(A,B,'linear');
h = plot(mdl);
hold on
% Get handles to plot components
dataHandle = findobj(h,'DisplayName','Data');
fitHandle = findobj(h,'DisplayName','Fit');
cbHandles = findobj(h,'DisplayName','Confidence bounds');
cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
dataHandle.Color = 'w'; 
fitHandle.Color = 'k';
fitHandle.LineWidth = 1;
fitHandle.LineStyle = '--';
sz=30;
set(cbHandles, 'Color', 'k', 'LineWidth', 1,'LineStyle','-')
scatter(SI_immobility_sus,immobility_sus_all.velocity,sz,'filled','MarkerEdgeColor',green,'MarkerEdgeAlpha',1,'MarkerFaceColor',green,'MarkerFaceAlpha',.2)
scatter(SI_immobility_res,immobility_res_all.velocity,sz,'filled','MarkerEdgeColor',purple,'MarkerEdgeAlpha',1,'MarkerFaceColor',purple,'MarkerFaceAlpha',.2)
xlabel('S       SI Time (%)     R')
ylabel({'% time spent';'immobile'})
title('Chamber Exploration Pre-CSDS') %,'FontSize',20)
xticks(0:40:100)
yticks(0:20:100)
xlim([-5 100])
%ylim([-5 50])
box off
end

%% Fig S1F: OFT control correlation
for xc = 1
clear
load('OFT_susVScontrol_stress.mat')
figure('DefaultAxesFontSize',12)
pos3 = [0.35 0.88 0.05 0.1];
subplot('Position',pos3)

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

% fit line to data
sz = 30;
A = OFT_control_SI;
B = OFT_control;
mdl = fitlm(A,B,'linear');
h = plot(mdl);
hold on
% Get handles to plot components
dataHandle = findobj(h,'DisplayName','Data');
fitHandle = findobj(h,'DisplayName','Fit');
cbHandles = findobj(h,'DisplayName','Confidence bounds');
cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
dataHandle.Color = 'w'; 
fitHandle.Color = 'k';
fitHandle.LineWidth = 1;
fitHandle.LineStyle = '--';
set(cbHandles, 'Color', 'k', 'LineWidth', 1,'LineStyle','-')
scatter(OFT_control_SI,OFT_control,sz,'filled','MarkerEdgeColor','k','MarkerEdgeAlpha',1,'MarkerFaceColor','k','MarkerFaceAlpha',.2)
xlabel('SI Time (%)')
ylabel('Time spent in inner zone')
title('OFT') %,'FontSize',20)
xlim([35 75])
%ylim([-5 40])
box off
end

%% Fig S2D: Plot some of the shuffles
for xc = 1
clear
load('example_shuffles_plot.mat')
% m569
figure('DefaultAxesFontSize',12)
pos8 = [0.2 0.54 0.3 0.1];
subplot('Position',pos8)
numcell = 5; 
time_len = size(C_raw_shuffles,2);
time_plot = 0:fps:time_len/fps;
imagesc(squeeze(C_raw_shuffles(numcell,:,1:10))')
colorbar
clim([-2 2])
a=colorbar;
xticks(0:fps*30:time_len)
xticklabels(time_plot)
ylabel(a,'\DeltaF/F (Z)','FontSize',12,'Rotation',90);
xlabel('Time (s)')
ylabel('Shuffles')
title('Example cell (shifted 5s)')
box off
end

%% Fig S2E: Plot timelocked traces
for xc = 1
clear
load('example_cell_sig_shuffles.mat')

numcell = 5; %m569
time_plot = -2:2:5;
figure('DefaultAxesFontSize',12)
pos8 = [0.2 0.54 0.05 0.1];
subplot('Position',pos8)
for s = 1:100%size(C_raw_shuffles,3)
    plot(squeeze(avg_start_shuffle_all(numcell,:,s)),'Color',[.7 .7 .7],'LineWidth',0.25)
    hold on
end
plot(avg_start_all(numcell,:),'r','LineWidth',1)
plot([60,60],[-2.5 4.5],'k--','linewidth',1)
xlim([0 210])
ylim([-2 2.5])
xticks(0:60:210)
%yticks(-0.5:.5:1.5)
xticklabels(time_plot)
xlabel('Time from social zone entry (s)')
ylabel('Mean \DeltaF/F (Z)')
title('Example Cell vs 100 shuffles')
box off
end

%% Fig S2F: Plot a histogram of shuffled peak averages vs real peak average
for xc = 1
clear
load('example_cell_histogram.mat')
numcell = 5; %4;
figure('DefaultAxesFontSize',12)
pos8 = [0.2 0.54 0.1 0.1];
subplot('Position',pos8)
y_val = 150;
histogram(pk_shuffles(numcell,:),'FaceColor',[.7 .7 .7])
hold on
plot([pk_real(numcell),pk_real(numcell)],[0,y_val],'r','LineWidth',1)
plot([for_up_plot,for_up_plot],[0,y_val])
plot([for_down_plot,for_down_plot],[0,y_val])
xlabel('Mean \DeltaF/F (Z)')
title('Example cell vs null distribution')
box off
end

%% Fig S2G: Avoidance vs magnitude: Imaging (all cells)
for xc = 1
clear
load('W:\Anna\Social_Defeat\CSDS_Imaging_Cohort1\SI_Post\SW\nose_sig_cells_res_sus_FIXED.mat')

% colors
green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

sz = 30;
figure('Renderer', 'painters', 'Position', [10 10 240 140],'DefaultAxesFontSize',12)
% sus
% sort mice by SI score
s_SI = cropped_SI_im_sorted( ismember(cropped_ID_im_sorted,sus_im.mice));
[val,idx] = sort(s_SI);
for m = 1:size(sus_im.mice,1)
    %swarmchart(m+zeros(size(susSWPost_im.nose_plot_start_Favg{idx(m)},1),1),susSWPost_im.nose_plot_start_Favg{idx(m)},sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',.2)
    swarmchart(m+zeros(size(sus_im.pk_real_sig_all{idx(m)},1),1),sus_im.pk_real_sig_all{idx(m)},sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor','w','MarkerFaceAlpha',.2)
    hold on
    if ~isnan(sus_im.sig_up_all{idx(m)})
        swarmchart(m+zeros(length(sus_im.pk_real_sig_all{idx(m)}(sus_im.sig_up_all{idx(m)})),1),sus_im.pk_real_sig_all{idx(m)}(sus_im.sig_up_all{idx(m)}),sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',.7)
        swarmchart(m+zeros(length(sus_im.pk_real_sig_all{idx(m)}(sus_im.sig_down_all{idx(m)})),1),sus_im.pk_real_sig_all{idx(m)}(sus_im.sig_down_all{idx(m)}),sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',.7)
    end
    %susSWPost_im.mice(idx(m))
end

ylabel('Mean \DeltaF/F (Z)')
%xticks(3:6)

% res
% sort mice by SI score
r_SI = cropped_SI_im_sorted(ismember(cropped_ID_im_sorted,res_im.mice));
[val,idx] = sort(r_SI);
for m = 1:size(res_im.mice,1)
   %swarmchart(6+m+zeros(size(resSWPost_im.nose_plot_start_Favg{idx(m)},1),1),resSWPost_im.nose_plot_start_Favg{idx(m)},sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',.2
   swarmchart(6+m+zeros(size(res_im.pk_real_sig_all{idx(m)},1),1),res_im.pk_real_sig_all{idx(m)},sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor','w','MarkerFaceAlpha',.2)
   hold on
   if ~isnan(res_im.sig_up_all{idx(m)})
        swarmchart(6+m+zeros(length(res_im.pk_real_sig_all{idx(m)}(res_im.sig_up_all{idx(m)})),1),res_im.pk_real_sig_all{idx(m)}(res_im.sig_up_all{idx(m)}),sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',.7)
        swarmchart(6+m+zeros(length(res_im.pk_real_sig_all{idx(m)}(res_im.sig_down_all{idx(m)})),1),res_im.pk_real_sig_all{idx(m)}(res_im.sig_down_all{idx(m)}),sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',.7)
    end
   %susSWPost_im.mice(idx(m))
end
title('SW SI Post-CSDS (Single-Cell)')
box off
end

%% Fig S2H (Left): CDF of inter-cell synchrony pre-CSDS (all mice)
for xc = 1
clear;
load('W:\Anna\Social_Defeat\CSDS_Imaging_Cohort1\SI_Pre\SW\correlations_res_sus_S.mat')
% exclude mice with too few cells
cutoff = 9;
sus_im.activity_avg_ex = sus_im.activity_avg(sus_im.total_cells>=cutoff,:);
res_im.activity_avg_ex = res_im.activity_avg(res_im.total_cells>=cutoff,:);

% plot CDF of each mouse
% colors
green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

figure('Renderer', 'painters', 'Position', [10 10 130 110],'DefaultAxesFontSize',12)

bin_size = 1;

% sus
for m = 1:size(sus_im.activity_avg_ex,1)
    bins_sus = min(sus_im.activity_avg_ex(m,:)):bin_size:max(sus_im.activity_avg_ex(m,:));
    cdf_sus = cumsum(histcounts(sus_im.activity_avg_ex(m,:),length(bins_sus)))/length(sus_im.activity_avg_ex(m,:))*100;

    sus_im.bin_sus_all{m,1} = bins_sus;
    sus_im.cdf_sus_all{m,1} = cdf_sus;

    sus_im.lengths(m,1) = length(bins_sus);
    
    plot(bins_sus,cdf_sus,'color',green,'linewidth',2)
    hold on
end

% res
for m = 1:size(res_im.activity_avg_ex,1)
    bins_res = min(res_im.activity_avg_ex(m,:)):bin_size:max(res_im.activity_avg_ex(m,:));
    cdf_res = cumsum(histcounts(res_im.activity_avg_ex(m,:),length(bins_res)))/length(res_im.activity_avg_ex(m,:))*100;

    res_im.bin_res_all{m,1} = bins_res;
    res_im.cdf_res_all{m,1} = cdf_res;

    res_im.lengths(m,1) = length(bins_res);
    
    plot(bins_res,cdf_res,'color',purple,'linewidth',2)
    hold on
end

xlabel('% Population active')
ylabel('Cumulative probability')
box off
end

%% Fig S2H (Right): CDF of inter-cell synchrony pre-CSDS (average)
for xc = 1
clear;
load('W:\Anna\Social_Defeat\CSDS_Imaging_Cohort1\SI_Pre\SW\correlations_res_sus_S.mat')
% exclude mice with too few cells
cutoff = 9;
sus_im.activity_avg_ex = sus_im.activity_avg(sus_im.total_cells>=cutoff,:);
res_im.activity_avg_ex = res_im.activity_avg(res_im.total_cells>=cutoff,:);

% calculate CDF for each mouse
bin_size = 1;

% sus
for m = 1:size(sus_im.activity_avg_ex,1)
    bins_sus = min(sus_im.activity_avg_ex(m,:)):bin_size:max(sus_im.activity_avg_ex(m,:));
    cdf_sus = cumsum(histcounts(sus_im.activity_avg_ex(m,:),length(bins_sus)))/length(sus_im.activity_avg_ex(m,:))*100;

    sus_im.bin_sus_all{m,1} = bins_sus;
    sus_im.cdf_sus_all{m,1} = cdf_sus;

    sus_im.lengths(m,1) = length(bins_sus);
end

% res
for m = 1:size(res_im.activity_avg_ex,1)
    bins_res = min(res_im.activity_avg_ex(m,:)):bin_size:max(res_im.activity_avg_ex(m,:));
    cdf_res = cumsum(histcounts(res_im.activity_avg_ex(m,:),length(bins_res)))/length(res_im.activity_avg_ex(m,:))*100;

    res_im.bin_res_all{m,1} = bins_res;
    res_im.cdf_res_all{m,1} = cdf_res;

    res_im.lengths(m,1) = length(bins_res);
end

% Calculate average 
% sus
max_bins = max(sus_im.lengths);
for m = 1:size(sus_im.activity_avg_ex,1)
    diff_length = max_bins - length(sus_im.cdf_sus_all{m});
    sus_im.cdf_sus_matrix(m,:) = [sus_im.cdf_sus_all{m},nan(1,diff_length)];
end
sus_im.cdf_sus_matrix_avg = nanmean(sus_im.cdf_sus_matrix,1);
sus_im.cdf_sus_matrix_std = nanstd(sus_im.cdf_sus_matrix,1)/sqrt(size(sus_im.cdf_sus_matrix,1));

% res
max_bins = max(res_im.lengths);
for m = 1:size(res_im.activity_avg_ex,1)
    diff_length = max_bins - length(res_im.cdf_res_all{m});
    res_im.cdf_res_matrix(m,:) = [res_im.cdf_res_all{m},nan(1,diff_length)];
end
res_im.cdf_res_matrix_avg = nanmean(res_im.cdf_res_matrix,1);
res_im.cdf_res_matrix_std = nanstd(res_im.cdf_res_matrix,1)/sqrt(size(res_im.cdf_res_matrix,1));

% colors
green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

figure('Renderer', 'painters', 'Position', [10 10 130 110],'DefaultAxesFontSize',12)

plot(sus_im.cdf_sus_matrix_avg,'color',green,'linewidth',2)
hold on
plot(res_im.cdf_res_matrix_avg,'color',purple,'linewidth',2)

xlim([0 100])
%xlim([-4 6])
xlabel('% Population active')
ylabel('Cumulative probability')
box off
end

%% EXTRA: Juvenile vs adult; Stress vs Control
for xc = 1
clear;
load('W:\Anna\Social_Defeat\juv_adult_BL6.mat')
figure('DefaultAxesFontSize',12)

pos3 = [0.35 0.88 0.05 0.07];
subplot('Position',pos3)

blue = [.301 .745 .933];

% datapoints
% stressed
sz=50;
swarmchart(zeros(size(juv_stressed,1),1),juv_stressed,sz,'filled','MarkerEdgeColor',blue,'MarkerFaceColor',blue,'MarkerFaceAlpha',0.2)
hold on
swarmchart(ones(size(adult_stressed,1),1),adult_stressed,sz,'filled','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerFaceAlpha',0.2)

% control
swarmchart(2+ones(size(juv_control,1),1),juv_control,sz,'filled','MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'MarkerFaceAlpha',0.2)
swarmchart(3+ones(size(adult_control,1),1),adult_control,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)

% errorbars
% stressed
errorbar(0,juv_stressed_avg,juv_stressed_std,'k');
errorbar(1,adult_stressed_avg,adult_stressed_std,'k');

% control
errorbar(3,juv_control_avg,juv_control_std,'k');
errorbar(4,adult_control_avg,adult_control_std,'k');

xticks(0:4)
xtickangle(45)
xticklabels({'Juvenile (stressed)';'Adult (stressed)';'';'Juvenile (control)';'Adult (control)'})
ylabel('% Time spent in social zone')
title('Post-CSDS') %,'FontSize',20)
xlim([-.5 4.5]);
box off
end

%% EXTRA: Immobility PRE bars
for xc = 1
clear
load('immobility_PRE.mat')
figure('DefaultAxesFontSize',12)
pos3 = [0.35 0.88 0.05 0.1];
subplot('Position',pos3)

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

% bars
bar(0,immobility_res_all.velocity_avg,'facecolor',purple,'EdgeColor',purple,'FaceAlpha',0.7)
hold on
bar(1,immobility_sus_all.velocity_avg,'facecolor',green,'EdgeColor',green,'FaceAlpha',0.7)
bar(2,immobility_ctrl_all.velocity_avg,'facecolor','k','EdgeColor','k','FaceAlpha',0.7)

% datapoints
sz=50;
swarmchart(zeros(size(immobility_res_all.velocity,1),1),immobility_res_all.velocity,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)
hold on
swarmchart(ones(size(immobility_sus_all.velocity,1),1),immobility_sus_all.velocity,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)
swarmchart(1+ones(size(immobility_ctrl_all.velocity,1),1),immobility_ctrl_all.velocity,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)

% errorbars
errorbar(0,immobility_res_all.velocity_avg,immobility_res_all.velocity_std,'k');
errorbar(1,immobility_sus_all.velocity_avg,immobility_sus_all.velocity_std,'k');
errorbar(2,immobility_ctrl_all.velocity_avg,immobility_ctrl_all.velocity_std,'k');

xticks(0:2)
xtickangle(45)
xticklabels({'Resilient';'Susceptible';'Control'})
ylabel('% Time spent immobile')
title('Pre-CSDS') %,'FontSize',20)
xlim([-.7 2.7]);
%ylim([0 50])
box off
end
