%% Load data
clear;

%% colors
for xc = 1
green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];
darkgreen = [.1 .4 .05];
darkpurple = [.2 0 .5];
blue = [.301 .745 .933];
gray = [.7 .7 .7];
darkblue = [0 .447 .741];
darkgray = [.15 .15 .15];

% create green-purple colormap
mycolormap = customcolormap(linspace(0,1,11), {'#410149','#762a84','#9b6fac','#c1a5cd','#e7d4e8','#faf6f7','#d7f1d6','#a6db9d','#5aae60','#1c7735','#014419'});
greenmap = customcolormap(linspace(0,1,2), {'#480091','#B66DFF'});
purplemap = customcolormap(linspace(0,1,2), {'#224624','#73B761'});
b = [.09,.08,.3];
w = [1,1,1];
r = [.5,.08,.09];
seismicmap = createcolormap(b,w,r); 
p = [.29,.11,.34];
g = [.04,.3,.15];
PRGn = createcolormap(p,w,g);
end

%% Calculate significance
for xc = 1
load('SI_behavior.mat')
load('homecage_behavior.mat')
load('EPM_susVScontrol_stress.mat')
load('NSF_susVScontrol_stress.mat')
load('immobility_post_FP_im_fos.mat')
load('OFT_susVScontrol_stress.mat')
load('Fig2_FP.mat')
load('Fig2_imaging.mat')


% 1B SI tests
[h,p_BL6res] = ttest(BL6_CSDS_Pre_res_wnan,BL6_CSDS_Post_res);
[h,p_BL6sus] = ttest(BL6_CSDS_Pre_sus_wnan,BL6_CSDS_Post_sus);
[h,p_BL6ctrl] = ttest(BL6_CTRL_Pre_wnan,BL6_CTRL_Post);
[h,p_AKRres] = ttest(AKR_CSDS_Pre_res_wnan,AKR_CSDS_Post_res);
[h,p_AKRsus] = ttest(AKR_CSDS_Pre_sus_wnan,AKR_CSDS_Post_sus);
[h,p_AKRctrl] = ttest(AKR_CTRL_Pre_wnan,AKR_CTRL_Post);

% correct for multiple comparisons
[corrected_p_SI, h]=bonf_holm([p_BL6res,p_BL6sus,p_BL6ctrl,p_AKRres,p_AKRsus,p_AKRctrl],0.05);
% [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh([p_BL6res,p_BL6sus,p_BL6ctrl,p_AKRres,p_AKRsus,p_AKRctrl],0.05);

% 1C homecage t-test
[h,p_HC_SW_C_R] = ttest2(HC_res_SW,HC_CTRL_SW);
[h,p_HC_SW_R_S] = ttest2(HC_res_SW,HC_sus_SW);
[h,p_HC_SW_C_S] = ttest2(HC_CTRL_SW,HC_sus_SW);
[h,p_HC_BL6_C_R] = ttest2(HC_res_BL6,HC_CTRL_BL6);
[h,p_HC_BL6_R_S] = ttest2(HC_res_BL6,HC_sus_BL6);
[h,p_HC_BL6_C_S] = ttest2(HC_CTRL_BL6,HC_sus_BL6);

% 1D: homecage correlation
[r,p_HC_corr] = corrcoef(HC_CSDS_SI,HC_CSDS_SW);

% 1E: EPM t-test
[h,p_EPM] = ttest2(opensum_res,opensum_sus);

% 1F: EPM correlations
[r,p_EPM_corr] = corrcoef([SI_res;SI_sus],[opensum_res;opensum_sus]);

% 1G: NSF t-test
[h,p_NSF] = ttest2(NSF_res,NSF_sus);

% 1H: NSF correlation
[r,p_NSF_corr] = corrcoef([SI_res;SI_sus(1:end-1)],[NSF_res;NSF_sus]);

% 1I: Chamber exploration t-test
[h,p_CE] = ttest2(immobility_res_all,immobility_sus_all);

% 1J: Chamber exploration correlation
[r,p_CE_corr] = corrcoef([SI_immobility_res;SI_immobility_sus],[immobility_res_all;immobility_sus_all]);

% 2F: FP SI
[h,p_SWres_fp] = ttest(SWresDist_all_pre.Favg,SWresDist_all_post.Favg);
[h,p_SWsus_fp] = ttest(SWsusDist_all_pre.Favg,SWsusDist_all_post.Favg);
[h,p_SWctrl_fp] = ttest(SWctrl_all_pre.Favg,SWctrl_all_post.Favg(2:end));
[h,p_BL6res_fp] = ttest(BL6resDist_all_pre.Favg,BL6resDist_all_post.Favg);
[h,p_BL6sus_fp] = ttest(BL6susDist_all_pre.Favg,BL6susDist_all_post.Favg);
[h,p_BL6ctrl_fp] = ttest(BL6ctrl_all_pre.Favg,BL6ctrl_all_post.Favg);
[h,p_AKRres_fp] = ttest(AKRresDist_all_pre.Favg,AKRresDist_all_post.Favg);
[h,p_AKRsus_fp] = ttest(AKRsusDist_all_pre.Favg,AKRsusDist_all_post.Favg);
[h,p_AKRctrl_fp] = ttest(AKRctrl_all_pre.Favg,AKRctrl_all_post.Favg);

% correct for multiple comparisons
[corrected_p_FP, h]=bonf_holm([p_SWres_fp,p_SWsus_fp,p_SWctrl_fp,p_BL6res_fp,p_BL6sus_fp,p_BL6ctrl_fp,p_AKRres_fp,p_AKRsus_fp,p_AKRctrl_fp],0.05);
% [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh([p_SWres_fp,p_SWsus_fp,p_SWctrl_fp,p_BL6res_fp,p_BL6sus_fp,p_BL6ctrl_fp,p_AKRres_fp,p_AKRsus_fp,p_AKRctrl_fp],0.05);

% 2O: imaging SI
[h,p_SWres_im] = ttest2(resSWPre_im.nose_plot_start_Favg_allcells,resSWPost_im.nose_plot_start_Favg_allcells);
[h,p_SWsus_im] = ttest2(susSWPre_im.nose_plot_start_Favg_allcells,susSWPost_im.nose_plot_start_Favg_allcells);
[h,p_BL6res_im] = ttest2(resBL6Pre_im.nose_plot_start_Favg_allcells,resBL6Post_im.nose_plot_start_Favg_allcells);
[h,p_BL6sus_im] = ttest2(susBL6Pre_im.nose_plot_start_Favg_allcells,susBL6Post_im.nose_plot_start_Favg_allcells);
[h,p_AKRres_im] = ttest2(resAKRPre_im.nose_plot_start_Favg_allcells,resAKRPost_im.nose_plot_start_Favg_allcells);
[h,p_AKRsus_im] = ttest2(susAKRPre_im.nose_plot_start_Favg_allcells,susAKRPost_im.nose_plot_start_Favg_allcells);

% correct for multiple comparisons
[corrected_p_imaging, h]=bonf_holm([p_SWres_im,p_SWsus_im,p_BL6res_im,p_BL6sus_im,p_AKRres_im,p_AKRsus_im],0.05);
% [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh([p_SWres_im,p_SWsus_im,p_BL6res_im,p_BL6sus_im,p_AKRres_im,p_AKRsus_im],0.05);

% 5G: ChR/ChRmine SI
[h,p_ChR_SI] = ttest2(ChR_SWpost,YFP_SWpost);

% 5H: ChR/ChRmine EPM
[h,p_EPM_SI] = ttest2(ChR_EPM,YFP_EPM);

% 5I ChR/ChRmine OFT
[h,p_OFT_SI] = ttest2(ChR_OFT,YFP_OFT);
end

%% Fig 1C (Left): SI behavior plot (SW Pre)
for xc = 1
clear;
%load('W:\Anna\Social_Defeat\paper\code_matlab\ALL_fig_data.mat')
load('SI_behavior.mat')
figure('DefaultAxesFontSize',12)
pos3 = [0.35 0.5 0.05 0.1];
subplot('Position',pos3)

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

% datapoints
% SW
sz=50;
swarmchart(zeros(size(SW_CTRL_Pre,1),1),SW_CTRL_Pre,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
hold on
swarmchart(.5+ones(size(SW_CSDS_Pre_sus,1),1),SW_CSDS_Pre_sus,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)
swarmchart(.5+ones(size(SW_CSDS_Pre_res,1),1),SW_CSDS_Pre_res,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)

% errorbars
errorbar(0,SW_CTRL_Pre_avg,nanstd(SW_CTRL_Pre)/sqrt(length(SW_CTRL_Pre)),'k');
errorbar(1.5,SW_CSDS_Pre_avg,nanstd(SW_CSDS_Pre)/sqrt(length(SW_CSDS_Pre)),'k');
%errorbar(0,SW_CSDS_Pre_res_avg,nanstd(SW_CSDS_Pre_res)/sqrt(length(SW_CSDS_Pre_res)),'k');
%errorbar(0,SW_CSDS_Pre_sus_avg,nanstd(SW_CSDS_Pre_sus)/sqrt(length(SW_CSDS_Pre_sus)),'k');

xlim([-.7 2.2])
xticks(0:1.5:2)
xtickangle(45)
xticklabels({'Control','Stressed'})
xlabel('Aggressor')
ylabel('Time spent in social zone')
title('Pre-CSDS') %,'FontSize',20)
ylim([0 110]);
box off
end

%% Fig 1C (Right): SI behavior plot (SW Post)
for xc = 1
clear;
%load('W:\Anna\Social_Defeat\paper\code_matlab\ALL_fig_data.mat')
load('SI_behavior.mat')
figure('DefaultAxesFontSize',12)
pos3 = [0.35 0.5 0.05 0.1];
subplot('Position',pos3)

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

% datapoints
% SW
sz=50;
swarmchart(zeros(size(SW_CTRL_Post,1),1),SW_CTRL_Post,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
hold on
swarmchart(.5+ones(size(SW_CSDS_Post_sus,1),1),SW_CSDS_Post_sus,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)
swarmchart(.5+ones(size(SW_CSDS_Post_res,1),1),SW_CSDS_Post_res,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)

% errorbars
errorbar(0,SW_CTRL_Post_avg,nanstd(SW_CTRL_Post)/sqrt(length(SW_CTRL_Post)),'k');
errorbar(1.5,SW_CSDS_Post_avg,nanstd(SW_CSDS_Post)/sqrt(length(SW_CSDS_Post)),'k');
%errorbar(0,SW_CSDS_Post_res_avg,nanstd(SW_CSDS_Post_res)/sqrt(length(SW_CSDS_Post_res)),'k');
%errorbar(0,SW_CSDS_Post_sus_avg,nanstd(SW_CSDS_Post_sus)/sqrt(length(SW_CSDS_Post_sus)),'k');

% add dashed line for cutoff
plot([-.7,2.2],[med,med],'k--','LineWidth',2)

text(.5,105,'***','FontSize',15)
plot([0,1.5],[99,99],'k')
xlim([-.7 2.2])
xticks(0:1.5:2)
xtickangle(45)
xticklabels({'Control','Stressed'})
xlabel('Aggressor')
ylabel({'% Time spent aggressor';'in Social Interaction test'})
title('Post-CSDS') %,'FontSize',20)
ylim([0 110]);
box off
end

%% Fig 1D (Left): SI test Pre vs Post; BL6
for xc = 1
figure('DefaultAxesFontSize',12)
pos8 = [0.2 0.54 0.05 0.1];
clear
load('SI_behavior.mat')
subplot('Position',pos8)
sz = 20;
% BL6
% sus
missing_val = length(BL6_CSDS_Post_sus) - length(BL6_CSDS_Pre_sus);
BL6_CSDS_Pre_sus_wnan = [BL6_CSDS_Pre_sus;nan(missing_val,1)];
xval_pre = 1.4*ones(size(BL6_CSDS_Pre_sus_wnan,1),1);
xval_post = 1.6*ones(size(BL6_CSDS_Post_sus,1),1);

for i = 1:size(BL6_CSDS_Post_sus,1)
    plot([xval_pre(i),xval_post(i)],[BL6_CSDS_Pre_sus_wnan(i),BL6_CSDS_Post_sus(i)],'k-')
    hold on
end
scatter([xval_pre;xval_post],[BL6_CSDS_Pre_sus_wnan;BL6_CSDS_Post_sus],sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)

% res
missing_val = length(BL6_CSDS_Post_res) - length(BL6_CSDS_Pre_res);
BL6_CSDS_Pre_res_wnan = [BL6_CSDS_Pre_res;nan(missing_val,1)];
xval_pre = 1.8*ones(size(BL6_CSDS_Pre_res_wnan,1),1);
xval_post = 2*ones(size(BL6_CSDS_Post_res,1),1);

for i = 1:size(BL6_CSDS_Post_res,1)
    plot([xval_pre(i),xval_post(i)],[BL6_CSDS_Pre_res_wnan(i),BL6_CSDS_Post_res(i)],'k-')
    hold on
end
scatter([xval_pre;xval_post],[BL6_CSDS_Pre_res_wnan;BL6_CSDS_Post_res],sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)

% ctrl
missing_val = length(BL6_CTRL_Post) - length(BL6_CTRL_Pre);
BL6_CTRL_Pre_wnan = [BL6_CTRL_Pre;nan(missing_val,1)];
xval_pre = ones(size(BL6_CTRL_Pre_wnan,1),1);
xval_post = 1.2*ones(size(BL6_CTRL_Post,1),1);

for i = 1:size(BL6_CTRL_Post,1)
    plot([xval_pre(i),xval_post(i)],[BL6_CTRL_Pre_wnan(i),BL6_CTRL_Post(i)],'k-')
    hold on
end
scatter([xval_pre;xval_post],[BL6_CTRL_Pre_wnan;BL6_CTRL_Post],sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)

xlim([.9 2.1])
ylim([-5 100])
xticks(1:.2:2)
xticklabels({'Pre','Post','Pre','Post','Pre','Post'})
ylabel({'% Time spent in social zone'})
title('Self Strain')
box off

% Add stars for significance 
% x_vals = [1.1,1.5,1.9];
% for pval = 4:6
%     if corrected_p_SI(pval) < 0.001
%         text(x_vals(pval-3),95,'***','FontSize',14)
%     elseif corrected_p_SI(pval) < 0.01 && corrected_p_SI(pval) > 0.001
%         text(x_vals(pval-3),95,'**','FontSize',14)
%     elseif corrected_p_SI(pval) < 0.05 && corrected_p_SI(pval) > 0.01
%         text(x_vals(pval-3),95,'*','FontSize',14)
%     end
% end

end

%% Fig 1D (Right): SI test Pre vs Post; AKR
for xc = 1
figure('DefaultAxesFontSize',12)
pos8 = [0.2 0.54 0.05 0.1];
clear
load('SI_behavior.mat')
subplot('Position',pos8)
sz = 20;
% AKR
% sus
missing_val = length(AKR_CSDS_Post_sus) - length(AKR_CSDS_Pre_sus);
AKR_CSDS_Pre_sus_wnan = [AKR_CSDS_Pre_sus;nan(missing_val,1)];
xval_pre = 1.4*ones(size(AKR_CSDS_Pre_sus_wnan,1),1);
xval_post = 1.6*ones(size(AKR_CSDS_Post_sus,1),1);

for i = 1:size(AKR_CSDS_Post_sus,1)
    plot([xval_pre(i),xval_post(i)],[AKR_CSDS_Pre_sus_wnan(i),AKR_CSDS_Post_sus(i)],'k-')
    hold on
end
scatter([xval_pre;xval_post],[AKR_CSDS_Pre_sus_wnan;AKR_CSDS_Post_sus],sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)

% res
missing_val = length(AKR_CSDS_Post_res) - length(AKR_CSDS_Pre_res);
AKR_CSDS_Pre_res_wnan = [AKR_CSDS_Pre_res;nan(missing_val,1)];
xval_pre = 1.8*ones(size(AKR_CSDS_Pre_res_wnan,1),1);
xval_post = 2*ones(size(AKR_CSDS_Post_res,1),1);

for i = 1:size(AKR_CSDS_Post_res,1)
    plot([xval_pre(i),xval_post(i)],[AKR_CSDS_Pre_res_wnan(i),AKR_CSDS_Post_res(i)],'k-')
    hold on
end
scatter([xval_pre;xval_post],[AKR_CSDS_Pre_res_wnan;AKR_CSDS_Post_res],sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)

% ctrl
missing_val = length(AKR_CTRL_Post) - length(AKR_CTRL_Pre);
AKR_CTRL_Pre_wnan = [AKR_CTRL_Pre;nan(missing_val,1)];
xval_pre = ones(size(AKR_CTRL_Pre_wnan,1),1);
xval_post = 1.2*ones(size(AKR_CTRL_Post,1),1);

for i = 1:size(AKR_CTRL_Post,1)
    plot([xval_pre(i),xval_post(i)],[AKR_CTRL_Pre_wnan(i),AKR_CTRL_Post(i)],'k-')
    hold on
end
scatter([xval_pre;xval_post],[AKR_CTRL_Pre_wnan;AKR_CTRL_Post],sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)

xlim([.9 2.1])
ylim([-5 100])
xticks(1:.2:2)
xticklabels({'Pre','Post','Pre','Post','Pre','Post'})
ylabel({'% Time spent in social zone'})
title('Other Strain')
box off

% Add stars for significance 
x_vals = [1.1,1.5,1.9];
for pval = 7:9
    if corrected_p_SI(pval) < 0.001
        text(x_vals(pval-6),95,'***','FontSize',14)
    elseif corrected_p_SI(pval) < 0.01 && corrected_p_SI(pval) > 0.001
        text(x_vals(pval-6),95,'**','FontSize',14)
    elseif corrected_p_SI(pval) < 0.05 && corrected_p_SI(pval) > 0.01
        text(x_vals(pval-6),95,'*','FontSize',14)
    end
end

end

%% Fig 1E: Homecage plots
for xc = 1
clear
load('homecage_behavior.mat')

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

figure('DefaultAxesFontSize',12)
pos3 = [0.35 0.88 0.07 0.1];
subplot('Position',pos3)

% datapoints
% SW
sz=10;
swarmchart(zeros(size(HC_CTRL_SW,1),1),HC_CTRL_SW,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
hold on
swarmchart(ones(size(HC_sus_SW,1),1),HC_sus_SW,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)
swarmchart(1+ones(size(HC_res_SW,1),1),HC_res_SW,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)

% BL6
swarmchart(4+zeros(size(HC_CTRL_BL6,1),1),HC_CTRL_BL6,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
swarmchart(4+ones(size(HC_sus_BL6,1),1),HC_sus_BL6,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)
swarmchart(4+1+ones(size(HC_res_BL6,1),1),HC_res_BL6,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)

% errorbars
% SW
errorbar(0,HC_CTRL_SW_avg,HC_CTRL_SW_std,'k');
hold on
errorbar(1,HC_sus_SW_avg,HC_sus_SW_std,'k');
errorbar(2,HC_res_SW_avg,HC_res_SW_std,'k');

% BL6
errorbar(4,HC_CTRL_BL6_avg,HC_CTRL_BL6_std,'k');
errorbar(5,HC_sus_BL6_avg,HC_sus_BL6_std,'k');
errorbar(6,HC_res_BL6_avg,HC_res_BL6_std,'k');

% text(1,105,'*','FontSize',15)
% text(5,115,'*','FontSize',15)
% plot([-.4,2.4],[95,95],'k','linewidth',.5)
% plot([3.6,6.6],[105,105],'k','linewidth',.5)
xlim([-.5 6.5])
xticks(0:1:6)
xtickangle(45)
xticks([])
ylabel({'% time spent investigating mouse';'in homecage assay'})
title('Post-CSDS') %,'FontSize',20)
ylim([0 100]);
box off
end

%% Fig 1F: Homecage correlations
for xc = 1
figure('DefaultAxesFontSize',12)
pos5 = [0.87 0.6 0.05 0.1];
subplot('Position',pos5)

clear
load('homecage_behavior.mat')

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

sz = 30;
% fit line to data
A = HC_CSDS_SI;
B = HC_CSDS_SW;
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
%scatter(HC_CSDS_SI,HC_CSDS_SW,sz,'filled','MarkerEdgeColor',green,'MarkerEdgeAlpha',1,'MarkerFaceColor',green,'MarkerFaceAlpha',.2)
scatter(HC_CSDS_SI_res,HC_CSDS_SW_res,sz,'filled','MarkerEdgeColor',purple,'MarkerEdgeAlpha',1,'MarkerFaceColor',purple,'MarkerFaceAlpha',.2)
scatter(HC_CSDS_SI_sus,HC_CSDS_SW_sus,sz,'filled','MarkerEdgeColor',green,'MarkerEdgeAlpha',1,'MarkerFaceColor',green,'MarkerFaceAlpha',.2)
xlabel('S       SI Time (%)     R')
ylabel('% time investigating')
title('HC SW') %,'FontSize',20)
xlim([-5 75])
ylim([-10 100])
box off
end

%% Fig 1G: EPM bars
for xc = 1
clear
load('EPM_susVScontrol_stress.mat')
figure('DefaultAxesFontSize',12)
pos3 = [0.35 0.88 0.05 0.1];
subplot('Position',pos3)

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

% bars
bar(0,opensum_control_avg,'facecolor','k','EdgeColor','k','FaceAlpha',0.7)
hold on
bar(1,opensum_sus_avg,'facecolor',green,'EdgeColor',green,'FaceAlpha',0.7)
bar(2,opensum_res_avg,'facecolor',purple,'EdgeColor',purple,'FaceAlpha',0.7)

% datapoints
sz=30;
swarmchart(zeros(size(opensum_control,1),1),opensum_control,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
swarmchart(ones(size(opensum_sus,1),1),opensum_sus,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)
swarmchart(1+ones(size(opensum_res,1),1),opensum_res,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)

% errorbars
errorbar(0,opensum_control_avg,opensum_control_std,'k');
errorbar(1,opensum_sus_avg,opensum_sus_std,'k');
errorbar(2,opensum_res_avg,opensum_res_std,'k');

text(1.3,40,'*','FontSize',15)
text(.2,40,'***','FontSize',15)
plot([0,.9],[38 38],'k','linewidth',.5)
plot([1,2],[38 38],'k','linewidth',.5)
xticks(0:2)
xtickangle(45)
xticklabels({'Control';'Susceptible';'Resilient'})
ylabel('% Time spent in open arms')
title('EPM') %,'FontSize',20)
xlim([-.7 2.7]);
ylim([0 45])
box off
end

%% Fig 1H: EPM correlation
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
A = [opencenter_res_SI;opencenter_sus_SI];
B = [opensum_res;opensum_sus];
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
%scatter([SI_res;SI_sus],[opencenter_res;opencenter_sus],sz,'filled','MarkerEdgeColor',rgb('DarkCyan'),'MarkerEdgeAlpha',1,'MarkerFaceColor',rgb('DarkCyan'),'MarkerFaceAlpha',.2)
scatter(opencenter_res_SI,opensum_res,sz,'filled','MarkerEdgeColor',purple,'MarkerEdgeAlpha',1,'MarkerFaceColor',purple,'MarkerFaceAlpha',.2)
scatter(opencenter_sus_SI,opensum_sus,sz,'filled','MarkerEdgeColor',green,'MarkerEdgeAlpha',1,'MarkerFaceColor',green,'MarkerFaceAlpha',.2)
xlabel('S       SI Time (%)     R')
ylabel({'% time spent';'in open arms'})
title('EPM') %,'FontSize',20)
xlim([-5 70])
ylim([-5 40])
box off
end

%% Fig 1I NSF bars
for xc = 1
clear
load('NSF_susVcontrol_stress.mat')
figure('DefaultAxesFontSize',12)
pos3 = [0.35 0.88 0.05 0.1];
subplot('Position',pos3)

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

% bars
bar(1,NSF_sus_avg,'facecolor',green,'EdgeColor',green,'FaceAlpha',0.7)
hold on
bar(0,NSF_control_avg,'facecolor','k','EdgeColor','k','FaceAlpha',0.7)
bar(2,NSF_res_avg,'facecolor',purple,'EdgeColor',purple,'FaceAlpha',0.7)

% datapoints
sz=30;
swarmchart(zeros(size(NSF_control,1),1),NSF_control,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
swarmchart(ones(size(NSF_sus,1),1),NSF_sus,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)
swarmchart(1+ones(size(NSF_res,1),1),NSF_res,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)

% errorbars
errorbar(0,NSF_control_avg,NSF_control_std,'k');
errorbar(1,NSF_sus_avg,NSF_sus_std,'k');
errorbar(2,NSF_res_avg,NSF_res_std,'k');

text(.2,685,'**','FontSize',15)
plot([0,1],[675,675],'k','linewidth',.5)
xticks(0:2)
xtickangle(45)
xticklabels({'Control';'Susceptible';'Resilient'})
ylabel('Latency to feed')
title('NSF') %,'FontSize',20)
xlim([-.7 2.7]);
ylim([0 750])
box off
end

%% Fig 1J: NSF correlations
for xc = 1
clear
load('NSF_susVcontrol_stress.mat')
figure('DefaultAxesFontSize',12)
pos5 = [0.87 0.6 0.05 0.1];
subplot('Position',pos5)
% fit line to data
sz = 30;
A = [NSF_res_SI;NSF_sus_SI];
B = [NSF_res;NSF_sus];
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
%scatter([SI_res;SI_sus(1:end-1)],[NSF_res;NSF_sus],sz,'filled','MarkerEdgeColor',rgb('DarkCyan'),'MarkerEdgeAlpha',1,'MarkerFaceColor',rgb('DarkCyan'),'MarkerFaceAlpha',.2)
scatter(NSF_res_SI,NSF_res,sz,'filled','MarkerEdgeColor',purple,'MarkerEdgeAlpha',1,'MarkerFaceColor',purple,'MarkerFaceAlpha',.2)
scatter(NSF_sus_SI,NSF_sus,sz,'filled','MarkerEdgeColor',green,'MarkerEdgeAlpha',1,'MarkerFaceColor',green,'MarkerFaceAlpha',.2)
xlabel('S       SI Time (%)     R')
ylabel({'latency to feed (s)'})
title('NSF') %,'FontSize',20)
xlim([-5 75])
ylim([-70 650])
box off
end

%% Fig 1K: Immobility bars
for xc=1
clear
load('immobility_post_FP_im_fos.mat')
figure('DefaultAxesFontSize',15)
pos2 = [0.7 0.87 0.05 0.1];
subplot('Position',pos2)

% bars
bar(0,immobility_control_all_avg,'facecolor','k','EdgeColor','k','FaceAlpha',0.7)
hold on
bar(1,immobility_sus_all_avg,'facecolor',green,'EdgeColor',green,'FaceAlpha',0.7)
bar(2,immobility_res_all_avg,'facecolor',purple,'EdgeColor',purple,'FaceAlpha',0.7)

% datapoints
sz=30;
swarmchart(zeros(size(immobility_control_all,1),1),immobility_control_all,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
swarmchart(ones(size(immobility_sus_all,1),1),immobility_sus_all,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)
swarmchart(1+ones(size(immobility_res_all,1),1),immobility_res_all,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)

% errorbars
errorbar(0,immobility_control_all_avg,immobility_control_all_std,'k');
errorbar(1,immobility_sus_all_avg,immobility_sus_all_std,'k');
errorbar(2,immobility_res_all_avg,immobility_res_all_std,'k');

text(1.3,57,'*','FontSize',15)
text(0.2,57,'***','FontSize',15)
text(.6,67,'**','FontSize',15)
plot([0,.9],[52,52],'k','linewidth',.5)
plot([1,2],[52,52],'k','linewidth',.5)
plot([0,2],[65,65],'k','linewidth',.5)
xlim([-.7 2.7])
xticks(0:2)
xtickangle(45)
xticklabels({'Control';'Susceptible';'Resilient'})%,'contrl'})
ylabel('% Time spent immobile')
title('Chamber exploration') %,'FontSize',20)
ylim([0 75]);
box off
end

%% Fig 1L: Immobility correlation
for xc = 1
clear
load('immobility_post_FP_im_fos.mat')
figure('DefaultAxesFontSize',15)
pos5 = [0.55 0.6 0.05 0.1];
subplot('Position',pos5)

% fit line to data
A = [SI_immobility_res;SI_immobility_sus];
B = [immobility_res_all;immobility_sus_all];
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
scatter(SI_immobility_sus,immobility_sus_all,sz,'filled','MarkerEdgeColor',green,'MarkerEdgeAlpha',1,'MarkerFaceColor',green,'MarkerFaceAlpha',.2)
scatter(SI_immobility_res,immobility_res_all,sz,'filled','MarkerEdgeColor',purple,'MarkerEdgeAlpha',1,'MarkerFaceColor',purple,'MarkerFaceAlpha',.2)
xlabel('S       SI Time (%)     R')
ylabel({'% time spent';'immobile'})
title('Chamber Exploration') %,'FontSize',20)
xticks(0:40:100)
yticks(0:20:100)
xlim([-5 100])
%ylim([-5 50])
box off
end

%% Fig 1M: OFT bars
for xc = 1
clear
load('OFT_susVScontrol_stress.mat')
figure('DefaultAxesFontSize',12)
pos3 = [0.35 0.88 0.05 0.1];
subplot('Position',pos3)

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

% bars
bar(0,OFT_control_avg,'facecolor','k','EdgeColor','k','FaceAlpha',0.7)
hold on
bar(1,OFT_sus_avg,'facecolor',green,'EdgeColor',green,'FaceAlpha',0.7)
bar(2,OFT_res_avg,'facecolor',purple,'EdgeColor',purple,'FaceAlpha',0.7)

% datapoints
sz=30;
swarmchart(zeros(size(OFT_control,1),1),OFT_control,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
hold on
swarmchart(ones(size(OFT_sus,1),1),OFT_sus,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)
swarmchart(1+ones(size(OFT_res,1),1),OFT_res,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)

% errorbars
errorbar(0,OFT_control_avg,OFT_control_std,'k');
errorbar(1,OFT_sus_avg,OFT_sus_std,'k');
errorbar(2,OFT_res_avg,OFT_res_std,'k');

text(.2,45,'**','FontSize',15)
text(1,55,'*','FontSize',15)
plot([0,1],[43 43],'k','linewidth',.5)
plot([0,2],[53 53],'k','linewidth',.5)
xticks(0:2)
xtickangle(45)
xticklabels({'Control';'Susceptible';'Resilient'})
ylabel('% Time spent in inner zone')
title('OFT') %,'FontSize',20)
xlim([-.7 2.7]);
ylim([0 60])
box off
end

%% Fig 1N: OFT correlation
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
A = [OFT_res_SI;OFT_sus_SI];
B = [OFT_res;OFT_sus];
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
%scatter([SI_res;SI_sus],[OFT_res;OFT_sus],sz,'filled','MarkerEdgeColor',rgb('DarkCyan'),'MarkerEdgeAlpha',1,'MarkerFaceColor',rgb('DarkCyan'),'MarkerFaceAlpha',.2)
scatter(OFT_res_SI,OFT_res,sz,'filled','MarkerEdgeColor',purple,'MarkerEdgeAlpha',1,'MarkerFaceColor',purple,'MarkerFaceAlpha',.2)
scatter(OFT_sus_SI,OFT_sus,sz,'filled','MarkerEdgeColor',green,'MarkerEdgeAlpha',1,'MarkerFaceColor',green,'MarkerFaceAlpha',.2)
xlabel('S       SI Time (%)     R')
ylabel('% Time spent in inner zone')
title('OFT') %,'FontSize',20)
xlim([-5 70])
ylim([-5 40])
box off
end

%% Fig  2B: Example Traces (FP)
for xc = 1
clear
load('example_traces_FP.mat')
figure('DefaultAxesFontSize',12)
pos4 = [0.5 0.88 0.1 0.09];
subplot('Position',pos4)
idx_stop = close_frames(find(diff(close_frames)>1));
idx_start = close_frames(find(diff(close_frames)>1)+1);
startstop(:,1) = [close_frames(1);idx_start];
startstop(:,2) = [idx_stop;close_frames(end)];

visit_lengths = startstop(:,2) - startstop(:,1);
[~,idx_visits] = sort(visit_lengths);

y = 0;
for i = 1:length(visit_lengths)
    start_plot = startstop(idx_visits(i),1)-round(steps_1)*10; %10s before start
    stop_plot = startstop(idx_visits(i),2)+round(steps_1)*10; % 10s before end
    
    y = y + 0.15;    % add y to stack plots 
    
    x_plot_full = 0:length(start_plot:stop_plot-1);
    x_plot_visit = 0:length(startstop(idx_visits(i),1):startstop(idx_visits(i),2)-1);
    
    plot(x_plot_full,F_final(start_plot:stop_plot) - y,'k','LineWidth',2)
    hold on
    plot(x_plot_visit+round(steps_1)*10,F_final(startstop(idx_visits(i),1):startstop(idx_visits(i),2)) - y,'color',[0 .447 .741],'LineWidth',2)
end

% add legend
x1 = 0;
x2 = round(steps_1)*10;
y1 = 0;
y2 = .1;
added_y = -.35;
added_x = 1200;
plot([x1+added_x, x2+added_x],[added_y, added_y],'k-')  % horizontal line: time
plot([x1+added_x, x1+added_x],[y1+added_y,y2+added_y],'k-') % vertical line: time

% add dashed line at start and end of visits
plot([round(steps_1)*10,round(steps_1)*10],[-.18,-.45],'--','LineWidth',2,'color','k')
x11 = visit_lengths(idx_visits(1))+round(steps_1)*10;%,round(steps_1)*10;
x22 = visit_lengths(idx_visits(2))+round(steps_1)*10;%,round(steps_1)*10;
plot([x11,x22],[-.18,-.32],'--','LineWidth',2,'color','k')
x33 = visit_lengths(idx_visits(2))+round(steps_1)*10;%,round(steps_1)*10;
x44 = visit_lengths(idx_visits(3))+round(steps_1)*10;%,round(steps_1)*10;
plot([x33,x44],[-.32,-.45],'--','LineWidth',2,'color','k')

set(gca,'visible','off')
set(gca,'xtick',[])
text(x1+added_x+60,added_y-.025,'10s','FontSize',10,'LineWidth',2)
h = text(x1+added_x-70,y1+added_y*1.1,'.1 \DeltaF/F','FontSize',10,'LineWidth',2);
text(0,-.55,'visit start','FontSize',10)
text(1000,-.55,'visit end','FontSize',10)
set(h,'Rotation',90);
title('example traces')
ax = gca;
ax.Title.Visible = 'on';
box off
end

%% Fig 2C: FP response across all mice
for xc = 1
clear 
load('Fig2_FP.mat')
figure('DefaultAxesFontSize',12)
pos5 = [0.3 0.73 0.11 0.1];
subplot('Position',pos5)
time_plot = -2:2:5;
[val,idx] = sort(cropped_SI_all);
imagesc(SWcsds_all_post.visits_plot_avg(idx,:))
caxis([-2.5 2.5])
a=colorbar;
ylabel(a,'z-scored \DeltaF/F','FontSize',10,'Rotation',90);
%ColorMap = get(gcf,'Colormap');
colormap(gca,PRGn)
hold on
plot([60,60],[0.5,20.5],'k--','linewidth',1) % red line at 0 sec
xticks(0:60:210)   % 30 = 1 sec
xticklabels(time_plot)
yticklabels({})
xlabel('time from social zone entry(s)')
ylabel({'mice';'res \leftrightarrow sus'})
title('Aggressor Strain Post-CSDS')
box off
hold off
end

%% Fig 2D: R vs S traces for SW; Pre-CSDS (FP)
for xc = 1
clear 
load('Fig2_FP.mat')

figure('DefaultAxesFontSize',15)
pos6 = [0.35 0.73 0.1 0.1];

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

subplot('Position',pos6)
time_plot = -2:2:5;
fill([60 60 120 120],[2 -.7 -.7 2],'k','FaceAlpha',0.07,'EdgeColor','none')
hold on
%plot(SWsusDist_all_pre.visits_trace_avg,'color',green,'linewidth',2)
shadedErrorBar(1:211,SWsusDist_all_pre.visits_trace_avg,SWsusDist_all_pre.visits_trace_std,'lineProps',{'color', green})
%plot(SWresDist_all_pre.visits_trace_avg,'color',purple,'linewidth',2)
shadedErrorBar(1:211,SWresDist_all_pre.visits_trace_avg,SWresDist_all_pre.visits_trace_std,'lineProps',{'color', purple})
plot([60,60],[-1 2],'k--','linewidth',1)
xlim([0 210])
ylim([-.7 2])
xticks(0:60:210)
yticks(-0.5:.5:2)
xticklabels(time_plot)
xlabel('time from social zone entry(s)')
ylabel({'LHb GCaMP';'z-scored \DeltaF/F'})
title('Pre-CSDS')
box off
end

%% Fig 2E: R vs S traces for SW; Post-CSDS (FP)
for xc = 1
clear 
load('Fig2_FP.mat')

figure('DefaultAxesFontSize',12)
pos7 = [0.5 0.73 0.1 0.1];

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

subplot('Position',pos7)
time_plot = -2:2:5;
fill([90 90 120 120],[2 -.7 -.7 2],'k','FaceAlpha',0.07,'EdgeColor','none')
hold on
%plot(SWsusDist_all_post.visits_trace_avg,'color',green,'linewidth',2)
shadedErrorBar(1:211,SWsusDist_all_post.visits_trace_avg,SWsusDist_all_post.visits_trace_std,'lineProps',{'color', green})
%plot(SWresDist_all_post.visits_trace_avg,'color',purple,'linewidth',2)
shadedErrorBar(1:211,SWresDist_all_post.visits_trace_avg,SWresDist_all_post.visits_trace_std,'lineProps',{'color', purple})
plot([60,60],[-1 2],'k--','linewidth',1)
xlim([0 210])
ylim([-.7 2])
xticks(0:60:210)
yticks(-0.5:.5:2)
xticklabels(time_plot)
%xlabel('time (s)')
%ylabel({'LHb GCaMP';'z-scored \DeltaF/F'})
title('Post-CSDS')
box off
end

%% Fig 2F: Favg responses for SW (FP)
for xc = 1
clear 
load('Fig2_FP.mat')

figure('DefaultAxesFontSize',15)
pos8 = [0.2 0.54 0.05 0.09];

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

subplot('Position',pos8)
sz = 20;
% SW
% res
xval = [1.8*ones(size(SWresDist_all_pre.Favg,1),1),2*ones(size(SWresDist_all_post.Favg,1),1)];
y = [SWresDist_all_pre.Favg,SWresDist_all_post.Favg];

for i = 1:size(SWresDist_all_post.Favg,1)
    plot([xval(i),xval(i+10)],[y(i),y(i+10)],'k-')
    hold on
end
scatter(xval(1:10),y(1:10),sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)
scatter(xval(11:20),y(11:20),sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)
%text(1.02,.36,'****','FontSize',20)

% sus
xval = [1.4*ones(size(SWsusDist_all_pre.Favg,1),1),1.6*ones(size(SWsusDist_all_post.Favg,1),1)];
y = [SWsusDist_all_pre.Favg,SWsusDist_all_post.Favg];

for i = 1:size(SWsusDist_all_post.Favg,1)
    plot([xval(i),xval(i+10)],[y(i),y(i+10)],'k-')
    hold on
end
scatter(xval(1:10),y(1:10),sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)
scatter(xval(11:20),y(11:20),sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)

% ctrl
xval = [ones(size(SWctrl_all_pre.Favg,1)+1,1),1.2*ones(size(SWctrl_all_post.Favg,1),1)];
y = [[nan;SWctrl_all_pre.Favg],SWctrl_all_post.Favg];

for i = 1:size(SWctrl_all_post.Favg,1)
    plot([xval(i),xval(i+10)],[y(i),y(i+10)],'k-')
    hold on
end
scatter(xval(1:10),y(1:10),sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
scatter(xval(11:20),y(11:20),sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)

xlim([.9 2.1])
ylim([-1 4])
xticks(1:.2:2)
xticklabels({'Pre','Post','Pre','Post','Pre','Post'})
%xlabel('aggressor             self           other              aggressor  self  other')
ylabel({'LHb GCaMP';'mean z-';'scored \DeltaF/F'})
text(1.5,3.8,'*','FontSize',14)
title('Aggressor Strain')
box off
end

%% Fig 2G: Favg responses BL6 (FP)
for xc = 1
clear 
load('Fig2_FP.mat')

figure('DefaultAxesFontSize',15)
pos9 = [0.35 0.54 0.05 0.09];

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

subplot('Position',pos9)
sz = 20;
% res
xval = [1.8*ones(size(BL6resDist_all_pre.Favg,1),1),2*ones(size(BL6resDist_all_post.Favg,1),1)];
y = [BL6resDist_all_pre.Favg,BL6resDist_all_post.Favg];

for i = 1:size(BL6resDist_all_post.Favg,1)
    plot([xval(i),xval(i+10)],[y(i),y(i+10)],'k-')
    hold on
end
scatter(xval(1:10),y(1:10),sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)
scatter(xval(11:20),y(11:20),sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)

% sus
xval = [1.4*ones(size(BL6susDist_all_pre.Favg,1),1),1.6*ones(size(BL6susDist_all_post.Favg,1),1)];
y = [BL6susDist_all_pre.Favg,BL6susDist_all_post.Favg];

for i = 1:size(BL6susDist_all_post.Favg,1)
    plot([xval(i),xval(i+10)],[y(i),y(i+10)],'k-')
    hold on
end
scatter(xval(1:10),y(1:10),sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)
scatter(xval(11:20),y(11:20),sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)

% ctrl
xval = [ones(size(BL6ctrl_all_pre.Favg,1),1),1.2*ones(size(BL6ctrl_all_post.Favg,1),1)];
y = [BL6ctrl_all_pre.Favg,BL6ctrl_all_post.Favg];

for i = 1:size(BL6ctrl_all_post.Favg,1)
    plot([xval(i),xval(i+10)],[y(i),y(i+10)],'k-')
    hold on
end
scatter(xval(1:10),y(1:10),sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
scatter(xval(11:20),y(11:20),sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)

xlim([.9 2.1])
ylim([-1 4])
xticks(1:.2:2.1)
xticklabels({'Pre','Post','Pre','Post','Pre','Post'})
%xlabel('aggressor             self           other              aggressor  self  other')
title('Self Strain')
box off
end

%% Fig 2H: Favg responses for AKR (FP)
for xc = 1
clear 
load('Fig2_FP.mat')

figure('DefaultAxesFontSize',15)
pos10 = [0.5 0.54 0.05 0.09];

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

subplot('Position',pos10)
sz = 20;
% res
xval = [1.8*ones(size(AKRresDist_all_pre.Favg,1),1),2*ones(size(AKRresDist_all_post.Favg,1),1)];
y = [AKRresDist_all_pre.Favg,AKRresDist_all_post.Favg];

for i = 1:size(AKRresDist_all_post.Favg,1)
    plot([xval(i),xval(i+5)],[y(i),y(i+5)],'k-')
    hold on
end
scatter(xval(1:5),y(1:5),sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)
scatter(xval(6:10),y(6:10),sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)

% sus
xval = [1.4*ones(size(AKRsusDist_all_pre.Favg,1),1),1.6*ones(size(AKRsusDist_all_post.Favg,1),1)];
y = [AKRsusDist_all_pre.Favg,AKRsusDist_all_post.Favg];

for i = 1:size(AKRsusDist_all_post.Favg,1)
    plot([xval(i),xval(i+5)],[y(i),y(i+5)],'k-')
    hold on
end
scatter(xval(1:5),y(1:5),sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)
scatter(xval(6:10),y(6:10),sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)

% ctrl
xval = [ones(size(AKRctrl_all_pre.Favg,1),1),1.2*ones(size(AKRctrl_all_post.Favg,1),1)];
y = [AKRctrl_all_pre.Favg,AKRctrl_all_post.Favg];

for i = 1:size(AKRctrl_all_post.Favg,1)
    plot([xval(i),xval(i+7)],[y(i),y(i+7)],'k-')
    hold on
end
scatter(xval(1:7),y(1:7),sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
scatter(xval(8:14),y(8:14),sz,'filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)

xlim([.9 2.1])
ylim([-1 4])
xticks(1:.2:2)
xticklabels({'Pre','Post','Pre','Post','Pre','Post'})
%xlabel('sus             res           control')
title('Other Strain')
box off
end

%% Fig 2I (Left): # Visit vs FP (Sus):
for xc = 1
clear;
load('W:\Anna\Social_Defeat\SW_BOTH_visits_pooled_peaks.mat')

% colors
green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];
gray = [.7 .7 .7];

figure('Renderer', 'painters', 'Position', [10 10 145 141],'DefaultAxesFontSize',12)

sz = 30;
scatter(1:length(SWsusDist_all_post.trial_avg),SWsusDist_all_post.trial_avg,sz,'filled','MarkerFaceColor',green,'MarkerEdgeColor',green,'MarkerFaceAlpha',.2)
hold on
errorbar(1:length(SWsusDist_all_post.trial_avg),SWsusDist_all_post.trial_avg,SWsusDist_all_post.trial_std,'k','LineStyle','none')

% fit line to data
A = 1:length(SWsusDist_all_post.trial_avg);
B = SWsusDist_all_post.trial_avg;
mdl = fitlm(A,B,'linear');
h = plot(mdl);

% Get handles to plot components
dataHandle = findobj(h,'DisplayName','Data');
fitHandle = findobj(h,'DisplayName','Fit');
cbHandles = findobj(h,'DisplayName','Confidence bounds');
cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
dataHandle.Color = 'w'; 
fitHandle.Color = green;
fitHandle.LineWidth = 1;
fitHandle.LineStyle = '--';
set(cbHandles, 'Color', 'k', 'LineWidth', 1,'LineStyle','-')

xlim([0 length(SWresDist_all_pre.trial_avg)+1])
ylim([-1 2])
xlabel('Visit number')
ylabel('Mean \DeltaF/F (Z)')

scatter(1:length(SWsusDist_all_pre.trial_avg),SWsusDist_all_pre.trial_avg,sz,'filled','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',.2)
hold on
errorbar(1:length(SWsusDist_all_pre.trial_avg),SWsusDist_all_pre.trial_avg,SWsusDist_all_pre.trial_std,'k','LineStyle','none')

% fit line to data
A = 1:length(SWsusDist_all_pre.trial_avg);
B = SWsusDist_all_pre.trial_avg;
mdl = fitlm(A,B,'linear');
h = plot(mdl);

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

xlim([0 length(SWsusDist_all_pre.trial_avg)+1])
ylim([-1 2])
xlabel('Visit number')
ylabel('Mean \DeltaF/F (Z)')
end

%% Fig 2I (Right): # Visit vs FP (Res)
for xc = 1
clear;
load('W:\Anna\Social_Defeat\SW_BOTH_visits_pooled_peaks.mat')

% colors
green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];
gray = [.7 .7 .7];

figure('Renderer', 'painters', 'Position', [10 10 145 141],'DefaultAxesFontSize',12)

sz = 30;

scatter(1:length(SWresDist_all_post.trial_avg),SWresDist_all_post.trial_avg,sz,'filled','MarkerFaceColor',purple,'MarkerEdgeColor',purple,'MarkerFaceAlpha',.2)
hold on
e = errorbar(1:length(SWresDist_all_post.trial_avg),SWresDist_all_post.trial_avg,SWresDist_all_post.trial_std,'k','LineStyle','none');

% Set transparency level (0:1)
alpha = 0.2;   
% Set transparency (undocumented)
set([e.Bar, e.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e.Line.ColorData(1:3); 255*alpha])

% fit line to data
A = 1:length(SWresDist_all_post.trial_avg);
B = SWresDist_all_post.trial_avg;
mdl = fitlm(A,B,'linear');
h = plot(mdl);

% Get handles to plot components
dataHandle = findobj(h,'DisplayName','Data');
fitHandle = findobj(h,'DisplayName','Fit');
cbHandles = findobj(h,'DisplayName','Confidence bounds');
cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
dataHandle.Color = 'w'; 
fitHandle.Color = purple;
fitHandle.LineWidth = 1;
fitHandle.LineStyle = '--';
set(cbHandles, 'Color', 'k', 'LineWidth', 1,'LineStyle','-')

xlim([0 length(SWresDist_all_pre.trial_avg)+1])
ylim([-1 2])
xlabel('Visit number')
ylabel('Mean \DeltaF/F (Z)')

scatter(1:length(SWresDist_all_pre.trial_avg),SWresDist_all_pre.trial_avg,sz,'filled','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',.2)
hold on
e = errorbar(1:length(SWresDist_all_pre.trial_avg),SWresDist_all_pre.trial_avg,SWresDist_all_pre.trial_std,'k','LineStyle','none');

% Set transparency level (0:1)
alpha = 0.2;   
% Set transparency (undocumented)
set([e.Bar, e.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e.Line.ColorData(1:3); 255*alpha])

% fit line to data
A = 1:length(SWresDist_all_pre.trial_avg);
B = SWresDist_all_pre.trial_avg;
mdl = fitlm(A,B,'linear');
h = plot(mdl);

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

xlim([0 length(SWresDist_all_post.trial_avg)+1])
ylim([-1 2])
xlabel('Visit number')
ylabel('Mean \DeltaF/F (Z)')
end

%% Fig 2J: FP Avoidance vs magnitude
for xc  =1
clear
load('FP_ALL_pval_05to25.mat')
figure('Renderer', 'painters', 'Position', [10 10 140 150],'DefaultAxesFontSize',12)
sz=30;

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

susSI = cropped_SI_all(ismember(cropped_ID_all,SWsusDist_all_post.mouse));
resSI = cropped_SI_all(ismember(cropped_ID_all,SWresDist_all_post.mouse));

scatter(resSI,SWresDist_all_post.Favg,sz,'filled','MarkerFaceColor',purple,'MarkerEdgeColor',purple,'MarkerFaceAlpha',.2)
hold on
scatter(susSI,SWsusDist_all_post.Favg,sz,'filled','MarkerFaceColor',green,'MarkerEdgeColor',green,'MarkerFaceAlpha',.2)

% fit line to data
A = [susSI;resSI];
B = [SWsusDist_all_post.Favg;SWresDist_all_post.Favg];
mdl = fitlm(A,B,'linear');
h = plot(mdl);
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

xlabel('% time spent near aggressor')
ylabel('Mean \Delta F/F (Z)')
title('SI test (post-CSDS)')
box off
end

%% Fig 2M: Example traces (single cell)
for xc = 1
clear
load('imaging_example_data.mat')
figure('DefaultAxesFontSize',12)
pos22 = [0.47 0.37 0.15 0.09];
subplot('Position',pos22)

fps = 30; %25;
start = 150*fps;
stop = 200*fps;
y_plot = 0;

%cells = [2,3,4,5,6]; % 29
cells = [10,41,17,20,8]; % 29

for cellnum = 1:length(cells)
    plot(start:stop,alldeltaF_social.C_zscore(cells(cellnum),start:stop) + y_plot)
    hold on
    
    y_plot = y_plot + 6;
end

% add legend
x1 = 0;
x2 = 10*fps;
y1 = 0;
y2 = 10;
added_y = 1;
added_x = 140*fps;
plot([x1+added_x, x2+added_x],[added_y, added_y],'k-')  % horizontal line: time
plot([x1+added_x, x1+added_x],[y1+added_y,y2+added_y],'k-') % vertical line: F
text(x1+added_x+50,added_y+0.02,'10s','FontSize',8,'LineWidth',2) % make sure the text = x2/fps
h = text(x1+added_x-fps,y1+added_y,'10 \DeltaF/F','FontSize',8,'LineWidth',2);% make sure the text = y2
set(h,'Rotation',90);

% get rid of axes
xlim([added_x stop])
%ylim([-1 y_plot])
set(gca,'visible','off')
set(gca,'xtick',[])
box off
end

%% Fig 2N: Single-cell response across all mice
for xc = 1
clear
load('Fig2_imaging.mat')

figure('DefaultAxesFontSize',12)
pos11 = [0.17 0.24 0.11 0.1];
subplot('Position',pos11)

time = -2:1/30:5;
time_plot = -2:2:5;

[val,idx] = sort(cropped_SI_im_sorted);

imagesc(SWPost.nose_plot_all_sorted(3:end,:)); caxis([-6 6])
hold on
plot([60,60],[0,size(SWPost.nose_plot_all_sorted,1)],'k--','linewidth',1)
plot([1,211],[size(susSWPost_im.nose_plot_start_Favg_allcells,1),size(susSWPost_im.nose_plot_start_Favg_allcells,1)],'k')

% %separate cells by mouse 
cell_num = [];
cellnum = 0;
% for mouse = 3:size(SWPost.mice,1)-1
%     plot([0,size(SWPost.nose_plot_all_sorted,2)],[size(SWPost.nose_plot_start_avg{idx(mouse)},1)+cellnum,size(SWPost.nose_plot_start_avg{idx(mouse)},1)+cellnum],'k-','linewidth',1.5)
%     cellnum = cellnum + size(SWPost.nose_plot_start_avg{idx(mouse)},1);
% end
xticks(0:60:210)   % 30 = 1 sec
xticklabels(time_plot)
colorbar;
caxis([-6 6])
%colormap(flipud(cool))
colormap(PRGn)
a=colorbar;
ylabel(a,{'Z-Scored';'\DeltaF/F'},'FontSize',10,'Rotation',90);
ylabel({'cell #';'res \leftrightarrow sus'})
xlabel('time from social zone entry(s)')
title('Aggressor Strain Post-CSDS')
box off
end

%% Fig 2O: R vs S traces for SW; Pre-CSDS (Single-cell)
for xc = 1
clear
load('Fig2_imaging.mat')

figure('DefaultAxesFontSize',12)
pos12= [0.35 0.24 0.1 0.1];
subplot('Position',pos12)

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

time_plot = -2:2:5;
fill([60 60 120 120],[1.5 -.3 -.3 1.5],'k','FaceAlpha',0.07,'EdgeColor','none')
hold on
plot(susSWPre_im.nose_plot_start_avg_FP_overall,'color',green,'linewidth',2)
shadedErrorBar(1:211,susSWPre_im.nose_plot_start_avg_FP_overall,susSWPre_im.nose_plot_start_std_FP_overall,'lineProps',{'color', green})
plot(resSWPre_im.nose_plot_start_avg_FP_overall,'color',purple,'linewidth',2)
shadedErrorBar(1:211,resSWPre_im.nose_plot_start_avg_FP_overall,resSWPre_im.nose_plot_start_std_FP_overall,'lineProps',{'color', purple})
plot([60,60],[-.3 1.5],'k--','linewidth',1)
xlim([0 210])
ylim([-.3 1.5])
xticks(0:60:210)
yticks(-0.5:.5:1.5)
xticklabels(time_plot)
xlabel('time from social zone entry(s)')
ylabel({'LHb GCaMP';'normalized F'})
title('Pre-CSDS')
box off
end

%% Fig 2P: R vs S traces for SW; Post-CSDS (Single-cell)
for xc = 1
clear
load('Fig2_imaging.mat')

figure('DefaultAxesFontSize',12)
pos13 = [0.5 0.24 0.1 0.1];
subplot('Position',pos13)

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

time_plot = -2:2:5;
fill([60 60 120 120],[1.5 -.3 -.3 1.5],'k','FaceAlpha',0.07,'EdgeColor','none')
hold on
plot(susSWPost_im.nose_plot_start_avg_FP_overall,'color',green,'linewidth',2)
shadedErrorBar(1:211,susSWPost_im.nose_plot_start_avg_FP_overall,susSWPost_im.nose_plot_start_std_FP_overall,'lineProps',{'color', green})
plot(resSWPost_im.nose_plot_start_avg_FP_overall,'color',purple,'linewidth',2)
shadedErrorBar(1:211,resSWPost_im.nose_plot_start_avg_FP_overall,resSWPost_im.nose_plot_start_std_FP_overall,'lineProps',{'color', purple})
plot([60,60],[-.3 1.5],'k--','linewidth',1)
xlim([0 210])
ylim([-.3 1.5])
xticks(0:60:210)
yticks(-0.5:.5:1.5)
xticklabels(time_plot)
%xlabel('time (s)')
%ylabel({'LHb GCaMP';'normalized \DeltaF/F'})
title('Post-CSDS')
box off
end

%% Fig 2Q: Favg responses for SW (Single-cell)
for xc = 1
clear
load('Fig2_imaging.mat')

figure('DefaultAxesFontSize',12)
pos14 = [0.2 0.05 0.05 0.09];
subplot('Position',pos14)

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

sz = 10;
% sus
swarmchart(ones(size(susSWPre_im.nose_plot_start_Favg_allcells,1),1),susSWPre_im.nose_plot_start_Favg_allcells,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',.2)
hold on
swarmchart(1+ones(size(susSWPost_im.nose_plot_start_Favg_allcells,1),1),susSWPost_im.nose_plot_start_Favg_allcells,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',.2)

% res
swarmchart(3+ones(size(resSWPre_im.nose_plot_start_Favg_allcells,1),1),resSWPre_im.nose_plot_start_Favg_allcells,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',.2)
hold on
swarmchart(4+ones(size(resSWPost_im.nose_plot_start_Favg_allcells,1),1),resSWPost_im.nose_plot_start_Favg_allcells,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',.2)

xlim([0 6])
ylim([-5 7.5])
xtickangle(45)
xticks(1:1:6)
xticklabels({'Pre','Post','','Pre','Post'})
%xlabel('aggressor             self           other              aggressor  self  other')
ylabel({'LHb GCaMP';'Normalized F'})
text(1,5,'***','FontSize',10)
title('Aggressor Strain')
box off
end

%% Fig 2R: Favg responses for BL6 (Single-cell)
for xc = 1
clear
load('Fig2_imaging.mat')

figure('DefaultAxesFontSize',12)
pos15 = [0.35 0.05 0.05 0.09];
subplot('Position',pos15)

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

sz = 10;
% sus
swarmchart(ones(size(susBL6Pre_im.nose_plot_start_Favg_allcells,1),1),susBL6Pre_im.nose_plot_start_Favg_allcells,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',.2)
hold on
swarmchart(1+ones(size(susBL6Post_im.nose_plot_start_Favg_allcells,1),1),susBL6Post_im.nose_plot_start_Favg_allcells,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',.2)

% res
swarmchart(3+ones(size(resBL6Pre_im.nose_plot_start_Favg_allcells,1),1),resBL6Pre_im.nose_plot_start_Favg_allcells,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',.2)
swarmchart(4++ones(size(resBL6Post_im.nose_plot_start_Favg_allcells,1),1),resBL6Post_im.nose_plot_start_Favg_allcells,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',.2)

xlim([0 6])
ylim([-5 7.5])
xtickangle(45)
xticks(1:1:6)
xticklabels({'Pre','Post','','Pre','Post'})
%ylabel({'LHb GCaMP';'Normalized F'})
title('Self Strain')
box off
end

%% Fig 2S: Favg responses for AKR (Single-cell)
for xc = 1
clear
load('Fig2_imaging.mat')

figure('DefaultAxesFontSize',12)
pos16 = [0.5 0.05 0.05 0.09];
subplot('Position',pos16)

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

sz = 10;
% sus
swarmchart(ones(size(susAKRPre_im.nose_plot_start_Favg_allcells,1),1),susAKRPre_im.nose_plot_start_Favg_allcells,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',.2)
hold on
swarmchart(1+ones(size(susAKRPost_im.nose_plot_start_Favg_allcells,1),1),susAKRPost_im.nose_plot_start_Favg_allcells,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',.2)

% res
swarmchart(3+ones(size(resAKRPre_im.nose_plot_start_Favg_allcells,1),1),resAKRPre_im.nose_plot_start_Favg_allcells,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',.2)
swarmchart(4+ones(size(resAKRPost_im.nose_plot_start_Favg_allcells,1),1),resAKRPost_im.nose_plot_start_Favg_allcells,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',.2)

xlim([0 6])
ylim([-5 7.5])
xtickangle(45)
xticks(1:1:6)
xticklabels({'Pre','Post','','Pre','Post'})
%ylabel({'LHb GCaMP';'Normalized F'})
title('Other Strain')
box off
end
 
%% Fig 2T: CDF of Pre-CSDS SW response
for xc = 1
clear
load('Fig2_imaging.mat')

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

bins_sus = min(susSWPre_im.nose_plot_start_Favg_allcells):.001:max(susSWPre_im.nose_plot_start_Favg_allcells);
cdf_sus = cumsum(histcounts(susSWPre_im.nose_plot_start_Favg_allcells,length(bins_sus)))/length(susSWPre_im.nose_plot_start_Favg_allcells)*100;
bins_res = min(resSWPre_im.nose_plot_start_Favg_allcells):.001:max(resSWPre_im.nose_plot_start_Favg_allcells);
cdf_res = cumsum(histcounts(resSWPre_im.nose_plot_start_Favg_allcells,length(bins_res)))/length(resSWPre_im.nose_plot_start_Favg_allcells)*100;

pos3 = [0.35 0.88 0.05 0.08];
subplot('Position',pos3)

plot(bins_sus,cdf_sus,'color',green,'linewidth',2)
hold on
plot(bins_res,cdf_res,'color',purple,'linewidth',2)
xlim([-4 6])
box off
xlabel('Mean Z-scored \DeltaF/F')
ylabel('% of neurons')
title('Pre-CSDS')
end

%% Fig 2U: CDF of Post-CSDS SW response
for xc = 1
clear
load('Fig2_imaging.mat')

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

bins_sus = min(susSWPost_im.nose_plot_start_Favg_allcells):.001:max(susSWPost_im.nose_plot_start_Favg_allcells);
cdf_sus = cumsum(histcounts(susSWPost_im.nose_plot_start_Favg_allcells,length(bins_sus)))/length(susSWPost_im.nose_plot_start_Favg_allcells)*100;
bins_res = min(resSWPost_im.nose_plot_start_Favg_allcells):.001:max(resSWPost_im.nose_plot_start_Favg_allcells);
cdf_res = cumsum(histcounts(resSWPost_im.nose_plot_start_Favg_allcells,length(bins_res)))/length(resSWPost_im.nose_plot_start_Favg_allcells)*100;

pos3 = [0.35 0.88 0.05 0.08];
subplot('Position',pos3)

plot(bins_sus,cdf_sus,'color',green,'linewidth',2)
hold on
plot(bins_res,cdf_res,'color',purple,'linewidth',2)
xlim([-4 6])
box off
xlabel('Mean Z-scored \DeltaF/F')
ylabel('% of neurons')
title('Post-CSDS')
end

%% Fig 2V: Significant cell percentages
for xc = 1
clear
load('nose_sig_cells_res_sus_FIXED.mat')

figure('Renderer', 'painters', 'Position', [10 10 130 130],'DefaultAxesFontSize',12)

labels = {'activated','inhibited','non-responsive'};

subplot(1,2,1)
piechart([res_im.sig_up_percent_allcells,res_im.sig_down_percent_allcells,100-(res_im.sig_up_percent_allcells+res_im.sig_down_percent_allcells)],labels)
title('Resilient')
colororder sail

subplot(1,2,2)
piechart([sus_im.sig_up_percent_allcells,sus_im.sig_down_percent_allcells,100-(sus_im.sig_up_percent_allcells+sus_im.sig_down_percent_allcells)],labels)
title('Susceptible')
colororder sail
end

%% Fig 2W: Avoidance vs magnitude: Imaging (activated cells)
for xc = 1
clear
load('nose_sig_cells_res_sus_FIXED.mat')

% colors
green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

s_SI = cropped_SI_im_sorted(ismember(cropped_ID_im_sorted,sus_im.mice));
r_SI = cropped_SI_im_sorted(ismember(cropped_ID_im_sorted,res_im.mice));

s_up_percent = sus_im.sig_up_percent(1:4);
r_up_percent = res_im.sig_up_percent;
r_up_percent(3:4) = [];
r_SI(3:4) = [];

figure('Renderer', 'painters', 'Position', [10 10 145 141],'DefaultAxesFontSize',12)
sz = 30;

% res
scatter(r_SI,r_up_percent,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',.2)
hold on

% sus
scatter(s_SI(1:4),s_up_percent,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',.2)

A = [s_SI(1:4);r_SI];
B = [s_up_percent;r_up_percent];
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
xlabel('S       SI Time (%)     R')
ylabel('Activated cells (%)')
ylim([-50 80])
title('Activated cells')
box off
end

%% Fig 2X: Avoidance vs magnitude: Imaging (inhibited cells)
for xc = 1
clear
load('nose_sig_cells_res_sus_FIXED.mat')

% colors
green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

s_SI = cropped_SI_im_sorted(ismember(cropped_ID_im_sorted,sus_im.mice));
r_SI = cropped_SI_im_sorted(ismember(cropped_ID_im_sorted,res_im.mice));

s_down_percent = sus_im.sig_down_percent(1:4);
r_down_percent = res_im.sig_down_percent;
r_down_percent(3:4) = [];
r_SI(3:4) = [];

figure('Renderer', 'painters', 'Position', [10 10 145 141],'DefaultAxesFontSize',12)
sz = 30;
% res
scatter(r_SI,r_down_percent,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',.2)
hold on

% sus
scatter(s_SI(1:4),s_down_percent,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',.2)

A = [s_SI(1:4);r_SI];
B = [s_down_percent;r_down_percent];
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
xlabel('S       SI Time (%)     R')
ylabel('Inhibited cells (%)')
title('SW SI Post-CSDS (single-cell)') %,'FontSize',20)
ylim([-50 80])
box off
end

%% Fig 2Y: Peaks Pre
for xc = 1
clear
load('res_sus_baseline_peaks_BOTH_cells.mat')
figure('Renderer', 'painters', 'Position', [10 10 134 128],'DefaultAxesFontSize',12)

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

total_time = 5*60; % 5min total recording
sz = 10;

% peak totals
swarmchart(zeros(size(Pre_sus_im.pk_totals_all_sorted,1),1),Pre_sus_im.pk_totals_all_sorted/total_time,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)
hold on
swarmchart(ones(size(Pre_res_im.pk_totals_all_sorted,1),1),Pre_res_im.pk_totals_all_sorted/total_time,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)

% errorbars
errorbar(0,Pre_sus_im.pk_totals_avg/total_time,Pre_sus_im.pk_totals_std/total_time,'k');
errorbar(1,Pre_res_im.pk_totals_avg/total_time,Pre_res_im.pk_totals_std/total_time,'k');

ylim([0 5])
title('Pre-CSDS baseline transient frequency')
ylabel('Peaks/s')
xticks(0:1)
xticklabels({'Sus';'Res'})
end

%% Fig 2Z: Peaks Post
for xc = 1
clear
load('res_sus_baseline_peaks_BOTH_cells.mat')
figure('Renderer', 'painters', 'Position', [10 10 134 128],'DefaultAxesFontSize',12)

green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

total_time = 5*60; % 5min total recording
sz = 10;

% peak totals
swarmchart(zeros(size(Post_sus_im.pk_totals_all_sorted,1),1),Post_sus_im.pk_totals_all_sorted/total_time,sz,'filled','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerFaceAlpha',0.2)
hold on
swarmchart(ones(size(Post_res_im.pk_totals_all_sorted,1),1),Post_res_im.pk_totals_all_sorted/total_time,sz,'filled','MarkerEdgeColor',purple,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.2)

% errorbars
errorbar(0,Post_sus_im.pk_totals_avg/total_time,Post_sus_im.pk_totals_std/total_time,'k');
errorbar(1,Post_res_im.pk_totals_avg/total_time,Post_res_im.pk_totals_std/total_time,'k');

ylim([0 5])
title('Post-CSDS baseline transient frequency')
text(.4,4.5,'**','FontSize',20)
ylabel('Peaks/s')
xticks(0:1)
xticklabels({'Sus';'Res'})
end

%% Fig 5C: Laser aligneed to attack
for xc = 1
clear
load('fig_5C_data.mat')
load('231002_anna_rf_annotations.mat')
start_stim = time_TD(stim_TD(1)-1);
start_RF = time_RF(stim_RF(1)-1)+camera_delay;

stop_stim = time_TD(stim_TD(end)+1);
stop_RF = time_RF(stim_RF(end)+1)+camera_delay;

start_time_idx = find(time_TD == time_TD(stim_TD(1)-1));
start_RF_idx = find(time_RF == time_RF(stim_RF(1)-1));

stop_time_idx = find(time_TD == time_TD(stim_TD(end)+1));
stop_RF_idx = find(time_RF == time_RF(stim_RF(end)+1));

figure('DefaultAxesFontSize',15);
% plot laser pulses
plot(time_TD(1:end-1),ard/5+1.2,'b')
hold on
% plot RF offline labels
for i = start_RF_idx:stop_RF_idx
    plot([time_RF(i)+camera_delay,time_RF(i)+camera_delay],[0,activityOI(i)],'r')
    hold on
    disp('running...')
end
ylim([-.01 2.02])
xlim([300 600])
xticks({})
yticks({})
axis off
box off
end

%% Fig 5F: Laser % of session
for xc = 1
clear
load('stim_info_both.mat')
figure('DefaultAxesFontSize',12)
pos16 = [0.5 0.05 0.1 0.1];
subplot('Position',pos16)

FR = 1000; % scan rate is 1000Hz for TD
fps = 100; % camera FR is 100Hz
edges = 0:4:56; %0:.5:10;
x_vals = 0:4:52; % 0:.5:9.5;

[N, edges] = histcounts(total_stim_all_trains,edges);

bar(x_vals,N./sum(N),'Facecolor',[.7 .7 .7])
hold on

% pdf plot
x_values = [-8,-4,x_vals,58];
pd = ksdensity(total_stim_all_trains,x_values);

scaled = 5.5; %.6 % scale kernel density to overlay histogram 

plot(x_values,pd*scaled,'LineWidth',2,'color','k')

xticks(0:16:56)
%xticklabels({'0';'50';'100';'150';'200'})
xlim([-8 58])
%ylim([0 .15])
xlabel('total time (%)')
ylabel('density')
box off
end

%% Fig 5G: ChR2 SI
for xc = 1
clear
load('all_post_tests_bothcohorts_opto.mat')
figure('DefaultAxesFontSize',15)
pos16 = [0.5 0.05 0.1 0.1];
subplot('Position',pos16)
% avgs
bar(1.2,YFP_SWpost_avg,'facecolor',rgb('Gray'));
hold on
bar(2.2,ChR_SWpost_avg,'facecolor',rgb('DodgerBlue'));

% errorbars
errorbar(1.2,YFP_SWpost_avg,YFP_SWpost_std,'k');
errorbar(2.2,ChR_SWpost_avg,ChR_SWpost_std,'k');

% datapoints
sz = 30;
swarmchart(1.2*ones(size(YFP_SWpost)),YFP_SWpost,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor',rgb('DarkGray'))
swarmchart(2.2*ones(size(ChR_SWpost)),ChR_SWpost,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor',rgb('Blue'))

xticks(1.2:2.2)
xticklabels({'YFP','ChR/ChRmine'})
xtickangle(45)
ylabel({'% time spent';'with aggressor'})
title('Social Interaction Test')
box off
end

%% Fig 5H: ChR2 EPM
for xc = 1
clear
load('all_post_tests_bothcohorts_opto.mat')
figure('DefaultAxesFontSize',12)
pos16 = [0.5 0.05 0.1 0.1];
subplot('Position',pos16)
% avgs
bar(1.2,YFP_EPM_avg,'facecolor',rgb('Gray'));
hold on
bar(2.2,ChR_EPM_avg,'facecolor',rgb('DodgerBlue'));

% datapoints
sz = 30;
swarmchart(1.2*ones(size(YFP_EPM)),YFP_EPM,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor',rgb('DarkGray'))
swarmchart(2.2*ones(size(ChR_EPM)),ChR_EPM,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor',rgb('Blue'))

% errorbars
errorbar(1.2,YFP_EPM_avg,YFP_EPM_std,'k');
errorbar(2.2,ChR_EPM_avg,ChR_EPM_std,'k');

xticks(1.2:2.2)
xticklabels({'YFP','ChR/ChRmine'})
xtickangle(45)
ylabel({'% time spent';'in open arms'})
title('Elevated Plus Maze')
box off
end

%% Fig 5I: ChR2 OFT
for xc = 1
clear
load('all_post_tests_bothcohorts_opto.mat') 
figure('DefaultAxesFontSize',12)
pos16 = [0.5 0.05 0.1 0.1];
subplot('Position',pos16)
% avgs
bar(1.2,YFP_OFT_avg,'facecolor',rgb('Gray'));
hold on
bar(2.2,ChR_OFT_avg,'facecolor',rgb('DodgerBlue'));

% datapoints
sz = 30;
swarmchart(1.2*ones(size(YFP_OFT)),YFP_OFT,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor',rgb('DarkGray'))
swarmchart(2.2*ones(size(ChR_OFT)),ChR_OFT,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor',rgb('Blue'))

% errorbars
errorbar(1.2,YFP_OFT_avg,YFP_OFT_std,'k');
errorbar(2.2,ChR_OFT_avg,ChR_OFT_std,'k');

xticks(1.2:2.2)
xticklabels({'YFP','ChR/ChRmine'})
xtickangle(45)
ylabel({'% time spent';'in center'})
title('Open Field Test')
box off
end

%% Fig 5J: ChR2 NSF
for xc = 1
clear
load('all_post_tests_bothcohorts_opto.mat')

figure('Renderer', 'painters', 'Position', [10 10 140 140],'DefaultAxesFontSize',12)
% avgs
bar(1.2,YFP_NSF_avg,'facecolor',rgb('Gray'));
hold on
bar(2.2,ChR_NSF_avg,'facecolor',rgb('DodgerBlue'));

% datapoints
sz = 30;
swarmchart(1.2*ones(size(YFP_NSF)),YFP_NSF,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor',rgb('DarkGray'))
swarmchart(2.2*ones(size(ChR_NSF)),ChR_NSF,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor',rgb('Blue'))

% errorbars
errorbar(1.2,YFP_NSF_avg,YFP_NSF_std,'k');
errorbar(2.2,ChR_NSF_avg,ChR_NSF_std,'k');

xticks(1.2:2.2)
xticklabels({'YFP','ChR/ChRmine'})
xtickangle(45)
ylabel({'latency to feed (s)'})
%title({'Novelty';'Suppressed Feeding'})
box off
end

%% Fig 5K: NpHr SI
for xc = 1
clear;
load('NpHr.mat')

figure('Renderer', 'painters', 'Position', [10 10 140 110],'DefaultAxesFontSize',12)
% avgs
bar(1.2,YFP_SWpost_avg,'facecolor',rgb('Gray'));
hold on
bar(2.2,NpHr_SWpost_avg,'facecolor',rgb('Gold'));

% errorbars
errorbar(1.2,YFP_SWpost_avg,YFP_SWpost_std,'k');
errorbar(2.2,NpHr_SWpost_avg,NpHr_SWpost_std,'k');

% datapoints
sz = 30;
swarmchart(1.2*ones(size(YFP_SWpost)),YFP_SWpost,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor',rgb('DarkGray'))
swarmchart(2.2*ones(size(NpHr_SWpost)),NpHr_SWpost,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor',rgb('DarkKhaki'))

xticks(1.2:2.2)
xticklabels({'YFP','NpHr'})
xtickangle(45)
xlim([.5 3])
ylim([0 75])
ylabel({'% time spent';'with aggressor'})
title('Social Interaction Test')
box off
end

%% Fig 5L: NpHr: EPM
for xc = 1
clear;
load('NpHr.mat')

figure('Renderer', 'painters', 'Position', [10 10 140 110],'DefaultAxesFontSize',12)
% avgs
bar(1.2,YFP_EPM_avg,'facecolor',rgb('Gray'));
hold on
bar(2.2,NpHr_EPM_avg,'facecolor',rgb('Gold'));

% errorbars
errorbar(1.2,YFP_EPM_avg,YFP_EPM_std,'k');
errorbar(2.2,NpHr_EPM_avg,NpHr_EPM_std,'k');

% datapoints
sz = 30;
swarmchart(1.2*ones(size(YFP_EPM)),YFP_EPM,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor',rgb('DarkGray'))
swarmchart(2.2*ones(size(NpHr_EPM)),NpHr_EPM,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor',rgb('DarkKhaki'))

% 
plot([1.2,2.2],[90 90], 'k-')

xticks(1.2:2.2)
xticklabels({'YFP','NpHr'})
xtickangle(45)
ylim([0 100])
xlim([.5 3])
text(1.6, 92, '**','FontSize',25)
ylabel({'% Time spent';'in open arms'})
title('Elevated Plus Maze')
box off
end

%% Fig 5M: NpHr: OFT
for xc = 1
clear;
load('NpHr.mat')

figure('Renderer', 'painters', 'Position', [10 10 140 110],'DefaultAxesFontSize',12)
% avgs
bar(1.2,YFP_OFT_avg,'facecolor',rgb('Gray'));
hold on
bar(2.2,NpHr_OFT_avg,'facecolor',rgb('Gold'));

% datapoints
sz = 30;
swarmchart(1.2*ones(size(YFP_OFT)),YFP_OFT,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor',rgb('DarkGray'))
swarmchart(2.2*ones(size(NpHr_OFT)),NpHr_OFT,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor',rgb('DarkKhaki'))

% errorbars
errorbar(1.2,YFP_OFT_avg,YFP_OFT_std,'k');
errorbar(2.2,NpHr_OFT_avg,NpHr_OFT_std,'k');

xticks(1.2:2.2)
xticklabels({'YFP','NpHr'})
xtickangle(45)
xlim([.5 3])
ylim([0 20])
ylabel({'% time spent';'in center'})
title('Open Field Test')
box off
end

%% Fig 5N: NpHr: NSF
for xc = 1
clear;
load('NpHr.mat')

figure('Renderer', 'painters', 'Position', [10 10 140 110],'DefaultAxesFontSize',12)
% avgs
bar(1.2,YFP_NSF_avg,'facecolor',rgb('Gray'));
hold on
bar(2.2,NpHr_NSF_avg,'facecolor',rgb('Gold'));

% errorbars
errorbar(1.2,YFP_NSF_avg,YFP_NSF_std,'k');
errorbar(2.2,NpHr_NSF_avg,NpHr_NSF_std,'k');

% datapoints
sz = 30;
swarmchart(1.2*ones(size(YFP_NSF)),YFP_NSF,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor',rgb('DarkGray'))
swarmchart(2.2*ones(size(NpHr_NSF)),NpHr_NSF,sz,'filled','MarkerEdgeColor','k','MarkerFaceColor',rgb('DarkKhaki'))

xticks(1.2:2.2)
xticklabels({'YFP','NpHr'})
xtickangle(45)
xlim([.5 3])
ylim([0 600])
ylabel({'latency to feed (s)'})
title({'Novelty Suppressed Feeding'})
box off
end
