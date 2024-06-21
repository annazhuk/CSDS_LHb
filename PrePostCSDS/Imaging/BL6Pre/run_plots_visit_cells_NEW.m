% run whole-field analysis on SI data
clear
load('mice.mat')

%% Collect data into a struct
BL6Pre = struct;
filename = '555_BL6pre_cells_visits_nose_NEW.mat';

for mouse = 1:size(mice,1)
    filename(1:3) = string(mice(mouse));
    load(filename)

    BL6Pre.mice(mouse,:) = mice(mouse);
    BL6Pre.nose_plot_start{mouse,1} = nose_plot_start;
    BL6Pre.nose_plot_stop{mouse,1} = nose_plot_stop;
    BL6Pre.nose_plot_start_avg_FP(mouse,:) = nose_plot_start_avg_FP;
    BL6Pre.nose_plot_stop_avg_FP(mouse,:) = nose_plot_stop_avg_FP;
    BL6Pre.nose_plot_start_avg{mouse,1} = nose_plot_start_avg;
    BL6Pre.nose_plot_stop_avg{mouse,1} = nose_plot_stop_avg;
end
save('BL6pre_SInose_cells_NEW.mat','BL6Pre')

%% Sort into res and sus
clear;
load('mice.mat')
load('SI_percent_ALL.mat')
load('BL6pre_SInose_cells_NEW.mat')
load('CSDS_infoID.mat')

SI_cells_nose_imaging_res_sus_NEW(BL6Pre,SW_CSDS_Post,mouse_CSDS,SW_CTRL_Post,CSDS_ID,CSDS_SI)

%% Get cell avgs etc
clear;
load('SI_cells_imaging_visits_nose_res_sus_NEW.mat')
load('mice.mat')

fps = 30;
time = -2:1/fps:5; 
time_plot = -2:1:5;
[val,idx] = sort(cropped_SI_im_sorted);

start_F = 2*fps; % start of proximity is at 3s 
stop_F = 5*fps; % end of proximity is at 4s

BL6Pre.nose_plot_all_sorted = [];
for mouse = 1:size(BL6Pre.mice,1)
    BL6Pre.nose_plot_all_sorted = [BL6Pre.nose_plot_all_sorted; BL6Pre.nose_plot_start_avg{idx(mouse)}];
end

susBL6Pre_im = sus_im;
resBL6Pre_im = res_im;

% get total cell # for res and sus
sus_cell_total = 0;
for m = 1:size(susBL6Pre_im.mice,1)
    sus_cell_total = sus_cell_total + size(susBL6Pre_im.nose_plot_start_avg{m,1},1);
end

% find cell avgs and stds
% sus
susBL6Pre_im.nose_plot_start_avg_FP_overall = nanmean(BL6Pre.nose_plot_all_sorted(1:sus_cell_total,:),1);
susBL6Pre_im.nose_plot_start_std_FP_overall = nanstd(BL6Pre.nose_plot_all_sorted(1:sus_cell_total,:),1)/sqrt(size(BL6Pre.nose_plot_all_sorted(1:sus_cell_total,:),1));
% res
resBL6Pre_im.nose_plot_start_avg_FP_overall = nanmean(BL6Pre.nose_plot_all_sorted(sus_cell_total+1:end,:),1);
resBL6Pre_im.nose_plot_start_std_FP_overall = nanstd(BL6Pre.nose_plot_all_sorted(sus_cell_total+1:end,:),1)/sqrt(size(BL6Pre.nose_plot_all_sorted(sus_cell_total+1:end,:),1));

% find avg activity 
susBL6Pre_im.nose_plot_start_Favg_allcells = nanmean(BL6Pre.nose_plot_all_sorted(1:sus_cell_total,start_F:stop_F),2);
resBL6Pre_im.nose_plot_start_Favg_allcells = nanmean(BL6Pre.nose_plot_all_sorted(sus_cell_total+1:end,start_F:stop_F),2);

% save
save('imaging_BL6Pre_nose_z-scored_NEW.mat','BL6Pre','susBL6Pre_im','resBL6Pre_im')