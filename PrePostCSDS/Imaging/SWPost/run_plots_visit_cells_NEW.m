% run whole-field analysis on SI data
clear
load('mice.mat')

%% Collect data into a struct
SWPost = struct;
filename = '555_SWpost_cells_visits_nose_NEW.mat';

for mouse = 1:size(mice,1)
    filename(1:3) = string(mice(mouse));
    load(filename)

    SWPost.mice(mouse,:) = mice(mouse);
    SWPost.nose_plot_start{mouse,1} = nose_plot_start;
    SWPost.nose_plot_stop{mouse,1} = nose_plot_stop;
    SWPost.nose_plot_start_avg_FP(mouse,:) = nose_plot_start_avg_FP;
    SWPost.nose_plot_stop_avg_FP(mouse,:) = nose_plot_stop_avg_FP;
    SWPost.nose_plot_start_avg{mouse,1} = nose_plot_start_avg;
    SWPost.nose_plot_stop_avg{mouse,1} = nose_plot_stop_avg;
end
save('SWpost_SInose_cells_NEW.mat','SWPost')

%% Sort into res and sus
clear;
load('mice.mat')
load('SI_percent_ALL.mat')
load('SWpost_SInose_cells_NEW.mat')
load('CSDS_infoID.mat')

SI_cells_nose_imaging_res_sus_NEW(SWPost,SW_CSDS_Post,mouse_CSDS,SW_CTRL_Post,CSDS_ID,CSDS_SI)

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

SWPost.nose_plot_all_sorted = [];
for mouse = 3:size(SWPost.mice,1)
    SWPost.nose_plot_all_sorted = [SWPost.nose_plot_all_sorted; SWPost.nose_plot_start_avg{idx(mouse)}];
end

susSWPost_im = sus_im;
resSWPost_im = res_im;

% get total cell # for res and sus
sus_cell_total = 0;
for m = 1:size(susSWPost_im.mice,1)-2 % last 2 mice did not visit
    sus_cell_total = sus_cell_total + size(susSWPost_im.nose_plot_start_avg{m,1},1);
end

% find cell avgs and stds
% sus
susSWPost_im.nose_plot_start_avg_FP_overall = nanmean(SWPost.nose_plot_all_sorted(1:sus_cell_total,:),1);
susSWPost_im.nose_plot_start_std_FP_overall = nanstd(SWPost.nose_plot_all_sorted(1:sus_cell_total,:),1)/sqrt(size(SWPost.nose_plot_all_sorted(1:sus_cell_total,:),1));
% res
resSWPost_im.nose_plot_start_avg_FP_overall = nanmean(SWPost.nose_plot_all_sorted(sus_cell_total+1:end,:),1);
resSWPost_im.nose_plot_start_std_FP_overall = nanstd(SWPost.nose_plot_all_sorted(sus_cell_total+1:end,:),1)/sqrt(size(SWPost.nose_plot_all_sorted(sus_cell_total+1:end,:),1));

% find avg activity 
susSWPost_im.nose_plot_start_Favg_allcells = nanmean(SWPost.nose_plot_all_sorted(1:sus_cell_total,start_F:stop_F),2);
resSWPost_im.nose_plot_start_Favg_allcells = nanmean(SWPost.nose_plot_all_sorted(sus_cell_total+1:end,start_F:stop_F),2);

% save
save('imaging_SWPost_nose_z-scored_NEW.mat','SWPost','susSWPost_im','resSWPost_im')