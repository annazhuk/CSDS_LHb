% run whole-field analysis on SI data
clear
load('mice.mat')

%% Collect data into a struct
BL6Post = struct;
filename = '555_BL6post_cells_visits_nose_NEW.mat';

for mouse = 1:size(mice,1)
    filename(1:3) = string(mice(mouse));
    load(filename)

    BL6Post.mice(mouse,:) = mice(mouse);
    BL6Post.nose_plot_start{mouse,1} = nose_plot_start;
    BL6Post.nose_plot_stop{mouse,1} = nose_plot_stop;
    BL6Post.nose_plot_start_avg_FP(mouse,:) = nose_plot_start_avg_FP;
    BL6Post.nose_plot_stop_avg_FP(mouse,:) = nose_plot_stop_avg_FP;
    BL6Post.nose_plot_start_avg{mouse,1} = nose_plot_start_avg;
    BL6Post.nose_plot_stop_avg{mouse,1} = nose_plot_stop_avg;
end
save('BL6post_SInose_cells_NEW.mat','BL6Post')

%% Sort into res and sus
clear;
load('mice.mat')
load('SI_percent_ALL.mat')
load('BL6post_SInose_cells_NEW.mat')
load('CSDS_infoID.mat')

SI_cells_nose_imaging_res_sus_NEW(BL6Post,SW_CSDS_Post,mouse_CSDS,SW_CTRL_Post,CSDS_ID,CSDS_SI)

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

BL6Post.nose_plot_all_sorted = [];
for mouse = 1:size(BL6Post.mice,1)
    BL6Post.nose_plot_all_sorted = [BL6Post.nose_plot_all_sorted; BL6Post.nose_plot_start_avg{idx(mouse)}];
end

susBL6Post_im = sus_im;
resBL6Post_im = res_im;

% get total cell # for res and sus
sus_cell_total = 0;
for m = 1:size(susBL6Post_im.mice,1)
    sus_cell_total = sus_cell_total + size(susBL6Post_im.nose_plot_start_avg{m,1},1);
end

% find cell avgs and stds
% sus
susBL6Post_im.nose_plot_start_avg_FP_overall = nanmean(BL6Post.nose_plot_all_sorted(1:sus_cell_total,:),1);
susBL6Post_im.nose_plot_start_std_FP_overall = nanstd(BL6Post.nose_plot_all_sorted(1:sus_cell_total,:),1)/sqrt(size(BL6Post.nose_plot_all_sorted(1:sus_cell_total,:),1));
% res
resBL6Post_im.nose_plot_start_avg_FP_overall = nanmean(BL6Post.nose_plot_all_sorted(sus_cell_total+1:end,:),1);
resBL6Post_im.nose_plot_start_std_FP_overall = nanstd(BL6Post.nose_plot_all_sorted(sus_cell_total+1:end,:),1)/sqrt(size(BL6Post.nose_plot_all_sorted(sus_cell_total+1:end,:),1));

% find avg activity 
susBL6Post_im.nose_plot_start_Favg_allcells = nanmean(BL6Post.nose_plot_all_sorted(1:sus_cell_total,start_F:stop_F),2);
resBL6Post_im.nose_plot_start_Favg_allcells = nanmean(BL6Post.nose_plot_all_sorted(sus_cell_total+1:end,start_F:stop_F),2);

% save
save('imaging_BL6Post_nose_z-scored_NEW.mat','BL6Post','susBL6Post_im','resBL6Post_im')