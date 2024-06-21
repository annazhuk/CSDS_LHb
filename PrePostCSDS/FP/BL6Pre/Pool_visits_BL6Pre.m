clear;
%% Pool visits data (BL6 Pre)
load('visits_plot_NEW_cohort1.mat');
load('SI_F_visits_res_sus_NEW_cohort1.mat');
load('visits_plot_NEW_cohort2.mat');
load('SI_F_visits_res_sus_NEW_cohort2.mat');

%% Pool data
% all mice
cropped_SI_all = [cropped_SI_sorted;cropped_SI_cb_sorted];
cropped_ID_all = [cropped_ID_sorted;cropped_ID_cb_sorted];

BL6csds_all_pre = struct;
BL6csds_all_pre.mouse = [csds.mouse;csds_cb.mouse];
BL6csds_all_pre.visits_plot_avg = [csds.visits_plot_avg;csds_cb.visits_plot_avg];
BL6csds_all_pre.num_visits = [csds.num_visits;csds_cb.num_visits];

[val,idx] = sort(cropped_SI_all);
BL6csds_all_pre.visits_plot_avg_sorted = BL6csds_all_pre.visits_plot_avg(idx,:);

% ctrl
BL6ctrl_all_pre = struct;
BL6ctrl_all_pre.mouse = [ctrl.mouse;ctrl_cb.mouse];
BL6ctrl_all_pre.visits_plot_avg = [ctrl.visits_plot_avg;ctrl_cb.visits_plot_avg];
BL6ctrl_all_pre.num_visits = [ctrl.num_visits;ctrl_cb.num_visits];

% res
BL6resDist_all_pre = struct;
BL6resDist_all_pre.mouse = [res.mouse;res_cb.mouse];
BL6resDist_all_pre.visits_trace = [res.visits_plot_avg;res_cb.visits_plot_avg];
BL6resDist_all_pre.num_visits = [res.num_visits;res_cb.num_visits];

% sus
BL6susDist_all_pre = struct;
BL6susDist_all_pre.mouse = [sus.mouse;sus_cb.mouse];
BL6susDist_all_pre.visits_trace = [sus.visits_plot_avg;sus_cb.visits_plot_avg];
BL6susDist_all_pre.num_visits = [sus.num_visits;sus_cb.num_visits];

%% Find avgs
fps= 30;
start_F = 3*fps; % start of proximity is at 3s 
stop_F = 4*fps; % end of proximity is at 4s

% all mice
% Favg 
BL6csds_all_pre.Favg = nanmean(BL6csds_all_pre.visits_plot_avg(:,start_F:stop_F),2);
BL6csds_all_pre.Fstd = nanstd(BL6csds_all_pre.visits_plot_avg(:,start_F:stop_F),0,2)/sqrt(size(BL6csds_all_pre.visits_plot_avg(:,start_F:stop_F),2));

% ctrl
% Favg 
BL6ctrl_all_pre.Favg = nanmean(BL6ctrl_all_pre.visits_plot_avg(:,start_F:stop_F),2);
BL6ctrl_all_pre.Fstd = nanstd(BL6ctrl_all_pre.visits_plot_avg(:,start_F:stop_F),0,2)/sqrt(size(BL6ctrl_all_pre.visits_plot_avg(:,start_F:stop_F),2));

% res
% trace
BL6resDist_all_pre.visits_trace_avg = nanmean(BL6resDist_all_pre.visits_trace,1);
BL6resDist_all_pre.visits_trace_std = nanstd(BL6resDist_all_pre.visits_trace,1)/sqrt(size(BL6resDist_all_pre.visits_trace,1));

% Favg
BL6resDist_all_pre.Favg = nanmean(BL6resDist_all_pre.visits_trace(:,start_F:stop_F),2);
BL6resDist_all_pre.Fstd = nanstd(BL6resDist_all_pre.visits_trace(:,start_F:stop_F),0,2)/sqrt(size(BL6resDist_all_pre.visits_trace(:,start_F:stop_F),2));

% sus
% trace
BL6susDist_all_pre.visits_trace_avg = nanmean(BL6susDist_all_pre.visits_trace,1);
BL6susDist_all_pre.visits_trace_std = nanstd(BL6susDist_all_pre.visits_trace,1)/sqrt(size(BL6susDist_all_pre.visits_trace,1));

% Favg
BL6susDist_all_pre.Favg = nanmean(BL6susDist_all_pre.visits_trace(:,start_F:stop_F),2);
BL6susDist_all_pre.Fstd = nanstd(BL6susDist_all_pre.visits_trace(:,start_F:stop_F),0,2)/sqrt(size(BL6susDist_all_pre.visits_trace(:,start_F:stop_F),2));

%% save
save('BL6Pre_visits_pooled_NEW.mat','cropped_ID_all','cropped_SI_all','BL6susDist_all_pre','BL6resDist_all_pre','BL6csds_all_pre','BL6ctrl_all_pre')