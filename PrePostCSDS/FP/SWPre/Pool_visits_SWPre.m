clear;
%% Pool visits data (SW Post)
load('visits_plot_NEW_cohort1.mat');
load('SI_F_visits_res_sus_NEW_cohort1.mat');
load('visits_plot_NEW_cohort2.mat');
load('SI_F_visits_res_sus_NEW_cohort2.mat');

%% Pool data
% all mice
cropped_SI_all = [cropped_SI_sorted;cropped_SI_cb_sorted];
cropped_ID_all = [cropped_ID_sorted;cropped_ID_cb_sorted];

SWcsds_all_pre = struct;
SWcsds_all_pre.mouse = [csds.mouse;csds_cb.mouse];
SWcsds_all_pre.visits_plot_avg = [csds.visits_plot_avg;csds_cb.visits_plot_avg];
SWcsds_all_pre.num_visits = [csds.num_visits;csds_cb.num_visits];

[val,idx] = sort(cropped_SI_all);
SWcsds_all_pre.visits_plot_avg_sorted = SWcsds_all_pre.visits_plot_avg(idx,:);

% ctrl
SWctrl_all_pre = struct;
SWctrl_all_pre.mouse = [ctrl.mouse;ctrl_cb.mouse];
SWctrl_all_pre.visits_plot_avg = [ctrl.visits_plot_avg;ctrl_cb.visits_plot_avg];
SWctrl_all_pre.num_visits = [ctrl.num_visits;ctrl_cb.num_visits];

% res
SWresDist_all_pre = struct;
SWresDist_all_pre.mouse = [res.mouse;res_cb.mouse];
SWresDist_all_pre.visits_trace = [res.visits_plot_avg;res_cb.visits_plot_avg];
SWresDist_all_pre.num_visits = [res.num_visits;res_cb.num_visits];

% sus
SWsusDist_all_pre = struct;
SWsusDist_all_pre.mouse = [sus.mouse;sus_cb.mouse];
SWsusDist_all_pre.visits_trace = [sus.visits_plot_avg;sus_cb.visits_plot_avg];
SWsusDist_all_pre.num_visits = [sus.num_visits;sus_cb.num_visits];

%% Find avgs
fps= 30;
start_F = 3*fps; % start of proximity is at 3s 
stop_F = 4*fps; % end of proximity is at 4s

% all mice
% Favg 
SWcsds_all_pre.Favg = nanmean(SWcsds_all_pre.visits_plot_avg(:,start_F:stop_F),2);
SWcsds_all_pre.Fstd = nanstd(SWcsds_all_pre.visits_plot_avg(:,start_F:stop_F),0,2)/sqrt(size(SWcsds_all_pre.visits_plot_avg(:,start_F:stop_F),2));

% ctrl
% Favg 
SWctrl_all_pre.Favg = nanmean(SWctrl_all_pre.visits_plot_avg(:,start_F:stop_F),2);
SWctrl_all_pre.Fstd = nanstd(SWctrl_all_pre.visits_plot_avg(:,start_F:stop_F),0,2)/sqrt(size(SWctrl_all_pre.visits_plot_avg(:,start_F:stop_F),2));

% res
% trace
SWresDist_all_pre.visits_trace_avg = nanmean(SWresDist_all_pre.visits_trace,1);
SWresDist_all_pre.visits_trace_std = nanstd(SWresDist_all_pre.visits_trace,1)/sqrt(size(SWresDist_all_pre.visits_trace,1));

% Favg
SWresDist_all_pre.Favg = nanmean(SWresDist_all_pre.visits_trace(:,start_F:stop_F),2);
SWresDist_all_pre.Fstd = nanstd(SWresDist_all_pre.visits_trace(:,start_F:stop_F),0,2)/sqrt(size(SWresDist_all_pre.visits_trace(:,start_F:stop_F),2));

% sus
% trace
SWsusDist_all_pre.visits_trace_avg = nanmean(SWsusDist_all_pre.visits_trace,1);
SWsusDist_all_pre.visits_trace_std = nanstd(SWsusDist_all_pre.visits_trace,1)/sqrt(size(SWsusDist_all_pre.visits_trace,1));

% Favg
SWsusDist_all_pre.Favg = nanmean(SWsusDist_all_pre.visits_trace(:,start_F:stop_F),2);
SWsusDist_all_pre.Fstd = nanstd(SWsusDist_all_pre.visits_trace(:,start_F:stop_F),0,2)/sqrt(size(SWsusDist_all_pre.visits_trace(:,start_F:stop_F),2));

%% save
save('SWPre_visits_pooled_NEW.mat','cropped_ID_all','cropped_SI_all','SWsusDist_all_pre','SWresDist_all_pre','SWcsds_all_pre','SWctrl_all_pre')