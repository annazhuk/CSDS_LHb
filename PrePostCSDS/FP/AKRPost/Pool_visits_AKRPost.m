clear;
%% Pool visits data (AKR Post)
% Only one cohort 
load('visits_plot_NEW_cohort2.mat');
load('SI_F_visits_res_sus_NEW_cohort2.mat');

%% Pool data
% all mice
cropped_SI_all = cropped_SI_cb_sorted;
cropped_ID_all = cropped_ID_cb_sorted;

AKRcsds_all_post = struct;
AKRcsds_all_post.mouse = csds_cb.mouse;
AKRcsds_all_post.visits_plot_avg = csds_cb.visits_plot_avg;
AKRcsds_all_post.num_visits = csds_cb.num_visits;

[val,idx] = sort(cropped_SI_all);
AKRcsds_all_post.visits_plot_avg_sorted = AKRcsds_all_post.visits_plot_avg(idx,:);

% ctrl
AKRctrl_all_post = struct;
AKRctrl_all_post.mouse = ctrl_cb.mouse;
AKRctrl_all_post.visits_plot_avg = ctrl_cb.visits_plot_avg;
AKRctrl_all_post.num_visits = ctrl_cb.num_visits;

% res
AKRresDist_all_post = struct;
AKRresDist_all_post.mouse = res_cb.mouse;
AKRresDist_all_post.visits_trace = res_cb.visits_plot_avg;
AKRresDist_all_post.num_visits = res_cb.num_visits;

% sus
AKRsusDist_all_post = struct;
AKRsusDist_all_post.mouse = sus_cb.mouse;
AKRsusDist_all_post.visits_trace = sus_cb.visits_plot_avg;
AKRsusDist_all_post.num_visits = sus_cb.num_visits;

%% Find avgs
fps= 30;
start_F = 3*fps; % start of proximity is at 3s 
stop_F = 4*fps; % end of proximity is at 4s

% all mice
% Favg 
AKRcsds_all_post.Favg = nanmean(AKRcsds_all_post.visits_plot_avg(:,start_F:stop_F),2);
AKRcsds_all_post.Fstd = nanstd(AKRcsds_all_post.visits_plot_avg(:,start_F:stop_F),0,2)/sqrt(size(AKRcsds_all_post.visits_plot_avg(:,start_F:stop_F),2));

% ctrl
% Favg 
AKRctrl_all_post.Favg = nanmean(AKRctrl_all_post.visits_plot_avg(:,start_F:stop_F),2);
AKRctrl_all_post.Fstd = nanstd(AKRctrl_all_post.visits_plot_avg(:,start_F:stop_F),0,2)/sqrt(size(AKRctrl_all_post.visits_plot_avg(:,start_F:stop_F),2));

% res
% trace
AKRresDist_all_post.visits_trace_avg = nanmean(AKRresDist_all_post.visits_trace,1);
AKRresDist_all_post.visits_trace_std = nanstd(AKRresDist_all_post.visits_trace,1)/sqrt(size(AKRresDist_all_post.visits_trace,1));

% Favg
AKRresDist_all_post.Favg = nanmean(AKRresDist_all_post.visits_trace(:,start_F:stop_F),2);
AKRresDist_all_post.Fstd = nanstd(AKRresDist_all_post.visits_trace(:,start_F:stop_F),0,2)/sqrt(size(AKRresDist_all_post.visits_trace(:,start_F:stop_F),2));

% sus
% trace
AKRsusDist_all_post.visits_trace_avg = nanmean(AKRsusDist_all_post.visits_trace,1);
AKRsusDist_all_post.visits_trace_std = nanstd(AKRsusDist_all_post.visits_trace,1)/sqrt(size(AKRsusDist_all_post.visits_trace,1));

% Favg
AKRsusDist_all_post.Favg = nanmean(AKRsusDist_all_post.visits_trace(:,start_F:stop_F),2);
AKRsusDist_all_post.Fstd = nanstd(AKRsusDist_all_post.visits_trace(:,start_F:stop_F),0,2)/sqrt(size(AKRsusDist_all_post.visits_trace(:,start_F:stop_F),2));

%% save
save('AKRPost_visits_pooled_NEW.mat','cropped_ID_all','cropped_SI_all','AKRsusDist_all_post','AKRresDist_all_post','AKRcsds_all_post','AKRctrl_all_post')