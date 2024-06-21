clear;
%% Pool visits data (SW Post)
load('visits_plot_NEW_cohort1.mat');
load('SI_F_visits_res_sus_trials_cohort1.mat');
load('visits_plot_NEW_cohort2.mat');
load('SI_F_visits_res_sus_trials_cohort2.mat');

%% Pool data
% all mice
cropped_SI_all = [cropped_SI_sorted;cropped_SI_cb_sorted];
cropped_ID_all = [cropped_ID_sorted;cropped_ID_cb_sorted];

SWcsds_all_post = struct;
SWcsds_all_post.mouse = [csds.mouse;csds_cb.mouse];
SWcsds_all_post.visits_plot = [csds.visits_plot;csds_cb.visits_plot];
SWcsds_all_post.visits_plot_avg = [csds.visits_plot_avg;csds_cb.visits_plot_avg];

[val,idx] = sort(cropped_SI_all);
SWcsds_all_post.visits_plot_avg_sorted = SWcsds_all_post.visits_plot_avg(idx,:);

% res
SWresDist_all_post = struct;
SWresDist_all_post.mouse = [res.mouse;res_cb.mouse];
SWresDist_all_post.visits_plot = [res.visits_plot;res_cb.visits_plot];
SWresDist_all_post.visits_trace = [res.visits_plot_avg;res_cb.visits_plot_avg];

% sus
SWsusDist_all_post = struct;
SWsusDist_all_post.mouse = [sus.mouse;sus_cb.mouse];
SWsusDist_all_post.visits_plot = [sus.visits_plot;sus_cb.visits_plot];
SWsusDist_all_post.visits_trace = [sus.visits_plot_avg;sus_cb.visits_plot_avg];

%% Calculate peak for each visit
fps = 30;
start = 3*fps;
stop = 4*fps;

% all mice
for m = 1:size(SWcsds_all_post.mouse,1)    
    SWcsds_all_post.peaks{m,1} = nanmean(SWcsds_all_post.visits_plot{m}(:,start:stop),2);
end

% res
for m = 1:size(SWresDist_all_post.mouse,1)    
    SWresDist_all_post.peaks{m,1} = nanmean(SWresDist_all_post.visits_plot{m}(:,start:stop),2);
end

% sus
for m = 1:size(SWsusDist_all_post.mouse,1)    
    SWsusDist_all_post.peaks{m,1} = nanmean(SWsusDist_all_post.visits_plot{m}(:,start:stop),2);
end

%% Find Average and Std of each trial
% Find max number of trials
for m = 1:size(SWresDist_all_post.mouse,1)
    SWresDist_all_post.num_visits(m,1) = length(SWresDist_all_post.peaks{m});
end
res_visit = max(SWresDist_all_post.num_visits);

for m = 1:size(SWsusDist_all_post.mouse,1)
    SWsusDist_all_post.num_visits(m,1) = length(SWsusDist_all_post.peaks{m});
end
sus_visit = max(SWsusDist_all_post.num_visits);

% Sort peaks by visit number instead of mouse
for v = 1:res_visit
    all_peaks = [];
    for m = 1:size(SWresDist_all_post.mouse,1) 
        if length(SWresDist_all_post.peaks{m}) >= v
            all_peaks = [all_peaks,SWresDist_all_post.peaks{m}(v)];
        end
    end
    SWresDist_all_post.peaks_by_trial{v,1} = all_peaks;
    SWresDist_all_post.trial_avg(v,1) = nanmean(SWresDist_all_post.peaks_by_trial{v,1});
    SWresDist_all_post.trial_std(v,1) = nanstd(SWresDist_all_post.peaks_by_trial{v,1})/sqrt(length(SWresDist_all_post.peaks_by_trial{v,1}));
end

% Sort peaks by visit number instead of mouse
for v = 1:sus_visit
    all_peaks = [];
    for m = 1:size(SWsusDist_all_post.mouse,1) 
        if length(SWsusDist_all_post.peaks{m}) >= v
            all_peaks = [all_peaks,SWsusDist_all_post.peaks{m}(v)];
        end
    end
    SWsusDist_all_post.peaks_by_trial{v,1} = all_peaks;
    SWsusDist_all_post.trial_avg(v,1) = nanmean(SWsusDist_all_post.peaks_by_trial{v,1});
    SWsusDist_all_post.trial_std(v,1) = nanstd(SWsusDist_all_post.peaks_by_trial{v,1})/sqrt(length(SWsusDist_all_post.peaks_by_trial{v,1}));
end

%% Find find decay for each mouse
for m = 1:size(SWresDist_all_post.mouse,1)
    num_trial = 1:length(SWresDist_all_post.peaks{m});
    model = fitlm(num_trial,SWresDist_all_post.peaks{m});
    SWresDist_all_post.model{m,1} = model;
    SWresDist_all_post.p_val(m,1) = model.Coefficients{2,4};
    SWresDist_all_post.decay(m,1) = model.Coefficients{2,1};
    SWresDist_all_post.sig(m,1) = model.Coefficients{2,4}<0.05;
end

for m = 1:size(SWsusDist_all_post.mouse,1)
    num_trial = 1:length(SWsusDist_all_post.peaks{m});
    model = fitlm(num_trial,SWsusDist_all_post.peaks{m});
    SWsusDist_all_post.model{m,1} = model;
    SWsusDist_all_post.p_val(m,1) = model.Coefficients{2,4};
    SWsusDist_all_post.decay(m,1) = model.Coefficients{2,1};
    SWsusDist_all_post.sig(m,1) = model.Coefficients{2,4}<0.05;
end

%% save
save('SWPost_visits_trials_pooled_NEW.mat','cropped_ID_all','cropped_SI_all','SWsusDist_all_post','SWresDist_all_post','SWcsds_all_post')