%% Collect all trials together 
clear;
load('mice_info_cohort2.mat')
total = size(mice_info.mice,1);

%% Save data into a struct
csds_cb = struct;
ctrl_cb = struct;

% csds mice
% cell bodies
k = 1;
for i = 1:total
    filename = 'm3801_SW_visitsplot_NEW.mat';
    filename(1:5) = mice_info.mice(i,:);
    load(filename)
    
    if mice_info.signal_cb(i) == 1
        if sum(mice_info.id(i,:) == 'csds') == 4
            csds_cb.mouse(k,1) = str2num(mice_info.mice(i,2:end));
            csds_cb.visits_plot_avg(k,:) = visits_plot_avg;
            csds_cb.visits_plot{k,1} = visits_plot;
            if isempty (csds_cb.visits_plot{k})
                csds_cb.visits_plot{k} = nan(size(csds_cb.visits_plot_avg(k,:)));
            end
            k = k + 1;
        end
    end
end

% ctrl mice
% cell bodies
k = 1;
for i = 1:total
    filename = 'm3801_SW_visitsplot_NEW.mat';
    filename(1:5) = mice_info.mice(i,:);
    load(filename)
    
    if mice_info.signal_cb(i) == 1
        if sum(mice_info.id(i,:) == 'ctrl') == 4
            ctrl_cb.mouse(k,1) = str2num(mice_info.mice(i,2:end));
            ctrl_cb.visits_plot_avg(k,:) = visits_plot_avg;
            ctrl_cb.visits_plot{k,1} = visits_plot;
            if isempty (ctrl_cb.visits_plot{k})
                ctrl_cb.visits_plot{k} = nan(size(ctrl_cb.visits_plot_avg(k,:)));
            end
            k = k + 1;
        end
    end
end

save('visits_plot_trials_cohort2.mat','csds_cb','ctrl_cb')

%% Split csds mice into res vs sus
clear;
neural_data = 'visits_plot_trials_cohort2.mat';
SI_data = 'SI_all_percent_fos.mat';
CSDSinfoID = 'CSDS_infoID_cohort2.mat';

load(neural_data)
load(SI_data)
load(CSDSinfoID)

SI_visits_res_sus_trials_cohort2(csds_cb,SW_CSDS_Post,mouse_CSDS,SW_CTRL_Post,CSDS_ID,CSDS_SI)