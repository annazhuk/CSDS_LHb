%% Collect all trials together 
clear;
load('mice_info_cohort1.mat')
total = size(mice_info.mice,1);

%% Save data into a struct
csds = struct;
ctrl = struct;

% csds mice
% cell bodies
k = 1;
for i = 2:total
    filename = 'm3800_SW_visitsplot_NEW.mat';
    filename(1:5) = mice_info.mice(i,:);
    load(filename)
    
    if mice_info.signal(i) == 1
        if sum(mice_info.id(i,:) == 'csds') == 4
            csds.mouse(k,1) = str2num(mice_info.mice(i,2:end));
            csds.visits_plot_avg(k,:) = visits_plot_avg;
            csds.visits_plot{k,1} = visits_plot;
            if isempty (csds.visits_plot{k})
                csds.visits_plot{k} = nan(size(csds.visits_plot_avg(k,:)));
            end
            k = k + 1;
        end
    end
end

% ctrl mice
% cell bodies
k = 1;
for i = 2:total
    filename = 'm3800_SW_visitsplot_NEW.mat';
    filename(1:5) = mice_info.mice(i,:);
    load(filename)
    
    if mice_info.signal(i) == 1
        if sum(mice_info.id(i,:) == 'ctrl') == 4
            ctrl.mouse(k,1) = str2num(mice_info.mice(i,2:end));
            ctrl.visits_plot_avg(k,:) = visits_plot_avg;
            ctrl.visits_plot{k,1} = visits_plot;
            if isempty (ctrl.visits_plot{k})
                ctrl.visits_plot{k} = nan(size(ctrl.visits_plot_avg(k,:)));
            end
            k = k + 1;
        end
    end
end

save('visits_plot_trials_cohort1.mat','csds','ctrl')

%% Split csds mice into res vs sus
clear;
neural_data = 'visits_plot_trials_cohort1.mat';
SI_data = 'SI_all_percent_fos.mat';
CSDSinfoID = 'CSDS_infoID_cohort1.mat';

load(neural_data)
load(SI_data)
load(CSDSinfoID)

SI_visits_res_sus_trials_cohort1(csds,SW_CSDS_Post,mouse_CSDS,SW_CTRL_Post,CSDS_ID,CSDS_SI)