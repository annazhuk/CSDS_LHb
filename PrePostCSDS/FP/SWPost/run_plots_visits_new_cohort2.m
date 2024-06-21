clear;
load('mice_info_cohort2.mat')
total = size(mice_info.mice,1);

%% Timelock F data to social target proximity
for i = 1:total
    F_data = 'data_foranalysis_m3801_SW.mat';
    dist_data = 'm3801_SW_distcent_data.mat';
    F_data(18:22) = mice_info.mice(i,:);
    dist_data(1:5) = mice_info.mice(i,:);
    
    load(F_data)
    load(dist_data)
    
    dist_visits_NEW_cohort2(dist_center_filtered,F_sync_cb,mouse,strain(1:3));
end

%% Save data into structs for CSDS and CTRL mice
clear;
load('mice_info_cohort2.mat')
total = size(mice_info.mice,1);

csds_cb = struct;
ctrl_cb = struct;

% csds mice
k = 1;
for i = 1:total
    filename = 'm3801_SW_visitsplot_NEW.mat';
    filename(1:5) = mice_info.mice(i,:);
    load(filename)
    
    if mice_info.signal_cb(i) == 1
        if sum(mice_info.id(i,:) == 'csds') == 4
            csds_cb.mouse(k,1) = str2num(mice_info.mice(i,2:end));
            csds_cb.visits_plot_avg(k,:) = visits_plot_avg;
            csds_cb.num_visits(k,1) = size(visits_plot,1);
            k = k + 1;
        end
    end
end

% ctrl mice
k = 1;
for i = 1:total
    filename = 'm3801_SW_visitsplot_NEW.mat';
    filename(1:5) = mice_info.mice(i,:);
    load(filename)
    
    if mice_info.signal_cb(i) == 1
        if sum(mice_info.id(i,:) == 'ctrl') == 4
            ctrl_cb.mouse(k,1) = str2num(mice_info.mice(i,2:end));
            ctrl_cb.visits_plot_avg(k,:) = visits_plot_avg;
            ctrl_cb.num_visits(k,1) = size(visits_plot,1);
            k = k + 1;
        end
    end
end

save('visits_plot_NEW_cohort2.mat','csds_cb','ctrl_cb')

%% Split csds mice into res vs sus
clear;
neural_data = 'visits_plot_NEW_cohort2.mat';
SI_data = 'SI_all_percent_fos.mat';
CSDSinfoID = 'CSDS_infoID_cohort2.mat';

load(neural_data)
load(SI_data)
load(CSDSinfoID)

SI_visits_res_sus_NEW_cohort2(csds_cb,SW_CSDS_Post,mouse_CSDS,SW_CTRL_Post,CSDS_ID,CSDS_SI)
