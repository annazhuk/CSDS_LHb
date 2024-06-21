%% 
clear;
load('mice_info.mat')

%% Find cells that have sig response to aggressor
shuffle_time = 5;  % in seconds
prct = 97.5;

filename = '555_SWpost_cells_visits_nose_NEW.mat';
% for example plots use m569 and numcell = 5 % for example plots use m577 and numcell = 4
for m = 1:size(mice_info.mice,1)
    filename(1:3) = mice_info.mice(m,2:end);
    load(filename)
    nose_imaging_cells_shuffles_FIXED(alldeltaF_social,time_stamps,mousename,shuffle_time,prct)
end

%% Collect data into struct
allmice = struct;
filename = '555_sig_nose_cells_FIXED.mat';

for m = 1:size(mice_info.mice,1)
    filename(1:3) = mice_info.mice(m,2:end);
    load(filename)

    allmice.mice(m,1) = str2double(mice_info.mice(m,2:end));
    allmice.sig_up_total(m,1) = sum(sig_up_start);
    allmice.sig_down_total(m,1) = sum(sig_down_start);
    allmice.sig_up_percent(m,1) = sum(sig_up_start)/length(sig_up_start)*100;
    allmice.sig_down_percent(m,1) = sum(sig_down_start)/length(sig_down_start)*100;
    allmice.pk_real_sig(m,1) = nanmean(pk_real); %nanmean(peak_max);
    %allmice.peak_min_sig(m,1) = nanmean(peak_min);
    allmice.pk_real_sig_all{m,1} = pk_real; %peak_max;
    %allmice.peak_min_sig_all{m,1} = peak_min;
    allmice.sig_up_all{m,1} = sig_up_start;
    allmice.sig_down_all{m,1} = sig_down_start;
end
% save
save('allmice_sig_nose_FIXED.mat','allmice')

%% Split into res and sus
clear; 
load('mice.mat')
load('SI_percent_ALL.mat')
load('allmice_sig_nose_FIXED.mat')
load('CSDS_infoID.mat')

nose_sig_cells_sus_res_FIXED(allmice,SW_CSDS_Post,mouse_CSDS,SW_CTRL_Post,CSDS_ID,CSDS_SI)