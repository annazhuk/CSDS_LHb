%% load data
clear;
load('mice.mat')

%% Sync neuronal data with camera frames
F_filename = 'm555_SI_neurondata_split.mat';

scope_FR = 25; 
cond = 'SWPost_baseline';

for mousenum = 1:size(mice,1)
    Fdata_file = 'm555_SW_Post_baseline.sch.csv';

    Fdata_file(2:4) = string(mice(mousenum));
    F_filename(2:4) = string(mice(mousenum));
    load(F_filename)
    mousename = string(mice(mousenum));

    % Get video data csv
    %[file,path] = uigetfile('*.mp4');
    videoname = 'Trial     4_xvid.mp4';%fullfile(path,file);

    Sync_scope_cell_etho(Fdata_file,alldeltaF_baseline,mousename,cond,scope_FR,videoname)
end

%% Find peaks in baseline pre-CSDS 
filename = 'celldata_555_SWPost_baseline.mat';
for m = 1:size(mice,1)

    filename(10:12) = string(mice(m,:));
    load(filename)
    mouse = filename(10:12);
    
    alldeltaF_baseline = alldeltaF_social;

    threshold = 1;

    findpeaks_baseline_cells_S(alldeltaF_baseline,mouse,threshold)
end

%% Create rasterplot of activity
filename = '555_baseline_pks_cells_S.mat';

for m = 1:size(mice,1)
    filename(1:3) = num2str(mice(m));
    load(filename)    
    threshold = 1;
    
    baseline_peaks_S_cosine(alldeltaF_baseline,alldeltaF_pkdata,threshold,mouse)
end

%% Collect data into struct
allmice = struct;
filename = '555_baseline_pks_cells_S.mat';
for m = 1:size(mice,1)
    filename(1:3) = string(mice(m,:));
    load(filename)

    allmice.mice(m,:) = str2num(mouse);
    allmice.pk_totals{m,1} = nansum(alldeltaF_pkdata.raster,2);
end
save('allmice_peaks_baseline_cells_S.mat','allmice')

%% Split into res and sus
clear;
peaks_data = 'allmice_peaks_baseline_cells_S.mat';
SI_data = 'SI_percent_ALL.mat';
CSDSinfoID = 'CSDS_infoID.mat';

load(peaks_data)
load(SI_data)
load(CSDSinfoID)

baseline_peaks_res_sus_cells_S(allmice,SW_CSDS_Post,mouse_CSDS,SW_CTRL_Post,CSDS_ID,CSDS_SI)

%% Create plots pooling cells
clear;
load('mice.mat')
load('baseline_peaks_res_sus_cells_S_NEW.mat')

[val,idx] = sort(cropped_SI_sorted_im);

allmice_im.pk_totals_all_sorted = [];
for mouse = 1:size(allmice_im.mice,1)
    allmice_im.pk_totals_all_sorted = [allmice_im.pk_totals_all_sorted; allmice_im.pk_totals{idx(mouse)}];
end

sus_cell_total = 0;
for m = 1:size(sus_im.mice,1)
    sus_cell_total = sus_cell_total + length(sus_im.pk_totals{m});
end

res_cell_total = 0;
for m = 1:size(res_im.mice,1)
    res_cell_total = res_cell_total + length(res_im.pk_totals{m});
end

% sus
sus_im.pk_totals_all_sorted = allmice_im.pk_totals_all_sorted(1:sus_cell_total);

% res
res_im.pk_totals_all_sorted = allmice_im.pk_totals_all_sorted(sus_cell_total+1:end);

% find avgs
% sus
sus_im.pk_totals_avg = nanmean(allmice_im.pk_totals_all_sorted(1:sus_cell_total));

% res
res_im.pk_totals_avg = nanmean(allmice_im.pk_totals_all_sorted(sus_cell_total+1:end));

% find std
% sus
sus_im.pk_totals_std = nanstd(allmice_im.pk_totals_all_sorted(1:sus_cell_total))/sqrt(size(allmice_im.pk_totals_all_sorted(1:sus_cell_total),1));

% res
res_im.pk_totals_std = nanstd(allmice_im.pk_totals_all_sorted(sus_cell_total+1:end))/sqrt(size(allmice_im.pk_totals_all_sorted(sus_cell_total+1:end),1));

%% Find avg of each mouse
% res
Post_res_im = res_im;
for m = 1:size(Post_res_im.mice,1)
    Post_res_im.peak_totals_mouse(m,1) = nanmean(Post_res_im.pk_totals{m});
end

% sus
Post_sus_im = sus_im;
for m = 1:size(Post_sus_im.mice,1)
    Post_sus_im.peak_totals_mouse(m,1) = nanmean(Post_sus_im.pk_totals{m});
end

%% save
sus_cell_total_post = sus_cell_total;
res_cell_total_post = res_cell_total;
save('baseline_peaks_POST_res_sus_NEW.mat','Post_sus_im','Post_res_im','sus_cell_total_post','res_cell_total_post','cropped_SI_sorted_im','cropped_ID_sorted_im')