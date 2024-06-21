%% load data
clear;
load('mice.mat')

%% get correlations 
filename = '555_baseline_pks_cells_S.mat';

for m = 1:size(mice,1)
    filename(1:3) = num2str(mice(m));
    load(filename)    

    threshold = 1;
    timebin = 0.5; %500ms time bins
    
    baseline_peaks_S_cosine(alldeltaF_baseline,alldeltaF_pkdata,threshold,timebin,mouse)
end

%% Plot a histogram of % of synchrony for each mouse
filename = '555_baseline_pks_cells_S.mat';

%edges = 0:5:100;
allmice = struct;
for m = 1:size(mice,1)
    filename(1:3) = num2str(mice(m));
    load(filename)    

    allmice.mice(m,1) = mice(m);
    allmice.activity_avg(m,:) = alldeltaF_pkdata.activity_avg;

    allmice.total_cells(m,1) = total_cell;
end
save('allmice_correlations_S.mat','allmice')

%% Split into res and sus
clear;
peaks_data = 'allmice_correlations_S.mat';
SI_data = 'SI_percent_ALL.mat';
CSDSinfoID = 'CSDS_infoID.mat';

load(peaks_data)
load(SI_data)
load(CSDSinfoID)

correlation_res_sus_S(allmice,SW_CSDS_Post,mouse_CSDS,SW_CTRL_Post,CSDS_ID,CSDS_SI)

%% Avg with excluding mice with < 10 cells
% exclude mice with < 9 cells
clear;
load('correlations_res_sus_S.mat')

cutoff = 9;
sus_im.activity_avg_ex = sus_im.activity_avg(sus_im.total_cells>=cutoff,:);
res_im.activity_avg_ex = res_im.activity_avg(res_im.total_cells>=cutoff,:);

% find avg
sus_im.activity_avg_ex_avg = nanmean(sus_im.activity_avg_ex,1);
res_im.activity_avg_ex_avg = nanmean(res_im.activity_avg_ex,1);

% find std
sus_im.activity_std_ex_std = nanmean(sus_im.activity_avg_ex,1)/sqrt(size(sus_im.activity_avg_ex,1));
res_im.activity_std_ex_std = nanmean(res_im.activity_avg_ex,1)/sqrt(size(res_im.activity_avg_ex,1));

%% Plot CDF (each mouse)
% colors
bin_size = 1;

% sus
for m = 1:size(sus_im.activity_avg_ex,1)
    bins_sus = min(sus_im.activity_avg_ex(m,:)):bin_size:max(sus_im.activity_avg_ex(m,:));
    cdf_sus = cumsum(histcounts(sus_im.activity_avg_ex(m,:),length(bins_sus)))/length(sus_im.activity_avg_ex(m,:))*100;

    sus_im.bin_sus_all{m,1} = bins_sus;
    sus_im.cdf_sus_all{m,1} = cdf_sus;

    sus_im.lengths(m,1) = length(bins_sus);
    
   
end

% res
for m = 1:size(res_im.activity_avg_ex,1)
    bins_res = min(res_im.activity_avg_ex(m,:)):bin_size:max(res_im.activity_avg_ex(m,:));
    cdf_res = cumsum(histcounts(res_im.activity_avg_ex(m,:),length(bins_res)))/length(res_im.activity_avg_ex(m,:))*100;

    res_im.bin_res_all{m,1} = bins_res;
    res_im.cdf_res_all{m,1} = cdf_res;

    res_im.lengths(m,1) = length(bins_res);
end

%% Calcuate CDF (avg)
% sus
max_bins = max(sus_im.lengths);
for m = 1:size(sus_im.activity_avg_ex,1)
    diff_length = max_bins - length(sus_im.cdf_sus_all{m});
    sus_im.cdf_sus_matrix(m,:) = [sus_im.cdf_sus_all{m},nan(1,diff_length)];
end
sus_im.cdf_sus_matrix_avg = nanmean(sus_im.cdf_sus_matrix,1);
sus_im.cdf_sus_matrix_std = nanstd(sus_im.cdf_sus_matrix,1)/sqrt(size(sus_im.cdf_sus_matrix,1));

% res
max_bins = max(res_im.lengths);
for m = 1:size(res_im.activity_avg_ex,1)
    diff_length = max_bins - length(res_im.cdf_res_all{m});
    res_im.cdf_res_matrix(m,:) = [res_im.cdf_res_all{m},nan(1,diff_length)];
end
res_im.cdf_res_matrix_avg = nanmean(res_im.cdf_res_matrix,1);
res_im.cdf_res_matrix_std = nanstd(res_im.cdf_res_matrix,1)/sqrt(size(res_im.cdf_res_matrix,1));

%%
save('correlations_res_sus_S_CDF.mat','res_im','sus_im','cropped_SI_sorted_im','cropped_ID_sorted_im',"cdf_res",'cdf_sus','bins_res','bins_sus','bin_size','allmice')

%% Plot new CDF avg and std
% colors
green = [.1367,.5625,.4102];%[.4 .7 .15];
purple = [.4921,.4063,.6797];%[.5 0 .9];

figure('Renderer', 'painters', 'Position', [10 10 130 110],'DefaultAxesFontSize',12)

%s = errorbar(sus_im.cdf_sus_matrix_avg,sus_im.cdf_sus_matrix_std,'k');
hold on
%r = errorbar(res_im.cdf_res_matrix_avg,res_im.cdf_res_matrix_std,'k');

plot(sus_im.cdf_sus_matrix_avg,'color',green,'linewidth',2)

plot(res_im.cdf_res_matrix_avg,'color',purple,'linewidth',2)

ylim([0 100])
xlabel('% Population active')
ylabel('Cumulative probability')
box off
