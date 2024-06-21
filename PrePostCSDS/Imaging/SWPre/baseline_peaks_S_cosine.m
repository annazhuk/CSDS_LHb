function baseline_peaks_S_cosine(alldeltaF_baseline,alldeltaF_pkdata,threshold,timebin,mouse)
% Find correlations using CNMFE identified peaks
F_mean = nanmean(alldeltaF_baseline.S_sync,2);
F_std = nanstd(alldeltaF_baseline.S_sync,0,2);

% threshold = 1;

alldeltaF_pkdata.raster = zeros(size(alldeltaF_baseline.S_sync));
for cellnum = 1:size(alldeltaF_baseline.S_sync,1)
    alldeltaF_baseline.S_zscore(cellnum,:) = (alldeltaF_baseline.S_sync(cellnum,:) - F_mean(cellnum))/F_std(cellnum);
    [val,loc] = findpeaks(alldeltaF_baseline.S_zscore(cellnum,:),'MinPeakHeight',threshold);
    
    alldeltaF_pkdata.raster(cellnum,loc) = 1;

    alldeltaF_pkdata.pks{cellnum,1} = val;
    alldeltaF_pkdata.locs{cellnum,1} = loc;
end

total_cell = size(alldeltaF_baseline.S_zscore,1);

if size(alldeltaF_pkdata.raster,2)>8992
    alldeltaF_pkdata.raster(:,8993:end) = [];
end



%% Bin peaks within a time window and get average activity across cells 
fps = 30;
bin_size = timebin*fps; %0.5*fps; % 5s
edges = 1:bin_size:size(alldeltaF_pkdata.raster,2);

edges = [edges,size(alldeltaF_pkdata.raster,2)];

alldeltaF_pkdata.binned = nan(size(alldeltaF_pkdata.raster,1),length(edges)-1);
for numcell = 1:size(alldeltaF_pkdata.raster,1)
    for bins = 1:length(edges)-1
        alldeltaF_pkdata.binned(numcell,bins) = nansum(alldeltaF_pkdata.raster(numcell,edges(bins):edges(bins+1)));
    end
end

% Convert into active vs inactive 
alldeltaF_pkdata.activity = alldeltaF_pkdata.binned>0;

total_cell = size(alldeltaF_pkdata.raster,1);

alldeltaF_pkdata.activity_avg = nansum(alldeltaF_pkdata.activity,1)/total_cell*100;
%edges_05s = edges;

%% Plot binned avg
% figure('DefaultAxesFontSize',12)
% subplot(2,1,1)
% imagesc(alldeltaF_pkdata.raster)
% ylabel('Cells')
% xlabel('Frames')
% box off
% 
% subplot(2,1,2)
% plot(edges(1:end-1),alldeltaF_pkdata.activity_avg,'LineWidth',1)
% xlim([1 9000])
% ylim([0 100])
% ylabel('% cells active')
% xlabel('Frames')
% title('500 ms bins')
% box off

%% Plot binned avg
% figure('DefaultAxesFontSize',12)
% subplot(3,1,1)
% imagesc(alldeltaF_pkdata.raster)
% ylabel('Cells')
% xlabel('Frames')
% box off
% 
% subplot(3,1,2)
% plot(edges_01s(1:end-1),alldeltaF_pkdata.activity_avg_01s,'LineWidth',1)
% xlim([1 9000])
% ylim([0 100])
% ylabel('% cells active')
% xlabel('Frames')
% title('100 ms bins')
% box off
% 
% subplot(3,1,3)
% plot(edges_05s(1:end-1),alldeltaF_pkdata.activity_avg_05s,'LineWidth',1)
% xlim([1 9000])
% ylim([0 100])
% ylabel('% cells active')
% xlabel('Frames')
% title('500 ms bins')
% box off

% subplot(5,1,4)
% plot(edges_10s(1:end-1),alldeltaF_pkdata.activity_avg_10s,'LineWidth',1)
% xlim([1 9000])
% ylim([0 100])
% ylabel('% cells active')
% xlabel('Frames')
% title('10s bins')
% box off
% 
% subplot(5,1,5)
% plot(edges_30s(1:end-1),alldeltaF_pkdata.activity_avg_30s,'LineWidth',1)
% xlim([1 9000])
% ylim([0 100])
% ylabel('% cells active')
% xlabel('Frames')
% title('30s bins')
% box off
%% Plot raster plot
% figure('DefaultAxesFontSize',12)
% imagesc(alldeltaF_pkdata.raster)
% ylabel('Cells')
% xlabel('Frames')
% box off

%% Get the correlation values for each cell
% cos_val = nan(size(alldeltaF_pkdata.raster,1),size(alldeltaF_pkdata.raster,1));
% for num_cell = 1:size(alldeltaF_pkdata.raster,1)
%     for cell_corr = 1:size(alldeltaF_pkdata.raster,1)
%         cos_val(num_cell,cell_corr) = dot(alldeltaF_pkdata.raster(num_cell,:),alldeltaF_pkdata.raster(cell_corr,:))/(norm(alldeltaF_pkdata.raster(num_cell,:))*norm(alldeltaF_pkdata.raster(cell_corr,:)));
%     end
% end
% alldeltaF_pkdata.cos_val = cos_val;

%% save
save(strcat(mouse,'_baseline_pks_cells_S'),'alldeltaF_baseline','alldeltaF_pkdata','mouse','total_cell')