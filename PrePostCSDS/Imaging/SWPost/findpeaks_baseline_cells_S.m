function findpeaks_baseline_cells_S(alldeltaF_baseline,mouse,threshold)
% Find the zscore 
F_mean = nanmean(alldeltaF_baseline.S_sync,2); % F_mean = movmean(F_final,1000,'omitnan');
F_std = nanstd(alldeltaF_baseline.S_sync,0,2);

alldeltaF_baseline.S_zscore = nan(size(alldeltaF_baseline.S_sync));

alldeltaF_pkdata = struct;
for numcell = 1:size(alldeltaF_baseline.S_sync,1)
    % z-score
    alldeltaF_baseline.S_zscore(numcell,:) = (alldeltaF_baseline.S_sync(numcell,:) - F_mean(numcell))/F_std(numcell); 

    % find peaks
    [pks,locs] = findpeaks(alldeltaF_baseline.S_zscore(numcell,:),'MinPeakHeight',threshold);

    alldeltaF_pkdata.Szscore(numcell,:) = alldeltaF_baseline.S_zscore(numcell,:);
    alldeltaF_pkdata.pks{numcell,1} = pks;
    alldeltaF_pkdata.locs{numcell,1} = locs;
end

%% Plot

%% save
save(strcat(mouse,'_baseline_pks_cells_S.mat'),'alldeltaF_baseline','alldeltaF_pkdata','mouse')

end