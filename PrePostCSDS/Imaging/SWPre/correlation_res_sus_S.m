function correlation_res_sus_S(allmice,SW_CSDS_Post,mouse_CSDS,SW_CTRL_Post,CSDS_ID,CSDS_SI)

%% Find where to split
avg = nanmean(SW_CTRL_Post);
std_val = nanstd(SW_CTRL_Post);
med = avg - std_val; %med = floor(avg - std_val);


%med = median(SW_CSDS_Post);

SI_CSDS = SW_CSDS_Post;

resID = mouse_CSDS(find(SI_CSDS>=med));
susID = mouse_CSDS(find(SI_CSDS<med));

res_im = struct;
sus_im = struct;

%% Split mice into res and sus
% cell bodies
% res
r = 1;
for i = 1:size(allmice.mice,1)
    if ismember(allmice.mice(i),resID)
        res_im.mice(r,1) = allmice.mice(i);
        % res_im.raster{r,1} = corr_val.raster{i,1};
        % res_im.cos_val{r,1} = corr_val.cos_val{i,1};
        % res_im.cos_val_mean(r,1) = corr_val.cos_val_mean(i,1);
        res_im.activity_avg(r,:) = allmice.activity_avg(i,:);
        res_im.total_cells(r,1) = allmice.total_cells(i,:);
        r = r + 1;
    end
end

% sus
s = 1;
for i = 1:size(allmice.mice,1)
    if ismember(allmice.mice(i),susID)
        sus_im.mice(s,1) = allmice.mice(i);
        % sus_im.raster{s,1} = allmice.raster{i,1};
        % sus_im.cos_val{s,1} = allmice.cos_val{i,1};
        % sus_im.cos_val_mean(s,1) = allmice.cos_val_mean(i,1);
        sus_im.activity_avg(s,:) = allmice.activity_avg(i,:);
        sus_im.total_cells(s,1) = allmice.total_cells(i,:);
        s = s + 1;
    end
end

%% Find avgs
% res 
% res_im.cos_val_overall_avg = nanmean(res_im.cos_val_mean);
% res_im.cos_val_overall_std = nanstd(res_im.cos_val_mean)/sqrt(length(res_im.cos_val_mean));
res_im.activity_avg_overall = nanmean(res_im.activity_avg,1);

% sus 
% sus_im.cos_val_overall_avg = nanmean(sus_im.cos_val_mean);
% sus_im.cos_val_overall_std = nanstd(sus_im.cos_val_mean)/sqrt(length(sus_im.cos_val_mean));
sus_im.activity_avg_overall = nanmean(sus_im.activity_avg,1);

%% Sort
% Crop out mice npt analyzed (bad signal, gfp, etc)
% load SI data and make sure the IDs and SIs are in the correct order!

CSDS_SI = CSDS_SI/300*100;

% cell bodies
cropped_ID = CSDS_ID(find(ismember(CSDS_ID,allmice.mice)));
cropped_SI = CSDS_SI(find(ismember(CSDS_ID,allmice.mice)));
% make sure croppedSIs are in the same order as struct with Favg data
[~,idx] = sort(allmice.mice);

% Favg data and IDs should be in order, but double check this 
if sum(idx == (1:length(idx))') == length(idx)
    disp('sorted order is correct')
else
    warning('Error in the mouse ID ordering!')
end

[cropped_ID_sorted_im,idx_sorted] = sort(cropped_ID);
cropped_SI_sorted_im = cropped_SI(idx_sorted);

%% save
save('correlations_res_sus_S.mat','sus_im','res_im','allmice','cropped_ID_sorted_im','cropped_SI_sorted_im')