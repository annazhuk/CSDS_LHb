function nose_sig_cells_sus_res_FIXED(allmice,SW_CSDS_Post,mouse_CSDS,SW_CTRL_Post,CSDS_ID,CSDS_SI)

%% Find where to split
avg = nanmean(SW_CTRL_Post);
std_val = nanstd(SW_CTRL_Post);
med = avg - std_val;

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
        res_im.sig_up_total(r,1) = allmice.sig_up_total(i,1);
        res_im.sig_up_percent(r,1) = allmice.sig_up_percent(i,1);
        res_im.sig_down_total(r,1) = allmice.sig_down_total(i,:);
        res_im.sig_down_percent(r,1) = allmice.sig_down_percent(i,:);
        res_im.sig_up_all{r,1} = allmice.sig_up_all{i,1};
        res_im.sig_down_all{r,1} = allmice.sig_down_all{i,1};
        %res_im.peak_max_sig(r,1) = allmice.peak_max_sig(i,1);
        res_im.pk_real_sig(r,1) = allmice.pk_real_sig(i,1);
        %res_im.peak_min_sig(r,1) = allmice.peak_min_sig(i,1);
        %res_im.peak_max_sig_all{r,1} = allmice.peak_max_sig_all{i,1};
        res_im.pk_real_sig_all{r,1} = allmice.pk_real_sig_all{i,1};
        %res_im.peak_min_sig_all{r,1} = allmice.peak_min_sig_all{i,1};
        r = r + 1;
    end
end

% sus
s = 1;
for i = 1:size(allmice.mice,1)
    if ismember(allmice.mice(i),susID)
        sus_im.mice(s,1) = allmice.mice(i);
        sus_im.sig_up_total(s,1) = allmice.sig_up_total(i,1);
        sus_im.sig_up_percent(s,1) = allmice.sig_up_percent(i,1);
        sus_im.sig_down_total(s,1) = allmice.sig_down_total(i,:);
        sus_im.sig_down_percent(s,1) = allmice.sig_down_percent(i,:);
        sus_im.sig_up_all{s,1} = allmice.sig_up_all{i,1};
        sus_im.sig_down_all{s,1} = allmice.sig_down_all{i,1};
        %sus_im.peak_max_sig(s,1) = allmice.peak_max_sig(i,1);
        sus_im.pk_real_sig(s,1) = allmice.pk_real_sig(i,1);
        %sus_im.peak_min_sig(s,1) = allmice.peak_min_sig(i,1);
        %sus_im.peak_max_sig_all{s,1} = allmice.peak_max_sig_all{i,1};
        sus_im.pk_real_sig_all{s,1} = allmice.pk_real_sig_all{i,1};
        %sus_im.peak_min_sig_all{s,1} = allmice.peak_min_sig_all{i,1};
        s = s + 1;
    end
end

%% Put all sig neurons for sus and res into a single vector
% res
all_neurons_up = [];
all_neurons_down = [];
for m = 1:size(res_im.mice,1)
    if ~isnan(res_im.sig_up_all{m})
        all_neurons_up = [all_neurons_up;res_im.pk_real_sig_all{m}(res_im.sig_up_all{m})];
        all_neurons_down = [all_neurons_down;res_im.pk_real_sig_all{m}(res_im.sig_down_all{m})];
    end
end
res_im.peak_max_sig_all_collected = all_neurons_up;
res_im.peak_min_sig_all_collected = all_neurons_down;

% sus
all_neurons_up = [];
all_neurons_down = [];
for m = 1:size(sus_im.mice,1)
    if ~isnan(sus_im.sig_up_all{m})
        all_neurons_up = [all_neurons_up;sus_im.pk_real_sig_all{m}(sus_im.sig_up_all{m})];
        all_neurons_down = [all_neurons_down;sus_im.pk_real_sig_all{m}(sus_im.sig_down_all{m})];
    end
end
sus_im.peak_max_sig_all_collected = all_neurons_up;
sus_im.peak_min_sig_all_collected = all_neurons_down;

%% Calculate percentages out of all cells
% res
res_cell_total = 0;
res_cell_up_sig = 0;
res_cell_down_sig = 0;
for m = 1:size(res_im.mice,1)
    res_cell_total = res_cell_total + size(res_im.sig_up_all{m},1);
    res_cell_up_sig = res_cell_up_sig + nansum(res_im.sig_up_all{m});
    res_cell_down_sig = res_cell_down_sig + nansum(res_im.sig_down_all{m});
end
res_im.sig_up_percent_allcells = res_cell_up_sig/res_cell_total*100;
res_im.sig_down_percent_allcells = res_cell_down_sig/res_cell_total*100;

% sus
sus_cell_total = 0;
sus_cell_up_sig = 0;
sus_cell_down_sig = 0;
for m = 1:size(sus_im.mice,1)
    sus_cell_total = sus_cell_total + size(sus_im.sig_up_all{m},1);
    sus_cell_up_sig = sus_cell_up_sig + nansum(sus_im.sig_up_all{m});
    sus_cell_down_sig = sus_cell_down_sig + nansum(sus_im.sig_down_all{m});
end
sus_im.sig_up_percent_allcells = sus_cell_up_sig/sus_cell_total*100;
sus_im.sig_down_percent_allcells = sus_cell_down_sig/sus_cell_total*100;

%% Find avg
% cell bodies
% res
% avg
res_im.sig_up_total_avg = nanmean(res_im.sig_up_total);
res_im.sig_up_percent_avg = nanmean(res_im.sig_up_percent);
res_im.sig_down_total_avg = nanmean(res_im.sig_down_total);
res_im.sig_down_percent_avg = nanmean(res_im.sig_down_percent);
%res_im.peak_max_sig_avg = nanmean(res_im.pk_real_sig);
res_im.pk_real_sig_avg = nanmean(res_im.pk_real_sig);
%res_im.peak_min_sig_avg = nanmean(res_im.peak_min_sig);

% std
res_im.sig_up_total_std = nanstd(res_im.sig_up_total)/sqrt(size(res_im.sig_up_total,1));
res_im.sig_up_percent_std = nanstd(res_im.sig_up_percent)/sqrt(size(res_im.sig_up_percent,1));
res_im.sig_down_total_std = nanstd(res_im.sig_down_total)/sqrt(size(res_im.sig_down_total,1));
res_im.sig_up_total_std = nanstd(res_im.sig_up_total)/sqrt(size(res_im.sig_down_percent,1));
%res_im.peak_max_sig_std = nanstd(res_im.pk_real_sig)/sqrt(size(res_im.pk_real_sig,1));
res_im.pk_real_sig_std = nanstd(res_im.pk_real_sig)/sqrt(size(res_im.pk_real_sig,1));
%res_im.peak_min_sig_std = nanstd(res_im.peak_min_sig)/sqrt(size(res_im.peak_min_sig,1));

% sus
% avg
sus_im.sig_up_total_avg = nanmean(sus_im.sig_up_total);
sus_im.sig_up_percent_avg = nanmean(sus_im.sig_up_percent);
sus_im.sig_down_total_avg = nanmean(sus_im.sig_down_total);
sus_im.sig_down_percent_avg = nanmean(sus_im.sig_down_percent);
%sus_im.peak_max_sig_avg = nanmean(sus_im.pk_real_sig);
sus_im.pk_real_sig_avg = nanmean(sus_im.pk_real_sig);
%sus_im.peak_min_sig_avg = nanmean(sus_im.peak_min_sig);

% std
sus_im.sig_up_total_std = nanstd(sus_im.sig_up_total)/sqrt(size(sus_im.sig_up_total,1));
sus_im.sig_up_percent_std = nanstd(sus_im.sig_up_percent)/sqrt(size(sus_im.sig_up_percent,1));
sus_im.sig_down_total_std = nanstd(sus_im.sig_down_total)/sqrt(size(sus_im.sig_down_total,1));
sus_im.sig_up_total_std = nanstd(sus_im.sig_up_total)/sqrt(size(sus_im.sig_down_percent,1));
%sus_im.peak_max_sig_std = nanstd(sus_im.pk_real_sig)/sqrt(size(sus_im.pk_real_sig,1));
sus_im.pk_real_sig_std = nanstd(sus_im.pk_real_sig)/sqrt(size(sus_im.pk_real_sig,1));
%sus_im.peak_min_sig_std = nanstd(sus_im.peak_min_sig)/sqrt(size(sus_im.peak_min_sig,1));

%% Sort
% Crop out mice npt analyzed (bad signal, gfp, etc)
% load SI data and make sure the IDs and SIs are in the correct order!

%CSDS_SI = CSDS_SI/300*100;

% cell bodies
cropped_ID_im = CSDS_ID(find(ismember(CSDS_ID,allmice.mice)));
cropped_SI_im = CSDS_SI(find(ismember(CSDS_ID,allmice.mice)));
% make sure croppedSIs are in the same order as struct with Favg data
[~,idx_im] = sort(allmice.mice);

% Favg data and IDs should be in order, but double check this 
if sum(idx_im == (1:length(idx_im))') == length(idx_im)
    disp('sorted order is correct')
else
    warning('Error in the mouse ID ordering!')
end

[cropped_ID_im_sorted,idx_im_sorted] = sort(cropped_ID_im);
cropped_SI_im_sorted = cropped_SI_im(idx_im_sorted);

%% save
save('nose_sig_cells_res_sus_FIXED.mat','sus_im','res_im','allmice','cropped_ID_im_sorted','cropped_SI_im_sorted')