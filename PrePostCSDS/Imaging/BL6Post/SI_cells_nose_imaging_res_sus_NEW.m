function SI_cells_nose_imaging_res_sus_NEW(BL6Post,SW_CSDS_Post,mouse_CSDS,SW_CTRL_Post,CSDS_ID,CSDS_SI)

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
for i = 1:size(BL6Post.mice,1)
    if ismember(BL6Post.mice(i),resID)
        res_im.mice(r,1) = BL6Post.mice(i);
        res_im.nose_plot_start{r,1} = BL6Post.nose_plot_start{i};
        res_im.nose_plot_stop{r,1} = BL6Post.nose_plot_stop{i};
        res_im.nose_plot_start_avg{r,1} = BL6Post.nose_plot_start_avg{i};
        res_im.nose_plot_stop_avg{r,1} = BL6Post.nose_plot_stop_avg{i};
        res_im.nose_plot_start_avg_FP(r,:) = BL6Post.nose_plot_start_avg_FP(i,:);
        res_im.nose_plot_stop_avg_FP(r,:) = BL6Post.nose_plot_stop_avg_FP(i,:);
        r = r + 1;
    end
end

% sus
s = 1;
for i = 1:size(BL6Post.mice,1)
    if ismember(BL6Post.mice(i),susID)
        sus_im.mice(s,1) = BL6Post.mice(i);
        sus_im.nose_plot_start{s,1} = BL6Post.nose_plot_start{i};
        sus_im.nose_plot_stop{s,1} = BL6Post.nose_plot_stop{i};
        sus_im.nose_plot_start_avg{s,1} = BL6Post.nose_plot_start_avg{i};
        sus_im.nose_plot_stop_avg{s,1} = BL6Post.nose_plot_stop_avg{i};
        sus_im.nose_plot_start_avg_FP(s,:)  = BL6Post.nose_plot_start_avg_FP(i,:);
        sus_im.nose_plot_stop_avg_FP(s,:)  = BL6Post.nose_plot_stop_avg_FP(i,:);
        s = s + 1;
    end
end

%% Find avg
% cell bodies
% % res
% res_im.visits_plot_avg_overall = nanmean(res_im.visits_plot_avg);
% res_im.visits_plot_std_overall = nanstd(res_im.visits_plot_avg)/sqrt(size(res_im.visits_plot_avg,1));
% res_im.sorted_center_avg_overall = nanmean(res_im.sorted_center_avg);
% res_im.sorted_center_std_overall = nanstd(res_im.sorted_center_avg)/sqrt(size(res_im.sorted_center_avg,1));
% 
% % sus
% sus_im.visits_plot_avg_overall = nanmean(sus_im.visits_plot_avg);
% sus_im.visits_plot_std_overall = nanstd(sus_im.visits_plot_avg)/sqrt(size(sus_im.visits_plot_avg,1));
% sus_im.sorted_center_avg_overall = nanmean(sus_im.sorted_center_avg);
% sus_im.sorted_center_std_overall = nanstd(sus_im.sorted_center_avg)/sqrt(size(sus_im.sorted_center_avg,1));

%% Sort
% Crop out mice npt analyzed (bad signal, gfp, etc)
% load SI data and make sure the IDs and SIs are in the correct order!

%CSDS_SI = CSDS_SI/300*100;

% cell bodies
cropped_ID_im = CSDS_ID(find(ismember(CSDS_ID,BL6Post.mice)));
cropped_SI_im = CSDS_SI(find(ismember(CSDS_ID,BL6Post.mice)));
% make sure croppedSIs are in the same order as struct with Favg data
[~,idx_im] = sort(BL6Post.mice);

% Favg data and IDs should be in order, but double check this 
if sum(idx_im == (1:length(idx_im))') == length(idx_im)
    disp('sorted order is correct')
else
    warning('Error in the mouse ID ordering!')
end

[cropped_ID_im_sorted,idx_im_sorted] = sort(cropped_ID_im);
cropped_SI_im_sorted = cropped_SI_im(idx_im_sorted);

%% save
save('SI_cells_imaging_visits_nose_res_sus_NEW.mat','sus_im','res_im','BL6Post','cropped_ID_im_sorted','cropped_SI_im_sorted')