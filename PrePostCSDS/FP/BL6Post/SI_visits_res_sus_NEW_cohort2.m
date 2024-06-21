function SI_visits_res_sus_NEW_cohort2(csds_cb,SW_CSDS_Post,mouse_CSDS,SW_CTRL_Post,CSDS_ID,CSDS_SI)

%% Find where to split
avg = nanmean(SW_CTRL_Post);
std_val = nanstd(SW_CTRL_Post);
med = avg - std_val;

SI_CSDS = SW_CSDS_Post;

resID = mouse_CSDS(find(SI_CSDS>=med));
susID = mouse_CSDS(find(SI_CSDS<med));

res_cb = struct;
sus_cb = struct;

%% Split mice into res and sus
% res
r = 1;
for i = 1:size(csds_cb.mouse,1)
    if ismember(csds_cb.mouse(i),resID)
        res_cb.mouse(r,1) = csds_cb.mouse(i);
        res_cb.visits_plot_avg(r,:) = csds_cb.visits_plot_avg(i,:);
        res_cb.num_visits(r,:) = csds_cb.num_visits(i,:);
        r = r + 1;
    end
end

% sus
s = 1;
for i = 1:size(csds_cb.mouse,1)
    if ismember(csds_cb.mouse(i),susID)
        sus_cb.mouse(s,1) = csds_cb.mouse(i);
        sus_cb.visits_plot_avg(s,:) = csds_cb.visits_plot_avg(i,:);
        sus_cb.num_visits(s,:) = csds_cb.num_visits(i,:);
        s = s + 1;
    end
end

%% Find avg
% res
res_cb.visits_plot_avg_overall = nanmean(res_cb.visits_plot_avg);
res_cb.visits_plot_std_overall = nanstd(res_cb.visits_plot_avg)/sqrt(size(res_cb.visits_plot_avg,1));
res_cb.num_visits_avg_overall = nanmean(res_cb.num_visits);
res_cb.num_visits_std_overall = nanstd(res_cb.num_visits)/sqrt(size(res_cb.num_visits,1));

% sus
sus_cb.visits_plot_avg_overall = nanmean(sus_cb.visits_plot_avg);
sus_cb.visits_plot_std_overall = nanstd(sus_cb.visits_plot_avg)/sqrt(size(sus_cb.visits_plot_avg,1));
sus_cb.num_visits_avg_overall = nanmean(sus_cb.num_visits);
sus_cb.num_visits_std_overall = nanstd(sus_cb.num_visits)/sqrt(size(sus_cb.num_visits,1));

%% Sort
% Crop out mice npt analyzed (bad signal, gfp, etc)
% load SI data and make sure the IDs and SIs are in the correct order!

CSDS_SI = CSDS_SI/300*100;

% cell bodies
cropped_ID_cb = CSDS_ID(find(ismember(CSDS_ID,csds_cb.mouse)));
cropped_SI_cb = CSDS_SI(find(ismember(CSDS_ID,csds_cb.mouse)));
% make sure croppedSIs are in the same order as struct with Favg data
[~,idx_cb] = sort(csds_cb.mouse);

% Favg data and IDs should be in order, but double check this 
if sum(idx_cb == (1:length(idx_cb))') == length(idx_cb)
    disp('sorted order is correct')
else
    warning('Error in the mouse ID ordering!')
end

[cropped_ID_cb_sorted,idx_cb_sorted] = sort(cropped_ID_cb);
cropped_SI_cb_sorted = cropped_SI_cb(idx_cb_sorted);

%% save
save('SI_F_visits_res_sus_NEW_cohort2.mat','sus_cb','res_cb','csds_cb','cropped_ID_cb_sorted','cropped_SI_cb_sorted')%,'sus_term','res_term','csds_term','cropped_ID_term_sorted','cropped_SI_term_sorted')