function SI_visits_res_sus_NEW_cohort1(csds,SW_CSDS_Post,mouse_CSDS,SW_CTRL_Post,CSDS_ID,CSDS_SI)

%% Find where to split
avg = nanmean(SW_CTRL_Post);
std_val = nanstd(SW_CTRL_Post);
med = avg - std_val;

SI_CSDS = SW_CSDS_Post;

resID = mouse_CSDS(find(SI_CSDS>=med));
susID = mouse_CSDS(find(SI_CSDS<med));

res = struct;
sus = struct;

%% Split mice into res and sus
% res
r = 1;
for i = 1:size(csds.mouse,1)
    if ismember(csds.mouse(i),resID)
        res.mouse(r,1) = csds.mouse(i);
        res.visits_plot_avg(r,:) = csds.visits_plot_avg(i,:);
        res.num_visits(r,:) = csds.num_visits(i,:);
        r = r + 1;
    end
end

% sus
s = 1;
for i = 1:size(csds.mouse,1)
    if ismember(csds.mouse(i),susID)
        sus.mouse(s,1) = csds.mouse(i);
        sus.visits_plot_avg(s,:) = csds.visits_plot_avg(i,:);
        sus.num_visits(s,:) = csds.num_visits(i,:);
        s = s + 1;
    end
end

%% Find avg
% res
res.visits_plot_avg_overall = nanmean(res.visits_plot_avg);
res.visits_plot_std_overall = nanstd(res.visits_plot_avg)/sqrt(size(res.visits_plot_avg,1));
res.num_visits_avg_overall = nanmean(res.num_visits);
res.num_visits_std_overall = nanstd(res.num_visits)/sqrt(size(res.num_visits,1));

% sus
sus.visits_plot_avg_overall = nanmean(sus.visits_plot_avg);
sus.visits_plot_std_overall = nanstd(sus.visits_plot_avg)/sqrt(size(sus.visits_plot_avg,1));
sus.num_visits_avg_overall = nanmean(sus.num_visits);
sus.num_visits_std_overall = nanstd(sus.num_visits)/sqrt(size(sus.num_visits,1));


%% Sort
% Crop out mice npt analyzed (bad signal, gfp, etc)
% load SI data and make sure the IDs and SIs are in the correct order!

CSDS_SI = CSDS_SI/300*100;

% cell bodies
cropped_ID = CSDS_ID(find(ismember(CSDS_ID,csds.mouse)));
cropped_SI = CSDS_SI(find(ismember(CSDS_ID,csds.mouse)));
% make sure croppedSIs are in the same order as struct with Favg data
[~,idx] = sort(csds.mouse);

% Favg data and IDs should be in order, but double check this 
if sum(idx == (1:length(idx))') == length(idx)
    disp('sorted order is correct')
else
    warning('Error in the mouse ID ordering!')
end

[cropped_ID_sorted,idx_sorted] = sort(cropped_ID);
cropped_SI_sorted = cropped_SI(idx_sorted);

%% save
save('SI_F_visits_res_sus_NEW_cohort1.mat','sus','res','csds','cropped_ID_sorted','cropped_SI_sorted')%,'sus_term','res_term','csds_term','cropped_ID_term_sorted','cropped_SI_term_sorted')