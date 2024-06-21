function SI_visits_res_sus_trials_cohort1(csds,SW_CSDS_Post,mouse_CSDS,SW_CTRL_Post,CSDS_ID,CSDS_SI)

%% Find where to split
avg = nanmean(SW_CTRL_Post);
std_val = nanstd(SW_CTRL_Post);
med = avg - std_val;

%med = median(SW_CSDS_Post);

SI_CSDS = SW_CSDS_Post;

resID = mouse_CSDS(find(SI_CSDS>=med));
susID = mouse_CSDS(find(SI_CSDS<med));

res = struct;
sus = struct;
% res_term = struct;
% sus_term = struct;

%% Split mice into res and sus
% cell bodies
% res
r = 1;
for i = 1:size(csds.mouse,1)
    if ismember(csds.mouse(i),resID)
        res.mouse(r,1) = csds.mouse(i);
        res.visits_plot_avg(r,:) = csds.visits_plot_avg(i,:);
        res.visits_plot{r,1} = csds.visits_plot{i,:};
        r = r + 1;
    end
end

% sus
s = 1;
for i = 1:size(csds.mouse,1)
    if ismember(csds.mouse(i),susID)
        sus.mouse(s,1) = csds.mouse(i);
        sus.visits_plot_avg(s,:) = csds.visits_plot_avg(i,:);
        sus.visits_plot{s,1} = csds.visits_plot{i,:};
        s = s + 1;
    end
end

%% Find avg
% cell bodies
% res
res.visits_plot_avg_overall = nanmean(res.visits_plot_avg);
res.visits_plot_std_overall = nanstd(res.visits_plot_avg)/sqrt(size(res.visits_plot_avg,1));

% sus
sus.visits_plot_avg_overall = nanmean(sus.visits_plot_avg);
sus.visits_plot_std_overall = nanstd(sus.visits_plot_avg)/sqrt(size(sus.visits_plot_avg,1));


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

% % terminals
% cropped_ID_term = CSDS_ID(find(ismember(CSDS_ID,csds_term.mouse)));
% cropped_SI_term = CSDS_SI(find(ismember(CSDS_ID,csds_term.mouse)));
% % make sure croppedSIs are in the same order as struct with Favg data
% [~,idx_term] = sort(csds_term.mouse);
% 
% % Favg data and IDs should be in order, but double check this 
% if sum(idx_term == (1:length(idx_term))') == length(idx_term)
%     disp('sorted order is correct')
% else
%     warning('Error in the mouse ID ordering!')
% end
% 
% [cropped_ID_term_sorted,idx_term_sorted] = sort(cropped_ID_term);
% cropped_SI_term_sorted = cropped_SI_term(idx_term_sorted);

%% save
save('SI_F_visits_res_sus_trials_cohort1.mat','sus','res','csds','cropped_ID_sorted','cropped_SI_sorted')%,'sus_term','res_term','csds_term','cropped_ID_term_sorted','cropped_SI_term_sorted')