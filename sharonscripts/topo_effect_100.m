%% This script is aimed for inspecting the spatial features of the effects.
%load the subjects, reduce frequency for a single chan
cd('C:\tvns_pilot\short_pulses_experiment')

folder='C:\tvns_pilot\short_pulses_experiment\eeg_all_chan_algo100\';
spect_ft_content = dir([folder '*_spect_ft_4parafac_algo.mat']);
model_folder = 'C:\tvns_pilot\short_pulses_experiment\eeg_all_chan_algo100\models\';
model_content=dir([model_folder '*.mat']);
temp= struct2cell(model_content);
allnames=temp(1,:);
allnames3=cellfun(@(x) x(1:3),allnames,'un',0);
whichfac_alpha= [3,2,1,4,1,4,3,3,2,2,1,2,3,2,1,2,1,2,3,4,1];
t_selection=10:83;
C=cell(length(model_content),2);
alpha_prctile = 1;

load('layout_193.mat');
all_labels = layout193.label(1:(end-2));

sens = ft_read_sens('C:\tvns_pilot\EGI_256.sfp'); %reading electrode general, relevant if we use the same coordinates for all subjects
sens.label{260}='VREF'; % because the data is recorded with Cz as VREF 

cfg=[];
cfg.channel=all_labels;
sens_ex=ft_selectdata(cfg,sens);

cfg=[];
cfg.method='distance';
cfg.neighbourdist= 4.5;
cfg.elec=sens_ex;
neighbours=ft_prepare_neighbours(cfg);
%% first let's look on the average spatial components themsleves:
% loading the models to look on the channel domain of the alpha-comp:
model_cell = cell(size(21,2),2);
% content_model = dir('C:\tvns_pilot\short_pulses_experiment\eeg_all_chan_algo100\models\*.mat');
content_model = dir([folder 'models\*.mat']);
content_w_trf = dir([folder 'weigted_tfr_percent\*.mat']);
sub_label= cell(length(content_w_trf),2);

for sub=1:length(content_model)
%     temp  = load([folder 'models\' content_model(sub).name]);
    temp  = load(['C:\tvns_pilot\short_pulses_experiment\eeg_all_chan_algo100\models\' content_model(sub).name]);
    temp1=temp.model;    
    temp_w = load([folder 'weigted_tfr_percent\' content_w_trf(sub).name]);
    for cond=1:2
    model_cell{sub,cond}=temp1{cond}; 
    sub_label{sub,cond} = temp_w.w_spect_ft{cond}.cfg.previous.channel;
    end
end
% small check that channels are the same for both conditions:
just_check= cellfun(@(x,y) setxor(x,y), sub_label(:,1),sub_label(:,2),'un',0);
% move the factor from cell to matrix: 
whichfac_alpha_cell=num2cell([whichfac_alpha'  whichfac_alpha']);
chan_model =cellfun(@(x,y) x{1}(:,y),model_cell,whichfac_alpha_cell,'un',0);
% problem cause each subject has other missing channels :(

% let's count the label in each:
[missing_labels xorA xorB]= cellfun(@(x) setxor(x,all_labels), sub_label,'un',0);
[intersect_labels interA interB]= cellfun(@(x) intersect(x,all_labels), sub_label,'un',0);

%% best is to add NaNs for the missing:
model_chan_mat= zeros(length(all_labels),size(sub_label,1));
for sub=1:size(sub_label,1)
model_chan_mat(interB{sub,1},sub) = chan_model{sub,cond}(interA{sub,cond}); %place intersection;
model_chan_mat(xorB{sub,1},sub) = NaN; %place intersection;
end
imagesc(isnan(model_chan_mat)) % let's see who is missing.

% median the data
med_chan_model = nanmedian(model_chan_mat,2);
nan_count_model =sum(isnan(model_chan_mat),2);
plot(med_chan_model)

%% plot the median component on topo:
load('C:\tvns_pilot\skelton.mat')
scale_map= sqrt(1:-0.01:0);
% red_map= [ones(1,length(scale_map)); scale_map; scale_map]';
pink_map  = [ones(1,length(scale_map)); scale_map; ones(1,length(scale_map))]';

skelton.label = all_labels;
skelton.powspctrm=med_chan_model'; % spatial
highlight_chan= skelton.label(med_chan_model'>390  & 1:192>80);

cfg=[];
cfg.colormap= pink_map;
cfg.comment='no';
cfg.layout=layout193;
cfg.marker = 'yes';
cfg.highlight          = 'on';
cfg.highlightchannel   =  highlight_chan;
cfg.highlightsymbol    = '.';
cfg.highlightcolor  = 'b';
cfg.highlightsize      = 20;

ft_topoplotER(cfg,skelton)
 % THE MEDIAN MODEL LOOKS REALY GOOD AND VISUAL, WITH NOT A LOT OF MISSING
 % CHANELS IN RELEVANT AREAS.
%% get the reduced time dimesion matrices, but leave the channel UN-reduced:
for sub=1:length(spect_ft_content) 
    load([folder spect_ft_content(sub).name])
    load([model_folder model_content(sub).name])
    comp=whichfac_alpha(sub); % select the alpha componenet
    for cond=1:2
        % reduce the freq domain first:
        freq_model =model{cond}{2};
        V_freq = freq_model(:,comp);
        X=permute(spect_ft{cond}.powspctrm,[1 2 4 3]); % move the freq domain to the end
        C_pluschan{sub,cond} = sum(bsxfun(@times,X,reshape(V_freq,1,1,1,numel(V_freq))),4);               
    end
end
time_axis = spect_ft{cond}.time; % extract time_axis for later
%% select trials according to order_ind (from "fix_weigting_select_trials")
load('C:\tvns_pilot\short_pulses_experiment\eeg_all_chan_algoclean\order_ind_1_1sec.mat')

% load('C:\tvns_pilot\short_pulses_experiment\order_ind.mat')
C_sorted =cellfun(@(x,y) x(y,:,:),C_pluschan,order_ind,'un',0);
trl_amount=cellfun(@(x) floor(size(x,1)*alpha_prctile),C_sorted,'un',0);
c_mean_up = cellfun(@(x,y) squeeze(mean(x(1:y,:,:),1)),C_sorted,trl_amount,'un',0);

% project the subjects into a shared matrix:
%%best is to add NaNs for the missing:
time_frame= 91;
data_chan_mat= zeros(size(sub_label,1),2,length(all_labels),time_frame);

for sub=1:size(sub_label,1)
    for cond=1:2
    data_chan_mat(sub,cond,interB{sub,1},:) = c_mean_up{sub,cond}(interA{sub,cond},:,:); %place intersection;
    data_chan_mat(sub,cond,xorB{sub,1},:) = NaN; %place intersection;
    end
end

%% baseline: -1 0: each subject
baseline_lim =[-1 0];
baseline_ind = time_axis>baseline_lim(1) & time_axis<baseline_lim(2);
P=permute(data_chan_mat,[1 2 4 3]);

baseline_val= nanmean(P(:,:,baseline_ind,:),3);
A=bsxfun(@(x,y) (minus(x,y)),P,baseline_val);
% A=bsxfun(@(x,y) (rdivide(x,y)*100),P,baseline_val);

Z=permute(A,[1 2 4 3]); 
%% Mean across subjects:
avg_topo_frame = squeeze(nanmean(Z,1)); % mean
F=permute(avg_topo_frame,[1 2 4 3]);

% create a multi-sub skelton for statistics:
for sub=1:21
    F_multisub{sub}= permute(Z(sub,:,:,:),[2 3 1 4]);
end
%% plot
for cond =1:2
    subplot(2,1,cond)
% imagesc(time_axis,1:192,squeeze(avg_topo_frame(cond,:,:)),[80 120])
imagesc(time_axis,1:192,squeeze(avg_topo_frame(cond,:,:)),[ -4 4])

xlim([-13 13])
% clim([ 0.8 1.2])
r=rectangle('Position',[ 0 1 3.4 192],'LineWidth',0.2,'LineStyle','--');
end
%% preare average skelton
load(['C:\tvns_pilot\short_pulses_experiment\eeg_all_chan_algoclean\' 'ES8_spect_ft_4parafac_algo.mat']) % load a random subject for skelton.

cfg=[];
cfg.avgoverfreq = 'yes';
cfg.avgfovertime = 'no';
cfg.avgoverrpt = 'yes';
skelton =cellfun(@(x) ft_selectdata(cfg,x),spect_ft,'un',0);
for cond=1:2
    skelton{cond}.previous= [];
    skelton{cond}.label =all_labels;
    skelton{cond}.powspctrm = permute(F(cond,:,:,:),[2 3 4 1]);
end

%% plotting topgographies
cntr=1;
tit={'tv','sh'};
% topo1=figure('units','normalized','outerposition',[0 0 1 1]);
% plustime= 3.4;
plustime = 1;
for cond=1:2
    for i=-1:6
%         topo= figure;
    cfg=[];
    cfg.comment= 'no';
    cfg.xlim = [i-plustime/2 i+plustime/2];
    cfg.zlim= [-4.5 4.5]; % absoulte values MEAN
%     cfg.zlim= [70 110]; % % values
    cfg.layout = layout193;
%     if ~isempty(highlight_chan{i+2})
    cfg.highlight          = 'on';
    cfg.highlightchannel   =  highlight_chan{i+2};
    cfg.highlightsymbol    = '.';
    cfg.highlightcolor  = 'y';  
    cfg.highlightsize  = 6;  
%     end
        subplot(2,8,cntr)
    set(gca,'fontsize',20)
    ft_topoplotTFR(cfg,skelton{cond})
     title([num2str(i) 's'])
%             h=colorbar;
%      ylabel(h, '\muV^2','fontSize',20,'fontName','Arial')

    hold on
    cntr=cntr+1;
%     saveas(topo,['C:\tvns_pilot\paper figures\topoNew\' 'cond=' num2str(cond) '_t=' num2str(i) '.png'])
    
    end
end

%% prepare a stat skelton:
load('C:\tvns_pilot\skelton.mat')
time_index=find(time_axis<4 & time_axis>0);

time_1_6=42:65;
% cfg=[];
% cfg.avgoverfreq = 'yes';
% cfg.avgfovertime = 'no';
% cfg.avgoverrpt = 'yes';
% skelton =cellfun(@(x) ft_selectdata(cfg,x),spect_ft,'un',0);
for cond=1:2
    for sub=1:21
    stat_skelton{sub,cond}=skelton;
%     stat_skelton{sub,cond}.time=time_axis(time_1_6);
    stat_skelton{sub,cond}.time=[-1:6];

    stat_skelton{sub,cond}.previous=[];
    stat_skelton{sub,cond}.label =all_labels;
%     stat_skelton{sub,cond}.powspctrm =permute(nanmean(F_multisub{sub}(cond,:,:,time_index),4),[2 3 4 1]);
    insert=permute(F_multisub{sub}(cond,:,:,time_1_6),[2 3 4 1]);
    %average on 3 values to get 8 values:
    insert_8=permute(mean(reshape(insert,192,1,3,8),3),[1 2 4 3]);
    stat_skelton{sub,cond}.powspctrm = insert_8;
    end
end

%% for topographical statistics:
% design=[ones(1,21);2*ones(1,21)];
% divide channels to 10 groups (19 channel in each):
stat={};
cfg=[];
% cfg.latency =[0 3.4];
cfg.correctm='cluster';
cfg.method ='montecarlo';
cfg.neighbours=neighbours;
cfg.avgovertime = 'no';
cfg.tail=0;
cfg.clusteralpha=0.05;
cfg.statistic = 'depsamplesT';
cfg.frequency = 'all';
cfg.numrandomization = 10000;
fg.clusterstatistic = 'maxsum';
% cfg.minnbchan=val;
% % cfg.clusterthreshold= 'nonparametric_individual';
% cfg.clustercritval = 4;
cfg.clustertail=0;
cfg.ivar = 1;
cfg.uvar = 2;
% cfg.clustercritval=val;
cfg.design = [ones(length(stat_skelton),1) (1:21)';
              2*ones(length(stat_skelton),1) (1:21)'];

stat= ft_freqstatistics(cfg,stat_skelton{:});

%% plot cluster:
load('C:\tvns_pilot\skelton.mat')
skelton.label = all_labels;
skelton.powspctrm=zeros(1,192); % spatial
sel_time=[-1:6];
cfg=[];
% cfg.colormap= pink_map;
cfg.comment='no';
cfg.layout=layout193;
cfg.marker = 'yes';
cfg.highlight          = 'on';
highlight_chan=cell(1,8);
for t=1:8
    subplot(1,8,t)
    title(num2str(sel_time(t)))
highlight_chan{t}= all_labels(stat.negclusterslabelmat(:,1,t)==1);
cfg.highlightchannel   =  highlight_chan{t};
cfg.highlightsymbol    = '.';
cfg.highlightcolor  = 'w';
cfg.highlightsize      = 5;
ft_topoplotER(cfg,skelton)
end

%% plot stuff
for val=1:6
subplot(6,1,val)
imagesc(stat{val}.time,1:192,squeeze(stat{val}.ref))
end
% colorbar
A=extractfield(stat.negclusters,'prob');

%%
for val=1
cfg = [];
cfg.alpha  = 0.2;
% cfg.parameter = 'chan_freq_time';
% cfg.zlim   = [-1e-27 1e-27];
cfg.layout = layout193;
% cfg.subplotsize=[1,1];
% cfg.highlightcolorneg=[1 0 0];
ft_clusterplot(cfg, stat);
end
%% let's do a priori comparision - accroding the 27 electrodes hilghted in the break comp topo:
time_index=find(time_axis<4 & time_axis>0);
[~,~,highlight_index]=intersect(highlight_chan, layout193.label);

area_values_cell=cellfun(@(x) squeeze(nanmean(nanmean(x(:,highlight_index,1,time_index),2),4)),F_multisub,'un',0);
area_values_mat=cell2mat(area_values)';

grand_avg=mean(area_values_mat);
real_result_value=grand_avg(2)-grand_avg(1);
grand_sem=std(area_values_mat)./sqrt(21);

[p,h]=signrank(area_values_mat(:,1),area_values_mat(:,2));

% this takes 192 continous areas in the same size:
cfg=[];
cfg.method='distance';
% cfg.neighbourdist= 6.6; % for 27
cfg.neighbourdist=4.5;
cfg.elec=sens_ex;
neighbours_wide1=ft_prepare_neighbours(cfg);

grand_avg_i=nan(192,1);
area_values_mat_contin=nan(21,2,192);
for i=1:length(neighbours_wide1)
    [~,~,neigh_index]=intersect(neighbours_wide1(i).neighblabel, layout193.label);
    area_values_cell_i=cellfun(@(x) squeeze(nanmean(nanmean(x(:,neigh_index,1,time_index),2),4)),F_multisub,'un',0);
    area_values_mat_contin(:,:,i)=cell2mat(area_values_cell_i)';
    grand_avg_temp=mean(area_values_mat_contin(:,:,i));
    grand_avg_i(i)=grand_avg_temp(2)-grand_avg_temp(1);
end
hist(grand_avg_i)
% since there are 192 values lower than the real result
p=1/193;
    
pvals=[];
for i=1:size(area_values_mat_contin,3)
    [pvals(i),h(i)]=signrank(area_values_mat_contin(:,1,i),area_values_mat_contin(:,2,i));
end
[p_fdr, p_masked] = fdr(pvals, 0.05);
sum(p_masked)    
%% plot the electrodes after FDR:



%% try 13+13 
cfg=[];
cfg.method='distance';
cfg.neighbourdist= 4.6; % for 27
cfg.elec=sens_ex;
neighbours_wide2=ft_prepare_neighbours(cfg);

combos_2 = combntns(1:192,2); % all possible paring of two groups;
grand_avg_2groups=nan(length(combos_2),1);
for i=1:length(combos_2)
    [~,~,neigh_index1]=intersect(neighbours_wide1(combos_2(i,1)).neighblabel, layout193.label); % group1
    [~,~,neigh_index2]=intersect(neighbours_wide1(combos_2(i,2)).neighblabel, layout193.label); % group2
    both_groups_index=unique([neigh_index1; neigh_index2]);
    area_values_cell_i=cellfun(@(x) squeeze(nanmean(nanmean(x(:,both_groups_index,1,time_index),2),4)),F_multisub,'un',0);
    area_values_mat_contin=cell2mat(area_values_cell_i)';
    grand_avg_temp=mean(area_values_mat_contin);
    grand_avg_2groups(i)=grand_avg_temp(2)-grand_avg_temp(1);
end
hist([grand_avg_2groups ;grand_avg_i])


%% or 9+9+9
cfg=[];
cfg.method='distance';
cfg.neighbourdist= 3.9; % for 27
cfg.elec=sens_ex;
neighbours_wide3=ft_prepare_neighbours(cfg);

combos_3 = combntns(1:192,3); % all possible paring of three groups;
grand_avg_3groups=nan(length(combos_3),1);
for i=1:length(combos_3)
    [~,~,neigh_index1]=intersect(neighbours_wide1(combos_3(i,1)).neighblabel, layout193.label); % group1
    [~,~,neigh_index2]=intersect(neighbours_wide1(combos_3(i,2)).neighblabel, layout193.label); % group2
    [~,~,neigh_index3]=intersect(neighbours_wide1(combos_3(i,3)).neighblabel, layout193.label); % group2
    both_groups_index=unique([neigh_index1; neigh_index2 ; neigh_index3]);
    area_values_cell_i=cellfun(@(x) squeeze(nanmean(nanmean(x(:,both_groups_index,1,time_index),2),4)),F_multisub,'un',0);
    area_values_mat_contin=cell2mat(area_values_cell_i)';
    grand_avg_temp=mean(area_values_mat_contin);
    grand_avg_3groups(i)=grand_avg_temp(2)-grand_avg_temp(1);
end
hist([grand_avg_3groups])


%% Altenrnatively pick random 27 channels:
permuations_num=10000;
grand_random=nan(1,permuations_num);
for i=1:size(area_values_mat_random,3)
    neigh_index=randperm(192,27);
    area_values_cell_i=cellfun(@(x) squeeze(nanmean(nanmean(x(:,neigh_index,1,time_index),2),4)),F_multisub,'un',0);
    area_values_mat_random=cell2mat(area_values_cell_i)';
    grand_avg_temp=nanmean(area_values_mat_random);
    grand_random(i)=grand_avg_temp(2)-grand_avg_temp(1);
end

hist(grand_random)



pvals=[];
for i=1:size(area_values_mat_contin,3)
    [pvals(i),h(i)]=signrank(area_values_mat_contin(:,1,i),area_values_mat_contin(:,2,i));
end
[p_fdr, p_masked] = fdr(pvals, 0.05);
sum(p_masked)

