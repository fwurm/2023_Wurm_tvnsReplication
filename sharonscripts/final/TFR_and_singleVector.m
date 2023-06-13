% alpha according to the componenet
% it produces the tfr and the single vector result

folder='C:\tvns_pilot\short_pulses_experiment\eeg_all_chan_algo100\';
spect_ft_content = dir([folder '*_spect_ft_4parafac_algo_tfr7cycle_100.mat']);
model_folder = 'C:\tvns_pilot\short_pulses_experiment\eeg_all_chan_algoclean\models\';
model_content=dir([model_folder '*.mat']);
temp= struct2cell(model_content);
allnames=temp(1,:);
allnames3=cellfun(@(x) x(1:3),allnames,'un',0);
whichfac_alpha= [3,2,1,4,1,4,3,3,2,2,1,2,3,2,1,2,1,2,3,4,1];
t_selection=10:83;
gold=(1+sqrt(5))/2;

%% get the time- only data using the parrafac model:
C=cell(length(spect_ft_content),2);
for sub=1:length(spect_ft_content)
    load([folder spect_ft_content(sub).name])
    load([model_folder model_content(sub).name])
    comp=whichfac_alpha(sub); % select the alpha componenet
    for cond=1:2
        % reduce the channels domain:
        chan_model=model{cond}{1}; % these are the channels from the model i.e from the break.
        V_chan= chan_model(:,comp); % get the right comp
        Y=permute(spect_ft{cond}.powspctrm,[1 3 4 2]); % move the channels dim to the end
        B = sum(bsxfun(@times,Y,reshape(V_chan,1,1,1,numel(V_chan))),4);

        B2 =reshape(B,size(B,1),1,size(B,2),size(B,3)); % trl, ~, freq, time
        % reduce the freq domain:
        freq_model =model{cond}{2};
        V_freq = freq_model(:,comp);
        X=permute(B2,[1 2 4 3]); % move the freq domain to the end
        C{sub,cond} = squeeze(sum(bsxfun(@times,X,reshape(V_freq,1,1,1,numel(V_freq))),4));
    end
end

%% compare # of trials:
trial_num=cellfun(@(x) size(x,1),C);

mean_trl_cond=mean(trial_num,1);
[p,h]= signrank(trial_num(:,1),trial_num(:,2));
sem_trl_cond=std(trial_num,[],1)/sqrt(21);
min_trl_cond=min(trial_num);
max_trl_cond=max(trial_num);

%% divide and then use perecent
time_axis = spect_ft{1}.time;
freq_axis = spect_ft{1}.freq;

c_mean=cellfun(@(x) mean(x,1),C,'un',0);
G=cell2mat(c_mean(:,1));
G(:,:,2) = cell2mat(c_mean(:,2));

grand_G=squeeze(nanmean(G,1));
baseline_lim =[-1 0];
baseline_ind = time_axis>baseline_lim(1) & time_axis<baseline_lim(2);

plot(time_axis,grand_G);legend('tv','sham');
hold on
plot(time_axis,baseline_ind*10^17)

baseline_val= mean(G(:,baseline_ind,:),2);
grand_base_values=mean(baseline_val(:,:,2)-baseline_val(:,:,1),1);

A=bsxfun(@(x,y) (rdivide(x,y)),G,baseline_val);

% average across subjects:
grand_avg_per=squeeze(nanmean(A,1));
std_per=squeeze(nanstd(A,[],1))/sqrt(21); % this is SEM
stim_lim =[0 4];


stim_ind = time_axis>stim_lim(1) & time_axis<stim_lim(2);
stim_values = nanmean(A(:,stim_ind,:),2);
values= stim_values*100;

report_mean=mean(values,1);
report_sem=std(values,1)./sqrt(21);
%% plot effect size figure
eff_size_fig=figure('Renderer', 'painters', 'Position', [10 10 900/gold^2 900/gold]);
color='kk';
kind={'-','--'};
for sub=1:20 
    col=(stim_values(sub,1,1)-stim_values(sub,1,2)>0)+1;
    plot([1 2],[values(sub,1,1) values(sub,1,2)],kind{col},'lineWidth',3,'Color',[0 0 0]+(col-1)*0.5)
    hold on
end
set(gca,'Fontsize',18,'fontName','Arial')

xlim([0.7 2.3])
ylim([75 130])

names = {'tVNS','Sham'};

set(gca,'xtick',[1:5],'xticklabel',names)

[p, h, stats] = signrank(values(:,1,1),values(:,1,2),'tail','left');

[h, p, stats] = ttest(values(:,1,1),values(:,1,2),'tail','left');

ylabel('% power')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NAMES=struct2cell(spect_ft_content);
sub_names=cellfun(@(x) x(1:3),NAMES(1,:),'un',0);
save('single_sub_eeg_eff','sub_names','values'); % percent per subject values

%% plot GA
plot(time_axis,grand_avg_per)
xlim([-1 8])
legend('tvns','sham')

%% plot the  grand average with std
colors='rk';
% for cond=1:2
% grand_avg_per_smooth(:,cond)=lowPassFilter(2,grand_avg_per(:,cond),20)+mean(grand_avg_per(:,cond));
% end

trace_fig=figure('Renderer', 'painters', 'Position', [10 10 910 910/gold^2]);
for cond=1:2
shadedErrorBar(time_axis,grand_avg_per(:,cond)*100,std_per(:,cond)*100,colors(cond),0.2);
hold on
xlim([-1 6])
ylim([0.81 1.16]*100)
end
set(gca,'Fontsize',18,'fontName','Arial')
stim_patch=patch([0,0,3.4,3.4],[0.8,1.2,1.2,0.8]*100,'k','LineStyle','none');
alpha(stim_patch,0.2)
colorbar
ylabel('Average alpha (%)')
xlabel('Time (sec)')
set(gca,'fontSize',20,'fontName','Arial')

%% statistics
time_lim=[-1 6];
time_ind=find(time_axis>time_lim(1) & time_axis<time_lim(2));

pvals=nan(1,size(time_ind,2)); h=nan(1,size(time_ind,2));
for i=1:length(time_ind)
[pvals(i),h(i)]=signrank(A(:,time_ind(i),1),A(:,time_ind(i),2));
end
plot(time_axis(time_ind),pvals<0.01)
[p_fdr, p_masked] = fdr(pvals, 0.05);
significant_ind=find(p_masked.*h==1);
significant_time=time_axis(significant_ind)./1000;
plot(time_axis,pvals.*p_masked)
plot(time_axis,h.*p_masked,'r.')

% using sign rank works as well:
[p, h, stats] = signrank(values(:,1,1),values(:,1,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  load the full spectrum trials to show the same effect in TFR:  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('order_ind_1_1sec')
for sub=1:length(spect_ft_content)
    load([folder spect_ft_content(sub).name])
    load([model_folder model_content(sub).name])
    comp=whichfac_alpha(sub); % select the alpha componenet
    for cond=1:2
        % reduce the channels domain:
        chan_model=model{cond}{1}; % these are the channels from the model i.e from the break.
        V_chan= chan_model(:,comp); % get the right comp
        Y=permute(spect_ft{cond}.powspctrm,[1 3 4 2]); % move the channels dim to the end
        B = sum(bsxfun(@times,Y,reshape(V_chan,1,1,1,numel(V_chan))),4);

        B2 =reshape(B,size(B,1),1,size(B,2),size(B,3)); % trl, ~, freq, time
        
        cfg=[];
        cfg.channel = 'E192'; %dummy
        cfg.trials = 1;
        w_spect_ft{cond}=ft_selectdata(cfg,spect_ft{cond});
        w_spect_ft{cond}.powspctrm = mean(B2,1); % median or mean on the trials, #how many trials persub?
        
        w_spect_ft{cond}.cfg.previous= []; % reduce size of the datasets

        save([folder 'weigthed_tfr_percent\' model_content(sub).name(1:3) '_w_spect_ft_1_1'],'w_spect_ft')         
    end
end   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loading all subject data and arranging it:
load('axis.mat')
content = dir([folder 'weigthed_tfr_percent\*_1_1.mat']);
sub_cell = cell(size(content,1),2);

for sub=1:length(content)
    temp  = load([folder 'weigthed_tfr_percent\' content(sub).name]);
    temp1=temp.w_spect_ft;
    for cond=1:2
    temp1{cond}.label = {'alpha'};
    sub_cell{sub,cond}=temp1{cond};  
    end
end

%% loading the models to look on the frequency of the alpha-comp:
model_content = dir([folder 'models\*.mat']);
model_cell = cell(size(model_content,1),2);

for sub=1:length(model_content)
    temp  = load([folder 'models\' model_content(sub).name]);
    temp1=temp.model;
    for cond=1:2
    model_cell{sub,cond}=temp1{cond};  
    end
end

% move the factor from cell to matrix: 
whichfac_alpha_cell=num2cell([whichfac_alpha'  whichfac_alpha']);

freq_model =cellfun(@(x,y) x{2}(:,y),model_cell,whichfac_alpha_cell,'un',0);

freq_model_mat=cell2mat(freq_model(:,1)');
freq_model_mat(:,:,2)=cell2mat(freq_model(:,2)');

avg_freq_model=squeeze(median(freq_model_mat,2));
sem_freq_model = squeeze(std(freq_model_mat,0,2))/sqrt(size(freq_model_mat,2)); % SEM of the model

%% plot median component on freq axis
figure('Renderer', 'painters', 'Position', [10 10 300 300]);
for cond=1
shadedErrorBar(freq_axis,avg_freq_model(:,cond),sem_freq_model(:,cond),'m')
hold on
end
ylim([0 0.25])
set(gca,'fontsize',18,'fontName','Arial')
xlabel('Frequency (Hz)')
ylabel('Weight')

% plot single subject components
figure
for sub=1:21
plot(freq_axis,freq_model{sub,1})
hold on
end

%% baseline each subject:
base_window=[-1 0];
cfg=[];
cfg.baseline = base_window;
cfg.baselinetype = 'relative'; % 

sub_baseline= cellfun(@(x) ft_freqbaseline(cfg,x),sub_cell,'un',0);
% sub_baseline=sub_cell;
% smoothing on the single subject level:
smoothingWindow = ones(3,3);
smoothingWindow = smoothingWindow./numel(smoothingWindow);

mat_smoothed = cellfun(@(x) permute(conv2(squeeze(x.powspctrm),smoothingWindow,'same'),[3 4 1 2]),sub_baseline,'un',0);
sub_smoothed=sub_baseline;
for i=1:numel(sub_baseline) % inject the smoothed data
    sub_smoothed{i}.powspctrm=mat_smoothed{i};
end
%% averaging in fieldtrip
cfg=[];
cfg.keeptrials='no';
sub_cell_desc=cellfun(@(x) ft_freqdescriptives(cfg,x),sub_baseline(:,:),'un',0);
ga=cell(2,1);

for cond=1:2
    cfg=[];    
    A=sub_cell_desc(:,cond)';
    ga{cond}=ft_freqgrandaverage(cfg,A{:});
end

% doing median outside of fieldtrip
temp_out=cellfun(@(x)  x.powspctrm,sub_cell_desc,'un',0);
data_out=cell2mat(temp_out(:,1));
data_out(:,:,:,2) =cell2mat(temp_out(:,2));

%%  is this alpha diffrence significant?
cfg= [];
ind_to_inc= 1:21;
cfg.latency  = [-1 6]; % this could be done using clustering;
cfg.frequency = [5 15];
cfg.alpha  =0.05;
cfg.parameter = 'powspctrm';
cfg.method = 'montecarlo';
cfg.avgovertime = 'no';
cfg.avgoverfreq = 'no';
cfg.correctm= 'cluster';
cfg.tail = -1;
cfg.clustertail=-1;
cfg.statistic = 'depsamplesT';
cfg.numrandomization = 10000;

D=sub_cell_desc(ind_to_inc,:); %excluding subject with no pupil response
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design = [ones(1,length(D)) ones(1,length(D))*2 ;1:21 1:21];
[stat] = ft_freqstatistics(cfg, D{:});

%% inject the mask
dif_tfr.mask=stat.mask;
median_data = mean(data_out,1);
ga_median= ga;
for cond=1:2
    ga_median{cond}.powspctrm = median_data(:,:,:,cond)*100;
    ga_median{cond}.mask = stat.mask;
%     ga_median{cond}.roi = logical(roi);
end
cfg=[];
cfg.latency =[-1 6];
ga_median612 = cellfun(@(x) ft_selectdata(cfg,x),ga_median,'un',0);
%% create the difference tfr:

cfg=[];
cfg.parameter='powspctrm';
cfg.operation='x1-x2';
dif_tfr=ft_math(cfg,ga_median612{1},ga_median612{2});
dif_tfr.mask=stat.mask;

%% plotttingtfr
% sub=3;
cfg=[];

titles= {'tVNS','sham'};
cfg.zlim= [0.9 1.2]*100; % rel baseline

cfg.ylim=[5.3 14.7];
cfg.xlim= [-1 6];
cfg.maskparameter ='mask';
cfg.maskstyle= 'outline';
cfg.maskalpha= 0.6;
cfg.fontsize  = 20;
for i=1:2
    pup_fig=figure('Renderer', 'painters', 'Position', [10 10 910 910/gold^2]);
%     ft_singleplotTFR(cfg,sub_cell_desc{sub,i});
    ft_singleplotTFR(cfg,ga_median612{i});
%     ft_singleplotTFR(cfg,D_lat{sub,i});
    title(titles{i})
    r=rectangle('Position',[ 0 5.2 3.4 9.5],'LineWidth',0.5,'LineStyle','--');
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    set(gca,'fontSize',20,'fontName','Arial')
    h=colorbar;
    ylabel(h, '% power','fontSize',20,'fontName','Arial')
end
%% plot the stats clusters:
figure
imagesc(stat.time,stat.freq,squeeze(stat.negclusters)) % relchange, relative and db works. but not absoulte
colorbar
set(gca,'YDir','normal')

%% plot the difference TFR
cfg=[];
cfg.zlim= [-20 20]; % rel baseline
cfg.ylim=[5.3 14.7];
cfg.xlim= [-1 6];
cfg.maskparameter ='mask';
cfg.maskstyle= 'outline';
cfg.maskalpha= 0.6;
cfg.fontsize  = 20;
pup_fig=figure('Renderer', 'painters', 'Position', [10 10 910 910/gold^2]);
ft_singleplotTFR(cfg,dif_tfr);
r=rectangle('Position',[ 0 5.2 3.4 9.5],'LineWidth',0.5,'LineStyle','--');
xlabel('Time (sec)')
title('difference TFR: tVNS-SHAM')
ylabel('Frequency (Hz)')
set(gca,'fontSize',20,'fontName','Arial')
h=colorbar;
ylabel(h, '% power','fontSize',20,'fontName','Arial')
%% a-priori  statistics:
%%getting single subject values in the cluster:
cfg=[];
cfg.latency  = [-1 8];
D_lat=cellfun(@(x) ft_selectdata(cfg,x),D,'un',0);

flatten=@(x) x(:);
roi= stat.mask; % using the cluster
roi_apriori_time = double(stat.time>0 & stat.time<4);
roi_apriori_freq = double(stat.freq>8 & stat.freq<12)';
roi_apriori=roi_apriori_freq*roi_apriori_time;
roi_apriori=shiftdim(roi_apriori,-1);
roi=roi_apriori;
% roi = ga_median612{2}.mymask;
cluster_value_persub=cellfun(@(x) nansum(flatten(x.powspctrm.*roi)./sum(roi(:))),D_lat);

avg=mean(cluster_value_persub*100)

sem_cluster=std(cluster_value_persub*100)./sqrt(21);

[m]=bootstrp(1000,@(x) median(x),cluster_value_persub(:,2));

[fi,xi] = ksdensity(m);
     plot(xi,fi);
     
     up=prctile(m,95);
     down=prctile(m,5);
%% Effect size figure for the cluster - defined by mask or a-priori above
colors='rk';

eff_size_fig=figure('Renderer', 'painters', 'Position', [10 10 900/gold^2 900/gold]);
for sub=1:21 
    if sub==3
    text(0.9,cluster_value_persub(sub,1,1),'c','HorizontalAlignment','left','fontSize',18,'Color',[1 0 0],'fontName','Arial')
    end
%     if sub==12
%     text(0.9,cluster_value_persub(sub,1,1),'d','HorizontalAlignment','left','fontSize',18,'Color',[1 0 0],'fontName','Arial')
%     end
    col=(cluster_value_persub(sub,1)-cluster_value_persub(sub,2)>0)+1;
%     if sum(sub==ind_to_exc)==1
%         col=3;
%     end
    plot([1 2],[cluster_value_persub(sub,1) cluster_value_persub(sub,2)],colors(col),'lineWidth',3)
    hold on
end
set(gca,'Fontsize',18,'fontName','Arial')
% title('stim-baseline')
xlim([0.7 2.3])
% ylim([-2 2])
names = {'tVNS','Sham'};
set(gca,'xtick',[1:5],'xticklabel',names)
%% Check the significance in the cluster
[p, h, stats] = signrank(cluster_value_persub(:,1),cluster_value_persub(:,2));
[h, p, stats] = ttest(cluster_value_persub(:,1),cluster_value_persub(:,2));

save('single_sub_eeg_eff_cluster','sub_names','cluster_value_persub'); % percent per subject values

%%  Check if the decrease in alpha is siginficant in each of the conditions alone:
    
cfg=[];
cfg.latency  = [-6 0]; % this could be done using clustering;
one_cond_base=cellfun(@(x) ft_selectdata(cfg,x),D_lat,'un',0);
cfg.latency  = [0.15 6];
one_cond_stim=cellfun(@(x) ft_selectdata(cfg,x),D_lat,'un',0);
for j=1:numel(one_cond_base)
    one_cond_base{j}.time = one_cond_stim{1}.time;
end

for cond=1:2
cfg= [];
ind_to_inc= 1:21;
cfg.frequency = [5 15];
cfg.alpha  =0.05;
cfg.parameter = 'powspctrm';
cfg.method = 'montecarlo';
cfg.avgovertime = 'no';
cfg.avgoverfreq = 'no';
cfg.correctm= 'cluster';
cfg.tail = 0;
cfg.clustertail=0;
cfg.statistic = 'depsamplesT';
cfg.numrandomization = 1000;

D=[one_cond_base(:,cond) one_cond_stim(:,cond)]; %excluding subject with no pupil response
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design = [ones(1,length(D)) ones(1,length(D))*2 ;1:21 1:21];

[stat_single{cond}] = ft_freqstatistics(cfg, D{:});
end
%%
gap=size(ga_median612{cond}.mask,3)- size(stat_single{1}.mask,3);

for cond=1:2
ga_median612{cond}.mymask = logical(cat(3,zeros(1,31,floor(gap/2)),stat_single{cond}.mask,zeros(1,31,ceil(gap/2))));
end
%% plot the tfr with mask that fits significance of each seperately
cfg=[];
titles= {'tVNS --','sham --'};
cfg.zlim= [0.9 1.2]*100; % rel baseline
cfg.ylim=[5.3 14.7];
cfg.xlim= [-6 12];
cfg.maskparameter ='mymask';
% cfg.maskparameter ='mask';

cfg.maskstyle= 'outline';
cfg.maskalpha= 0.9;
cfg.fontsize  = 18;
sub=20;
for i=1:2
    pup_fig=figure('Renderer', 'painters', 'Position', [10 10 910 910/gold^2]);
    ft_singleplotTFR(cfg,ga_median612{i});
    title(titles{i})
    r=rectangle('Position',[ 0 5 3.4 10],'LineWidth',0.5,'LineStyle','--');
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    set(gca,'fontSize',18,'fontName','Arial')
    
    hold on    
end

