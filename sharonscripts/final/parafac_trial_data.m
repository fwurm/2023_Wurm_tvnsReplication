% this script uses data from the real expriment and then applies components
% achieved from the break data

folder='C:\tvns_pilot\short_pulses_experiment\';
load([folder 'layout_193.mat']);
legit_labels=layout193.label(1:end-3); % skull channels

% trl_data_folder='eeg_all_chan_algoclean\';
trl_data_folder='eeg_all_chan_algo100\';
break_data_folder='break_data\';

trl_datasets=dir([folder  trl_data_folder '*eeg_all_cond_avref*']); %datasets with the trialed stimulation data

%% run on subjects:
for sub=1:length(trl_datasets)
    trl_datasets(sub).name(1:3)
    %load the break model:
    % load the channels that were excluded in the break model:
    load([folder break_data_folder trl_datasets(sub).name(1:3) '_breakdata_cut_clean.mat'],'exc_chan') 
    include_chan=setxor(1:size(legit_labels),exc_chan); % these are the channels to include
    include_label=legit_labels(include_chan);

    load([folder trl_data_folder trl_datasets(sub).name])     %load the trialed data
    % select only the relevant channels:
    cfg=[];
    cfg.channel=include_label;
    data_cond_sel=cellfun(@(x) ft_selectdata(cfg,x),data_cond,'un',0);
    clear data_cond % clear some memory
    %  frequency decomp following the same resoultion that used in the model:
    cfg=[];
    cfg.method='tfr'; % 
    cfg.foi=[5:0.33:15];
    cfg.tapsmofrq=ones(size(cfg.foi))*0.1; % smoothing - i.e +-0.1Hz
    cfg.toi=[-6:0.05:12];
    cfg.width = 7;
    cfg.gwidth=3;
    cfg.keeptrials='yes';
    spect_ft=cellfun(@(x) ft_freqanalysis(cfg,x),data_cond_sel,'un',0); %results in win, chan, freq, original_trl  %takes 14 minutes
    save([folder trl_data_folder trl_datasets(sub).name(1:3) '_spect_ft_4parafac_algo_tfr7cycle_100'],'spect_ft');
end
%% load in case it was premade:
spect_ft_content = dir([folder trl_data_folder '*_4parafac_algo_tfr7cycle_100.mat']);

% load([folder trl_data_folder trl_datasets(sub).name(1:3) '_spect_ft_4parafac_algo'],'spect_ft');

%% Topogrpahy in  7-12 Hz - naively
naivetopo=figure;
cfg=[];
cfg.channel= {'all'};
cfg.parameter = 'powspctrm';
cfg.xlim  =[-15 15];
cfg.ylim  =[7 12];
%     cfg.zlim = [0 6];
cfg.layout= layout193;
cfg.marker='labels';
%     cfg.highlight='labels';
%     cfg.highlightchannel=chan_labels;
%     cfg.highlightcolor=[1 0 0];
ft_topoplotER(cfg,spect_ft{2})
% saveas(naivetopo,[folder trl_data_folder trl_datasets(sub).name(1:3) '_trl_naiveTopo7_12'],'fig');

%% let's look again on the tfr a of a specific channel:
figure
cfg=[];
cfg.channel = 'E126'; % or E107
cfg.trial = 'all';
% cfg.zlim= [0 1*10^17];
for i=1:2
subplot(2,1,i)
% ft_singleplotTFR(cfg,w_data_tfr_cond{i})
ft_singleplotTFR(cfg,spect_ft{i})
% title([num2str(i) ' comp2'])
end
%% let's try to fit a prafac model to these data:
Nk=2;
factorsc={};
spect_eeg=permute(spect_ft{2}.powspctrm,[2 3 4 1]);
free_model = parafac(spect_eeg,Nk,[1e-8 0 0],[0 0 0 0]); %convergence criterion, 2- for nonnegativity

factosFig=figure;
disp(['try ' num2str(Nk) ' components:'])
factorsFig=figure;
disp(['try ' num2str(Nk) ' components:'])
mode_strings ={'Channel';'Frequency (Hz)';'time (sec)';'trial'};
plotfac(free_model,[],NaN,spect_ft{1}.freq,spect_ft{1}.time,NaN,mode_strings)

%% Use componenets from break to fit the data:
% pup_reponders=[3 5 7 8 10 11 12 14 15 17 18 19 22 23 24];
% pup_non_responders=   [4 6 16 21];
% sub_select=all_sub_with_alpha;
% multisub_trlmodel={};
how_many_fac  = [3,2,2,4,2,4,4,3,2,2,2,2,3,2,2,2,3,2,3,4,2];
whichfac_alpha= [3,2,1,4,1,3,3,3,2,2,1,2,3,2,1,2,1,2,3,4,1];

% sub_par=10; % temporary
% disp ( ['# of comp is: ' num2str(how_many_fac(sub_par))] )
% disp ( ['alpha comp is:' num2str(whichfac_alpha{sub_par})])

for sub=1:length(spect_ft_content)
    % load the trial data ready for parafac: - it does not fit
    load([folder  trl_data_folder spect_ft_content(sub).name]);
    sub_folder=[spect_ft_content(sub).name(1:3) '_plots_fitted\'];
    load([folder break_data_folder sub_folder spect_ft_content(sub).name(1:3) '_' num2str(how_many_fac(sub)) '_sel_factors.mat'])
    
    % now use the component in every trial:
    spect_trl=cellfun(@(x) x.powspctrm,spect_ft,'un',0); % cell is for tvns/sham then inside: trl, chan, freq, time
    
    % run on trials independently: // WHY? beacause it allows averaging.
    t_selection=10:83;

    fac=size(sel_factors{1},2);
    
    for cond=1:size(spect_trl,1)
            % permute the trl data to fit the parafac as the breakmodel data:
            % spect_eeg=permute(one_trial,[2,3,1]);% move to: chan, freq, win - this is 1 trial
            trl_data= permute(spect_trl{cond}(:,:,:,t_selection),[2 3 4 1]);
            %  2. Make a cell array that holds initial values for the model for the test data. The loadings should equal the calibration model loadings and the scores should be random of appropriate size (number of samples in test data times number of factors)	Itest = size(Xtest,1); % Number of samples in test set
            Itest3=size(trl_data,3);% time, of the samples
            Itest4=size(trl_data,4);% time, of the samples

            Loads = sel_factors; %chan, spect, window
            Loads{3} = rand(Itest3,fac);
            Loads{4} = rand(Itest4,fac);
            
            % 3. Fit a PARAFAC model to the test data but make sure that the initial loadings are input and that the loadings in the variable modes are fixed.
            
            opt = [1e-8 0 0]; % same as model, I enter only convergence parametes
            const = [0 0 2 2]; % Change to [2 0 0] if you want nonnegative scores. In the second and third mode, nothing can be changed in prediction because B and C are fixed.
            Fix = [1 1 0 0]; % Defines that A (channels) and B (frequency) should maintain their initial values (hence no estimation of these and the output A and B will equal the ones from the calibration model)
            
            model{cond}=parafac(trl_data,fac,opt,const,Loads,Fix);
            
            %let's look
            factorsFig=figure;
            disp(['try ' num2str(fac) ' components:'])
            mode_strings ={'Channel';'Frequency (Hz)';'Time (sec)';'trial'};
            plotfac(model{cond},[],NaN,spect_ft{1}.freq,spect_ft{1}.time(t_selection),NaN,mode_strings)
    end
    save([folder 'eeg_all_chan_algoclean\models\' spect_ft_content(sub).name(1:3) '_model'],'model')
end
%% let's look on the single subject TFR on a specifc channels:
% weight according to the comp
comp=whichfac_alpha{sub_par}; % PAY ATTENTION TO PUT THE RIGHT COMPONENT MANUALLY
for cond=1:2
    A=model{cond}{1}; % these are the channels from the model i.e from the break.
    V= A(:,comp)/sum(A(:,comp)); % + A(:,comp(2)); % take compenets from comp

    Y=permute(spect_ft{cond}.powspctrm,[1 3 4 2]);
    B = sum(bsxfun(@times,Y,reshape(V,1,1,1,numel(V))),4);

    B2 =reshape(B,size(B,1),1,size(B,2),size(B,3));
    % replace the weigted average with whole data
    cfg=[];
    cfg.channel = 'E192';
    cfg.trials = 1;
    w_spect_ft{cond}=ft_selectdata(cfg,spect_ft{cond});
    w_spect_ft{cond}.powspctrm = mean(B2,1); % median or mean on the trials, #how many trials persub?
end
save([folder trl_data_folder 'weigthed_trf\' trl_datasets(sub).name(1:3) '_w_spect_ft'],'w_spect_ft')
%% let's look on the w TFR
trf_figure =figure;
cfg=[];
cfg.channel = 'E192';
% cfg.channels = 'VREF'; % or E139
cfg.trial = 'all';
% cfg.zlim= [0 8*10^6];
for i=1:2
subplot(2,1,i)
% ft_singleplotTFR(cfg,w_data_tfr_cond{i})
ft_singleplotTFR(cfg,w_spect_ft{i})
title(['cond' num2str(i)])
end
saveas(trf_figure,[folder trl_data_folder 'figures_time\' trl_datasets(sub).name(1:3) '_trf_parafac_comp.fig'])
saveas(trf_figure,[folder trl_data_folder 'figures_time\' trl_datasets(sub).name(1:3) '_trf_parafac_comp.png'])

%% let's look on a-priori single chan TFR
trf_figure = figure;
cfg=[];
cfg.channel = 'E192';
% cfg.channels = 'VREF'; % or E139
cfg.trial = 'all';
cfg.zlim= [0 20];
for i=1:2
subplot(2,1,i)
% ft_singleplotTFR(cfg,w_data_tfr_cond{i})
ft_singleplotTFR(cfg,spect_ft{i})
title(['cond' num2str(i)])
end
saveas(trf_figure,'trf_figure.png')
%% plot
sel_sub=multisub_trlmodel(sub_select);

sel_sub_avg=cellfun(@(x) cellfun(@(y) squeeze(nanmean(y,1)),x,'un',0),sel_sub,'un',0);
% avg_trl_model{sub}=cellfun(@(x) squeeze(nanmedian(x,1)),trial_model,'un',0);

sel_sub_avg_alpha=cellfun(@(x,y) cellfun(@(z) sum(z(:,y),2),x,'un',0),sel_sub_avg,whichfac_alpha(sub_select),'un',0);

temp=cellfun(@(x) cell2mat(x'),sel_sub_avg_alpha,'un',0);

sel_sub_avg_alpha_mat=cat(3,temp{:});

% bc correction
bc_persub=repmat((mean(sel_sub_avg_alpha_mat,1)),[28 1]);
sel_sub_avg_alpha_mat_bc=arrayfun(@(x,y) x-y,sel_sub_avg_alpha_mat,bc_persub);

cfg.toi=[-14.5:1:16.5]; %3sec windows


sel_sub_avg_alpha_mat_bc_exc=sel_sub_avg_alpha_mat_bc(:,:,[1:4 6 8:19]);
avg_alpha_timecourse=nanmedian(sel_sub_avg_alpha_mat_bc,3);

std_alpha_timecourse=std(sel_sub_avg_alpha_mat_bc,0,3)/sqrt(size(sel_sub_avg_alpha_mat_bc,3));

%% plotting

for cond=1:2
plot(cfg.toi(t_selection),avg_alpha_timecourse(:,cond)-mean(avg_alpha_timecourse(:,cond)))
hold on
end
legend('tvns','sham')
xlim([-15 15])


%%
color='rb';
for cond=1:2
shadedErrorBar(cfg.toi(t_selection),avg_alpha_timecourse(:,cond),std_alpha_timecourse(:,cond),color(cond))
hold on
end

%% singlesub plotting
for cond=1:2
subplot(2,1,cond)
plot(cfg.toi(t_selection),squeeze(sel_sub_avg_alpha_mat_bc_exc(:,cond,:)))
xlim([-15 15])

end

%%
imagesc(squeeze(sel_sub_avg_alpha_mat_bc(:,1,:))')