function batch_peeg_sharon
% Batch file for the replication of Sharon et al (2021) "Transcutaneous Vagus 
% Nerve Stimulation in Humans Induces Pupil Dilation and Attenuates Alpha 
% Oscillations" - pipeline scripts directly from Sharon

%%% Directories %%%
fdir.dir = 'C:\Users\wurmf\Dropbox\GitHub\2023_Wurm_tvnsReplication';
fdir.dir_data = 'C:\Users\wurmf\surfdrive\Shared\raw_data\EEG_data\'; %cloud
fdir.EEGdir = fdir.dir_data;
fdir.TFAdir = strcat(fdir.dir_data, 'TimeFreq\');


%%% Add specific paths %%%
addpath('Routines')
addpath('Routines/nway331-toolbox\')
addpath('Routines/eeg_2d_minimal_rejection\')

% Add fieltrip
ftpath = 'C:\Users\wurmf\Documents\Github\fieldtrip\';
% enableFieldTrip(1,ftpath) %TODO before using fieltrip routines


%%% Load participant structure %%%
% SCode = [1:4 6:27 29:39]; %full datasets
SCode = [6:11 13:17 19 20 24:27 29 31:33 35 36 39]; %%final participants
[S] = getS(SCode,fdir);


%%% load electrode layout for this study (standard 64 electrode)
fid = fopen(fullfile(fdir.dir,'elecLabels_64_standard.txt'));
elecinfo = textscan(fid,'%d %s','delimiter','\t');
fclose(fid);

% load channel maps and create neighbours map accordingly
% prepare_neighbours
cfg_layout = [];
cfg_layout.layout = 'biosemi64.lay';
cfg_neighb           = [];
cfg_neighb.layout = ft_prepare_layout(cfg_layout);
cfg_neighb.method    = 'triangulation';
cfg_neighb.feedback  = 'no';
neighbours       = ft_prepare_neighbours(cfg_neighb);


% eeg1_import(S,chanFile,'highpass',0.1,'in','\','relabel',elecinfo,'referenceElectrodes',[48]) %reref to Cz
% eeg2_epoch(S,33,[-10 10],'out','stimLocked\e')

%%%%%%%%%%%%%%%%%%
%% sharon approach
%%%%%%%%%%%%%%%%%%

%load info for stimulation artefacts
load('ExperimentInfo_final.mat','INFO')


outpath = fullfile(fdir.EEGdir, 'stimLocked\e\sharon');
if ~exist(outpath, 'dir')
    mkdir(outpath)
end

%% Step1: preprocessing: experiment data
if 0
    % remove artefacts, filter, interpolation
    for isub = 1:length(S)
        
        %load EEGLAB file (lp filter & epoched only)
        inpath = fullfile(S(isub).EEGdir, 'stimLocked\e');
        raw = pop_loadset('filename',[S(isub).EEGfn '.set'],'filepath',inpath);
        
        %export to fieldtrip format
        FFT = eeglab2fieldtrip(raw,'preprocessing','none');
        
        % manually exclude trials without stimulation artefact
        cfg=[];
        if ~isnan(INFO(S(isub).index).noArtefact)
            cfg.trials= setxor(INFO(S(isub).index).noArtefact,1:length(FFT.trial));
        else
            cfg.trials = 1:length(FFT.trial)
        end
        FFT_trl=ft_selectdata(cfg,FFT);
        clear FFT
        
        % reduce data to skull channels only
        cfg=[];
        cfg.channel = ft_channelselection({'eeg'}, FFT_trl.label);
        FFT = ft_selectdata(cfg,FFT_trl);
        clear FFT_trl
        
        % filter out the stimulation frequency in 25Hz and its harmonics:
        cfg=[];
        cfg.bsfilter      = 'yes';
        notch_f=25:25:250;
        cfg.bsfreq=[notch_f'-3 notch_f'+3];
        cfg.bsfiltord = 3;
        cfg.bpinstabilityfix = 'reduce';
        cfg.bpfilter= 'yes';
        cfg.bpfreq = [5 15];
        FFT_dft=ft_preprocessing(cfg,FFT);
        clear FFT
        
        % use a absloute crietria to reject trials/channels:
        MAXABS_THRESHOLD=100;
        max_abs_values=cell2mat(cellfun(@(x) max(abs(x),[],2),FFT_dft.trial,'un',0)); % max_ abs values chan*trl
        unvalid=max_abs_values>MAXABS_THRESHOLD; % this a boolean matrix of the currenly unvalid chan*trl combinations
        
        view_on=0;
        % [exc_chan, exc_trl] = minimal_automatic_rejection(unvalid,view_on,'trl');
        [exc_chan, exc_trl] = minimal_automatic_rejection(unvalid,view_on,'chan');
        
        allchan{isub} = exc_chan;
        alltrl{isub} = exc_trl;
        
        % interpolate bad channels and reject bad _trials:
        if ~isempty(exc_chan)
            cfg=[];
            cfg.trials=setxor(1:1:length(FFT_dft.trial),exc_trl);
            cfg.neighbours=neighbours;
            cfg.method='spline';
            cfg.badchannel=FFT_dft.label(exc_chan);
            cfg.layout = 'biosemi64.lay';
            FFT_intrp=ft_channelrepair(cfg,FFT_dft);
        else
            cfg=[];
            cfg.trials=setxor(1:1:length(FFT_dft.trial),exc_trl);
            FFT_intrp=ft_selectdata(cfg,FFT_dft);
        end
        
        
        % save file
        fn = fullfile(outpath,[S(isub).EEGfn '.mat']);
        parsave(fn,FFT_intrp,'FFT')
        
        gnu = 1;
    end
end
%% Step 2: splitting into condition and average reference
if 0
    cond = [];
    cond(1).name = 'sham';
    cond(1).num = 0;
    cond(2).name = 'stim';
    cond(2).num = 1;
    
    warning('changed for loop index to start with 5')
    for isub = 5:length(S)
        fn = fullfile(outpath,[S(isub).EEGfn '.mat']);
        load(fn,'FFT')
        
        % split sham and stim
        b = cell(1,2);
        for icond=1:2
            cfg=[];
            b{icond}=find(FFT.trialinfo.condition==cond(icond).num);
            cfg.trials=b{icond};
            FFT_cond{icond}=ft_selectdata(cfg,FFT);
        end
        
        %average reference
        cfg=[];
        cfg.reref='yes';
        cfg.refchannel = 'all';
        cfg.refmethod=   'avg';
        FFT_cond_avref=cellfun(@(x) ft_preprocessing(cfg,x),FFT_cond,'un',0);
        
        %time-frequency
        %  frequency decomp following the same resoultion that used in the model:
        cfg=[];
        cfg.method='tfr'; %
        cfg.foi=[5:0.33:15];
        cfg.tapsmofrq=ones(size(cfg.foi))*0.1; % smoothing - i.e +-0.1Hz
        %     cfg.toi=[-3:0.05:8];
        timewin = [-8 8];
        timeidx = dsearchn(FFT_cond_avref{1}.time{1,1}',timewin')
        cfg.toi=FFT_cond_avref{1}.time{1,1}(timeidx(1):16:timeidx(2));
        cfg.width = 7;
        cfg.gwidth=3;
        cfg.keeptrials='yes';
        FFT_ft=cellfun(@(x) ft_freqanalysis(cfg,x),FFT_cond_avref,'un',0); %results in win, chan, freq, original_trl  %takes 14 minutes
        
        fn = fullfile(outpath,[S(isub).EEGfn '_fft.mat']);
        parsave(fn,FFT_ft)
        
        
        
    end
end

%% Step3: preprocessing: break data
outpath_break = fullfile(fdir.EEGdir, 'breakLocked\e\sharon');
if ~exist(outpath_break, 'dir')
    mkdir(outpath_break)
end
if 0
    % remove artefacts, filter, interpolation
    for isub = 1:length(S)
        %load EEGLAB file (lp filter & epoched only)
        inpath = fullfile(S(isub).EEGdir, 'breakLocked\e');
        raw = pop_loadset('filename',[S(isub).EEGfn '.set'],'filepath',inpath);
        
        %export to fieldtrip format
        FFT = eeglab2fieldtrip(raw,'preprocessing','none');
        
        % reduce data to skull channels only
        cfg=[];
        cfg.channel = ft_channelselection({'eeg'}, FFT.label);
        FFT = ft_selectdata(cfg,FFT);
        
        % filter out the stimulation frequency in 25Hz and its harmonics:
        cfg=[];
        cfg.bpinstabilityfix = 'reduce';
        cfg.bpfilter= 'yes';
        cfg.bpfreq = [5 15];
        FFT=ft_preprocessing(cfg,FFT);
        
        % use a absloute crietria to reject trials/channels:
        MAXABS_THRESHOLD=100;
        max_abs_values=cell2mat(cellfun(@(x) max(abs(x),[],2),FFT.trial,'un',0)); % max_ abs values chan*trl
        unvalid=max_abs_values>MAXABS_THRESHOLD; % this a boolean matrix of the currenly unvalid chan*trl combinations
        
        view_on=0;
        % [exc_chan, exc_trl] = minimal_automatic_rejection(unvalid,view_on,'trl');
        [exc_chan, exc_trl] = minimal_automatic_rejection(unvalid,view_on,'chan');
        
        fprintf('   bad channels: %s\n',strjoin(FFT.label(exc_chan)))
        fprintf('   bad trials: %s\n',num2str(exc_trl))
        badstuff(isub,:) = [length(exc_chan) length(exc_trl)];
        
        % interpolate bad channels and reject bad _trials:
        if ~isempty(exc_chan)
            cfg=[];
            cfg.trials=setxor(1:1:length(FFT.trial),exc_trl);
            cfg.neighbours=neighbours;
            cfg.method='spline';
            cfg.badchannel=FFT.label(exc_chan);
            cfg.layout = 'biosemi64.lay';
            FFT=ft_channelrepair(cfg,FFT);
        else
            cfg=[];
            cfg.trials=setxor(1:1:length(FFT.trial),exc_trl);
            FFT=ft_selectdata(cfg,FFT);
        end
        
        %average reference
        cfg=[];
        cfg.reref='yes';
        cfg.refchannel = 'all';
        cfg.refmethod=   'avg';
        FFT_avref=ft_preprocessing(cfg,FFT);
        
        %time-frequency
        %  frequency decomp following the same resoultion that used in the model:
        cfg=[];
        cfg.method='tfr'; %
        cfg.foi=[5:0.33:15];
        cfg.tapsmofrq=ones(size(cfg.foi))*0.1; % smoothing - i.e +-0.1Hz
        %     cfg.toi=[-3:0.05:8];
        timewin = [0 3];
        timeidx = dsearchn(FFT_avref.time{1,1}',timewin');
        cfg.toi=FFT_avref.time{1,1}(timeidx(1):16:timeidx(2));
        cfg.width = 7;
        cfg.gwidth=3;
        cfg.keeptrials='yes';
        FFT_ft=ft_freqanalysis(cfg,FFT_avref); %results in win, chan, freq, original_trl  %takes 14 minutes
        
        
        % save file
        fn = fullfile(outpath_break,[S(isub).EEGfn '_fft.mat']);
        parsave(fn,FFT_ft,'FFT')
        
        gnu = 1;
    end
end

%% Step4: run parafac on break data
if 0
    for isub = 1:length(S)
        fn = fullfile(outpath_break,[S(isub).EEGfn '_fft.mat']);
        load(fn,'FFT')
        
        %     break_data= permute(FFT.powspctrm,[2 3 4 1]);
        break_data = FFT.powspctrm;
        
        opt = [1e-8 0 0]; % same as model, I enter only convergence parametes
        const = [2 2 2 2];
        
        [ssX,Corco,It] = pftest(3,break_data,4,opt,const);
        
        out.ssX = ssX;
        out.Corco = Corco;
        out.It = It;
        out.nComp = find(mean(Corco,2)>55,1,'last');
        
        nFactors = out.nComp;
        Factors = parafac(break_data,nFactors,opt,const);
        
        out.Factors = Factors;
        
        save(fullfile(outpath_break,[S(isub).EEGfn '_pftest.mat']),'out');
        
        % % plot some of the results
        %     factorsFig=figure;
        %     disp(['try ' num2str(nFactors) ' components:'])
        %     mode_strings ={'trial';'Channel';'Frequency (Hz)';'Time (sec)'};
        % %     plotfac(Factors,[],NaN,NaN,FFT.freq,FFT.time,mode_strings)
        % %     plotfac(Factors)
        %
        %     load('chanlocs.mat','chanlocs')
        %     figure
        %     topoplotMS(Factors{1,2}(:,1),chanlocs(1:64),'maplimits',[-0.4 0.4]);
        %     colormap('redblue')
        %
        %     gnu =1;
        
    end
    
    load('chanlocs.mat','chanlocs')
    for isub = 1:length(S)
        load(fullfile(outpath_break,[S(isub).EEGfn '_pftest.mat']),'out');
        
        f = figure
        
        plotfac(out.Factors)
        
        for i = 1:out.nComp
            ax = axes('parent',f)
            ax.Position = [0.55+(i-1)*0.11 0.7 0.1 0.15];
            topoplotMS(out.Factors{1,2}(:,i),chanlocs(1:64),'maplimits',[-0.4 0.4]);
            colormap('redblue')
        end
        
        ax = axes('parent',f)
        ax.Position = [0.03 0.8 0.06 0.04];
        text(0.5,0.5,S(isub).EEGfn)
        axis off
        
        
        close all
        
    end
end

sub_number = [6,7,8,9,10,11,13,14,15,16,17,19,20,24,25,26,27,29,31,32,33,35,36,39];
how_many_fac  = [1,2,2,2,2,2,1,1,1,2,2,3,2,1,3,1,2,2,1,1,2,2,1,1];
whichfac_alpha= [1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1];

Sidx = ismember(sub_number,[S.index]);
how_many_fac(~Sidx) = [];
whichfac_alpha(~Sidx) = [];

%% Step5: apply parafac to trial data
if 0
    for isub = 1:length(S)
        %load trial data
        fn = fullfile(outpath,[S(isub).EEGfn '_fft.mat']);
        load(fn,'data');
        
        %load parafac break results
        load(fullfile(outpath_break,[S(isub).EEGfn '_pftest.mat']),'out');
        
        % now use the component in every trial:
        spect_trl=cellfun(@(x) x.powspctrm,data,'un',0); % cell is for tvns/sham then inside: trl, chan, freq, time
        
        fac = how_many_fac(find(sub_number==S(isub).index));
        
        for cond=1:size(spect_trl,2)
            
            % 1. permute the trl data to fit the parafac as the breakmodel data:
            trl_data= permute(spect_trl{cond},[2 3 4 1]);
            
            %  2. Make a cell array that holds initial values for the model for the test data. The loadings should equal the calibration model loadings and the scores should be random of appropriate size (number of samples in test data times number of factors)	Itest = size(Xtest,1); % Number of samples in test set
            Itest3=size(trl_data,3);% time, of the samples
            Itest4=size(trl_data,4);% time, of the samples
            
            Loads{1} = out.Factors{2}; %chan, spect, window
            Loads{2} = out.Factors{3}; %chan, spect, window
            Loads{3} = rand(Itest3,fac);
            Loads{4} = rand(Itest4,fac);
            
            % 3. Fit a PARAFAC model to the test data but make sure that the initial loadings are input and that the loadings in the variable modes are fixed.
            opt = [1e-8 0 0]; % same as model, I enter only convergence parametes
            const = [0 0 2 2]; % Change to [2 0 0] if you want nonnegative scores. In the second and third mode, nothing can be changed in prediction because B and C are fixed.
            Fix = [1 1 0 0]; % Defines that A (channels) and B (frequency) should maintain their initial values (hence no estimation of these and the output A and B will equal the ones from the calibration model)
            
            model{cond}=parafac(trl_data,fac,opt,const,Loads,Fix);
            
            %let's plot the results
            factorsFig=figure;
            disp(['try ' num2str(fac) ' components:'])
            mode_strings ={'Channel';'Frequency (Hz)';'Time (sec)';'trial'};
            plotfac(model{cond})
            
        end
        
        save(fullfile(outpath,[S(isub).EEGfn '_model.mat']),'model');
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting the PARAFAC results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    model_content = dir([outpath '\*model.mat']);
    model_cell = cell(size(model_content,1),2);
    
    for isub=1:length(model_content)
        temp  = load([outpath '\' model_content(isub).name]);
        temp1=temp.model;
        for cond=1:2
            model_cell{isub,cond}=temp1{cond};
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
    ylim([0 0.4])
    set(gca,'fontsize',18,'fontName','Arial')
    xlabel('Frequency (Hz)')
    ylabel('Weight')
    
    % plot single subject components
    figure
    for isub=1:length(S)
        plot(freq_axis,freq_model{isub,1})
        hold on
    end
end

%%%%%%%%%%%%%%%%%%%%
%% Data analysis
%%%%%%%%%%%%%%%%%%%


%% get the time- only data using the parafac model:

gold=(1+sqrt(5))/2;

load_only = 1;

if ~load_only
    C=cell(length(S),2);
    for isub=1:length(S)
        fprintf('#########\n')
        fprintf('## %s \n',S(isub).EEGfn)
        fprintf('#########\n')
        
        %load trial data
        fn = fullfile(outpath,[S(isub).EEGfn '_fft.mat']);
        load(fn,'data');
        
        %load trial data with parafac
        load(fullfile(outpath,[S(isub).EEGfn '_model.mat']),'model');
        
        comp=whichfac_alpha(find(sub_number==S(isub).index)); % select the alpha componenet
        for cond=1:2
            % reduce the channels domain:
            chan_model=model{cond}{1}; % these are the channels from the model i.e from the break.
            V_chan= chan_model(:,comp); % get the right comp
            Y=permute(data{cond}.powspctrm,[1 3 4 2]); % move the channels dim to the end
            B = sum(bsxfun(@times,Y,reshape(V_chan,1,1,1,numel(V_chan))),4);
            
            B2 =reshape(B,size(B,1),1,size(B,2),size(B,3)); % trl, ~, freq, time
            % reduce the freq domain:
            freq_model =model{cond}{2};
            V_freq = freq_model(:,comp);
            X=permute(B2,[1 2 4 3]); % move the freq domain to the end
            C{isub,cond} = squeeze(sum(bsxfun(@times,X,reshape(V_freq,1,1,1,numel(V_freq))),4));
        end
    end
    
    time_axis = data{1}.time;
    freq_axis = data{1}.freq;
    save(fullfile(outpath,'Trialdata.mat'),'C','time_axis','freq_axis');
else
    %load pre-saved data
    load(fullfile(outpath,'Trialdata.mat'),'C','time_axis','freq_axis');
end



%% compare # of trials:
trial_num=cellfun(@(x) size(x,1),C);

mean_trl_cond=mean(trial_num,1);
[p,h]= signrank(trial_num(:,1),trial_num(:,2));
sem_trl_cond=std(trial_num,[],1)/sqrt(length(S));
min_trl_cond=min(trial_num);
max_trl_cond=max(trial_num);

%% Baselining according to sharon (average level)

% take mean over trials
c_mean=cellfun(@(x) mean(x,1),C,'un',0);

%cenvert cell to matrix (subj,time,condition)
G=cell2mat(c_mean(:,1));
G(:,:,2) = cell2mat(c_mean(:,2));

% get grandAVG
grand_G=squeeze(nanmean(G,1));

%identify baseline time points, plot and extract baseline mean
baseline_lim =[-1 0];
baseline_ind = time_axis>baseline_lim(1) & time_axis<baseline_lim(2);
baseline_val= mean(G(:,baseline_ind,:),2);
grand_base_values=mean(baseline_val(:,:,2)-baseline_val(:,:,1),1);

%apply subject baseline
A=bsxfun(@(x,y) (rdivide(x,y)),G,baseline_val);

%compare baseline values
[p,h]= signrank(baseline_val(:,1,1),baseline_val(:,1,2));


%% Compare grandaverages.

% average across subjects:
grand_avg_per=squeeze(nanmean(A,1));
std_per=squeeze(nanstd(A,[],1))/sqrt(length(S)); % this is SEM
stim_lim =[0 4];


stim_ind = time_axis>stim_lim(1) & time_axis<stim_lim(2);
stim_values = nanmean(A(:,stim_ind,:),2);
values= stim_values*100;

save(fullfile(outpath,'single_sub_eeg_eff'),'S','values'); % percent per subject values

%directed tests
fprintf('Compare directed (tvns<sham)\n')
[p, h, stats] = signrank(values(:,1,1),values(:,1,2),'tail','left');
fprintf('\tWilcoxon: z = %.2f, p = %.2f\n',stats.zval,p)
[h, p, ~, stats] = ttest(values(:,1,1),values(:,1,2),'tail','left');
fprintf('\tt-test: t(%d) = %.2f, p = %.2f\n',stats.df,stats.tstat,p)

% undirected test
fprintf('Compare uncirecteddirected (tvns~=sham)\n')
[p, h, stats] = signrank(values(:,1,1),values(:,1,2));
fprintf('\tWilcoxon: z = %.2f, p = %.2f\n',stats.zval,p)
[h, p, ~, stats] = ttest(values(:,1,1),values(:,1,2));
fprintf('\tt-test: t(%d) = %.2f, p = %.2f\n',stats.df,stats.tstat,p)

%against zero
fprintf('Compare sham against zero\n')
[p, h, stats] = signrank(values(:,1,1),100);
fprintf('\tWilcoxon: z = %.2f, p = %.2f\n',stats.zval,p)

fprintf('Compare tvns against zero\n')
[p, h, stats] = signrank(values(:,1,2),100);
fprintf('\tWilcoxon: z = %.2f, p = %.2f\n',stats.zval,p)

report_mean=squeeze(mean(values,1));
report_sem=std(squeeze(values),1)./sqrt(length(S));

% plot effect size figure (3B)
eff_size_fig=figure('Renderer', 'painters', 'Position', [10 10 900/gold^2 900/gold]);
color='kk';
kind={'-','--'};
for isub=1:length(S) 
    col=(stim_values(isub,1,1)-stim_values(isub,1,2)>0)+1;
    plot([1 2],[values(isub,1,1) values(isub,1,2)],kind{col},'lineWidth',3,'Color',[0 0 0]+(col-1)*0.5)
    hold on
end
set(gca,'Fontsize',18,'fontName','Arial')

xlim([0.7 2.3])
ylim([70 140])

names = {'Sham','tVNS'};

set(gca,'xtick',[1:5],'xticklabel',names)
ylabel('% power')


%% plot grandaverages
figure
plot(time_axis,grand_avg_per)
xlim([-1 8])
legend('tvns','sham')

% plot the  grand average with std
colors='kr';
colors = [49, 116, 161; 225, 129, 43];

% %apply smoothing
% for cond=1:2
% grand_avg_per_smooth(:,cond)=lowPassFilter(2,grand_avg_per(:,cond),20)+mean(grand_avg_per(:,cond));
% end

trace_fig=figure('Renderer', 'painters', 'Position', [10 10 910 910/gold^2]);
for cond=1:2
ax(cond) = shadedErrorBar(time_axis,grand_avg_per(:,cond)*100,std_per(:,cond)*100,{'color',colors(cond,:)./255,'markerfacecolor',colors(cond,:)./255},1);
% ,colors(cond,:),0.2);
hold on
xlim([-1 6])
ylim([0.85 1.25]*100)
% ylim([0.9 1.9]*100)
end

set(gcf,'Color','w')
set(gca,'Fontsize',18,'fontName','Arial')
stim_patch=patch([0,0,3.4,3.4],[0.8,1.3,1.3,0.8]*100,'k','LineStyle','none');
alpha(stim_patch,0.2)
% colorbar
ylabel('Average alpha (%)')
xlabel('Time (sec)')
set(gca,'fontSize',20,'fontName','Arial')
line([-1 6],[100 100],'Color','k','LineWidth',2,'Linestyle','--')
line([0 0],[50 150],'Color','k','LineWidth',2,'Linestyle','--')
line([4 4],[50 150],'Color','k','LineWidth',2,'Linestyle','--')
legend([ax(1).mainLine ax(2).mainLine],{'sham' 'tvns'},'Box','off')

%% compare and plot each time point for difference
figure
hold on
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
plot(time_axis(time_ind),pvals.*p_masked)
plot(time_axis(time_ind),h.*p_masked,'r.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  load the full spectrum trials to show the same effect in TFR:  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('order_ind_1_1sec')

%% reduce to time x frequency
if 0
    for isub=1:length(S)
        fprintf('#########\n')
        fprintf('## %s \n',S(isub).EEGfn)
        fprintf('#########\n')
        
        %load trial data
        fn = fullfile(outpath,[S(isub).EEGfn '_fft.mat']);
        load(fn,'data');
        
        %load trial data with parafac
        load(fullfile(outpath,[S(isub).EEGfn '_model.mat']),'model');
        comp=whichfac_alpha(isub); % select the alpha componenet
        
        
        for cond=1:2
            % reduce the channels domain:
            chan_model=model{cond}{1}; % these are the channels from the model i.e from the break.
            V_chan= chan_model(:,comp); % get the right comp
            Y=permute(data{cond}.powspctrm,[1 3 4 2]); % move the channels dim to the end
            B = sum(bsxfun(@times,Y,reshape(V_chan,1,1,1,numel(V_chan))),4);
            
            B2 =reshape(B,size(B,1),1,size(B,2),size(B,3)); % trl, ~, freq, time
            
            cfg=[];
            cfg.channel = 'Cz'; %dummy
            cfg.trials = 1;
            w_spect_ft{cond}=ft_selectdata(cfg,data{cond});
            w_spect_ft{cond}.powspctrm = mean(B2,1); % median or mean on the trials, #how many trials persub?
            
            w_spect_ft{cond}.cfg.previous= []; % reduce size of the datasets
            
            %         save([folder 'weigthed_tfr_percent\' model_content(isub).name(1:3) '_w_spect_ft_1_1'],'w_spect_ft')
            
        end
        save(fullfile(outpath,[S(isub).EEGfn '_w_spect_ft_1_1.mat']),'w_spect_ft')
    end
end

%% reduce to time x channel
if 0
    for isub=1:length(S)
        fprintf('#########\n')
        fprintf('## %s \n',S(isub).EEGfn)
        fprintf('#########\n')
        
        %load trial data
        fn = fullfile(outpath,[S(isub).EEGfn '_fft.mat']);
        load(fn,'data');
        
        %load trial data with parafac
        load(fullfile(outpath,[S(isub).EEGfn '_model.mat']),'model');
        comp=whichfac_alpha(isub); % select the alpha componenet
        
        
        for cond=1:2
            % reduce the frequency domain:
            freq_model=model{cond}{2}; % these are the channels from the model i.e from the break.
            V_freq= freq_model(:,comp); % get the right comp
            Y=permute(data{cond}.powspctrm,[1 2 4 3]); % move the channels dim to the end
            B = sum(bsxfun(@times,Y,reshape(V_freq,1,1,1,numel(V_freq))),4);
            
            B2 =reshape(B,size(B,1),size(B,2),1,size(B,3)); % trl, chan, ~, time
            
            cfg=[];
            cfg.frequency = 10; %dummy
            cfg.trials = 'all'; %new
            cfg.avgoverrpt = 'no'; %new
            w_chan_ft{cond}=ft_selectdata(cfg,data{cond});
            w_chan_ft{cond}.powspctrm = mean(B2,1); % median or mean on the trials, #how many trials persub?
            
            w_chan_ft{cond}.cfg.previous= []; % reduce size of the datasets
            
            %         save([folder 'weigthed_tfr_percent\' model_content(isub).name(1:3) '_w_spect_ft_1_1'],'w_spect_ft')
            
        end
        save(fullfile(outpath,[S(isub).EEGfn '_w_chan_ft_1_1.mat']),'w_chan_ft')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Follow up analysis 1 (time x channel data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    % load('axis.mat')
    content = dir([outpath '\*_w_chan_ft_1_1.mat']);
    sub_cell = cell(size(content,1),2);
    
    for isub=1:length(content)
        temp  = load([outpath '\' content(isub).name]);
        temp1=temp.w_chan_ft;
        for cond=1:2
            %     temp1{cond}.label = {'alpha'};
            sub_cell{isub,cond}=temp1{cond};
        end
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
    
    mat_smoothed = cellfun(@(x) permute(conv2(squeeze(x.powspctrm),smoothingWindow,'same'),[3 1 4 2]),sub_baseline,'un',0);
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
    temp_out=cellfun(@(x)  permute(x.powspctrm,[2 1 3]),sub_cell_desc,'un',0);
    data_out=cell2mat(temp_out(:,1));
    data_out(:,:,:,2) =cell2mat(temp_out(:,2));
    
    %%  clusterbased permutation test
    cfg= [];
    ind_to_inc= 1:length(S);
    cfg.latency  = [-1 6]; % this could be done using clustering;
    cfg.frequency = [10];
    cfg.alpha  =0.05;
    cfg.parameter = 'powspctrm';
    cfg.method = 'montecarlo';
    cfg.avgovertime = 'no';
    cfg.avgoverfreq = 'no';
    cfg.neighbours = neighbours;
    cfg.correctm= 'cluster';
    % cfg.tail = -1;
    % cfg.clustertail=-1;
    cfg.tail = 0;
    cfg.clustertail=0;
    cfg.statistic = 'depsamplesT';
    cfg.numrandomization = 10000;
    
    D=sub_cell_desc(ind_to_inc,:); %excluding subject with no pupil response
    cfg.ivar = 1;
    cfg.uvar = 2;
    cfg.design = [ones(1,length(D)) ones(1,length(D))*2 ;1:length(S) 1:length(S)];
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
    cfg.operation='x2-x1';
    dif_tfr=ft_math(cfg,ga_median612{1},ga_median612{2});
    dif_tfr.mask=stat.mask;
    
    %% plotting
    cfg = [];
    cfg.xlim = [-1:1:6];
    cfg.ylim = [15 20];
    cfg.zlim = [80 120];
    cfg.layout = 'biosemi64.lay';
    
    
    cnames = {'sham' 'tvns'};
    tpnts = [0:1:6];
    for iC =1:2
        figure
        ft_topoplotTFR(cfg,ga_median612{iC});
        colorbar
        colormap('redblue')
    end
    
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Follow up analysis 2 (time x frequency data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1
    % load('axis.mat')
    content = dir([outpath '\*_w_spect_ft_1_1.mat']);
    sub_cell = cell(size(content,1),2);
    
    for isub=1:length(content)
        temp  = load([outpath '\' content(isub).name]);
        temp1=temp.w_spect_ft;
        for cond=1:2
            temp1{cond}.label = {'alpha'};
            sub_cell{isub,cond}=temp1{cond};
        end
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
    
    %%  clusterbased permutation test
    cfg= [];
    ind_to_inc= 1:length(S);
    cfg.latency  = [-1 6]; % this could be done using clustering;
    cfg.frequency = [5 15];
    cfg.alpha  =0.05;
    cfg.parameter = 'powspctrm';
    cfg.method = 'montecarlo';
    cfg.avgovertime = 'no';
    cfg.avgoverfreq = 'no';
    cfg.correctm= 'cluster';
    % cfg.tail = -1;
    % cfg.clustertail=-1;
    cfg.tail = 0;
    cfg.clustertail=0;
    cfg.statistic = 'depsamplesT';
    cfg.numrandomization = 10000;
    
    D=sub_cell_desc(ind_to_inc,:); %excluding subject with no pupil response
    cfg.ivar = 1;
    cfg.uvar = 2;
    cfg.design = [ones(1,length(D)) ones(1,length(D))*2 ;1:length(S) 1:length(S)];
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
    cfg.operation='x2-x1';
    dif_tfr=ft_math(cfg,ga_median612{1},ga_median612{2});
    dif_tfr.mask=stat.mask;
    
    %% plotttingtfr
    % sub=3;
    cfg=[];
    
    titles= {'sham','tVNS'};
    cfg.zlim= [0.8 1.2]*100; % rel baseline
    
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
        set(gcf,'Color','white')
        colormap('redblue')
        h=colorbar;
        ylabel(h, '% power','fontSize',20,'fontName','Arial')
    end
    
    %% plot the stats clusters:
    figure
    imagesc(stat.time,stat.freq,squeeze(stat.negclusterslabelmat-stat.posclusterslabelmat)) % relchange, relative and db works. but not absoulte
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
    colormap('redblue')
    
    
    %% a-priori  statistics:
    %%getting single subject values in the cluster:
    cfg=[];
    cfg.latency  = [-1 6];
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
    
    sem_cluster=std(cluster_value_persub*100)./sqrt(length(S));
    
    [m]=bootstrp(1000,@(x) median(x),cluster_value_persub(:,2));
    
    figure
    [fi,xi] = ksdensity(m);
    plot(xi,fi);
    
    up=prctile(m,95);
    down=prctile(m,5);
    
    %% Effect size figure for the cluster - defined by mask or a-priori above
    colors='kr';
    
    eff_size_fig=figure('Renderer', 'painters', 'Position', [10 10 900/gold^2 900/gold]);
    for isub=1:length(S)
        if isub==3
            text(0.9,cluster_value_persub(isub,1,1),'c','HorizontalAlignment','left','fontSize',18,'Color',[1 0 0],'fontName','Arial')
        end
        %     if sub==12
        %     text(0.9,cluster_value_persub(sub,1,1),'d','HorizontalAlignment','left','fontSize',18,'Color',[1 0 0],'fontName','Arial')
        %     end
        col=(cluster_value_persub(isub,1)-cluster_value_persub(isub,2)>0)+1;
        %     if sum(sub==ind_to_exc)==1
        %         col=3;
        %     end
        plot([1 2],[cluster_value_persub(isub,1) cluster_value_persub(isub,2)],colors(col),'lineWidth',3)
        hold on
    end
    set(gca,'Fontsize',18,'fontName','Arial')
    % title('stim-baseline')
    xlim([0.7 2.3])
    % ylim([-2 2])
    names = {'Sham','tVNS'};
    set(gca,'xtick',[1:5],'xticklabel',names)
    
    %% Check the significance in the cluster
    fprintf('Significance in cluster\n')
    [p, h, stats] = signrank(cluster_value_persub(:,1),cluster_value_persub(:,2));
    fprintf('\tWilcoxon: z = %.2f, p = %.2f\n',stats.zval,p)

    [h, p, ~, stats] = ttest(cluster_value_persub(:,1),cluster_value_persub(:,2));
    fprintf('\tt-test: t(%d) = %.2f, p = %.2f\n',stats.df,stats.tstat,p)
    
%     save('single_sub_eeg_eff_cluster','cluster_value_persub'); % percent per subject values
    
end


gnu = 1;


    function [S] = getS(SCode,dir)
        nS = length(SCode);
        for iVP = 1:nS
            
            S(iVP).index = SCode(iVP);
            S(iVP).dir = dir.dir;
            S(iVP).dir_data = dir.dir_data;
            S(iVP).EEGdir = dir.EEGdir;
            S(iVP).TFAdir = dir.TFAdir;
            
            S(iVP).code = sprintf('P%03d',SCode(iVP));
            S(iVP).EEGfn = sprintf('P%03d',SCode(iVP));
            
            S(iVP).suffix = '';
        end
        
    end




end