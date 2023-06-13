function batch_peeg_sharon
% Batch file for the replication of Sharon et al (2021) "Transcutaneous Vagus 
% Nerve Stimulation in Humans Induces Pupil Dilation and Attenuates Alpha 
% Oscillations" - pipeline scripts directly from Sharon
%
% adapted for revision

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


%%%%%%%%%%%%%%%%%%
%% sharon approach
%%%%%%%%%%%%%%%%%%

load('PEEG_ExperimentInfo_39participants.mat','INFO')
outpath = fullfile(fdir.EEGdir, 'stimLocked\e\sharon');
if ~exist(outpath, 'dir')
    mkdir(outpath)
end

sub_number = [6,7,8,9,10,11,13,14,15,16,17,19,20,24,25,26,27,29,31,32,33,35,36,39];
how_many_fac  = [1,2,2,2,2,2,1,1,1,2,2,3,2,1,3,1,2,2,1,1,2,2,1,1];
whichfac_alpha= [1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1];

Sidx = ismember(sub_number,[S.index]);
sub_number(~Sidx) = [];
how_many_fac(~Sidx) = [];
whichfac_alpha(~Sidx) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting the PARAFAC results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(outpath,'Trialdata_mean.mat'),'freq_axis');
load('chanlocs.mat','chanlocs')

model_content = dir([outpath '\*model.mat']);
model_cell = cell(size(model_content,1),2);

for isub=1:length(model_content)
    temp  = load([outpath '\' model_content(isub).name]);
    temp1=temp.model;
    for cond=1:2
    model_cell{isub,cond}=temp1{cond};  
    end
end

%%%%%
%% Parafac from Sharon
%%%%%

folder = 'C:\Users\wurmf\Documents\Github\PEEG_tVNS\sharonscripts\models';
model_content=dir(fullfile(folder, '*.mat'));
temp= struct2cell(model_content);
allnames=temp(1,:);
allnames3=cellfun(@(x) x(1:3),allnames,'un',0);

whichfac_alpha_sharon= [3,2,1,4,1,4,3,3,2,2,1,2,3,2,1,2,1,2,3,4,1];


for iS = 1:21
    load(fullfile(folder,allnames{iS}),'model');    
    freqcomps(iS,:) = model{2}{2}(:,whichfac_alpha_sharon(iS));   
end
avg_freq_model_sharon=squeeze(median(freqcomps,1));
sem_freq_model_sharon = squeeze(std(freqcomps,0,1))/sqrt(size(freqcomps,1)); % SEM of the model

avg_freq_model_sharon_norm = avg_freq_model_sharon./sum(avg_freq_model_sharon);
sem_freq_model_sharon_norm = squeeze(std(freqcomps./sum(avg_freq_model_sharon),0,1))/sqrt(size(freqcomps,1)); % SEM of the model


%% Mean alpha component (standard PARAFAC)

gold=(1+sqrt(5))/2;

method = 'standard'; %'cutoff', '', '', ''
calculate_results = 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% different methods to calculate grandaverage

switch method
    case 'standard'
        %% option 1
        % get the time- only data using the standard parafac model
        % This is the approach by Sharon
        if calculate_results            
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
            load(fullfile(outpath,'Trialdata.mat'),'C','time_axis','freq_axis');
        end
    case 'cutoff'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% option 2 (cutoff model)
        % get the time- only data using the cut-off parafac model
        % based on reviews 
        % apply threshold (mean) and consider only channels above threshold
        if calculate_results
            
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
                
                % select the alpha componenet
                comp=whichfac_alpha(find(sub_number==S(isub).index));
                
                
                for cond=1:2
                    
                    V_chan_cutoff = zeros(64,1);
                    V_freq_cutoff = zeros(31,1);
                    
                    fprintf('Condition %d\n',cond)
                    % reduce the channels domain:
                    chan_model=model{cond}{1}; % these are the channels from the model i.e from the break.
                    V_chan= chan_model(:,comp); % get the right comp
                    V_chan_cutoff(V_chan>mean(V_chan)) = 1;
                    fprintf('\tCHANNELS\n')
                    fprintf('\t mean PARAFAC weight: %.2f\n',mean(V_chan))
                    fprintf('\t channels above cutoff: %d\n',sum(V_chan_cutoff))
                    Y=permute(data{cond}.powspctrm,[1 3 4 2]); % move the channels dim to the end
                    B = sum(bsxfun(@times,Y,reshape(V_chan_cutoff,1,1,1,numel(V_chan))),4);
                    
                    B2 =reshape(B,size(B,1),1,size(B,2),size(B,3)); % trl, ~, freq, time
                    % reduce the freq domain:
                    freq_model =model{cond}{2};
                    V_freq = freq_model(:,comp);
                    V_freq_cutoff(V_freq>mean(V_freq)) = 1;
                    fprintf('\tFREQUENCIES\n')
                    fprintf('\t mean PARAFAC weight: %.2f\n',mean(V_freq))
                    fprintf('\t frequencies above cutoff: %d\n',sum(V_freq_cutoff))
                    X=permute(B2,[1 2 4 3]); % move the freq domain to the end
                    C{isub,cond} = squeeze(sum(bsxfun(@times,X,reshape(V_freq_cutoff,1,1,1,numel(V_freq))),4));
                    
    
                end
                
                avgchans(isub,:) = V_chan_cutoff;
                avgfreqs(isub,:) = V_freq_cutoff;
                
                gnu = 1 ;
            end
            
            time_axis = data{1}.time;
            freq_axis = data{1}.freq;
            save(fullfile(outpath,'Trialdata_cutoff.mat'),'C','time_axis','freq_axis','avgchans','avgfreqs');
        else
            load(fullfile(outpath,'Trialdata_cutoff.mat'),'C','time_axis','freq_axis','avgchans','avgfreqs');
        end
        
        %     figure
        allavgchans = mean(avgchans,1);
        chanidx = find(allavgchans>0.6)
        %     topoplotMS(allavgchans, chanlocs,'emarker2',{chanidx,'d','k',3},'maplimits',[-1 1],'whitebk','on','electrodes','off','headrad','rim');
        %     colormap('redblue')
        
        titletxt = 'Average alpha - PARAFAC cutoff';
        
        avg_freq_model=squeeze(median(avgfreqs,1));
        sem_freq_model = squeeze(std(avgfreqs,0,1))/sqrt(size(avgfreqs,1)); % SEM of the model
        
    case 'alternative_components'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% option 3 (alternative components model)
        % get the time- only data using the cut-off parafac model
        % based in reviews, use alternative components for ambiguous participants
                
        whichfac_alpha_alternat= [1,1,1,1,2,1,1,1,1,2,1,3,1,1,1,1,1,2,1,1,2,1,1,1,1];\
        
        if calculate_results
            
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
                
                % select the alpha componenet
                comp=whichfac_alpha_alternat(find(sub_number==S(isub).index));
                
                
                for cond=1:2
                    
                    
                    fprintf('Condition %d\n',cond)
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
                
                avgchans(isub,:) = V_chan;
                avgfreqs(isub,:) = V_freq;
                
                gnu = 1 ;
            end
            
            time_axis = data{1}.time;
            freq_axis = data{1}.freq;
            save(fullfile(outpath,'Trialdata_alternat.mat'),'C','time_axis','freq_axis','avgchans','avgfreqs');
        else
            load(fullfile(outpath,'Trialdata_alternat.mat'),'C','time_axis','freq_axis','avgchans','avgfreqs');
        end
        
        %     figure
        allavgchans = mean(avgchans,1);
        chanidx = find(allavgchans>mean(allavgchans))
        %     topoplotMS(allavgchans, chanlocs,'emarker2',{chanidx,'d','k',3},'maplimits',[-1 1],'whitebk','on','electrodes','off','headrad','rim');
        %     colormap('redblue')
        
        titletxt = 'Average alpha - alternative comps';
        
        avg_freq_model=squeeze(median(avgfreqs,1));
        sem_freq_model = squeeze(std(avgfreqs,0,1))/sqrt(size(avgfreqs,1)); % SEM of the model
        
    case 'alpha_predefined'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% option 4 (predefined alpha model)
        % get the time- only data using predefined ROIs for channels and freqs
        % based on reviews, ALPHA (8-12)
        
        
        channels = {'Oz','PO7','PO3','O1','O2','PO4','PO8', 'P9','P7','P5','P3','P1','Pz','P2','P4','P6','P8','P10','POz'};
        chanweight = double(ismember({chanlocs(1:64).labels},channels));
        
        freq_lims = [8 12];
        freqnums = dsearchn(freqs',freq_lims');
        freqweight = zeros(31,1);
        freqweight(freqnums(1):freqnums(2)) = 1;
        if calculate_results
            
            C=cell(length(S),2);
            for isub=1:length(S)
                fprintf('#########\n')
                fprintf('## %s \n',S(isub).EEGfn)
                fprintf('#########\n')
                
                %load trial data
                fn = fullfile(outpath,[S(isub).EEGfn '_fft.mat']);
                load(fn,'data');
                
                for cond=1:2
                    fprintf('Condition %d\n',cond)
                    % reduce the channels domain:
                    Y=permute(data{cond}.powspctrm,[1 3 4 2]); % move the channels dim to the end
                    B = sum(bsxfun(@times,Y,reshape(chanweight,1,1,1,numel(chanweight))),4);
                    B2 =reshape(B,size(B,1),1,size(B,2),size(B,3)); % trl, ~, freq, time
                    
                    % reduce the freq domain
                    X=permute(B2,[1 2 4 3]); % move the freq domain to the end
                    C{isub,cond} = squeeze(sum(bsxfun(@times,X,reshape(freqweight,1,1,1,numel(freqweight))),4));
                    
                    
                end
                
                gnu = 1 ;
            end
            
            time_axis = data{1}.time;
            freq_axis = data{1}.freq;
            save(fullfile(outpath,'Trialdata_predef_alpha.mat'),'C','time_axis','freq_axis');
        else
            load(fullfile(outpath,'Trialdata_predef_alpha.mat'),'C','time_axis','freq_axis');
        end
        
        allavgchans = chanweight;
        chanidx = find(chanweight==1);
        
        avg_freq_model=freqweight;
        
        titletxt = 'Average alpha - predefined ROI';

        
    case 'theta_predefined'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% option 5 (predefined theta model)
        % get the time- only data using predefined ROIs for channels and freqs
        % based on reviews, THETA (4-8)
        
        channels = {'Oz','PO7','PO3','O1','O2','PO4','PO8', 'P9','P7','P5','P3','P1','Pz','P2','P4','P6','P8','P10','POz'};
        chanweight = double(ismember({chanlocs(1:64).labels},channels));
        
        freq_lims = [4 8];
        freqnums = dsearchn(freqs',freq_lims');
        freqweight = zeros(31,1);
        freqweight(freqnums(1):freqnums(2)) = 1;
        
        if calculate_results
            C=cell(length(S),2);
            for isub=1:length(S)
                fprintf('#########\n')
                fprintf('## %s \n',S(isub).EEGfn)
                fprintf('#########\n')
                
                %load trial data
                fn = fullfile(outpath,[S(isub).EEGfn '_fft.mat']);
                load(fn,'data');
                
                for cond=1:2
                    fprintf('Condition %d\n',cond)
                    % reduce the channels domain:
                    Y=permute(data{cond}.powspctrm,[1 3 4 2]); % move the channels dim to the end
                    B = sum(bsxfun(@times,Y,reshape(chanweight,1,1,1,numel(chanweight))),4);
                    B2 =reshape(B,size(B,1),1,size(B,2),size(B,3)); % trl, ~, freq, time
                    
                    % reduce the freq domain
                    X=permute(B2,[1 2 4 3]); % move the freq domain to the end
                    C{isub,cond} = squeeze(sum(bsxfun(@times,X,reshape(freqweight,1,1,1,numel(freqweight))),4));
                    
                    
                end
                
                gnu = 1 ;
            end
            
            time_axis = data{1}.time;
            freq_axis = data{1}.freq;
            save(fullfile(outpath,'Trialdata_predef_theta.mat'),'C','time_axis','freq_axis');
        else
            
            load(fullfile(outpath,'Trialdata_predef_theta.mat'),'C','time_axis','freq_axis');
        end
        
        allavgchans = chanweight;
        chanidx = find(chanweight==1);
        
        avg_freq_model = freqweight;
        
        titletxt = 'Average theta - predefined ROI';
        
    otherwise
        error('no method selected')
end
        



%% Baselining (average level)

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
plot(time_axis,grand_G);legend('sham','tvns');
hold on
plot(time_axis,baseline_ind*10^17)
baseline_val= mean(G(:,baseline_ind,:),2);
grand_base_values=mean(baseline_val(:,:,2)-baseline_val(:,:,1),1);

%apply subject baseline
A=bsxfun(@(x,y) (rdivide(x,y)),G,baseline_val);


%% average across subjects:
grand_avg_per=squeeze(nanmean(A,1));
std_per=squeeze(nanstd(A,[],1))/sqrt(length(S)); % this is SEM
stim_lim =[0 4];


stim_ind = time_axis>stim_lim(1) & time_axis<stim_lim(2);
stim_values = nanmean(A(:,stim_ind,:),2);
values= stim_values*100;

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
report_sem=squeeze(std(values,1)./sqrt(length(S)));

%%%%%%%%%%%%%%%%%%%%%%%
%% plot results
%%%%%%%%%%%%%%%%%%%%%%

f2 = figure('Renderer','painters','Units','normalized')
f2.Position = [0.1 0.1 2/gold^2 1/gold];

%% plot grandaverage

ax = axes;
ax.Position = [0.1 0.2 1.3/gold^2 0.8/gold];

% plot the  grand average with std
colors='kr';
colors = [49, 116, 161; 225, 129, 43];
colors = [0, 0, 0; 255, 0, 0];
% for cond=1:2
% grand_avg_per_smooth(:,cond)=lowPassFilter(2,grand_avg_per(:,cond),20)+mean(grand_avg_per(:,cond));
% end

% trace_fig=figure('Renderer', 'painters', 'Position', [10 10 910 910/gold^2]);
for cond=1:2
    lx(cond) = shadedErrorBar(time_axis,grand_avg_per(:,cond)*100,std_per(:,cond)*100,{'color',colors(cond,:)./255,'markerfacecolor',colors(cond,:)./255,'linewidth',2},1);
    % ,colors(cond,:),0.2);
    hold on
    xlim([-1 6])
    ylim([0.85 1.25]*100)
    % ylim([0.9 1.9]*100)
end

set(gcf,'Color','w')
set(gca,'Fontsize',18,'fontName','Arial')
stim_patch=patch([0,0,3.4,3.4],[0.8,1.3,1.3,0.8]*100,'k','LineStyle','none');
alpha(stim_patch,0.1)
% colorbar
ylabel('Average alpha (%)')
xlabel('Time (sec)')
set(gca,'fontSize',18,'fontName','Arial')
line([-1 6],[100 100],'Color','k','LineWidth',2,'Linestyle','--')
line([0 0],[50 150],'Color','k','LineWidth',2,'Linestyle','--')
line([4 4],[50 150],'Color','k','LineWidth',2,'Linestyle','--')
legend([lx(1).mainLine lx(2).mainLine],{'sham' 'tvns'},'Box','off')
title(ax,titletxt)


%% plot channel weights
ax2 = axes;
ax2.Position = [0.85 0.46 0.3/gold^2 0.3/gold]; 
topoplotMS(allavgchans, chanlocs,'emarker2',{chanidx,'d','k',3},'maplimits',[-max(allavgchans) max(allavgchans)],'whitebk','on','electrodes','off','headrad','rim');
colormap('redblue')
title(ax2,'Topo mask','FontSize',18)


%% plot frequency weights
ax4 = axes;
ax4.Position = [0.85 0.2 0.3/gold^2 0.3/gold]; 
ax4.FontSize = 18;

if ~exist('sem_freq_model','var')
    plot(freq_axis,avg_freq_model,'k','LineWidth',2)
    ax4.YLim = [0 1.2];
    title('Freq mask','FontSize',18)
else
    shadedErrorBar(freq_axis,avg_freq_model,sem_freq_model,{'color','k','markerfacecolor','k','linewidth',1},1)
    title('Mean freq mask','FontSize',18)
end

    
%% plot single-subject data
ax3 = axes;
ax3.Position = [0.68 0.25 0.3/gold^2 0.6/gold]; 
jit = normrnd(0,0.1,[length(values),1])
    
color='kk';
kind={'-','-'};
for isub=1:length(values)
    col=(values(isub,1,1)-values(isub,1,2)>0)+1;
    plot([2 1]+jit(isub,1),[values(isub,1,1) values(isub,1,2)],kind{col},'lineWidth',1.5,'Color',[0 0 0]+(col-1)*0.5)
    hold on
    allcol(isub,1) = col-1;
end

s = scatter(2*ones(24,1)+jit,values(:,1,1))
s.MarkerEdgeColor = colors(1,:)./255;
s.MarkerFaceColor = colors(1,:)./255;

s2 = scatter(ones(24,1)+jit,values(:,1,2))
s2.MarkerEdgeColor = colors(2,:)./255;
s2.MarkerFaceColor = colors(2,:)./255;

set(gca,'Fontsize',18,'fontName','Arial')

xlim([0.6 2.4])
ylim(ax3,[70 150])

names = {'tVNS','sham'};

set(gca,'xtick',[1:5],'xticklabel',names,'ytick',[80 100 120 140])
ylabel('Power [%]')
% title('Mean alpha','Fontsize',18)
ax3.LineWidth = 2;
ax3.Box = 'off'

%% plot significance
ax5 = axes;
ax5.Position = [0.68 0.65 0.3/gold^2 0.2/gold]; 
[p, h, stats] = signrank(values(:,1,1),values(:,1,2));
text(0.5,0.5,sprintf('h = %d, p = %.3f',h,p),'HorizontalAlignment','center','FontSize',14)
ax5.XColor = 'none';
ax5.YColor = 'none';

f2.Color = 'w';



    function [S] = getS(SCode,dir)
        nS = length(SCode);
        for iVP = 1:nS
            
            S(iVP).index = SCode(iVP);
            S(iVP).dir = dir.dir;
            S(iVP).dir_data = dir.dir_data;
            S(iVP).EEGdir = dir.EEGdir;
            S(iVP).ERPdir = dir.ERPdir;
            S(iVP).TFAdir = dir.TFAdir;
            
            %             S(iVP).code = sprintf('Pilot%d',SCode(iVP));
            %             S(iVP).EEGfn = sprintf('Pilot%d',SCode(iVP));
            S(iVP).code = sprintf('P%03d',SCode(iVP));
            S(iVP).EEGfn = sprintf('P%03d',SCode(iVP));
            
            S(iVP).historyFile = dir.historyFile;
            
            S(iVP).suffix = '';
        end
        
    end




end