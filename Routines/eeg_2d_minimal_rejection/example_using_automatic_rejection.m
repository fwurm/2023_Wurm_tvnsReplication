
% This is an example sciprt for using minimal automatic rejection script
% I load an example dataset with the fieldtrip toolbox:
trigger_value=65402;
file= 'F:\MEGA\EEG course\fieldtrip Examples\ex1\hemiTrans1\001.bdf';

cfg=[];
cfg.dataset=file;
cfg.trialdef.eventtype = 'STATUS';
cfg.trialdef.eventvalue= 65402;
cfg.trialdef.prestim= 1; % in seconds
cfg.trialdef.poststim= 1;
cfg=ft_definetrial(cfg);
cfg.demean = 'yes'; % DC shift
cfg.bpfilter='yes'; % bandpass filter
cfg.bpfreq=[0.5 30];
cfg.bpinstabilityfix = 'reduce';
data=ft_preprocessing(cfg);

% use a absloute crietria to reject trials/channels:
MAXABS_THRESHOLD=100;
max_abs_values=cell2mat(cellfun(@(x) max(abs(x),[],2),data.trial,'un',0)); % max_ abs values chan*trl
unvalid=max_abs_values>MAXABS_THRESHOLD; % this a boolean matrix of the currenly unvalid chan*trl combinations

view_on=0;
[exc_chan, exc_trl] = minimal_automatic_rejection(unvalid,view_on,'trl');
    
%remove channels
cfg=[];
cfg.channel=setxor(1:size(data.label,1),exc_chan);  % channels to exclude
temp= ft_selectdata(cfg,data);

%remove trials
cfg=[];
cfg.trials=setxor(1:size(data.label,1),exc_trl);  % channels to exclude
data_clean= ft_selectdata(cfg,temp);


%% alternatively you can load the boolean matrix directly if you don't use fieltrip:
load('unvalid_example')
view_on=1; % this allows to see the iterative process
[exc_chan, exc_trl] = minimal_automatic_rejection(unvalid,view_on, 'trl');

%%  plot unvalid matrix and mark the removed channels and trials
imagesc(unvalid)
ylabel('Channels')
xlabel('Trials')
b=get(gca,'XTickLabel');
set(gca,'XTickLabel',b,'fontSize',14)
for chan=exc_chan
text(1,chan,'*','Color','r','fontSize',30)
end
for trl=exc_trl
text(trl,1,'*','Color','r','fontSize',30)
end