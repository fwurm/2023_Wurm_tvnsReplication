%% this is an improved final prepocessing for the short pulses tvns data
home_folder='C:\tvns_pilot\short_pulses_experiment';

cd(home_folder)
content=dir('*18');
sub_names={content.name};

% occipital_chan={'E82','E74','E119','E118','E117','E125','E126','E127','E137','E138','E139'}; % currently occiptal electrodes +2 for tvns
% central_chan= {'E190','E189','E178','E81','E45','VREF','E9','E186','E132','E131','E44','E17','E8','E198','E185','E144','E131','E90','E80','E53'};

nt=11; % trl in one block

% load excel from google drive
DOCID='1ux1ZIx27yot0tFWSuqnJmQztaWqhSrT2te1Au59f3AE';
subjects_data = GetGoogleSpreadsheet(DOCID);

sub_code_ind=find(cellfun(@(x) strcmp(x,'subject code'),subjects_data(:,1))==1);
which_is_first_ind = find(cellfun(@(x) strcmp(x,'which is first'),subjects_data(:,1))==1);
exclude_trl_ind = find(cellfun(@(x) strcmp(x,'exclude trials'),subjects_data(:,1))==1);
tvns_trl_ind = find(cellfun(@(x) strcmp(x,'tvns trials'),subjects_data(:,1))==1);
sham_trl_ind = find(cellfun(@(x) strcmp(x,'sham trials'),subjects_data(:,1))==1);

which_is_first=subjects_data(which_is_first_ind,:);
which_is_first{4}='TVNS'; which_is_first{15}='SHAM';
sub_code=subjects_data(sub_code_ind,:);
sheet_sub=cellfun(@(x) x(1:3),sub_names,'un',0);

% load channel maps and create neighbours map accordingly
cfg=[];
cfg.method='distance';
cfg.neighbourdis=3;
cfg.elecfile='C:\tvns_pilot\EGI_256.sfp';
neighbours=ft_prepare_neighbours(cfg);

load(['C:\tvns_pilot\' 'layout_193.mat']) % skull layout, var name is layout193

%% these paramas are specific for this run
sub_nums=[1:19 21:25];  %20 is RS9- missing EEG file
target_folder='eeg_all_chan_algoclean\';

target_folder_new='eeg_all_chan_algo100\';
%% this part is intented to create eeg_all_trl/mat per subject
errors={};
all_chan = 'all';
ref_chan = 'all'; % av refrence
resample_on=0;
for sub=1:length(sub_nums)
    try
        cd(home_folder)
        % extract the trials in occiptal channels and save them
        single_sub_eeg_tvns_pulses(target_folder,sub_names{sub_nums(sub)},resample_on); 
    catch ME
        errors{sub}=ME.message;
    end
end
% this pae
%% select relevant trials
for sub = 1:length(sub_nums)
    [x1, x2, where_is_in_excel]=intersect(sheet_sub{sub_nums(sub)},sub_code);
    exclude_trl{sub_nums(sub)} = sscanf(subjects_data{exclude_trl_ind,where_is_in_excel},'%f');
    sham_trl{sub_nums(sub)}    = sscanf(subjects_data{sham_trl_ind,where_is_in_excel},'%f');
    tvns_trl{sub_nums(sub)}    = sscanf(subjects_data{tvns_trl_ind,where_is_in_excel},'%f');
    % verfications:
    check_len=length(exclude_trl{sub_nums(sub)}) + length(sham_trl{sub_nums(sub)}) + length(tvns_trl{sub_nums(sub)});
    missing_trl=setxor([exclude_trl{sub_nums(sub)}; sham_trl{sub_nums(sub)}; tvns_trl{sub_nums(sub)}],1:88);
    intersect(exclude_trl{sub_nums(sub)},tvns_trl{sub_nums(sub)})
    intersect(exclude_trl{sub_nums(sub)},sham_trl{sub_nums(sub)})
    intersect(tvns_trl{sub_nums(sub)},sham_trl{sub_nums(sub)})
end

%%  exclude bad trials - this is according to the excel (manual inspection chekcing if there was stimulation).
sub_nums=[3]; %2:5 8:19 22:24]; %[1:19 21:26];  %20 is RS9- missing EEG file

for sub =1:length(sub_nums)
    cd('C:\tvns_pilot\short_pulses_experiment\eeg_allchan_trials')
    load([sub_names{sub_nums(sub)} '_eeg_all_trl.mat'])
    load('tvnstrlinfo.mat')
    load('shamtrlinfo.mat')
    % before removing trials and changing the order, get the sampleinfo of
    % the sham/tvns trials
    sham_trlinfo{sub_nums(sub)}= data_trl.sampleinfo(sham_trl{sub_nums(sub)},:);
    tvns_trlinfo{sub_nums(sub)} =  data_trl.sampleinfo(tvns_trl{sub_nums(sub)},:);
    
    cfg=[];
    cfg.trials= setxor(exclude_trl{sub_nums(sub)},1:length(data_trl.trial)); % manually select tvns/sham trials
    data_trl=ft_selectdata(cfg,data_trl);
    save([sub_names{sub_nums(sub)} '_eeg_all_trl_cln'],'data_trl','-v7.3')
    save('shamtrlinfo','sham_trlinfo')
    save('tvnstrlinfo','tvns_trlinfo')
end

%% extra preprocessing steps I added in July 2019
for sub =2:length(sub_nums)
    disp(sub_names{sub_nums(sub)})
    load(['C:\tvns_pilot\short_pulses_experiment\eeg_allchan_trials\' sub_names{sub_nums(sub)} '_eeg_all_trl_cln'],'data_trl')
    % reduce data to skull channels only :
    cfg=[];
    cfg.channel = layout193.label;
    data = ft_selectdata(cfg,data_trl);
    clear data_trl

    %  filter out the stimulation frequency in 25Hz and its harmonics:
    cfg=[];
    cfg.bsfilter      = 'yes';
    notch_f=25:25:475;
    cfg.bsfreq=[notch_f'-3 notch_f'+3];
    cfg.bsfiltord = 3;
    cfg.bpinstabilityfix = 'reduce';
    cfg.bpfilter= 'yes';
    cfg.bpfreq = [5 15]; 
    data_dft=ft_preprocessing(cfg,data);
    
    % use  to clean using naive cleaning used in break data
%     MAXABS_THRESHOLD=150; % uV
    MAXABS_THRESHOLD=100; % uV % reduced for the review purposes
    max_abs_values=cell2mat(cellfun(@(x) max(abs(x),[],2),data_dft.trial,'un',0)); % max_ abs values chan*trl
    unvalid=max_abs_values>MAXABS_THRESHOLD; % this a 1/0 matrix of the currenly unvalid chan*trl combinations
    view_on=0;
    [exc_chan, exc_trl] = naive_unvalid_removal(unvalid,view_on,'chan');
        
    % interpolate bad channels and reject bad _trials:
    cfg=[];
    cfg.trials=setxor(1:1:length(data_dft.trial),exc_trl);
    cfg.neighbours=neighbours;
    cfg.method='spline';
    cfg.badchannel=data_dft.label(exc_chan);
    cfg.layout = layout193;
    data_intrp=ft_channelrepair(cfg,data_dft);
    save(['C:\tvns_pilot\short_pulses_experiment\' target_folder_new  sub_names{sub_nums(sub)} '_eeg_all_trl_intrp_wide_100'],'data_intrp','exc_chan','exc_trl','-v7.3')
end

%% divide to sham and real using trlinfo:
data_cond=cell(2,1);
load(['C:\tvns_pilot\short_pulses_experiment\eeg_allchan_trials\shamtrlinfo.mat'])
load(['C:\tvns_pilot\short_pulses_experiment\eeg_allchan_trials\tvnstrlinfo.mat'])

b=cell(2,1);
for sub=1:length(sub_nums)
    load(['C:\tvns_pilot\short_pulses_experiment\' target_folder_new sub_names{sub_nums(sub)} '_eeg_all_trl_intrp_wide_100'])
    trlinfos={tvns_trlinfo{sub_nums(sub)}; % tvns
              sham_trlinfo{sub_nums(sub)}};  %sham
    for cond=1:length(trlinfos) % (1) tvns (2) sham
        cfg=[];
        [x,a,b{cond}]=intersect(trlinfos{cond}(:,1), data_intrp.sampleinfo(:,1));
        cfg.trials=b{cond};
        data_cond{cond}=ft_selectdata(cfg,data_intrp);
    end
    save(['C:\tvns_pilot\short_pulses_experiment\' target_folder_new sub_names{sub_nums(sub)} '_eeg_all_cond_wide_100.mat'],'data_cond','-v7.3')
end
%% average ref
% now that data_cond is clean I can make average ref
for sub=1:length(sub_nums)
   load([target_folder_new sub_names{sub_nums(sub)} '_eeg_all_cond_wide_100.mat'])
   cfg=[];
   cfg.reref='yes';
   cfg.refchannel = 'all';
   cfg.refmethod=   'avg';
   data_cond_avref=cellfun(@(x) ft_preprocessing(cfg,x),data_cond);
   save([target_folder_new sub_names{sub_nums(sub)} '_eeg_all_cond_avref_wide_100.mat'],'data_cond','-v7.3')
end
