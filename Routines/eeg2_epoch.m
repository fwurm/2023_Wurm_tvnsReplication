function eeg2_epoch_peeg(S,targetEventI,epoch1,varargin)

warning on

options.epoch = epoch1;
options.baseline = [];
options.targetEvent = targetEventI;

options.inprefix = '';
options.outprefix = '';

for i = 1:length(varargin)
    if strcmp(varargin{i}, 'in')
        if length(varargin)>i
            options.inprefix = varargin{i+1};
        else
            disp('ERROR: Input for parameter ''in'' is not valid!');
        end
    end
    if strcmp(varargin{i}, 'out')
        if length(varargin)>i
            options.outprefix = varargin{i+1};
        else
            disp('ERROR: Input for parameter ''out'' is not valid!');
        end
    end
    if strcmp(varargin{i}, 'baseline')
        if length(varargin)>i
            options.baseline = varargin{i+1};
            options.baseline_ind = 1;
        else
            disp('ERROR: Input for parameter ''baseline'' is not valid!');
        end
    end
end

%%%%%%%%%%%%%
% Portcodes %
%%%%%%%%%%%%%
trigger.stimon = [33];
trigger.stimoff = [34];
trigger.condtype1 = [31]; %condition types
trigger.condtype2 = [32];
trigger.breakon = [45];
trigger.breakoff = [46];
options.trigger = trigger;

if length(targetEventI)==1 && ~isempty(epoch1)
    options.epochtype = 'stimLocked';
elseif length(targetEventI)==2 && ~isempty(epoch1)
    options.epochtype = 'breakLocked';
else
    error('no correct epoching type can be specified\ncheck input parameters')
end
    

% information about the experiment
experiment.nBlock = 8;
experiment.allTrial = 88;
experiment.nTrial = 11;
options.experiment = experiment;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%      Step 2: epoching      %');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

fprintf('   Time-Locking Event: %s\n', options.targetEvent);
fprintf('   options.baseline: %s\n', mat2str(options.baseline));
fprintf('   options.epoch: %s\n', mat2str(options.epoch));


for i = 1:length(S)
    subfun(S(i),options);
end

function subfun(S,options)

if exist(S.historyFile)
    load(S.historyFile,'INFO'); %load file to note down rejected trials (visual inspection of tvns artefact)
else
    INFO = struct();
    save(S.historyFile,'INFO')
    fprintf('creating history file...')
end

cS = S.index; %current subject
fprintf('=> Subject: %d\n', cS);
INFO(cS).subject = cS;

%% epoching

% import EEG
inpath = [S.EEGdir options.inprefix ''];
raw = pop_loadset('filename',[S.EEGfn '.set'],'filepath',inpath);

events = [raw.event];

if isstr(events(1).type)
    eventnums = cell2mat(cellfun(@(x) str2num(erase(x,'condition ')),{events.type},'UniformOutput',false));
    
else
    eventnums = [events.type];
end

INFO(S.index).imported = length(eventnums);

    
%% removing outliers (trigger offset due to external buttons)

external_trigger = 8192;

newenums = rem(eventnums,external_trigger);
if ~isequal(rem(newenums,external_trigger),eventnums)
    warning('new outlier detected')
end
eventnums = newenums;


%remove distracting triggers (from the triggerbox on the seat)
events(find(eventnums==0)) = [];
eventnums(find(eventnums==0)) = [];

triggerevents = find(eventnums==33);
lats = [events(triggerevents).latency];
difflats = diff(lats);
latsident = find(difflats < 10);
events(triggerevents(latsident)) = [];
eventnums(triggerevents(latsident)) = [];

INFO(S.index).outlierRem = length(eventnums);

if (S.index == 8)
    %insert missing condition trigger
    events(47:end+1) = events(46:end);
    events(46) = struct('type','32','edftype',[],'latency',350000,'urevent',16000);
    
    eventnums(47:end+1) = eventnums(46:end);
    eventnums(46) = 32;
end

if (S.index == 16)
   % delete unreliable events
   events(50:57) = []; 
   eventnums(50:57) = []; 
end

if (S.index == 18)
    %insert missing condition trigger
    events(2:end+1) = events(1:end);
    events(1) = struct('type',31,'latency',1,'urevent',1);
    
    eventnums(2:end+1) = eventnums(1:end);
    eventnums(1) = 31;
end






%% find blockstarts
blockstart = find(ismember(eventnums,[31 32])==1);
trialidx = 0;
if (length(blockstart)==8)
    blockur = [events(blockstart).urevent];
    trialidx = 1;
elseif (length(blockstart)>8)
    blockdiff = diff(blockstart);
    blockident = find(blockdiff < 10);
    blockstart(blockident) = [];
    blockur = [events(blockstart).urevent];
    trialidx = 1;
else
    blocklat = [events(blockstart).latency];
    blockdiff = diff(blocklat);
    warning('fewere blocks detected!!!')
end

if (S.index == 35)
    blockur = [events(blockstart).urevent];
end

%% recoding of triggers (each experiment starts with 31 first)
trigger = options.trigger;
vptype = rem(S.index,2);
if (vptype == 1)
    blocktypes = [1 1 0 0 1 1 0 0];
    trigger_tvns = trigger.condtype1;
    trigger_sham = trigger.condtype2;
elseif (vptype == 0)
    blocktypes = [0 0 1 1 0 0 1 1];
    trigger_tvns = trigger.condtype2;
    trigger_sham = trigger.condtype1;
end
    
if (S.index == 9) || (S.index == 19) %reversed conditions
    blocktypes = [0 0 1 1 0 0 1 1];
    trigger_tvns = trigger.condtype2;
    trigger_sham = trigger.condtype1;
end

if (S.index == 22) || (S.index == 24) %reversed conditions
    blocktypes = [1 1 0 0 1 1 0 0];
    trigger_tvns = trigger.condtype1;
    trigger_sham = trigger.condtype2;
end

if (S.index == 8)
    blocktypes = [0 0 1 1 1 1 0 0];
end

%% get trials (stimlocked or breaklocked)
stimtrigger = [];
condtype = nan;
condcount = 0;
trialcount = 0;

if strcmp(options.epochtype,'stimLocked')
    for iE = 1:length(events)-1 
        
        
        if ismember(events(iE).urevent,blockur)
            blocknum = find(blockur==events(iE).urevent);
            condtype = blocktypes(blocknum);
        end
          
        if eventnums(iE) == 33 %& eventnums(iE+1) == 34
            trialcount = trialcount+1;
            stimtrigger(trialcount).type = eventnums(iE);
            stimtrigger(trialcount).urevent = events(iE).urevent;
            stimtrigger(trialcount).latency =events(iE).latency;
            stimtrigger(trialcount).blocknum = blocknum; %condcount;
            stimtrigger(trialcount).trialnum = trialcount;
            stimtrigger(trialcount).condition = condtype;
        end
    end    

    INFO(S.index).B_counted = length(blockstart);       
    INFO(S.index).E_counted = length(stimtrigger);
       
    % identify events that are temporally (too) close
    timediff = diff([stimtrigger.latency]);
    tooshort = find(timediff<100); %sampling points?
    if ~isempty(tooshort)
        stimtrigger(tooshort+1) = []; %remove duplicate trials
    end
    
    
    
    raw.event = stimtrigger;
    INFO(S.index).E_duplicateRem = length(stimtrigger);
    
    save(S.historyFile,'INFO')
    
elseif strcmp(options.epochtype,'breakLocked')

    for iE = 1:length(events)-1        
        if eventnums(iE) == options.targetEvent(1)
            trialcount = trialcount+1;
            breakLatency(trialcount,1) = events(iE).latency;
        end        
    end
    breakLatency(find(diff(breakLatency)<10*raw.srate))=[]; %clear duplicate events (closer than 10 seconds) 
    
    breakends = find(eventnums==options.targetEvent(2));
    breakends(find(diff(breakends)<10))=[]; %clear duplicate events (more than 10 events inbetween)
    if length(breakends)==length(breakLatency)
        breakLatency(:,2) = [events(breakends).latency];
    else
        error('break starts and ends still mismatch')
    end
    
    fprintf('  total break duration: %.2fs\n',sum(diff(breakLatency,[],2)./512))
    
    load('PEEG_BreakInfo.mat','breakduration','breaknum')
    breakduration(S.index) = sum(diff(breakLatency,[],2)./512);
    breaknum(S.index) = size(breakLatency,1);
    save('PEEG_BreakInfo.mat','breakduration','breaknum')
    
    breaktrigger = [];
    breaktrial = 0;
    for iB = 1:size(breakLatency,1)
        currLat = breakLatency(iB,1);
        while (currLat < breakLatency(iB,2)-options.epoch(2)*raw.srate)
            breaktrial = breaktrial + 1;
            
            breaktrigger(breaktrial).type = 1;
            breaktrigger(breaktrial).urevent = breaktrial;
            breaktrigger(breaktrial).latency = currLat + (-1*options.epoch(1)*raw.srate); %
            breaktrigger(breaktrial).breaknum = iB;
            breaktrigger(breaktrial).trialnum = breaktrial;
            
            
            currLat = currLat+(sum(abs(options.epoch(1:2)))-options.epoch(3))*raw.srate;
        end
    end
    
    raw.event = breaktrigger;
    
%     EEG = eeg_regepochs(raw, 'recurrence',4, 'eventtype','3','extractepochs','off');
    
else
    error()
    
end


data = pop_epoch(raw,{},options.epoch(1:2), 'newname', [options.outprefix S.EEGfn '_epochs'], 'epochinfo', 'yes');

if strcmp(options.epochtype,'breakLocked')
    %delete duplicate events (end of each trial as they are overlapping)
    [~,ia,~] = unique([data.event.urevent],'last');
    data.event = data.event(ia);
end


%% save data
outpath = [S.EEGdir options.outprefix '\'];
if ~exist(outpath, 'dir')
    mkdir(outpath)
end
finaldata = pop_saveset( data, 'filename' , [S.EEGfn '.set'],'filepath',outpath);



























