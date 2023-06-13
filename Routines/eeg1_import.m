function eeg1_import(S,chanFileFolder,varargin)
% eeg1_import() - Imports raw eeg data and does some basic preparatory stuff, 
%                 such as filtering. 
%
% Usage: 
%   >> eeg1_import(S,channelFileFolder,optionalParameters)
% 
%   Mandatory Inputs:
%           S           - Information on the subject sample (a structure
%                         that has to be created in the batch file)
%           channelFileFolder - A folder has to be specified that contains
%                         a copy of the .elp channel file for each of the
%                         subjects (chanfile - Kopie (1).elp, 
%                                   chanfile - Kopie (2).elp,
%                                   chanfile - Kopie (3).elp, ...)                    
%
%   Optional Inputs:
%           'lowpass'   - Applies a lowpass filter to the eeg data that
%                         cuts off all frequencies above the specified
%                         value (e.g., 40 Hz)
%           'highpass'  - Applies a highpass filter to the eeg data that
%                         cuts off all frequencies below the specified
%                         value (e.g., 0.1 Hz)
%           'filterMethod' - The filter methods for lowpass and highpass
%                         filters can be changed from the default value
%                         'firls' to 'fir1'. EEGLAB help recommendation: 
%                         "We recommend that you use fir1, which yields
%                         larger attenuation"
%           'referenceElectrodes'  - You can change the electrodes for 
%                         offline rereferencing (Default: [69 70], linked
%                         mastoids in the BioSemi ActiveTwo system)
%           'resample'  - Change the sampling rate of the eeg data (e.g. to
%                         500 Hz)
%           'relable'   - Relable individual electrodes; channel numbers 
%                         and (new) electrode names are specified as a
%                         cell with two rows and n columns (n=number of 
%                         electrodes): 
%                         {6,15,41,52; 'PO9','I1','PO10','I2'} 
%                         => Labels channel 6 as PO9, channel 15 as I1,...
%           'in'        - folder name and filename prefix for input eeg
%                         files to be loaded. Low dashes ("_") indicate 
%                         new folders ("\")
%           'out'       - folder name and filename prefix for output eeg
%                         files to be saved. Low dashes ("_") indicate 
%                         new folders ("\")
%
%   Example: 
%   >> eeg1_import(S(cS), 'C:\EEG\Channelfiles\', 'lowpass',40, 'highpass', 0.1,...
%           'filterMethod', 'fir1', 'referenceElectrodes', [65 66], 'resample', 500,...
%           'in', 'RawData_', 'out', 'Filtered01-40_' );
%                                         
% (c) 2017 by Robert Steinhauser
% (c) 2023 adapted by Franz Wurm

%% settings

inprefix = '';
outprefix = '';
lp = 0;
hp = 0;
resample = 0; % 0 = no resampling;
filterMethod = 'firls'; % default: fir1
referenceElectrodes = [33:34]; % default: reref to linked mastoids
relabel = {};
nochan = [];

for i = 1:length(varargin)
    if strcmp(varargin{i}, 'in') 
       if length(varargin)>i  
           inprefix = varargin{i+1}; 
       else
           disp('ERROR: Input for parameter ''in'' is not valid!');
       end
    end   
    if strcmp(varargin{i}, 'out') 
       if length(varargin)>i  
           outprefix = varargin{i+1}; 
       else
           disp('ERROR: Input for parameter ''out'' is not valid!');
       end
    end       
   if strcmp(varargin{i}, 'lowpass') 
       if length(varargin)>i && isnumeric(varargin{i+1}) 
           lp = varargin{i+1}; 
       else
           disp('ERROR: Input for parameter ''lowpass'' is not valid!');
       end
    end
    if strcmp(varargin{i}, 'highpass')
        if length(varargin)>i && isnumeric(varargin{i+1}) 
           hp = varargin{i+1};
        else
           disp('ERROR: Input for parameter ''highpass'' is not valid!'); 
        end
   
    end   
     if strcmp(varargin{i}, 'filterMethod')
           filterMethod = varargin{i+1};     
     end     
     if strcmp(varargin{i}, 'referenceElectrodes')
           referenceElectrodes = varargin{i+1};     
     end   
     
    if strcmp(varargin{i}, 'resample')
        if length(varargin)>i && isnumeric(varargin{i+1}) 
           resample = varargin{i+1};
        else
           disp('ERROR: Input for parameter ''resample'' is not valid!');
        end
    
    end   
    if strcmp(varargin{i}, 'relabel')
        if length(varargin)>i
            relabel = varargin{i+1};
        else
            disp('ERROR: Input for parameter ''relabel'' is not valid!');
        end
    end
    
    if strcmp(varargin{i}, 'nochan')
        if length(varargin)>i
            nochan = varargin{i+1};
        else
            disp('ERROR: Input for parameter ''nochan'' is not valid!');
        end
    end
    
end

%%%%%%%%%
% START %
%%%%%%%%%

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Step 1: Data Import and Filtering %');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');


aTot = cputime();

% parfor i = 1:length(S)
for i = 1:length(S)
    subfun(S(i),i,chanFileFolder,lp,hp,inprefix,outprefix,resample,filterMethod, referenceElectrodes,nochan,relabel);
end

diffTot = cputime()-aTot;

fprintf('   Total Time: %f ms\n',diffTot);



function subfun(S,index,chanFileFolder,lp,hp,inprefix,outprefix,resample,filterMethod,referenceElectrodes,nochan,relabel)


% Set of chan files necessary for parallel computing
currentChanFile = [chanFileFolder 'chanfile - Kopie (' num2str(index) ').elp'];

fprintf('=> Subject: %d START (LP= %d HP=%d)\n',S.index,lp,hp)

aVP = cputime();

%% STEP 1: import dataset
fn = [S.EEGdir inprefix S.EEGfn '.bdf'];
disp(['Load bdf file: ' fn]);
raw = pop_biosig([fn], 'ref', referenceElectrodes); 

if ~isempty(nochan)
    raw = pop_select(raw, 'nochannel',nochan);
end

if ~isempty(relabel)
    for i = 1:length(relabel{1})
        raw.chanlocs(relabel{1}(i)).labels = relabel{2}{i};
    end
end

raw = pop_chanedit(raw,  'lookup', currentChanFile);


%%  STEP 2: resample 

if resample > 0     
    raw = pop_resample(raw,resample);
end

%% STEP 3: filter 

if strcmp(filterMethod,'firls')
    if lp>0
        raw = pop_eegfilt( raw, 0, lp,[],0,0,0,'firls',0);        
    end
    if hp>0
        raw = pop_eegfilt( raw, hp, 0,[],0,0,0,'firls',0);
    end
    
elseif strcmp(filterMethod,'fir1')
    if lp>0
        raw = pop_eegfilt( raw, 0, lp,[],0,0,0,'fir1',0);
        
    end
    if hp>0
        raw = pop_eegfilt( raw, hp, 0,[],0,0,0,'fir1',0);
    end   
else
    error('No proper filter method selected!');
    
end


    
%% STEP 4: Save
outpath = [S.EEGdir outprefix '\'];
if ~exist(outpath, 'dir')
    mkdir(outpath)
end
raw = pop_saveset( raw, 'filename' , [S.EEGfn '.set'],'filepath',outpath);

    
diffVP = cputime()-aVP;

fprintf('         VP %d END (%s.set)=> %f min\n', S.index ,outpath,floor(((diffVP/1000)/60)))


