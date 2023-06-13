function [exc_chan, exc_trl] = minimal_automatic_rejection(unvalid,view_on, prefer)
% This a simpe script to itertively removes channels and trials from an EEG
% dataset. In each iteration a the most noisy channel/trial is removed.
%input:
%unvalid - a matrix of booolens in a size of channels (row) and trial(col)
% view on - 0/1 if equals 1 shows the iterative process visually using imagesc.
% prefer - 'chan; /'trl'indicates which dimension channel/trials to remove in case
% maximam number of noisy pieces is eqaul in a channel and a trial
% output:
% exc_chan: indexes of the channels to remove/interpolate.
% exce_trl: indexes of the trials to reject
% @Omer Sharon omerxsharon@gmail.com
exc_chan=[];
exc_trl=[];
iteration=0;
% prefer: indicates if when chan removal and trial removal are qual which will be removed
while sum(sum(unvalid))>0 % remove until there no more 
    % sum the unvalid points in each dimesion
    sum_trl=sum(unvalid,1); 
    sum_chan=sum(unvalid,2);
    if strcmp(prefer,'chan') % prefer to remove channels and spare trials - this is better for trialed data
       condition =  max(sum_chan)>= max(sum_trl);
    else if strcmp(prefer,'trl') % prefer to remove trials and spare channels
       condition =  max(sum_chan)> max(sum_trl);
        else
            disp('decide: chan or trl')
        end
    end
    % each iteration choose the best dimesnsion to remove from:
    if condition % remove a channel is prefferred     
        [val, this_exc_chan]=max(sum_chan);
        exc_chan=[exc_chan this_exc_chan]; % select the channel to output
        unvalid(this_exc_chan,:)=0; % zeros this channel
    else % remove a trl
        [val, this_exc_trl]=max(sum_trl); % where is the max
        exc_trl = [exc_trl this_exc_trl]; % add the selected trl  to output
        unvalid(:,this_exc_trl)=0; % zero this trial
    end
    if view_on
    imagesc(unvalid)
    title(iteration)
    pause(1)
    iteration=iteration+1;
    end
end
