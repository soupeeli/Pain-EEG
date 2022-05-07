function [plv] = eegPLV(eegData,toi)
% Computes the Phase Locking Value (PLV) for an EEG dataset.
% 
%
% Input parameters:
%   eegData is a 3D matrix numChannels x numTimePoints x numTrials
%   srate is the sampling rate of the EEG data
%   toi is the time of interest for analysis. form e.g.: [300 100] for
%   300-1000 points of the eeg sequences.
%
% Output parameters:
%   plv is a 2D matrix - 
%     numChannels x numChannels
%--------------------------------------------------------------------------
% Example: Consider a 28 channel EEG data sampled @ 1000 Hz with 231 trials,
% where each trial lasts for 1 seconds. Below is an example of how to
% use this function.
%
%   eegData = rand(28, 1000, 231);
%   toi = [500 1000]; % only calculating the 500-1000 points
%   [plv] = eegPLV(eegData, toi);
% 
%   NOTE:
% 
% Also note that in order to extract the PLV between channels 17 and 20, 
% use plv(:, 17, 20, :) and NOT plv(:, 20, 17, :). The smaller channel 
% number is to be used first.
%--------------------------------------------------------------------------
% 
% Reference:
%   Lachaux, J P, E Rodriguez, J Martinerie, and F J Varela. 
%   Measuring phase synchrony in brain signals.? 
%   Human brain mapping 8, no. 4 (January 1999): 194-208. 
%   http://www.ncbi.nlm.nih.gov/pubmed/10619414.
% 
%--------------------------------------------------------------------------
% Written by: 
% Yamin Li
% Neuroengineering Lab, SJTU
% 01 Dec 2020



numChannels = size(eegData, 1);
numTrials = size(eegData, 3);
filteredData = eegData;

% disp(['Calculating PLV for ' mat2str(numTrials) ' trials...']);
for channelCount = 1:numChannels
    filteredData_com(channelCount, :, :) = angle(hilbert(squeeze(filteredData(channelCount, :, :))));
end

filteredData_com = filteredData_com(:,toi,:);


% plv = zeros(size(filteredData, 2), numChannels, numChannels, numConditions);
for channelCount = 1:numChannels
    channelData = squeeze(filteredData_com(channelCount, :, :));
    for compareChannelCount = 1:numChannels
        compareChannelData = squeeze(filteredData_com(compareChannelCount, :, :));
        phasediff = channelData - compareChannelData;
        plv(channelCount, compareChannelCount) = sum(abs(mean(exp(1i*phasediff), 1))) /numTrials; % checked: no problem 2021/10/13
    end
end
plv = squeeze(plv);

return;