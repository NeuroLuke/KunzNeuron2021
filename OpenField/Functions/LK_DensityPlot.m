function [out] = LK_DensityPlot(time, spikes, idx, donorm)
%
% LK_DensityPlot plots spike waveforms as a density plot.
%
% Input:
%   time                --> time points of the waveform data points (x-axis)
%   spikes              --> spike waveforms
%   idx (optional)      --> which spikes to plot
%   donorm (optional)   --> whether to normalize the density plot
%
% Adapted from Reber et al., PLOS Biology, 2019.
%
% Lukas Kunz, 2021

% waveforms for plotting
if ~exist('idx', 'var')
    wavs = spikes;
else
    % if you want to plot specific waveforms
    wavs = spikes(idx, :);
end

% whether to do a normalization
if ~exist('donorm', 'var')
    donorm = true;
end

% figure bounds
lbound  = floor(min(min(wavs))); % lower bound of the figure
ubound  = ceil(max(max(wavs))); % upper bound of the figure

% settings for the 2D histogram
numXbins    = size(spikes, 2);
ybins       = linspace(lbound, ubound, 150);
ybinSize    = ybins(2) - ybins(1);
numYbins    = length(ybins);

% preallocate 2D histogram that contains bin counts
n   = zeros(numXbins, numYbins); % spike-time * microvolt
% loop through spike-time bins
for k = 1:numXbins
    % loop through voltage bins
    for j = 1:numYbins
        n(k, j) = sum(wavs(:, k) <= ybins(j) + ybinSize / 2 & wavs(:, k) > ybins(j) - ybinSize / 2);
    end
end

% normalization by highest count
if donorm
    n = n ./ max(n(:));
end

% clip extreme outliers
cutoff          = 5 * std(n(:));
n(n > cutoff)   = cutoff;

% plot the 2D histogram
pcolor(time, ybins, n');
shading interp;

%% output
out             = [];
out.lbound      = lbound;
out.ubound      = ubound;
