%% Peak latency topographic plot
% In this exercise we produce the topographic plot of an ERP's peak 
% latency (the time elapsed to reach the peak) on a specific time 
% interval, for both filtered and unfiltered erp
%% prepare data

addpath('utils');
load data/sampleEEGdata.mat

time_period = [100, 400];

%% make a low-pass filter (windowed sinc function)
lowcut = 15;
filttime = -.3:1/EEG.srate:.3;
filtkern = sin(2*pi*lowcut*filttime) ./ filttime;
% adjust NaN and normalize filter to unit-gain
filtkern(~isfinite(filtkern)) = max(filtkern);
filtkern = filtkern./sum(filtkern);
% windowed sinc filter
filtkern = filtkern .* hann(length(filttime))';

%% make plot

% plot the unfiltered peak times
figure(1), clf
subplot(121)
make_peak_times_topoplot(EEG, time_period);
title({'Unfiltered ERP peak times';[' (' num2str(time_period(1)) '-' num2str(time_period(2)) ' ms)' ]})
set(gca,'clim',time_period)
colormap hot
colorbar

% plot the filtered peak times
subplot(122)
make_peak_times_topoplot(EEG, time_period, filtkern);
title({'Filtered ERP peak times';[ ' (' num2str(time_period(1)) '-' num2str(time_period(2)) ' ms)' ]})
set(gca,'clim',time_period)
colormap hot
colorbar

%% function that do the actual work
function make_peak_times_topoplot(EEG, time_period, varargin)
    % process optional variables
    opt_vars = inputParser;
    addOptional(opt_vars, 'filt_kern', []); 
    parse(opt_vars, varargin{:});

    peak_times = zeros(1, EEG.nbchan);
    % convert time period to indices
    time_idx = dsearchn(EEG.times', time_period')';
    
    % loop over channels
    for chani=1:EEG.nbchan
        % for each channel compute erp averaging over trials
        erp = double( mean(EEG.data(chani, :, :), 3) );
        % if a filter has been provided filter the erp
        if ~isempty(opt_vars.Results.filt_kern)
            erp = filtfilt(opt_vars.Results.filt_kern, 1, erp);
        end
        % find erp peak time index in specific time interval 
        [~,peak_time_idx] = max(erp(time_idx(1):time_idx(2)));
        % adjust the index from smaller erp selection to full erp
        peak_time_idx = peak_time_idx + time_idx(1) - 1;
        % save it in peak times list
        peak_times(chani) = EEG.times(peak_time_idx);
    end
    
    % plot peak times
    topoplotIndie(peak_times, EEG.chanlocs, 'numcontour', 4, 'electrodes', 'numbers');
end
%%
% Conclusion
% You may notice some differences between the same channel on the filtered 
% and unfiltered version. That is because the filtering operation may
% change the global maxima of the erp thus changing the peak time, so 
% changing the peak delay and the plot. This is not desirable for your
% analysis (operation like filtering should not affect your analysis that
% much) so it may be helpful in this cases to pick another time interval
%%



