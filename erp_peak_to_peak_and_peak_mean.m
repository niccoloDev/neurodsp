%%
% This file contains an exercise to analyze an ERP (Event Related Potential)
% computing its peak-to-peak both using only the peak and the peak-mean.
% peak-to-peak: it's the voltage difference between the peak and the trough
%               in given intevals
% peak-mean: instead of taking the peak as a single value we compute it
% as an average over a window of time around the actual peak.
% We compute these values for both the unfiltered and filtered erp
%%

load data/sampleEEGdata.mat

% channel to pick
chan2use = 'o1';

% time window for negative peak (expressed as indices)
negpeaktime = dsearchn(EEG.times',[  50 110 ]')';
pospeaktime = dsearchn(EEG.times',[ 110 170 ]')';

%%% compute ERP
chanidx = strcmpi({EEG.chanlocs.labels}, chan2use);
erp = double( mean(EEG.data(chanidx,:,:),3) );

% plot ERP
figure(1), clf
plot(EEG.times,erp,'k','linew',1)
set(gca,'xlim',[-300 1000])

% plot patches over areas
ylim = get(gca,'ylim');
ph = patch(EEG.times(negpeaktime([1 1 2 2])),ylim([1 2 2 1]),'y');
set(ph,'facealpha',.8,'edgecolor','none')
ph2 = patch(EEG.times(pospeaktime([1 1 2 2])),ylim([1 2 2 1]),'g');
set(ph2,'facealpha',.8,'edgecolor','none')

% move the patches to the background
set(gca,'Children',flipud( get(gca,'Children') ))

%% make a low-pass filter (windowed sinc function)

lowcut = 15;
filttime = -.3:1/EEG.srate:.3;
filtkern = sin(2*pi*lowcut*filttime) ./ filttime;

% adjust NaN and normalize filter to unit-gain
filtkern(~isfinite(filtkern)) = max(filtkern);
filtkern = filtkern./sum(filtkern);

% windowed sinc filter
filtkern = filtkern .* hann(length(filttime))';

% inspect the filter kernel
figure(2), clf
subplot(211)
plot(filttime,filtkern,'k','linew',2)
xlabel('Time (s)')
title('Time domain')

subplot(212)
hz = linspace(0,EEG.srate,length(filtkern));
plot(hz,abs(fft(filtkern)),'ks-','linew',2)
set(gca,'xlim',[0 lowcut*3])
xlabel('Frequency (Hz)'), ylabel('Gain')
title('Frequency domain')

%% filter the ERP and replot

% apply filter
filt_erp = filtfilt(filtkern, 1, erp);

% plot on top of unfiltered ERP
hold on;
plot(EEG.times,filt_erp,'r','linew',1)

%% peak-to-peak voltages and timings

%%%% first for unfiltered ERP

% find minimum/maximum peak values and peak times
[min_volt,minidx] = min(erp(negpeaktime(1):negpeaktime(end)));
[max_volt,maxidx] = max(erp(pospeaktime(1):pospeaktime(end)));

% ERP timings
mint = EEG.times(negpeaktime(1)+minidx-1);
maxt = EEG.times(pospeaktime(1)+maxidx-1);

% get results (peak-to-peak voltage and latency)
erpP2P = max_volt - min_volt;
erpP2Plat = maxt - mint;


%%%% then for low-pass filtered ERP
% find minimum/maximum peak values and peak times
[min_volt,minidx] = min(filt_erp(negpeaktime(1):negpeaktime(end)));
[max_volt,maxidx] = max(filt_erp(pospeaktime(1):pospeaktime(end)));

% ERP timings
mint = EEG.times(negpeaktime(1)+minidx-1);
maxt = EEG.times(pospeaktime(1)+maxidx-1);

% get results (peak-to-peak voltage and latency)
erpFP2P = max_volt - min_volt;
erpFP2Plat = maxt - mint;

%% Report the results in the command window

% clear the screen
clc

fprintf('\nRESULTS FOR PEAK POINT:')
fprintf('\n   Peak-to-peak on unfiltered ERP: %5.4g muV, %4.3g ms span.',erpP2P,erpP2Plat)
fprintf('\n   Peak-to-peak on filtered ERP:   %5.4g muV, %4.3g ms span.\n\n',erpFP2P,erpFP2Plat)

%%

%% repeat for mean around the peak

% time window for averaging (one-sided!!)
win = 10; % in ms
% now convert to indices
win_idx = round( EEG.srate * (win / 1000) );

%%%% first for unfiltered ERP

% find minimum/maximum peak times
[~,minidx] = min(erp(negpeaktime(1):negpeaktime(end)));
[~,maxidx] = max(erp(pospeaktime(1):pospeaktime(end)));

% adjust ERP timings (readjust indices to be used on the original erp array)
minidx = minidx + (negpeaktime(1)-1);
maxidx = maxidx + (pospeaktime(1)-1);

% now find average values around the peak time
avg_min = mean(erp(minidx-win_idx:minidx + win_idx));
avg_max = mean(erp(maxidx-win_idx:maxidx + win_idx));

% ERP timings
mint = EEG.times(minidx); 
maxt = EEG.times(maxidx);

% get results (peak-to-peak voltage and latency)
erpP2P = avg_max - avg_min;
erpP2Plat = maxt - mint;

%%%% then for low-pass filtered ERP

% find minimum/maximum peak values and peak times
[min_volt,minidx] = min(filt_erp(negpeaktime(1):negpeaktime(end)));
[max_volt,maxidx] = max(filt_erp(pospeaktime(1):pospeaktime(end)));


% adjust ERP timings
minidx = minidx + (negpeaktime(1)-1);
maxidx = maxidx + (pospeaktime(1)-1);

% now find average values around the peak time
avg_min = mean(filt_erp(minidx-win_idx:minidx + win_idx));
avg_max = mean(filt_erp(maxidx-win_idx:maxidx + win_idx));

% adjust ERP timings
mint = EEG.times(minidx); 
maxt = EEG.times(maxidx);

% get results (peak-to-peak voltage and latency)
erpFP2P = avg_max - avg_min;
erpFP2Plat = maxt - mint;

%% Report the results in the command window

fprintf('\nRESULTS FOR WINDOW AROUND PEAK:')
fprintf('\n   Peak-to-peak using mean-peak on unfiltered ERP: %5.4g muV, %4.3g ms span.',erpP2P,erpP2Plat)
fprintf('\n   Peak-to-peak using mean-peak on filtered ERP:   %5.4g muV, %4.3g ms span.\n\n',erpFP2P,erpFP2Plat)

%% done.
% Conclusions
% We can clearly see how in the mean-peak case the difference (in voltage values)
% between the filtered and unfiltered case is smaller than the peak-to-peak case
% This is because averaging over a window of time instead of taking a
% single value make the analysis more robust and less vulnerable to high
% peaks caused by random noise
%%
