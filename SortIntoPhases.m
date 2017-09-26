function [resf, E, E2, Phases1, Phases2] = SortIntoPhases(signal, nBins);
    
% Read the scan
D = ReadPhilipsScanPhysLog(signal);

%% Plot the respiratory motion
[x,y] = size(D.C);
SR = 500;
t = (1:x)*(1/SR);
[one,L] = size(t);
col = D.C(1:x,6);
    % breath = plot(t, col);

%% Plot the frecuency of the respiratory motion
FT = fft(col);
P2 = abs(FT/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = SR*(0:(L/2))/L;
%  plot(f,P1);
%      xlabel('f (Hz)');
%      axis([0 10 0 inf]);
% Respiratory frequency > output breaths per minute
resf = f(find(P1 == max(P1)))*60;

%% Find peaks and minimals in the signal
dcol = double(col);
minbreath = 25;
minpeakdistance = SR*(60/minbreath);
Prom = max(col)*0.1;
[pks,maxlocs, w1, p1] = findpeaks(double(col), 'MinPeakProminence', Prom, 'MinPeakDistance', minpeakdistance);
[mins,minlocs, w2, p2] = findpeaks(-double(col), 'MinPeakProminence', Prom, 'MinPeakDistance', minpeakdistance);
mxlcs = maxlocs*(1/SR);
mnlcs = minlocs*(1/SR);
mns = -mins;

%% My own try: Phases 1

%nBins = 10;
A = nan(length(dcol),2);
A(:,1) = dcol;

for i=2:length(maxlocs)
    lengthphase = maxlocs(i) - maxlocs(i-1); % define length of the phase
        empty = nan(lengthphase,3);
    empty(:,1) = dcol(maxlocs(i-1):maxlocs(i)-1);
    empty(:,2) = maxlocs(i-1):maxlocs(i)-1;
    s2l = sortrows(empty,1); % sort the points in the phase in ascending order
    bindex = round(linspace(0.5,nBins+0.5,lengthphase)); % evenly spaced bins
    bindex(bindex>nBins) = nBins; %Bin 11 doesn't exist
    s2l(:,3)=bindex;
    bininfo = sortrows(s2l,2); %sort the rows back in normal order
    A(maxlocs(i-1):maxlocs(i)-1,2) = bininfo(:,3);
end
Phases1 = figure; hold on;
set(gca,'Color','k');
C = hsv(nBins);
E = zeros(nBins,1);
for ii=1:nBins
    B = A(:,2)==ii;
    X = A(B,1);
    T = t(B);
    scatter(T,X,[],C(ii,:),'.');
        axis([max(t)/4 max(t)/2 min(dcol) max(dcol)]);
        xlabel('Phases 1');
    E(ii,1)=length(find(A(:,2)==ii));
end
hold off

%% My own try: Phases 2

%nBins2 = 5;
A2 = nan(length(dcol),2);
A2(:,1) = dcol;

for i=1:length(maxlocs)-1
    lengthphase2 = minlocs(i+1) - maxlocs(i); % define length of the phase
    empty2 = nan(lengthphase2,3);
    empty2(:,1) = dcol(maxlocs(i):minlocs(i+1)-1);
    empty2(:,2) = maxlocs(i):minlocs(i+1)-1;
    s2l2 = sortrows(empty2,1); % sort the points in the phase in ascending order
    bindex2 = round(linspace(0.5,(nBins/2)+0.5,lengthphase2)); % evenly spaced bins
    bindex2(bindex2>nBins2) = nBins2;
    s2l2(:,3)=bindex2;
    bininfo2 = sortrows(s2l2,2); %sort the rows back in normal order
    A2(maxlocs(i):minlocs(i+1)-1,2) = bininfo2(:,3);
end
for i=1:length(minlocs)
    lengthphase2B = maxlocs(i) - minlocs(i); % define length of the phase
    empty2B = nan(lengthphase2B,3);
    empty2B(:,1) = dcol(minlocs(i):maxlocs(i)-1);
    empty2B(:,2) = minlocs(i):maxlocs(i)-1;
    s2l2B = sortrows(empty2B,1); % sort the points in the phase in ascending order
    bindex2B = round(linspace((nBins/2)+0.5,nBins+0.5,lengthphase2B)); % evenly spaced bins
    bindex2B(bindex2B>2*nBins2) = 2*nBins2;
    s2l2B(:,3)=bindex2B;
    bininfo2B = sortrows(s2l2B,2); %sort the rows back in normal order
    A2(minlocs(i):maxlocs(i)-1,2) = bininfo2B(:,3);
end
Phases2 = figure; hold on;
set(gca,'Color','k');
C2 = hsv(2*nBins2);
E2 = zeros(2*nBins2,1);
for ii=1:2*nBins2
    B2 = A2(:,2)==ii;
    X2 = A2(B2,1);
    T2 = t(B2);
    scatter(T2,X2,[],C2(ii,:),'.');
        axis([max(t)/4 max(t)/2 min(dcol) max(dcol)]);
        xlabel('Phases2');
        colorbar
    E2(ii,1)=length(find(A2(:,2)==ii));
end
end

