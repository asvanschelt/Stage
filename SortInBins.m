function [F, V1, freq] = SortInBins(dcol,t,nBins,SR)

%% Find peaks and minimals in the signal
L = length(dcol);
Prom = max(dcol)*0.1;
minpeakdistance = SR*(60/30);
[peaks,maxlocs, ~, ~] = findpeaks(dcol, 'MinPeakProminence', Prom, 'MinPeakDistance', minpeakdistance);
[mins,minlocs, ~, ~] = findpeaks(-dcol, 'MinPeakProminence', Prom, 'MinPeakDistance', minpeakdistance);
mxlcs = maxlocs*(1/SR);
mnlcs = minlocs*(1/SR);
mns = -mins;
npeaks = length(maxlocs);
nmins = length(minlocs);
freq = 60/((t(maxlocs(npeaks))-t(maxlocs(1)))/npeaks);

if SR ==1
    for m=1:npeaks
    mxlcs(m) = t(maxlocs(m),1);
    mnlcs(m) = t(minlocs(m),1);
    end
    dots = 'o';
else
    dots = '.';
end

%% Phases

F = nan(L,3);
F(:,1) = dcol;

tops = nan(npeaks+nmins);
tops(1:npeaks,1) = maxlocs;
tops(npeaks+1:length(tops),1) = minlocs;
sort_tops = sortrows(tops,1);
ntops = length(sort_tops);

up_down = nan(L,4);
up_down(:,1) = dcol;
up_down(:,2) = (1:L);

for i=2:npeaks-1
   phase_length_1 = maxlocs(i-1):maxlocs(i);
   phase_signal_1 = dcol(phase_length_1);
   empty_1=nan(length(phase_length_1),3);
   empty_1(:,1)=phase_signal_1;
   empty_1(:,2)=phase_length_1;
   bindex_1 = round(linspace(0.5,nBins+0.5,length(phase_length_1)));
   empty_1_sort = sortrows(empty_1,1);
   empty_1_sort(:,3) = bindex_1;
   empty_1_sortback = sortrows(empty_1_sort,2);
   
   F(phase_length_1,2) = empty_1_sortback(:,3);
end

for i=2:ntops
    phase_length = sort_tops(i-1):sort_tops(i)-1;
    phase_signal = dcol(phase_length);
    
    empty=nan(length(phase_length),3);
    empty(:,1)=phase_signal;
    empty(:,2)=phase_length;
    
    if dcol(sort_tops(i))>dcol(sort_tops(i-1)) % phase up
        bindex = round(linspace(0.5,(nBins/2)+0.5,length(phase_length)));
        bindex(bindex>nBins/2) = nBins;
        up_down(phase_length,4)=1;
    else % phase down
        bindex = round(linspace((nBins/2)+0.5,nBins+0.5,length(phase_length)));
        bindex(bindex>nBins) = nBins;
        up_down(phase_length,4)=0;
    end
    
    empty_sort = sortrows(empty,1);
    empty_sort(:,3) = bindex;
    empty_sortback = sortrows(empty_sort,2);
    F(phase_length,3)=empty_sortback(:,3);
end

% Phases1 = figure; hold on;
% set(gca,'Color','k');
% C = hsv(nBins);
% E = zeros(nBins,1);
% plot(mnlcs, mns, 'w*'); hold on;
% plot(mxlcs, peaks, 'w*'); hold on;
% plot(t,dcol,'w'); hold on;
% for ii=1:nBins
%     B = F(:,2)==ii;
%     X = F(B,1);
%     T = t(B);
%     scatter(T,X,[],C(ii,:),dots);
%         axis([min(t) max(t) min(dcol) max(dcol)]);
%         xlabel('Phases1');
%     E(ii,1)=length(find(F(:,2)==ii));
% end
% hold off

% Phases2 = figure; hold on;
% set(gca,'Color','k');
% E2 = zeros(nBins,1);
% plot(mnlcs, mns, 'w*'); hold on;
% plot(mxlcs, peaks, 'w*'); hold on;
% plot(t,dcol,'w'); hold on;
% for ii=1:nBins
%     B2 = F(:,3)==ii;
%     X2 = F(B2,1);
%     T2 = t(B2);
%     scatter(T2,X2,[],C(ii,:),dots);
%         axis([min(t) max(t) min(dcol) max(dcol)]);
%         xlabel('Phases2');
%     E2(ii,1)=length(find(F(:,3)==ii));
% end
% hold off

%% Value binning

% Value 1
V1 = nan(L,5);
V1(:,1) = dcol;
V1(:,2) = t;
bindex_v = round(linspace(0.5,nBins+0.5,L));
sort_A = sortrows(V1,1);
sort_A(:,3) = bindex_v;
V1 = sortrows(sort_A,2);

% Value1 = figure; hold on;
% set(gca,'Color','k');
% E3 = zeros(nBins,1);
% plot(mnlcs, mns, 'w*'); hold on;
% plot(mxlcs, peaks, 'w*'); hold on;
% plot(t,dcol,'w'); hold on;
% for ii=1:nBins
%     B3 = V1(:,3)==ii;
%     X3 = V1(B3,1);
%     T3 = t(B3);
%     scatter(T3,X3,[],C(ii,:),dots);
%         axis([min(t) max(t) min(dcol) max(dcol)]);
%         xlabel('Value');
%     E3(ii,1)=length(find(V1(:,2)==ii));
% end
% hold off

% Value 2
up = up_down(up_down(:,4)==1,:);
down = up_down(up_down(:,4)==0,:);

sort_up = sortrows(up,1);
sort_down = sortrows(down,1);

bindex_up = round(linspace(0.5,nBins/2+0.5,length(up(:,1))));
bindex_down = round(linspace((nBins/2)+0.5,nBins+0.5,length(down(:,1))));
    bindex_up(bindex_up>nBins/2) = nBins/2;
    bindex_down(bindex_down>nBins) = nBins;
    
sort_up(:,3) = bindex_up;
sort_down(:,3) = bindex_down;

sortback_up = sortrows(sort_up,2);
sortback_down = sortrows(sort_down,2);

A2 = nan(length(sort_up)+length(sort_down),3);
    A2(1:length(sort_up),1) = sortback_up(:,1);
    A2(length(sort_up)+1:end,1) = sortback_down(:,1);
    A2(1:length(sort_up),2) = sortback_up(:,2);
    A2(length(sort_up)+1:end,2) = sortback_down(:,2);
    A2(1:length(sort_up),3) = sortback_up(:,3);
    A2(length(sort_up)+1:end,3) = sortback_down(:,3);
V2 = sortrows(A2,2);
V1((V2(1,2):V2(length(V2),2)),4) = V2(:,3);

% Value2 = figure; hold on;
% set(gca,'Color','k');
% E4 = zeros(nBins,1);
% plot(mnlcs, mns, 'w*'); hold on;
% plot(mxlcs, peaks, 'w*'); hold on;
% plot(t,dcol,'w'); hold on;
% for ii=1:nBins
%     B4 = V1(:,4)==ii;
%     X4 = V1(B4,1);
%     T4 = t(B4);
%     scatter(T4,X4,[],C(ii,:),dots);
%         axis([min(t) max(t) min(dcol) max(dcol)]);
%         xlabel('value2');
%     E4(ii,1)=length(find(V1(:,4)==ii));
% end
% hold off

Fdcol = gradient(dcol);
A3 = nan(L,4);
A3(:,1) = dcol;
A3(:,2) = (1:L);
A3(:,3) = Fdcol;

up3 = A3(A3(:,3)>=0,:);
down3 = A3(A3(:,3)<0,:);

sort_up3 = sortrows(up3,1);
sort_down3 = sortrows(down3,1);

bindex_up3 = round(linspace(0.5,nBins/2+0.5,length(up3(:,1))));
bindex_down3 = round(linspace((nBins/2)+0.5,nBins+0.5,length(down3(:,1))));
    bindex_up3(bindex_up3>nBins/2) = nBins/2;
    bindex_down3(bindex_down3>nBins) = nBins;
    
sort_up3(:,4) = bindex_up3;
sort_down3(:,4) = bindex_down3;

sortback_up3 = sortrows(sort_up3,2);
sortback_down3 = sortrows(sort_down3,2);

A4 = nan(length(sort_up3)+length(sort_down3),3);
    A4(1:length(sort_up3),1) = sortback_up3(:,1);
    A4(length(sort_up3)+1:end,1) = sortback_down3(:,1);
    A4(1:length(sort_up3),2) = sortback_up3(:,2);
    A4(length(sort_up3)+1:end,2) = sortback_down3(:,2);
    A4(1:length(sort_up3),3) = sortback_up3(:,4);
    A4(length(sort_up3)+1:end,3) = sortback_down3(:,4);
sort_A4 = sortrows(A4,2);
V1(:,5)=sort_A4(:,3);

% Value3 = figure; hold on;
% set(gca,'Color','k');
% E5 = zeros(nBins,1);
% plot(mnlcs, mns, 'w*'); hold on;
% plot(mxlcs, peaks, 'w*'); hold on;
% plot(t,dcol,'w'); hold on;
% for ii=1:nBins
%     B5 = sort_A4(:,3)==ii;
%     X5 = sort_A4(B5,1);
%     T5 = t(B5);
%     scatter(T5,X5,[],C(ii,:),dots);
%         axis([min(t) max(t) min(dcol) max(dcol)]);
%         xlabel('value3');
%     E5(ii,1)=length(find(sort_A4(:,3)==ii));
% end
% hold off

end