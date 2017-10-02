clear all
clc

nBins = 10;

[dcol_nav, t_nav] = ReadNav ('SCANPHYSLOG_dc_27092017_1839033_19_1_dce_dix_radial_navV4.log');
[dcol_fys, t_fys] = ReadFysLog('devlogcurrent_1757-1900.log');

[F_fys, V_fys, freq_fys] = SortInBins(dcol_fys,t_fys,nBins,1);
[F_nav, V_nav, freq_nav] = SortInBins(dcol_nav,t_nav,nBins,500);

for p=2:3
figure; hold on;
set(gca,'Color','k');
E = zeros(nBins,1);
E1 = zeros(nBins,1);
C = hsv(nBins);
h(1) = subplot(1,2,1);
    plot(t_fys, dcol_fys); hold on;
    for ii=1:nBins
        B = F_fys(:,p)==ii;
        X = F_fys(B,1);
        T = t_fys(B);
        scatter(T,X,[],C(ii,:),'o');
        axis([min(t_fys) max(t_fys) min(dcol_fys) max(dcol_fys)]);
        xlabel('Phases fys');
        E(ii,1)=length(find(F_fys(:,p)==ii));
    end
    hold off
h(2) = subplot(1,2,2);
    plot(t_nav, dcol_nav); hold on;
    for ii=1:nBins
        B1 = F_nav(:,p)==ii;
        X1 = F_nav(B1,1);
        T1 = t_nav(B1);
        scatter(T1,X1,[],C(ii,:),'.');
        axis([min(t_nav)+15 max(t_nav) min(dcol_nav) max(dcol_nav)]);
        xlabel('Phases nav');
        E1(ii,1)=length(find(F_nav(:,p)==ii));
    end
    hold off
end

for p=3:5
figure; hold on;
set(gca,'Color','k');
E2 = zeros(nBins,1);
E3 = zeros(nBins,1);
C = hsv(nBins);
h(1) = subplot(1,2,1);
    plot(t_fys, dcol_fys); hold on;
    for ii=1:nBins
        B = V_fys(:,p)==ii;
        X = V_fys(B,1);
        T = t_fys(B);
        scatter(T,X,[],C(ii,:),'o');
        axis([min(t_fys) max(t_fys) min(dcol_fys) max(dcol_fys)]);
        xlabel('Value fys');
        E(ii,1)=length(find(V_fys(:,p)==ii));
    end
    hold off
h(2) = subplot(1,2,2);
    plot(t_nav, dcol_nav); hold on;
    for ii=1:nBins
        B1 = V_nav(:,p)==ii;
        X1 = V_nav(B1,1);
        T1 = t_nav(B1);
        scatter(T1,X1,[],C(ii,:),'.');
        axis([min(t_nav)+15 max(t_nav) min(dcol_nav) max(dcol_nav)]);
        xlabel('Value nav');
        E1(ii,1)=length(find(V_nav(:,p)==ii));
    end
    hold off
end

figure; hold on;
plot(t_fys,dcol_fys,'r'); hold on;
plot(t_nav,dcol_nav/500,'b'); hold off;