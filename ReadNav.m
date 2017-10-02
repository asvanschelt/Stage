function [dcol_nav, t_nav] = ReadNav(signal)

%% Read the scan of signal 1
D = ReadPhilipsScanPhysLog(signal);

% Plot the respiratory motion
L = length(D.C);
SR = 500;
t_nav = (1:L)*(1/SR);
col = D.C(1:L,6);
dcol_nav_1 = double(col);

fileID= fopen(signal); 
text = textscan(fileID,'%s');
fclose(fileID); 

search_markers = 'mark2';
search_start = '0010';
search_end = '0020';

ntitel = size(D.C,2);
nmark = 10;

loc_markers= find(strcmp(text{1},search_markers));
loc_start= find(strcmp(text{1},search_start));
loc_end_all= find(strcmp(text{1},search_end));

loc_start = ((loc_start-loc_markers)-nmark) / ntitel;
loc_end = ((loc_end_all(:)-loc_markers)-nmark) / ntitel;

dcol_nav = dcol_nav_1(loc_start:loc_end(3));
dcol_nav = -dcol_nav;
t_nav = t_nav(1:length(dcol_nav)).';
end
