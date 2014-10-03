function PVTTrigger(Triggerlogname,Endlabtimesync)

% % Load database file
% SubID=num2str(SubID);
% PVTlogname=strcat(SubID,'GX_PVT LOG.txt');
%Triggerlogname='TSTMarks_3319GX.txt';
Triggerlog=importdata(Triggerlogname);
Triggerlog.textdata(1:5,:)
Triggerlog.data(1:5,:)

TriggerLabtimes=(Triggerlog.textdata);
TriggerLabtimes(1:20,:)

% Create matrix of PVT tests 
 
PVTtimes=strcmp(TriggerLabtimes(:,5),'PVT');
PVTsessions=TriggerLabtimes(PVTtimes,:);

% Identify Labtime from TSTMarks log closest to input Labtime. Outputs labtime PVT starts of 7 cogs of
% recording

% % NEED TO CODE TO ENABLE INPUT OF END LABTIME SYNC --> manually entered
% as function
% Forced Endlabtimesync
%Endlabtimesync='2898:52';

% Remove colon from endlabtime sync 
[Endlabtime,b] = strtok(Endlabtimesync,':');
[b,c]=strtok(b,':');
Endlabtime=strcat(Endlabtime,b);

% Remove colon from PVT sessions time
% [PVT, c] = strtok(PVTsessions,':');
PVT = PVTsessions(:,2);
[PVT,b]=strtok(PVT,':');
[b,c]=strtok(b,':');
PVT=strcat(PVT,b);

% Convert Endlabtimesync from string to number
Endlabtime = str2num(Endlabtime);
PVT = str2double(PVT);
% PVT=unique(PVT)

% Find closest PVT time to End Recording labtime sync
% % NEED TO CODE METHOD TO CIRCUMVENT SITUATION WHERE ENDLABTIME IS CLOSER
% % TO A LATER PVT SESSION THAN A PREVIOUS SESSION

%[~, h] = min(abs(PVT - (Endlabtime)));
for i = 1:length(PVT)
    if PVT(i) < Endlabtime
        PVTnew(i) = PVT(i);
    end
end
[~,h] = max(PVTnew);  
PVTWP=PVT(h);


% Identify all PVT session times that occurred during same WP (i.e., within
% 2 hours of the previous)

for i=0:30
    if PVT(h-i)-PVT(h-i-1)<300
    PVTWP(end+1,:)=PVT(h-i);
    else
    break
    end
end

% Displays Labtimes of start of PVT sessions of WP
PVTsessionlabtimes=PVTsessions(h-(length(PVTWP)-1):2:h,:);
PVT1=strcat(PVTsessionlabtimes(:,2),':');
PVTsessionlabtimes=strcat(PVT1,PVTsessionlabtimes(:,4))

% PVTsessionlabtimes=strcat(PVTsessionlabtimes(:,2),(PVTsessionlabtimes(:,5)

% fileName = 'PVTtrigger.dat'
% fileID = fopen('PVTtrigger.txt','w');
fopen('PVTtrigger.txt','w');
formatspec='%s\n';
dlmwrite('PVTtrigger.txt',PVTsessionlabtimes,'');


% printf(fileID,formatspec,PVTsessionlabtimes);
% fprintf(fileID,formatspec,PVTsessionlabtimes(:,:));
% fclose(fileID);


% nrows=length(PVTsessionlabtimes)
% for i=1:nrows
% fprintf(fileID,'%s\n',PVTsessionlabtimes(i,:));
% end



