function PVTTriggerML(Triggerlogname,PVTlogname,Endlabtimesync)

% % Load database file
% SubID=num2str(SubID);
% PVTlogname=strcat(SubID,'GX_PVT LOG.txt');
%Triggerlogname='TSTMarks_3319GX.txt';
Triggerlog=importdata(Triggerlogname);
TriggerLabtimes=(Triggerlog.textdata);


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

% Remove colons from PVT sessions time
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
PVTsessionlabtimes(:,2);
% Extracts seconds from "Date/time" column 
[a,b]=strtok(PVTsessionlabtimes(:,3),' ');
[a,b]=strtok(b,':');
[a,b]=strtok(b,':');
[seconds,b]=strtok(b,':');

% Creates PVT sessions in labtime in format hhhh:mm:ss:msec
PVTsessiontimes=strcat(PVTsessionlabtimes(:,2),':',seconds,':',PVTsessionlabtimes(:,4));

%Import PVT data from PVT database file

PVTdata = importdata(PVTlogname);
% PVTdata.textdata(1:10,:)

% Create array of PVT presentation times (Sum trial start and trial delay)

PVTpresent= str2double(PVTdata.textdata(:,9))+  str2double(PVTdata.textdata(:,10));

% Separate PVT session start time by colon delimiter
[PVTHour,PVTSec] = strtok(PVTdata.textdata(:,3),':');
[PVTSec,~] = strtok(PVTSec,':');
PVTHour = str2double(PVTHour);
PVTSec = str2double(PVTSec);

% 
[PVTsessionstimesdelim,b] = strtok(PVTsessiontimes,':');
[PVTsessionstimesdelim(:,2),b]=strtok(b,':');
[PVTsessionstimesdelim(:,3),b]=strtok(b,':');
[PVTsessionstimesdelim(:,4),~]=strtok(b,':')

% Search PVT database for PVT session times for WP calculated from end
% recording sync time. When matched, add PVT trial elasped time PVT session
% start time (milliseconds)

for i = 1:length(PVTsessiontimes)
% for i = 1
h=0;
    for j=1:length(PVTHour)
      if PVTHour(j,1)== str2double(cell2mat(PVTsessionstimesdelim(i,1)))
        h=h+1;
          PVTtrials(h,i)=str2double(cell2mat(PVTsessionstimesdelim(i,4)))+PVTpresent(j,:);
      end
    end
end

% Remove 1st trial due to error in recording time
PVTtrials=PVTtrials(2:end,:)
format long

% Convert PVT session time to millseconds and add to PVT trial elasped time
for i = 1:length(PVTsessiontimes)
   PVTtrials(:,i)=PVTtrials(:,i)+60*60*1000*(str2double(cell2mat(PVTsessionstimesdelim(i,1))))+60*1000*(str2double(cell2mat(PVTsessionstimesdelim(i,2))))+1000*(str2double(cell2mat(PVTsessionstimesdelim(i,3))))+(str2double(cell2mat(PVTsessionstimesdelim(i,4))));
end

PVTtrials


% fopen('PVTtrigger.txt','w');
% formatspec='%s\n';
% dlmwrite('PVTtrigger.txt',PVTsessionlabtimes,'');