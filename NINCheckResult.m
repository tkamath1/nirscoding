function [OHRaw, OHFil, ECG, Resp, tVec, nData, event2Index, event1Index,WP,cog,SubID]= NINCheckResult(dataFileNm, targetDataSec, baselineSec, adaptiveFilterFlag, secLabel,baselineFileNm,cog,EndLabtimeSync)
% TO DO
%     i) Need to split files that are greater than 12 hours due to
%     "out-of-memory" error

%[OHRaw, OHFil, ECG, Resp, tVec, nData, event2Index, event1Index]= NINCheckResult(dataFileNm, targetDataSec, baselineSec, adaptiveFilterFlag, secLabel, baselineFileNm)

%cog=input('Enter Cog #:')
%Labtime=input('Enter labtime of beginning of NIRS recording#:')
set(0,'DefaultFigureVisible','off');
if nargin == 1
%Preview
[OHRaw, ECG, Resp, tVec, nData, event2Index, event1Index] = OH3Input(dataFileNm); 

[b,a]=butter(3,5/125);
OHFil=filtfilt(b, a, OHRaw);
% Fig=figure; plot(tVec,OHFil(:,1:4))
% [event1Marks, event2Marks ] = OHMarkEventNum( Fig, tVec, event1Index, event2Index);
% 
% set(gcf,'position',get(0,'screensize'));
% print(figure(Fig),'-dmeta',[dataFileNm, secLabel, 'Result']);

return;

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           1. CONFIGURATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% dataFileNm = 'cog1';
% baselineSec = [1291 1294];%in seconds
% targetDataSec = [2300 3200];%PVT section; in seconds 


sampleRate=250;
halfSampleRate=125;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           2. LOAD DATA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NIN channel configuration:
% channel 1, .1k1, near detector, 650Hz, 775nm
% channel 2, .8k1, near detector, 6KHz, 826nm
% channel 3, .1k2, far detector, 650Hz, 775nm
% channel 4, .8k2, far detector, 6KHz, 826nm
% channel 5, ECG
% channel 6, Accelerometer axis 1
% channel 7, Accelerometer axis 2
% channel 8, Accelerometer axis 3


[OHRaw, ECG, Resp, tVec, nData, event2Index, event1Index] = OH3Input(dataFileNm); 

[b,a]=butter(3,5/halfSampleRate);
OHFil=filtfilt(b, a, OHRaw);
% Fig=figure; plot(tVec,OHFil(:,1:4))
% [event1Marks, event2Marks ] = OHMarkEventNum( Fig, tVec, event1Index, event2Index);
% grid;
% set(gcf,'position',get(0,'screensize'))
% print(figure(Fig),'-dmeta',[dataFileNm, secLabel,  'Result']);

% % TUSHAR : Ignore the Baselinesec and baselineFilenm for now. We can use
% the same method we used before; use the baseline data from our favorite
% recording

if nargin == 8
    [OHRawTemp, ECGBase, RespBase, tVecBase, dataLengthBase, event2IndexBase, event1IndexBase] = OH3Input(baselineFileNm); 
    OHRawBase = OHRawTemp((baselineSec(1)*sampleRate):(baselineSec(2)*sampleRate), :);
else
    OHRawBase = OHRaw((baselineSec(1)*sampleRate):(baselineSec(2)*sampleRate), :); 
end;

offsetMean=median(OHRawBase);
optOffset=offsetMean([1:4]);
optOffset(4)=0.7*optOffset(4); %too large some times. 

% % TUSHAR : Create time matrix. 1st column in NIRS elapsed time, 2nd
% column in labtime. That way tVec and labtime can be retrieved. The
% targetdatasec input will be in labtime, and converted to NIRS time for
% tVec. 


    
    
 
optData=OHFil((targetDataSec(1)*sampleRate):(targetDataSec(2)*sampleRate),1:4);
tVecData=tVec((targetDataSec(1)*sampleRate):(targetDataSec(2)*sampleRate));
ECGData=OHFil((targetDataSec(1)*sampleRate):(targetDataSec(2)*sampleRate), 5);
AccData=OHFil((targetDataSec(1)*sampleRate):(targetDataSec(2)*sampleRate),6:8);

dataLength=length(optData);
optData=optData-ones(dataLength,1)*optOffset;


event1Time=(event1Index-1)/sampleRate; %convert to seconds
event2Time=(event2Index-1)/sampleRate;

event1SecIdx = event1Time>targetDataSec(1) & event1Time<targetDataSec(2);
event1SecTime=event1Time(event1SecIdx);
event2SecIdx = event2Time>targetDataSec(1) & event2Time<targetDataSec(2);
event2SecTime=event2Time(event2SecIdx);

% Fig=figure; plot(tVecData, optData);
% grid;
% xlabel('time(sec)');
% ylabel('AD Unit');
% title([dataFileNm, secLabel,'rawData']);
% NinMarkEvent(event1SecTime, event2SecTime, max(optData(:)), min(optData(:)));
% 
% set(gcf,'position',get(0,'screensize'))
% print(figure(Fig),'-dmeta',[dataFileNm, secLabel,  'rawData']);

% 
% Fig=figure; plot(tVecData, AccData);
% grid;
% xlabel('time(sec)');
% ylabel('AD Unit');
% title([dataFileNm, secLabel,' accelerometer']);
% NinMarkEvent(event1SecTime, event2SecTime, max(AccData(:)), min(AccData(:)));
% 
% set(gcf,'position',get(0,'screensize'))
% print(figure(Fig),'-dmeta',[dataFileNm, secLabel,  'accData']);

%xlswrite('Result_rawDataOutput.xls', [tVecData(1:10:end)' optData(1:10:end, :)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           3. CONVERSION TO HEMOGLOBINE CONTENT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dpf=[6 6];
sep(1)=1.5;%1.1;
sep(2)=4;%3.8;
path_factor = [ dpf(2)*sep(1)
    dpf(1)*sep(1)
    dpf(2)*sep(2)
    dpf(1)*sep(2)]';


normBase=median(optData(1:1000,:));
optData=log(ones(dataLength,1)*normBase./optData);
optData=optData./(ones(dataLength,1)*path_factor);



disp(sprintf('Converting to concentrations.'));
e_hb826 = 778.5;
e_hb775 = 1216.3;
e_hbo826 = 990.5;	
e_hbo775 = 733.7;

BT = [e_hb775 e_hb826
    e_hbo775 e_hbo826];
BTI=inv(BT);

for chanIdx=1:size(optData,2)/2
    AT(:,2)=optData(:,chanIdx*2);%Mua826
    AT(:,1)=optData(:,chanIdx*2-1);%Mua775
    CT(:,:,chanIdx)=AT*BTI;
end;


Hb=squeeze(CT(:,1,:));
HbO2=squeeze(CT(:,2,:));

% 
% Fig=figure; 
% subplot(4,1,1); plot(tVecData, HbO2(:,1), 'r'); title('HbO2, Near'); grid;
% NinMarkEvent(event1SecTime, event2SecTime, max(HbO2(:,1)), min(HbO2(:,1)));
% subplot(4,1,2); plot(tVecData, Hb(:,1)); title('Hb, Near');grid;
% NinMarkEvent(event1SecTime, event2SecTime, max(Hb(:,1)), min(Hb(:,1)));
% subplot(4,1,3); plot(tVecData, HbO2(:,2), 'r'); title('HbO2, Far'); xlabel('time(sec)');grid;
% NinMarkEvent(event1SecTime, event2SecTime, max(HbO2(:,2)), min(HbO2(:,2)));
% subplot(4,1,4); plot(tVecData, Hb(:,2)); title('Hb, Far');grid;
% NinMarkEvent(event1SecTime, event2SecTime, max(Hb(:,2)), min(Hb(:,2)));
% 
% set(gcf,'position',get(0,'screensize'))
% print(figure(Fig),'-dmeta',[dataFileNm, secLabel,  'BloodDetails']);

% Fig=figure; 
% subplot(2,1,1); plot(tVecData, HbO2(:,1), 'r'); 
% hold on; plot(tVecData, Hb(:,1)); title('Hb & HbO2, Near');grid;
% NinMarkEvent(event1SecTime, event2SecTime, max([HbO2(:,1); Hb(:,1)]), min([HbO2(:,1); Hb(:,1)]));
% subplot(2,1,2); plot(tVecData, HbO2(:,2), 'r'); xlabel('time(sec)');
% hold on; plot(tVecData, Hb(:,2)); title('Hb & HbO2, Far');grid;
% NinMarkEvent(event1SecTime, event2SecTime, max([HbO2(:,2); Hb(:,2)]), min([HbO2(:,2); Hb(:,2)]));
% 
% set(gcf,'position',get(0,'screensize'))
% print(figure(Fig),'-dmeta',[dataFileNm, secLabel,  'HbHbO2']);

xlswrite('Result_HbHbO2Output.xls', [tVecData(1:10:end)' Hb(1:10:end,:) HbO2(1:10:end,:)]);

xlswrite('Result_ECGAccOutput.xls', [tVecData(1:10:end)' ECGData(1:10:end, :) AccData(1:10:end, :)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           4. OUTPUT CODE TO COMPARE NIRS ABSORPTION WITH PVT TRIAL PRESENTATION. MLL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % % 
% % % 
format long e
% Subject ID
delimiter='GX';
SubID = strtok(dataFileNm,delimiter);
SubID=str2num(SubID);
SubjectID(1:length(tVecData),1)=SubID;

% WP

delimiter='WP';
[x,WP] = strtok(dataFileNm,delimiter);
WP=strtok(WP,delimiter);
delimiter='_';
WP=strtok(WP,delimiter);
WP=str2num(WP);
WakePeriod(1:length(tVecData),1)=WP;

% COG number

cogmatrix(1:length(tVecData),1)=cog;

% Week in FD

if WP(1,1)<11
    FD=0;
end
if WP(1,1)>10 && WP(1,1)<17
    FD=1;
end
if WP(1,1)>26 && WP(1,1)<34
    FD=2;
end
if WP(1,1)>34 
    FD=3;
end

FD
FD(1:length(tVecData),1)=FD;

% Identify beginning of PVT sessions

fileID = strcat(dataFileNm,'.event.txt');
A=importdata(fileID);
tmp = regexp(A,'([^ ,:]*)','tokens');
D = cellfun(@(c) c{4},tmp);
S = sprintf('%s', D{:});
z=strread(S,'%s','delimiter','(');
y=str2mat(z);
eventmarks=str2num(y);


%Scan event marks for beginning of PVT  sessinos

eventmarkersize=length(eventmarks);
pvtbegin=0;
pvtend=0;
for i=1:eventmarkersize-6
    if eventmarks(i+1,1)-eventmarks(i,1)>248 && eventmarks(i+1,1)-eventmarks(i,1)<251
        if eventmarks(i+2,1)-eventmarks(i+1,1)>60 && eventmarks(i+2,1)-eventmarks(i+1,1)<65
            if eventmarks(i+3,1)-eventmarks(i+2,1)>60 && eventmarks(i+3,1)-eventmarks(i+2,1)<65
                if eventmarks(i+4,1)-eventmarks(i+3,1)>60 && eventmarks(i+4,1)-eventmarks(i+3,1)<65
                    if eventmarks(i+5,1)-eventmarks(i+4,1)>60 && eventmarks(i+5,1)-eventmarks(i+4,1)>150000
                        if eventmarks(i+6,1)-eventmarks(i+5,1)>150 && eventmarks(i+6,1)-eventmarks(i+5,1)<160
                            if eventmarks(i+7,1)-eventmarks(i+6,1)>60 && eventmarks(i+7,1)-eventmarks(i+6,1)<65
                                if eventmarks(i+8,1)-eventmarks(i+7,1)>60 && eventmarks(i+8,1)-eventmarks(i+7,1)<65
                                    if eventmarks(i+9,1)-eventmarks(i+8,1)>60 && eventmarks(i+9,1)-eventmarks(i+8,1)<65
                                        if eventmarks(i+10,1)-eventmarks(i+9,1)>65 
                                            if pvtbegin(1,1)==0
                                                pvtbegin(length(pvtbegin),1)= eventmarks(i,1);
                                                pvtend(length(pvtbegin),1) =eventmarks(i+5,1);
                                            else
                                                pvtbegin(1+length(pvtbegin),1)=eventmarks(i,1);
                                                pvtend(1+length(pvtend),1) =eventmarks(i+5,1);
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

pvtbeginsize=size(pvtbegin);
pvtbegintime=pvtbegin/250
pvtendtime=pvtend/250
tVecsize=size(tVecData);
 
pvtbegin(1:length(tVecData),1)=0;
pvtend(1:length(tVecData),1)=0;

% Load PVT database file
SubID=num2str(SubID);
PVTlogname=strcat('VPVTONE_',SubID,'GX.txt');
PVTlog=importdata(PVTlogname);

% Identify Labtime from PVT log closest to input Labtime. Outputs labtime PVT starts of 7 cogs of
% recording

PVTLabtimes=(PVTlog.textdata(:,3));
[b,c]=strtok(PVTLabtimes,':');
[d,e]=strtok(c,':');
f=strcat(b,d);
g=str2double(f(:,:));
i=unique(g);
[~, h] = min(abs(i - (EndLabtimeSync)));

PVTlabtime = i(h:h+6);

% Assumes 7 cogs. If 7th cog is for next recording, it is deleted
if PVTlabtime(end)-PVTlabtime(1)>1999
PVTlabtime=PVTlabtime(1:end-1);
end

% Removes first 2 columns of PVTlog.data and adds column summing trial
% start and delay to return trial presentation

PVTlogdata=PVTlog.data(:,3:6);
PVTlogdata(:,5)=PVTlogdata(:,2)+PVTlogdata(:,3);

% Adds Labtime column to PVTlogdata

PVTlogdata(:,6)=g;

% Create Cog number column
PVTlogdata(:,7)=zeros;

for n=1:length(PVTlabtime)
    for i=1:length(PVTlogdata)
        if PVTlabtime(n)==PVTlogdata(i,6);
        PVTlogdata(i,7)=n;
        end
    end
end

% Removes all other cogs except current recording

% % NEED TO DEBUG. Have to manually input COG number. 
PVTlogdata=PVTlogdata(PVTlogdata(:,7)==cog,:);

% PVTlogdata=PVTlogdata(PVTlogdata(:,7)==1,:);

% Create a column for converting PVT time to NIRS time. Code uses 5x Event
% mark timepoint as beginning of PVT session.
PVTlogdata(:,8)=zeros;
pvtbegintime=pvtbegintime*1000;
PVTlogdata(1,8)=pvtbegintime(cog);        
for i=2:length(PVTlogdata)
    PVTlogdata(i,8)=PVTlogdata(1,8)+PVTlogdata(i,5);
end
        
% Individual PVT trial start times
% % NEED PVT TRIAL START TO FIND CLOSEST VALUE OF ABSORB DATA

PVTtrialtime=targetDataSec(1)+PVTlogdata(:,5)/1000;

for i=1:length(PVTtrialtime)
    
[~, h] = min(abs(tVecData-PVTtrialtime(i)));
PVTtrialtimeconvert(i)=tVecData(h);

end
% Export PVTtrialtime (individual PVT trial times to txt file)
%cog=input('Enter Cog #:')
fileName3 = strcat('PVTtrial',SubID,'_WP',num2str(WP),'_Cog',num2str(cog));
fileID = fopen(fileName3,'w');
fprintf(fileID,'%d \r\n',PVTtrialtime');
fclose(fileID);

%Creates fileName for Combine and AbsData files
fileName = strcat('CombineData_',SubID,'.txt');
fileName2 = strcat('AbsData_',SubID,'.txt');

% Creates matrix with NIRS time (tVec) and Hb, Hbo2
Absdata=horzcat(tVecData',Hb,HbO2);
fileID = fopen(fileName2,'w');
fprintf(fileID,'%d %d %d %d %d \r\n',Absdata');



% Merge SubjectID, WakePeriod, tVecData, Hb, and HbO2 matrices
AllData=horzcat(WakePeriod, cogmatrix, tVecData', Hb, HbO2);
AllDataSize=size(AllData);


% fileID = strcat(dataFileNm,'.event.txt');
fileID = fopen(fileName,'w');
fprintf(fileID,'%6s %6s %6s %6s %6s %6s %6s\r\n','AllData');
fprintf(fileID,'%d %d %d %d %d %d %d \r\n',AllData');

fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           5. SPECTRAL ANALYSIS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HbO2Near = HbO2(:,1);
HbNear = Hb(:,1);
HbO2Far = HbO2(:,2);
HbFar = Hb(:,2);

% Fig=figure;
% subplot(3,1,1); psd(HbO2Far, dataLength, sampleRate); 
% subplot(3,1,2); psd(HbO2Far, dataLength, sampleRate); axis([0 5 -160 -70])
% subplot(3,1,3); psd(HbO2Far, dataLength, sampleRate); axis([0 1 -160 -70])
% set(gcf,'position',get(0,'screensize'))
% print(figure(Fig),'-dmeta',[dataFileNm ,secLabel, 'HbO2FarPsd']);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           6. ADAPTIVE FILTERING
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%further filering and downsampling before adaptive filtering
if adaptiveFilterFlag == 1
disp('Adaptive Filter is on')
bandwidthLow=3;
bandwidthHigh=0.005;
[b, a]=butter(3,bandwidthLow/halfSampleRate);
[bH, aH]=butter(3,0.005/halfSampleRate,'high');

HbO2Fil=filtfilt(b, a, HbO2);
HbO2Fil=filtfilt(bH,aH,HbO2Fil);
HbFil=filtfilt(b, a, Hb);
HbFil=filtfilt(bH,aH,HbFil);

% Possible issue - we don't downsample...this is downsampling
%tVecDown=tVecData(1:10:end);
%HbO2Fil=HbO2Fil(1:10:end, :);
%HbFil=HbFil(1:10:end, :);

HbO2Near = HbO2Fil(:,1);
HbNear = HbFil(:,1);
HbO2Far = HbO2Fil(:,2);
HbFar = HbFil(:,2);

% Fig=figure; 
% subplot(4,1,1); plot(tVecDown, HbO2Near, 'r'); title('HbO2, Near');grid;
% subplot(4,1,2); plot(tVecDown, HbNear); title('Hb , Near');grid;
% subplot(4,1,3); plot(tVecDown, HbO2Far, 'r'); title('HbO2, Far');grid;
% subplot(4,1,4); plot(tVecDown, HbFar); title('Hb, Far');xlabel('time(sec)');grid;
% 
% set(gcf,'position',get(0,'screensize'))
% print(figure(Fig),'-dmeta',[dataFileNm, secLabel, 'HbHbO2DownFil']);

%adaptive filtering for HbO2

[HbO2Adap, hRecHbO2]=adaptiveQuan((HbO2Far)/std(HbO2Far),(HbO2Near)/std(HbO2Near),[1 zeros(1, 99)]);
[HbO2Adap, hRecHbO2]=adaptiveQuan((HbO2Far)/std(HbO2Far),(HbO2Near)/std(HbO2Near),hRecHbO2(end,:));

%%%%%%%%%%%%%%Low-pass filter on Hb data to remove cardiac rhythm%%%%%%%
[b2, a2] = butter(3,0.5/halfSampleRate);

HbO2Adap = filtfilt(b2,a2,HbO2Adap);


HbO2Adapdata = vertcat(tVecData, HbO2Adap');
IDoffile2 = fopen(strcat(SubID,'HbO2Adap.txt'),'w');
fprintf(IDoffile2,'%d %d \r\n',HbO2Adapdata);
fclose(IDoffile2);
%[HbO2Adap, hRecHbO2]=adaptiveQuan((HbO2Far)/std(HbO2Far),(HbO2Near)/std(HbO2Near),hRecHbO2(end,:));%removing respirati
% 
% Fig=figure; 
% subplot(4,1,1); plot(tVecDown, HbO2Far, 'r'); title('HbO2, Far');grid;
% subplot(4,1,2); plot(tVecDown, HbO2Near, 'r'); title('HbO2, Near');grid;
% subplot(4,1,3); plot(tVecDown, HbO2Adap, 'r'); title('HbO2, adaptive filtered');grid;
% subplot(4,1,4); plot(tVecDown, hRecHbO2); title('adaptive filtering node update');xlabel('time(sec)');grid;
% set(gcf,'position',get(0,'screensize'))
% print(figure(Fig),'-dmeta',[dataFileNm, secLabel,  'HbO2AdaptiveFiltering']);


%adaptive filtering for Hb
[HbAdap, hRecHb]=adaptiveQuan((HbFar)/std(HbFar),(HbNear)/std(HbNear),[1 zeros(1, 99)]);
[HbAdap, hRecHb]=adaptiveQuan((HbFar)/std(HbFar),(HbNear)/std(HbNear),hRecHb(end,:));


HbAdap = filtfilt(b2,a2,HbAdap);

HbAdapdata = vertcat(tVecData, HbAdap');
IDoffile = fopen(strcat(SubID,'HbAdap.txt'),'w');
fprintf(IDoffile,'%d %d \r\n',HbAdapdata);
fclose(IDoffile);

Adapdata= vertcat(tVecData,HbO2Adap',HbAdap');
IDoffile3 = fopen(strcat(SubID,'_WP',num2str(WP),'_cogno',num2str(cog),'AdapData.txt'),'w');
fprintf(IDoffile3,'%d %d %d \r\n',Adapdata);
fclose(IDoffile3);
%[HbAdap, hRecHb]=adaptiveQuan((HbFar)/std(HbFar),(HbNear)/std(HbNear),hRecHb(end,:));%removing respirati

% Fig=figure; 
% subplot(4,1,1); plot(tVecDown, HbFar); title('Hb, Far');grid;
% subplot(4,1,2); plot(tVecDown, HbNear); title('Hb, Near');grid;
% subplot(4,1,3); plot(tVecDown, HbAdap); title('Hb, adaptive filtered');grid;
% subplot(4,1,4); plot(tVecDown, hRecHb); title('adaptive filtering node update');xlabel('time(sec)');grid;
% set(gcf,'position',get(0,'screensize'))
% print(figure(Fig),'-dmeta',[dataFileNm, secLabel,  'HbAdaptiveFiltering']);
h = msgbox('Operation Completed');
set(0,'DefaultFigureVisible','on');
else
    h = msgbox('Operation Completed');
    set(0,'DefaultFigureVisible','on');
end;




