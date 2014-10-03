function [OHData, ECG, Resp, tVec, nData, event1Index, event2Index] = OH3Input(file)
%[OHData, ECG, Resp, tVec, nData] = OH2Input(file)
% this program performs 5 point median filtering to get rid of the snow and
% pepper noise in the digital board AD channels.

[OHData(:,1), tVec, nData]=readOH([file '.1k1']);
[OHData(:,2), tVec, nData]=readOH([file '.8k1']);
[OHData(:,3), tVec, nData]=readOH([file '.1k2']);
[OHData(:,4), tVec, nData]=readOH([file '.8k2']);
[OHData(:,5), tVec, nData]=readOH([file '.ecg']);
[OHData(:,6), tVec, nData]=readOH([file '.amx']);
[OHData(:,7), tVec, nData]=readOH([file '.amy']);
[OHData(:,8), tVec, nData]=readOH([file '.amz']);
[event1Index, event2Index]=readEvent2([file]);
%[event1Index, event2Index]=readEvent([file]);

%[b, a]=butter(3,10/125);
%OHData=filtfilt(b, a,OHData);

%figure; plot(tVec, OHData(:,7:8));

figure;
for ii=1:8
    subplot(2,4,ii); plot(tVec, OHData(:,ii));
    xlabel('time(sec)');
    title(['chan ' num2str(ii)]);
    %axis([1 tVec(end) 0 4096]);
end;

OHData(:,1:4)=medfilt2(OHData(:,1:4),[5 1]);

figure;
for ii=1:8
    subplot(2,4,ii); plot(tVec, OHData(:,ii)/4096*5);
    xlabel('time(sec)');
    title(['chan ' num2str(ii)]);
    %axis([1 tVec(end) 0 5]);
end;

ECG=OHData(:,5);
% [b, a]=butter(3,50/125);
% [bH, aH]=butter(3,0.02/125,'high');
% ECGRaw=OHData(:,5);
% ECGFil=filtfilt(bH, aH, ECGRaw);
% ECG=-filtfilt(b,a,ECGFil);
%figure; plot(tVec, ECG);

RespRaw=ECG-median(ECG);
[b,a]=butter(3,0.5/125);
Resp=filtfilt(b,a,RespRaw);
%figure; plot(tVec, Resp);
figure; plot(tVec, ECG-mean(ECG));
hold on; plot(tVec, Resp, 'r')

function [data, tVec, nData]=readOH(fileNm)
sampRate=250;
fid = fopen(fileNm, 'rb');
[data, nData] = fread(fid,'short');
fclose(fid);
tVec=(1:nData)/250;
tVec=tVec-tVec(1);



