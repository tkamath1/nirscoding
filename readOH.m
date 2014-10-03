function [data, tVec, nData]=readOH(fileNm)
sampRate=250;
fid = fopen(fileNm, 'rb');
[data, nData] = fread(fid,'short');
fclose(fid);

tVec=(1:nData)/250;
tVec=tVec-tVec(1);
