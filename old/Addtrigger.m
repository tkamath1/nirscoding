function [x] = Addtrigger(Triggerlogname, PVTlogname)
   Triggerlogdata = importdata(Triggerlogname);
   PVTdata = importdata(PVTlogname);
   
   
   for i = 1:size(Triggerlogdata,1)
       for j = 2:length(PVTdata.textdata(:,3))
           [Hour,Sec] = strtok(PVTdata.textdata(j,3),':');
           z = char(Sec);
           y = z(2:length(z));
           
           if str2num(char(Hour)) == Triggerlogdata(i,1)
            if str2num(y) == Triggerlogdata(i,2)
              PVTdata.textdata(j,3) = cellstr(strcat(char(PVTdata.textdata(j,3)),':',num2str(Triggerlogdata(i,3))))
               end
           end
       end
   end
   x = PVTdata.textdata(:,3);
end
   