function saveCircos(T,varNames,Du,D,fname)

T = round(T*1000);
T(isnan(T)) = 0;

Tstring = cell(D+1,D-Du+1);
for x = 1:D
 for y = Du+1:D
  Tstring{x+1,y-Du+1} = num2str(T(x,y));
 end
end

Tstring(1,2:end) = varNames(Du+1:end);
Tstring(2:end,1) = varNames;
Tstring(1,1) = {'labels'};

fid = fopen(fname,'w');
format=repmat('%s ',[1,D-Du+1]);
format = strcat(format,'\n');
for x = 1:D+1
 wstring = Tstring(x,:);
 fprintf(fid,format,wstring{1:end});
end
fclose(fid);


