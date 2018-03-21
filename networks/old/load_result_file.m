function T = load_result_file(fname)

T = zeros(10,5);

fid = fopen(fname,'r');
d = fgets(fid);

for i = 1:10
 d = fgets(fid);
 a = strsplit(d); a(1) = [];
 for j = 1:5
  T(i,j) = str2num(a{j});
 end
end

T = T/1000;

fclose(fid);



