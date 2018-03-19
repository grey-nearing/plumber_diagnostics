function [y,d,h] = convert_time(time,y1)
 
for y = 0:15

 ndays = 365;
 if rem(y+y1,4)==0; ndays = 366; end;
 
 if time > ndays*24*60^2
  time = time - ndays*24*60^2;
  continue
 end
 
 d = floor(time/24/60^2)+1;
 if d < 1; keyboard; end;
 time = time - (d-1)*(24*60^2);

 h = round(time/60^2 * 2)/2;
 
 break
end

if d > ndays; y = y+1; d = 1; h = 0; end;
 
y = y + y1;
