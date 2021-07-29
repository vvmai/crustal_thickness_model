function [x,y] = calc_err(a,da,b,db)

x = [];
y = [];
%y error bars
for i=1:length(a)
  x = [x a(i) a(i) NaN];
  y = [y b(i)-db(i) b(i)+db(i) NaN];
end
%x error bars
for i=1:length(a)
  x = [x a(i)-da(i) a(i)+da(i) NaN];
  y = [y b(i) b(i) NaN];
end
