function exp = getExponent(value)
%
%
%
if ~isfloat(value)
   exp = 1;
   return;
end

string = num2str(value);
e = 1;
for i = 1:size(string,2)
    if strcmp(string(1,i),'e')
        e = string(1,(i+1):end);
        break;
    end
end
if e==1
    return;
end
exp = str2double(['1e' e]);

end