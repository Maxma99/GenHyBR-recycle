function out=inverse_x(x)
len = size(x,1);
for i = 1:len
    out(i,:,:) = x(len+1-i,:,:);
end
end