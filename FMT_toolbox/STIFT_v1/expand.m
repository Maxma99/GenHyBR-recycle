function x_expand = expand(x_display)
grd = [55,55,15];
x_expand = zeros(grd);
for i=1:grd(1)
    for j = 1:grd(2)
        for k = 1:grd(3)
            if(floor(i/6)>0 && floor(j/6)>0 && floor(k/6)>0)
            x_expand(i,j,k) = (x_display(ceil(i/6),ceil(j/6),ceil(k/6))+ x_display(floor(i/6),floor(j/6),floor(k/6)))/2;
        end
    end
end
end