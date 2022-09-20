function x_true = exp_dynamic_2D_FMT(T, n)
grd = [n,n,T];
fltgt_g = zeros(grd);
nblobs = 2;
for l = 1:T
    fcontrast(l,:) = [1-0.01*l 1.15-0.015*l];
    blbrad(l,:) = grd(1)*[0.03+l*0.008 0.03+l*0.005]; % radius of blobs
end

blbcx = grd(1)*[0.7 0.3]; % xcentre of blobs
blbcy = grd(2)*[0.3 0.7]; % ycentre of blobs

for l = 1:T
    for i = 1:grd(1)
      for j = 1:grd(2)
          for k = 1:nblobs
            if( (i-blbcx(k))^2 + (j-blbcy(k))^2 < blbrad(l,k)^2)
                fltgt_g(i,j,l) = fltgt_g(i,j,l)+fcontrast(l,k);
            end
          end
      end
    end
end

x_true = fltgt_g;
end