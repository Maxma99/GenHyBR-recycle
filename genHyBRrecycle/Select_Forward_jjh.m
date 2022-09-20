function [A,bn,R,nlevel] = Select_Forward_jjh(prob,n,x_true,nlevel,n_t)

switch prob
    case 1
        % Sheprical tomography from IRtool
        %     angles = linspace(0,90,n/2+1);
        %     angles = angles(2:end);
        %     opt = PRset('angles',angles,'numCircles',90);
        %     [A,~,~,ProbInfo] = PRspherical(n,opt);
        angles = linspace(0,360,n_t+1);
        %         angles = linspace(0,90,n/2+1);
        angles = angles(2:end);
        opt = PRset('angles',angles,'numCircles',90);
        [A,~,~,ProbInfo] = PRspherical(n,opt);
    case 2
        % raytracing
        nx = n; ny = n;
%    ns =round(n/2);    
     ns =round(floor(n/2));    
    nr = round(floor(n/2));
    [~,~,A] = raytracing(nx,ny,ns, nr);
    case 3
    % PSF
    [PSF, center] = getPSF(6, n);
    A = psfMatrix(PSF,center,'reflexive');
end

b = A*x_true(:);
b = b(:);

%Noise covariance
m = size(b,1);
R = speye(m,m);

%Adding noise to observation (As nlevle->0 , approx sol get closer to nonapprox sol)
% nlevel = 0.01;
[N,sigma] = WhiteNoise(b(:),nlevel);
bn = b + N;