% The best codes handle N = 20,000 as long as the matrix is very sparse.
% N   = 3000; M = 4000; % Large scale. Things start to get interesting

weighting_Matrix_s = load(ff.weighting_Matrix);
A = weighting_Matrix_s.Jnew;
clear weighting_Matrix_s
[M, N] = size (A);
    [pixel_n laser_n] = size(ff.det_Em);
    measure_array = zeros(pixel_n*laser_n,1);
    nsol = size(A,2);
    for i= 1:laser_n
        measure_array(((i-1)*pixel_n+1):(pixel_n*i))=ff.det_Em(:,i)./ff.det_Ex(:,i);
    end
    if ff.noise.option == 1
        mean_m_array = mean(measure_array);
        measure_noise = rand(pixel_n*laser_n,1).*ff.noise.level.*mean_m_array;
        measure_array = measure_array + measure_noise;
    end
    if ~isempty(ff.recon_mask)
        tmp_measure_array = measure_array(ff.recon_mask);
        clear measure_array
        measure_array = tmp_measure_array;
        clear tmp_measure_array
    end
b = measure_array ;
%N   = 1000; M = 1500;     % at this size, some algo take a long time!
% N   = 100; M = 150;     % at this size, all algorithms take < 14 seconds
%A   = randn(M,N);
%b   = randn(M,1);

fcn     = @(x) norm( A*x - b)^2;
% here are two equivalent ways to make the gradient. grad2 is sometimes faster
grad1    = @(x) 2*A'*(A*x-b);
AtA     = A'*A; Ab = A'*b;
grad2    = @(x) 2*( AtA*x - Ab );

grad    = grad2;


x = [];
time = [];

%% Solve NNLS with L-BFGS-B

l  = zeros(N,1);    % lower bound
u  = inf(N,1);      % there is no upper bound
tstart=tic;
fun     = @(x)fminunc_wrapper( x, fcn, grad); 
% Request very high accuracy for this test:
opts    = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 10);
opts.printEvery     = 5;
if N > 10000
    opts.m  = 50;
end
% Run the algorithm:
[xk, ~, info] = lbfgsb(fun, l, u, opts );
t=toc(tstart)
% Record results
x.lbfgsb    = xk;
time.lbfgsb = t;