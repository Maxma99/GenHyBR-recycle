% v1
function obj = inversSTIFT_1(obj)    
    fprintf(1,'Applying inversion ');
    tic
    %% loading weighting matrix
    weighting_Matrix_s = load([obj.data_buffer_directory '/weighting_Matrix.mat']);
    weighting_Matrix = weighting_Matrix_s.Jnew;
    clear weighting_Matrix_s
    %%%%%
    % reconstruction parameter setting
    if ~isempty(obj.fluorecon)
        clear obj.fluorecon;
    end
    iterationTime = obj.reconSetting.iterTime;
    % regularisation proportional to trace
    lambda = obj.reconSetting.lambdaFactor*sum(weighting_Matrix(:).^2); 
    gamma = 0.3; %need to be specified beforehand for case 11!
    lambda_TV=1e-6;
    tolerror = obj.reconSetting.tolerror;
    %%%%%%%%%%%%%%%%%%% phamton

    [pixel_n,laser_n] = size(obj.det_Em);
    measure_array = zeros(pixel_n*laser_n,1);
    nsol = size(weighting_Matrix,2);
    for i= 1:laser_n
        measure_array(((i-1)*pixel_n+1):(pixel_n*i))=obj.det_Em(:,i)./obj.det_Ex(:,i);
        %measure_array(((i-1)*pixel_n+1):(pixel_n*i))=obj.det_Em(:,i)./obj.det_Ex_recon(:,i);
    end
    if ~isempty(obj.mea_mask)
        tmp_measure_array = measure_array(obj.mea_mask);
        clear measure_array
        measure_array = tmp_measure_array;
        clear tmp_measure_array
    end
    save([obj.data_buffer_directory '/measure_array.mat'],'measure_array','-v7.3');
    
    obj.measure_array = [obj.data_buffer_directory '/measure_array.mat'];
    obj.geneDisposeMeaMaskSTIFT;
    if ~isempty(obj.dispose_mea_mask)
        tmp_measure_array = measure_array(obj.dispose_mea_mask);
        tmp_weighting_Matrix = weighting_Matrix(obj.dispose_mea_mask,:);
        clear measure_array
        clear weighting_Matrix
        measure_array = tmp_measure_array;
        weighting_Matrix = tmp_weighting_Matrix;
        clear tmp_measure_array
        clear tmp_weighting_Matrix
    end
    save([obj.data_buffer_directory '/measure_array.mat'],'measure_array','-v7.3');
    save([obj.data_buffer_directory '/weighting_Matrix.mat'],'weighting_Matrix','-v7.3');
    for i = 1:size(weighting_Matrix,2)
        if(norm(weighting_Matrix(:,i))~=0)
            weighting_Matrix(:,i) = weighting_Matrix(:,i)/norm(weighting_Matrix(:,i));
        end
    end
    switch obj.reconSetting.Method
        case 1 % ART without regularization
            fprintf(1,'using ART without regularization...\n');
            obj.fluorecon = kaczmarz(weighting_Matrix,measure_array,iterationTime);
        case 2 % ART with regularization
            fprintf(1,'using ART with regularization...\n');
            tr = weighting_Matrix'*measure_array;
            tJ = weighting_Matrix'*weighting_Matrix+lambda*eye(nsol);
            obj.fluorecon = kaczmarz(tJ,tr,iterationTime);
        case 3 % MATLAB simple backslash without regularization
            fprintf(1,'using backslash without regularization...\n');
            obj.fluorecon = weighting_Matrix\measure_array;
        case 4 % MATLAB simple backslash with regularization
            fprintf(1,'using backslash with regularization...\n');
            tr = weighting_Matrix'*measure_array;
            tJ = weighting_Matrix'*weighting_Matrix+lambda*eye(nsol);
            obj.fluorecon = tJ\tr;
        case 5 % CG method: pcg without regularization
            fprintf(1,'CG method: pcg without regularization...\n');
            nr = size(weighting_Matrix,1);
            RR = spdiags(1./(sum(weighting_Matrix,2)),0:0,nr,nr);
            RJr = RR*weighting_Matrix;  
            tr = RJr'*RR*measure_array;
            %%%%%%% simon's code %%%%
            %tJ = RJr'*RJr+lambda*eye(nsol);
            %Rlambda = 1e-4*(trace(tJ));
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            tJ = RJr'*RJr;
            [obj.fluorecon,~,~,niter] = pcg(tJ,tr,tolerror,iterationTime);
            disp(['CG iterations: ',num2str(niter)]);
        case 6 % CG method #1: pcg with regularization
            fprintf(1,'CG method #1: pcg with regularization...\n');
            nr = size(weighting_Matrix,1);
            RR = spdiags(1./(sum(weighting_Matrix,2)),0:0,nr,nr);
            RJr = RR*weighting_Matrix;     
            Rlambda =  obj.reconSetting.lambdaFactor*sum( RJr(:).^2);
            tr = RJr'*RR*measure_array;          
            tJ = RJr'*RJr+Rlambda*eye(nsol);
            [obj.fluorecon,~,~,niter] = pcg(tJ,tr,tolerror,iterationTime);
            disp(['CG iterations: ',num2str(niter)]);
         case 7 % lsqr without regularization lsqr:least spuare QR-factorization
            fprintf(1,'lsqr without regularization...\n');
            nr = size(weighting_Matrix,1);
            RR = spdiags(1./(sum(weighting_Matrix,2)),0:0,nr,nr); %normalize
            [obj.fluorecon,~,~,niter] = lsqr(RR*weighting_Matrix, RR*measure_array,tolerror,iterationTime);
            disp(['lsqr iterations: ',num2str(niter)]);
         case 8 % lsqr lsqr with regularization
            fprintf(1,'lsqr with regularization...\n');
            nr = size(weighting_Matrix,1);
            RR = spdiags(1./(sum(weighting_Matrix,2)),0:0,nr,nr);
            RJr = RR*weighting_Matrix;
            Rlambda = obj.reconSetting.lambdaFactor*sum( RJr(:).^2);
            [obj.fluorecon,~,~,niter] = lsqr([RJr; sqrt(Rlambda)*speye(nsol)], [RR*measure_array;zeros(nsol,1)],tolerror,iterationTime);
            disp(['lsqr iterations: ',num2str(niter)]);
            %obj.fluorecon =  lsqnonneg([RJr; sqrt(Rlambda)*speye(nsol)], [RR*measure_array;zeros(nsol,1)]);
         case 9 % l1 regularization with ISTA
             fprintf(1,'l1 regularization with ISTA...\n');
             obj.fluorecon = l1I(weighting_Matrix, measure_array, lambda, iterationTime);
             disp(['l1 using ISTA iterations:', num2str(iterationTime)]);
         case 10 % l1 regularization with FISTA
             fprintf(1,'l1 regularization with FISTA...\n');
             obj.fluorecon = l1F(weighting_Matrix,measure_array,lambda,iterationTime);
             disp(['l1 using FISTA iterations:', num2str(iterationTime)]);
         case 11 %joint l1 and l2 regularization(lasso)
             fprintf(1,'joint l1 and l2 regularization...\n');
             obj.fluorecon = l1_l2(weighting_Matrix, measure_array, lambda, iterationTime, gamma);
             disp(['joint l1 and l2 regularization:', num2str(iterationTime)]);             
         case 12 %joint l1 and TV reularization
             fprintf(1, 'joint l1 and TV regularization...\n');
             obj.fluorecon = l1_TV(weighting_Matrix, measure_array, lambda, lambda_TV, iterationTime, obj);
             disp(['joint l1 and TV regularization:', num2str(iterationTime)]);
        case 13 % non negative least square without regularization
            fprintf(1, 'non negative least square without regularization...\n');
            A = weighting_Matrix;
            b = measure_array;
            [~, N] = size (A);
            fcn     = @(x) norm( A*x - b)^2;
            % here are two equivalent ways to make the gradient. grad2 is sometimes faster
            %grad1    = @(x) 2*A'*(A*x-b);
            AtA     = A'*A; Ab = A'*b;
            grad2    = @(x) 2*( AtA*x - Ab );
            grad    = grad2;
            %% Solve NNLS with L-BFGS-B
            l  = zeros(N,1);    % lower bound
            u  = inf(N,1);      % there is no upper bound
            %tstart=tic; 
            fun = @(x)fminunc_wrapper( x, fcn, grad); 
            % Request very high accuracy for this test:
            opts = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 10); %
            % defalt
            opts.printEvery = 5;
            if N > 10000
                opts.m  = 50;
            end
            % Run the algorithm:
            [obj.fluorecon, ~, ~] = lbfgsb(fun, l, u, opts );
            %t=toc(tstart)
        case 14 % non negative least square with regularization
            fprintf(1, 'non negative least square with regularization\n');
            nr = size(weighting_Matrix,1);
            RR = spdiags(1./(sum(weighting_Matrix,2)),0:0,nr,nr);
            RJr = RR*weighting_Matrix; 
            Rlambda = obj.reconSetting.lambdaFactor*sum( RJr(:).^2);
            A = [RJr; sqrt(Rlambda)*speye(nsol)];
            b = [RR*measure_array;zeros(nsol,1)];
            [~, N] = size (A);
            fcn     = @(x) norm( A*x - b)^2;
            % here are two equivalent ways to make the gradient. grad2 is sometimes faster
            %grad1    = @(x) 2*A'*(A*x-b);
            AtA     = A'*A; Ab = A'*b;
            grad2    = @(x) 2*( AtA*x - Ab );
            grad    = grad2;
            %% Solve NNLS with L-BFGS-B
            l  = zeros(N,1);    % lower bound
            u  = inf(N,1);      % there is no upper bound
            %tstart=tic; 
            fun = @(x)fminunc_wrapper( x, fcn, grad); 
            % Request very high accuracy for this test:
            opts = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 10);
            opts.printEvery = 5;
            if N > 10000
                opts.m  = 50;
            end
            % Run the algorithm:
            [obj.fluorecon, ~, ~] = lbfgsb(fun, l, u, opts );
            %t=toc(tstart)
        case 15
            fprintf(1, 'genHyBR recycle...\n');
            A = weighting_Matrix;
            b = measure_array;
            solver = 'tikhonov';
            % Set trunction options and matrices
            trunc_options.nOuter = 15; % number of outer iterations
            trunc_options.nInner = 25; % maximum storage of solution vector space
            trunc_options.max_mm = 20;  % maximum number of vectors to save at compression
            trunc_options.compress = 'SVD'; 
            trunc_mats = [];
%             sigma = 0.1;
            input = HyBRset('InSolv', solver, 'Iter', iterationTime, 'ResTol', [tolerror, tolerror]); 
            % xmin = [0 0];         %Coordinates of left corner
            % xmax = [1 1];         %Coordinates of right corner
            % nvec = [n, n];        %Number of points in grid
            % 
            % nu = 1.5; ell = .001;
            % k = @(r) matern(r,nu,ell);
            % %Additional parameters governing length scales. Can be used to model
            % %nonisotropic covariance kernels
            % theta = [1.0 1.0];      %For now set them as isotropic
            % 
            % % Build row/column of the Toeplitz matrix
            % Qr = createrow(xmin,xmax,nvec,k,theta);
            % Qfun = @(x)toeplitzproduct(x, Qr, nvec);
            % Q = funMat(Qfun,Qfun, nvec.^2);

            % % Set Q as identity matrix
            Q = speye(size(A,2));
            
            % (b) Set R as identity matrix
            R = speye(size(A,1));
            [obj.fluorecon,~,~] = genHyBRrecycle(A,b,[],Q,R,input,trunc_options,trunc_mats); 

    end
    fprintf(1,'done\n');
    toc
end