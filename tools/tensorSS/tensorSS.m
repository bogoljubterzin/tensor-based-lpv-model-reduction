classdef tensorSS
    % Tensor representation of LPV models (affine) using Tensor Toolbox 
    % for MATLAB by Sandia National Laboratories.
    % Tensor Decompositions: HOSVD - Tensor Toolbox for MATLAB
    % TSVD and succR1 - by prof. Siep Weiland

    properties
        A, B, C, D, eta_map, Ts
        domain, Nx, Ny, Nu, Np
    end
    methods(Static)

        % Constructor
        % Returns a tensorSS object that is constructed either from:
        % 1) Directly from tensors A,B,C,D, eta_map, Ts
        % obj = tensorSS(A,B,C,D,eta_map,Ts);
        % 2) From LPV-SS model and scheduling map 
        % obj = tensorSS(lpvss,eta_map)
        function obj = tensorSS(varargin)
            if nargin == 0
                return;
            elseif nargin == 6 || nargin == 5
                obj.A = tensor(varargin{1});
                obj.B = tensor(varargin{2});
                obj.C = tensor(varargin{3});
                obj.D = tensor(varargin{4});
                obj.eta_map = varargin{5};
                if nargin == 5 || varargin{6} == 0
                    obj.domain = 'ct';
                    obj.Ts = 0;
                else
                    obj.domain = 'dt';
                    obj.Ts = varargin{6};
                end    

            elseif nargin == 2 % tensorSS object from lpvss object
                if ~isa(varargin{1},'LPVcore.lpvss')
                    error('First argument must be an LPVSS object');
                end
                lpvssModel = varargin{1};
                eta_map = varargin{2};
                obj.A = tensor(lpvssModel.A.matrices);
                obj.B = tensor(lpvssModel.B.matrices);
                obj.C = tensor(lpvssModel.C.matrices);
                obj.D = tensor(lpvssModel.D.matrices);
                obj.Ts = lpvssModel.Ts;
                obj.eta_map = eta_map.map; % eta_map.map !!!
                % obj.eta_map = eta_map; % for PCA i used this
                if obj.Ts == 0
                    obj.domain = 'ct';
                else
                    obj.domain = 'dt';
                end
            end

            obj.Nx = size(obj.A,1);
            if numel(size(obj.A))==2
                obj.Np = 0;
            else
                obj.Np = size(obj.A,3)-1;
            end
            obj.Nu = size(obj.B,2);
            obj.Ny = size(obj.C,1);
        end
        
        function tensorObj = fromLPVSS(lpvssModel,eta_map)
            % Construct tensorSS object from affine lpvss object
            A_data = lpvssModel.A.matrices;
            B_data = lpvssModel.B.matrices;
            C_data = lpvssModel.C.matrices;
            D_data = lpvssModel.D.matrices;
            Ts = lpvssModel.Ts;
            % Create tensorSS object
            tensorObj = tensorSS(A_data, B_data, C_data, D_data, eta_map,Ts);
        end
        
    end
    methods
        function [tout,xout,yout,pout] = simulateSS(obj,u,tin,x0,constantTerm)
            % Tensor-based simulation - however - more efficient if we
            % convert tensorSS into LPVSS and then simulate LPVSS
            
            % Simulation of tensorSS given the input trajectory
            % A(p)*x is calculated as Atens * (x,p)
            % constantTerm is important to differentiate the cases when
            % the  model contains the constant part, and when it doesn't
            % Careful about constantTerm! if constantTerm = 1, p0 = 1; else p0 ~= 1
            
            if nargin < 5
                constantTerm = 1;
            end
            if constantTerm ~= 0
                constantTerm = 1;
            end
            if ismatrix(obj.A) 
                % if np = 1, A,B,C,D are transformed to be tensors instead
                % of matrices
                obj.A = reshape(obj.A,[size(obj.A),1]);
                obj.B = reshape(obj.B,[size(obj.B),1]);
                obj.C = reshape(obj.C,[size(obj.C),1]);
                obj.D = reshape(obj.D,[size(obj.D),1]);
            end
            if length(u) ~= length(tin)
                error("u and tin must be of same length");
            elseif size(x0,1) ~= size(obj.A,1)
                error("x0 should be of dimension (Nx,1)");
            end

            if obj.domain == "dt"
                N = length(tin);
                xout = [x0, zeros(size(x0,1),N)];
                yout = zeros(size(obj.C,1),N);
                pout = zeros(obj.Np+double(constantTerm==0),N);
                for i = 1:N
                    if i == 55
                        fprintf("Come here")
                    end
                    if constantTerm ~= 1
                        pi = obj.eta_map(xout(:,i),u(:,i));
                        pout(:,i) = pi;
                    else
                        if size(obj.eta_map(xout(:,i),u(:,i)),1) == 1
                            pt = obj.eta_map(xout(:,i),u(:,i)).';
                        else
                            pt = obj.eta_map(xout(:,i),u(:,i));
                        end
                        pout(:,i) = pt;
                        pi = [1; pt];
                    end
                    UA{1,1} = xout(:,i); UA{1,2} = pi;
                    UB{1,1} = u(:,i); UB{1,2} = pi;
                    UC = UA; UD = UB;
                    xout(:,i+1) = ttm(obj.A,UA,[2,3],'t') + ttm(obj.B,UB,[2,3],'t');
                    yout(:,i) = ttm(obj.C, UC,[2,3], 't') + ttm(obj.D, UD,[2,3],'t');
                end
                tout = tin;
                xout = xout(:,1:end-1);

            elseif obj.domain == "ct"
                ufun = @(tau) interp1(tin,u,tau,'previous');
                xdot = @(t,x) xfun(x,ufun(t),obj);
                [tout,xout] = ode45(xdot,tin,x0);
                [yout, pout] = hfun(xout,ufun(tout),obj);
                

            end

            function xdot = xfun(x,u,obj)
                p = [1; obj.eta_map(x,u)];
                mode = [2,3];
                UA{1,1} = x'; UA{1,2} = p';
                UB{1,1} = u'; UB{1,2} = p';
                xdot = double(ttm(obj.A,UA,mode) + ttm(obj.B,UB,mode));
            end
            
            function [yout,pout] = hfun(x,u,obj)
                yout = zeros(obj.Ny,length(x));
                pout = zeros(obj.Np,length(x));
                mode = [2,3];
                for i = 1:size(x,1)
                    p = [1; obj.eta_map(x(i,:)',u(i,:)')];
                    UC{1,1} = x(i,:); UC{1,2} = p'; 
                    UD{1,1} = u(i,:); UD{1,2} = p';
                    yout(:,i) = double(ttm(obj.C,UC,mode) + ttm(obj.D,UD,mode));
                    pout(:,i) = p(2:end);
                end
            end

            
        end
        
        function lpvss = tensSS2lpvss(obj)
            % tensorSS to LPVSS converter
            Ap = double(obj.A(:,:,1));
            Bp = double(obj.B(:,:,1));
            Cp = double(obj.C(:,:,1));
            Dp = double(obj.D(:,:,1));

            for i = 1:obj.Np
                pi = preal(sprintf('pr(%i)', i), obj.domain);
                Ap = Ap + double(obj.A(:,:,i+1))*pi;
                Bp = Bp + double(obj.B(:,:,i+1))*pi;
                Cp = Cp + double(obj.C(:,:,i+1))*pi;
                Dp = Dp + double(obj.D(:,:,i+1))*pi;
            end

            lpvss = LPVcore.lpvss(Ap,Bp,Cp,Dp,obj.Ts);
        end
        
        function [tensSSr, V,sv] = PCA_coef_matrices(obj, train_data, decomp, rx, rp)
            % Tensor-based PCA approach for model reduction: States are
            % reduced the same as in "POD_reduction" method, while the
            % scheduling variables are kept the same. After that, PCA
            % scheduling reduction method is used to find a reduced affine
            % dependancy based on the coefficients of the system matrices
            % throughout time
            % 
            % Returns a reduced-order model "sysr" as tensorSS
            % object
            % Inputs:
            %   u_train: input signal that is used for gathering data
            %            size: Nu x N, type: double
            %   decomp: string that determines the decomposition
            %            "TSVD","succR1","HOSVD"
            %   rx: ROM's number of state (rx <= nx) type: integer
            %   rp: ROM's number of scheduling (rp <= np) type: integer
            O = obj.Nx + obj.Ny;
            I = obj.Nx + obj.Nu;
            
            u_train = train_data.input;
            x0 = train_data.x0;
            tin = 0:obj.Ts:(length(u_train)-1)*obj.Ts;
            N = size(u_train,2);

            [~,~,~,p1_train] = obj.simulateSS(u_train,tin,x0);

            data = zeros(O,I,N);
            for i = 1:N
                p = [1; p1_train(:,i)];
                Ap = double(ttv(obj.A,p,3));
                Bp = double(ttv(obj.B,p,3));
                Cp = double(ttv(obj.C,p,3));
                Dp = double(ttv(obj.D,p,3));
                Sp = [Ap,Bp;Cp,Dp];
                data(:,:,i) = double(Sp);
            end
            R = max(rx,rp);
            sv = zeros(R,1);

            decomp = lower(decomp);
            switch decomp
                case "tsvd"
                    [~,x2,~,~,~,~]=TSVD3D_SW(data,R);
                    V = x2(1:obj.Nx,1:rx);
                    W = V;
                    Z = eye(obj.Np);

                    S = TSVcore(data,R);                
                    for i = 1:R
                        sv(i) = S(i,i,i);
                    end
                case "succr1"
                    [~,x2,~,~,~,~]=succR1_SW(data,R);
                    V = x2(1:obj.Nx,1:rx);
                    W = V;
                    Z = eye(obj.Np);
                case "hosvd"
                    T = hosvd(tensor(data),2*sqrt(eps));
                    V = T{2}(1:obj.Nx,1:rx);
                    W = V;
                    Z = eye(obj.Np);
                otherwise
                    error("Decomposition unsupperted. Choose TSVD, succR1, or HOSVD (case insensitive)!");
            end
            tensSSr1 = obj.PetrovGalerkinLPV(W,V,Z,1); % tensSS model with reduced state order (rx) and full scheduling order (np)
            lpvssr1 = tensSSr1.tensSS2lpvss(); % lpvss model with reduced state order (rx) and full scheduling order (np)
            etar = tensSSr1.eta_map;

            [lpvssr2, mapping_fnc2,~] = lpvpcared(lpvssr1, rp, p1_train.','trajectory');
            etar2 = @(x,u) mapping_fnc2(etar(x,u).'); 

            tensSSr = tensorSS.fromLPVSS(lpvssr2,etar2);
        end

        function [tensSSr, V, info] = POD4D(obj, train_data, decomp, rx, rp)
            decomp = lower(decomp);
            u_train = train_data.input;
            x0 = train_data.x0;
            tin = 0:obj.Ts:(length(u_train)-1)*obj.Ts;
            [~,x,~,p] = obj.simulateSS(u_train,tin,x0);
            xdot = x(:,2:end) - x(:,1:end-1);

            data = zeros(obj.Nx, obj.Nx, obj.Np, length(tin)-1);
            for i = 1:length(tin)-1
                t1 = ttt(tensor(x(:,i)),tensor(xdot(:,i)));
                t1 = squeeze(t1);
                t2 = ttt(t1,tensor(p(:,i)));
                t2 = squeeze(t2);
                data(:,:,:,i) = t2;
            end

            R = max(rx,rp);
            switch decomp
                case "tsvd"
                    tic
                    [x,gains,iterations,Xa]=TSVDND_SW(data,R);
                    tcomp = toc;
                    V = x{1}(:,1:rx);
                    W = x{2}(:,1:rx);
                    Z = x{3}(:,1:rp);

                    decomp_info.gains = gains;
                    decomp_info.iterations = iterations;
                    decomp_info.low_rank_app = Xa;
                case "succr1"
                    tic
                    [x1,x2,x3,gains,iterations,Xa]=succR1_SW(data,R);
                    tcomp = toc;
                    V = x1(:,1:rx);
                    W = V;
                    Z = x2(:,1:rp);

                    decomp_info.gains = gains;
                    decomp_info.iterations = iterations;
                    decomp_info.low_rank_app = Xa;
                case "hosvd"
                    tic
                    T = hosvd(tensor(data),2*sqrt(eps));
                    tcomp = toc;
                    V = T{1}(:,1:rx);
                    W = T{2}(:,1:rx);
                    Z = T{3}(:,1:rp);
                    
                    % info
        
                    decomp_info.S = T.core;
                    decomp_info.T = T;
                    
            end

            info.tcomp = tcomp;
            info.decomp = decomp_info;

            tensSSr = obj.PetrovGalerkinLPV(W,V,Z,1);
        end
        function [tensSSr,V, info] = POD_reduction(obj, train_data, decomp, rx, rp)
            % Tensor-based POD approach for joint reduction of states and
            % scheduling. Returns a reduced-order model "sysr" as tensorSS
            % object
            % Inputs:
            %   u_train: input signal that is used for gathering data
            %            size: Nu x N, type: double
            %   decomp: string that determines the decomposition
            %            "TSVD","succR1","HOSVD"
            %   rx: ROM's number of state (rx <= nx) type: integer
            %   rp: ROM's number of scheduling (rp <= np) type: integer
            decomp = lower(decomp);

            u_train = train_data.input;
            x0 = train_data.x0;
            tin = 0:obj.Ts:(length(u_train)-1)*obj.Ts;

            [~,x,~,p] = obj.simulateSS(u_train,tin,x0);

            data = zeros(obj.Nx,obj.Np,length(tin));
            for i = 1:length(tin)
                data(:,:,i) = x(:,i)*p(:,i).';
            end
            R = max(rx,rp);
            switch decomp
                case "tsvd"
                    tic
                    [x1,x2,x3,gains,iterations,Xa]=TSVD3D_SW(data,R);
                    tcomp = toc;
                    V = x1(:,1:rx);
                    W = V;
                    Z = x2(:,1:rp);
                    
                    decomp_info.x1 = x1;
                    decomp_info.x2 = x2;
                    decomp_info.x3 = x3;
                    decomp_info.gains = gains;
                    decomp_info.iterations = iterations;
                    decomp_info.low_rank_app = Xa;
                case "succr1"
                    tic
                    [x1,x2,x3,gains,iterations,Xa]=succR1_SW(data,R);
                    tcomp = toc;
                    V = x1(:,1:rx);
                    W = V;
                    Z = x2(:,1:rp);


                    decomp_info.x1 = x1;
                    decomp_info.x2 = x2;
                    decomp_info.x3 = x3;
                    decomp_info.gains = gains;
                    decomp_info.iterations = iterations;
                    decomp_info.low_rank_app = Xa;
                case "hosvd"
                    tic
                    T = hosvd(tensor(data),2*sqrt(eps));
                    tcomp = toc;
                    V = T{1}(:,1:rx);
                    W = V;
                    Z = T{2}(:,1:rp);
                    
                    % info
                    
                    decomp_info.S = T.core;
                    decomp_info.x1 = T{1};
                    decomp_info.x2 = T{2};
                    decomp_info.x3 = T{3};

                    
            end

            info.tcomp = tcomp;
            info.decomp = decomp_info;

            tensSSr = obj.PetrovGalerkinLPV(W,V,Z,1);
        end
        function [sysr,V] = coef_level_reduction(obj, decomp, type, rx, rp)
            % Tensor-based coefficient level joint reduction of LPV
            % models. Returns a reduced-order model "sysr" as tensorSS
            % object
            % Inputs:
            %   decomp - string that determines the decomposition
            %            "TSVD","succR1","HOSVD",'CP'
            %   type - indicates what is used for the data (all values or
            %   just values from tensor A)
            %       type = 1 - use only tensor A
            %       type = 2 - use [A B; C D]
            %   rx: ROM's number of state (rx <= nx) type: integer
            %   rp: ROM's number of scheduling (rp <= np) type: integer
            
            decomp = lower(decomp);
            switch type
                case 1
                    coef = double(obj.A);
                    R = max(rx,rp);
                    if decomp == "tsvd"
                        [~,x2,x3,~,~,~]=TSVD3D_SW(coef,R);
                        V = x2(:,1:rx); W = V; Z = x3(:,1:rp);
                    elseif decomp == "succr1"
                        [~,x2,x3,~,~,~]=succR1_SW(coef,R);
                        V = x2(:,1:rx); W = V; Z = x3(:,1:rp);
                    elseif decomp == "hosvd"
                        T = hosvd(tensor(coef),2*sqrt(eps));
                        V = T{2}(:,1:rx);
                        W = V;
                        Z = T{3}(:,1:rp);
                    elseif decomp == "cp"
                        T = cp_als(tensor(coef),R);
                        V = T{2}(:,1:rx); W = V; Z = T{3}(:,1:rp);
                    end
                case 2
                    AB = cat(2,double(obj.A),double(obj.B));
                    CD = cat(2,double(obj.C),double(obj.D));
                    coef = cat(1,AB,CD);
                    R = max(rx,rp);
                    if decomp == "tsvd"
                        [~,x2,x3,~,~,~]=TSVD3D_SW(coef,R);
                        V = x2(1:obj.Nx,1:rx); W = V; Z = x3(1:obj.Np+1,1:rp);
                    elseif decomp == "succr1"
                        [~,x2,x3,~,~,~]=succR1_SW(coef,R);
                        V = x2(1:obj.Nx,1:rx); W = V; Z = x3(:,1:rp);
                    elseif decomp == "hosvd"
                        T = hosvd(tensor(coef),2*sqrt(eps));
                        V = T{2}(1:obj.Nx,1:rx);
                        W = V;
                        Z = T{3}(:,1:rp);
                    elseif decomp == "cp"
                        T = cp_als(tensor(coef),R);
                        V = T{2}(1:obj.Nx,1:rx); W = V; Z = T{3}(:,1:rp);
                    end
            end
            sysr = PetrovGalerkinLPV(obj,W,V,Z,0);
        end

        

        function tensSSr = PetrovGalerkinLPV(obj,W,V,Z,constantTerm)
            % LPV Petrov Galerkin method with projection matrices W (for
            % the residual), V (for states) and Z (for scheduling); 
            % x_hat = Vxr; p_hat = Zpr
            
            if constantTerm == 1
                if size(Z,1) ~= obj.Np
                    error("Dimension of Z needs to be Np x Rp");
                end
                A0 = obj.A(:,:,1); At = obj.A(:,:,2:end);
                B0 = obj.B(:,:,1); Bt = obj.B(:,:,2:end);
                C0 = obj.C(:,:,1); Ct = obj.C(:,:,2:end);
                D0 = obj.D(:,:,1); Dt = obj.D(:,:,2:end);
                
                A0r = double(ttm(A0,{(W'*V)^-1*W', V'},[1,2]));
                B0r = double(ttm(B0,{(W'*V)^-1*W'},1));
                C0r = double(ttm(C0,{V'},2));
                D0r = double(D0);
                
                Ar = double(ttm(At,{(W'*V)^-1*W', V',Z'},[1,2,3]));
                Br = double(ttm(Bt,{(W'*V)^-1*W',Z'},[1,3]));
                Cr = double(ttm(Ct,{V',Z'},[2,3]));
                Dr = double(ttm(Dt,{Z'},3)); 
                
                Zpi = (Z'*Z)^-1*Z';
                etar = @(xr,u) Zpi*obj.eta_map(V*xr,u);

                Ared = cat(3,A0r,Ar);
                Bred = cat(3,B0r,Br);
                Cred = cat(3,C0r,Cr);
                Dred = cat(3,D0r,Dr);
                tensSSr = tensorSS(Ared,Bred,Cred,Dred,etar,obj.Ts);

            else
                if size(Z,1) ~= obj.Np+1
                    error("Dimension of Z needs to be (Np+1) x Rp");
                end
                Ar = double(ttm(obj.A,{(W'*V)^-1*W', V', Z'},[1,2,3]));
                Br = double(ttm(obj.B,{(W'*V)^-1*W',Z'},[1,3]));
                Cr = double(ttm(obj.C,{V',Z'},[2,3]));
                Dr = double(ttm(obj.D,{Z'},3));
                Zpi = (Z'*Z)^-1*Z';
                eta_red = @(xr,u) Zpi*[1; obj.eta_map(V*xr,u)];
                % eta_red = @(xr,u) Zpi*[obj.eta_map(V*xr,u)];

                tensSSr = tensorSS(Ar,Br,Cr,Dr,eta_red,obj.Ts);
            end
            
        end

        function V = reachabilityN(obj,N)
            % Construction of N-partial reachability matrix (im V =
            % \mathcal(R)_N
            R0 = tens2mat(double(obj.B),1,[2,3]);
            V = orth(R0);
            for k = 1:N
                V = orth([V, tens2mat(double(ttm((obj.A),V',2)),1,[2,3])]);
            end
        end
        function W = observabilityN(obj,N)
            % Construction of N-partial observability matrix (ker W =
            % \mathcal(O)_N
            O0 = tens2mat(permute(double(obj.C),[2,1,3]),1);
            W = orth(O0);
            At = permute(double(obj.A),[2,1,3]);
            for k = 1:N
                W = orth([W, tens2mat(double(ttm(tensor(At), W',2)),1,[2,3])]);
            end
        end

        function Qn = ObservabilityTensors(obj,N)
            % N-partial observability tensors stored in cell structure:
            % Q1 = Qn{1},..., Qn = Qn{n}
            A = obj.A;
            C = obj.C;
            Qn = cell(N,1);
            Qn{1} = permute(C,[1,3,2]);
            Ap = permute(A,[1,3,2]);
            for i = 2:N
                Qn{i} = ttt(Qn{i-1},Ap,i+1,1);
            end
        end
        function Wn = ReachabilityTensors(obj, n)
            % n-partial reachability tensor stored in cell structure
            % W1 = Wn{1},..., Wn = Wn{n}
            A = obj.A;
            B = obj.B;

            Wn = cell(n,1);
            Wn{1} = permute(B,[1,3,2]);
            Ap = permute(A,[1,3,2]);
            for i = 2:n
                Wn{i} = ttt(Ap,Wn{i-1},3,1);
            end
        
        end

        function [tensSSr, Tx, Tp, info] = lpvTensMM(obj, N, decomp, mode)
            % Tensor-based Moment Matching: Joint Reduction

            % N - horizon: (there are certain restrictions regarding the
            % horizon, as the tensor orders grow with N - it cant be very
            % large)
            % mode = 'T' - optimal; mode = 'R' - based on reachability;
            % mode = 'O' - based on observability;
            % decomp - choice of tensor decomposition
            % ('hosvd','tsvd','succr1')

            if nargin <= 3
                mode = 'T';
            end
            if nargin <= 2
                mode = 'T';
                decomp = 'hosvd';
            end

            assert(ismember(mode, {'T', 'O', 'R'}), 'Unsupported mode %s', mode);

            if strcmp(mode,'R') || strcmp(mode,'T')
                Wn = ReachabilityTensors(obj,N+1);

                Tstates = [];
                Tsched = [];
    
                for n = 1:N+1
                    R = modrank(Wn{n});
                    K = min(R);
    
                    if lower(decomp) == "tsvd"
                        % for every time horizon perform TSVD
                        tic
                        [x,gains,iterations,Xa] = TSVDND_SW(Wn{n},K);
                        tsvd_info{n}.CompTime = toc;
                        % separate singular vectors: 
                    
                        % for the state projection V, we gather the
                        % vectors from first dimension x{1}
                        Tstates = [Tstates, x{1}];
                    
                        % for the scheduling projection Z, we gather vectors from dimensions
                        % x{2},..., x{n}
                        for j = 2:n+1
                            Tsched = [Tsched, x{j}];
                            % Tsched{j-1} = null(x{j}.')
                        end
                
                        % save info from decomposition from every level
                        tsvd_info{n}.x = x;
                        tsvd_info{n}.gains = gains;
                        tsvd_info{n}.iterations = iterations;
                
                    elseif lower(decomp) == 'hosvd'
                        tic
                        T = hosvd(tensor(Wn{n}),2*sqrt(eps));
                        hosvd_info{n}.CompTime = toc;
                        hosvd_info{n}.T = T;

                        Tstates = [Tstates, T{1}];
                        for j = 2:n+1
                            Tsched = [Tsched, T{j}];
                        end
                    end
                end
                if strcmp(lower(decomp),'tsvd')
                    info = tsvd_info;
                elseif strcmp(lower(decomp),'hosvd')
                    info = hosvd_info;
                end
            end
            if strcmp(mode,'O') || strcmp(mode,'T')
                Qn = ObservabilityTensors(obj,N+1);
                Tstatesq = [];
                Tschedq = [];

                for n = 1:N+1
                    R = modrank(Qn{n});
                    K = min(R);
                    if lower(decomp) == "tsvd"
                        % for every time horizon perform TSVD
                        [x,gains,iterations,Xa] = TSVDND_SW(Qn{n},K);
                    
                        % separate singular vectors: 
                    
                        % for the state projection V, we gather the
                        % vectors from first dimension x{1}
                        Tstatesq = [Tstatesq, x{n+2}];
                    
                        % for the scheduling projection Z, we gather vectors from dimensions
                        % x{2},..., x{n}
                        for j = 2:n+1
                            Tschedq = [Tschedq, x{j}];
                            % Tsched{j-1} = null(x{j}.')
                        end
                
                        % save info from decomposition from every level
                        tsvd_info{n}.x = x;
                        tsvd_info{n}.gains = gains;
                        tsvd_info{n}.iterations = iterations;
                        
                    elseif lower(decomp) == 'hosvd'
                        tic
                        Tq = hosvd(tensor(Qn{n}),2*sqrt(eps));
                        hosvd_info{n}.CompTime = toc;
                        hosvd_info{n}.T = Tq;
                        Tstatesq = [Tstatesq, Tq{n+2}];
                        for j = 2:n+1
                            Tschedq = [Tschedq, Tq{j}];
                        end
                    end
                end
                if strcmp(lower(decomp),'tsvd')
                    info = tsvd_info;
                elseif strcmp(lower(decomp),'hosvd')
                    info = hosvd_info;
                end
            end

            if strcmp(mode,'T')
                Rn = orth(Tstates);
                Nn = orth(Tstatesq);
                
                if rank(Nn'*Rn,eps^2) == size(Rn,2)
                    [N,R] = qr(Nn'*Rn);

                    V = Rn*R;
                    W = Nn/N;
                    Z = orth(Tsched);
                else
                    mode = 'R';
                    fprintf("Rank condition not satisfied: switching to mode = 'R' \n");
                end
            end
            if strcmp(mode,'R')
                V = orth(Tstates);
                W = V;
                Z = orth(Tsched);
            elseif strcmp(mode,'O')
                V = orth(Tstatesq);
                W = V;
                Z = orth(Tschedq);
            end

            Tx = V;
            Tp = Z;
            tensSSr = PetrovGalerkinLPV(obj,W,V,Z,0);
            

        end

        function [tensSSr, Tx] = lpvTensMMStates(obj, N, mode)
            % Moment Matching: State Reduction (like lpvmmred)

            if nargin <= 2
                mode = 'T';
            end

            assert(ismember(mode, {'T', 'O', 'R'}), 'Unsupported mode %s', mode);
            V = reachabilityN(obj, N);
            W = observabilityN(obj,N).';
            Vinv = pinv(V);
            Winv = pinv(W);
            if (rank(V) == rank(W)) && (rank(V) == rank(W * V)) && strcmp(mode, 'T')
                RH = V / (W * V);
                tempA = ttm(obj.A,RH',2);
                Ar = ttm(tempA,W,1);
                Br = ttm(obj.B,W,1);
                Cr = ttm(obj.C,RH',2);
                Dr = obj.D;

                eta_map_r = @(x,u) obj.eta_map(RH*x,u);
                tensSSr = tensorSS(Ar,Br,Cr,Dr,eta_map_r,obj.Ts);
                Tx = RH;
            elseif strcmp(mode, 'T')
                warning('Rank condition not satisfied for ''mode=T'', switching to ''mode=R''');
                mode = 'R';
            end

            if strcmp(mode, 'R')
                tempA = ttm(obj.A,V',2);
                Ar = ttm(tempA,Vinv,1);
                Br = ttm(obj.B,Vinv,1);
                Cr = ttm(obj.C,V',2);
                Dr = obj.D;

                eta_map_r = @(x,u) obj.eta_map(V*x,u);
                tensSSr = tensorSS(Ar,Br,Cr,Dr,eta_map_r,obj.Ts);
                Tx = V;
            elseif strcmp(mode, 'O')
                tempA = tmprod(obj.A,Winv,2);
                Ar = ttm(tempA,W',1);
                Br = ttm(obj.B,W',1);
                Cr = ttm(obj.C,Winv,2);
                Dr = obj.D;

                eta_map_r = @(x,u) obj.eta_map(W*x,u);
                tensSSr = tensorSS(Ar,Br,Cr,Dr,eta_map_r,obj.Ts);
                Tx = W;
            end            
        end
           
    end
end

