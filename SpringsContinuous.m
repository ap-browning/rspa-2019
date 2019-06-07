function varargout = SpringsContinuous(pparams,nparams,varargin)
%SpringsContinuous  Solve the continuous metamaterial model described by
%Eq 2.6 and Eq 2.7 in the main document:
%   
%       0 = u_tt - u_xx + u_t + nu * (u^2 - 1)*(u/eta - a)      (2.6)
%       0 = a_t  - kappa * (u/eta - a)                          (2.7)
%   
%       u(x,0) = +1, x <  1
%       u(x,0) = -1, x >= 1
%       a(x,0) = -1/eta
%
%   Usage options: 
%       [U,W,A]     = SpringsContinuous(pparams,nparams)
%       [T,F]       = SpringsContinuous(pparams,nparams)
%       [T,F,U,W,A] = SpringsContinuous(pparams,nparams)
%       [U,W,A]     = SpringsContinuous(pparams,nparams,picard)
%
%   Inputs:
%       pparams     a structure array containing the physical params:
%           pparams.nu
%           pparams.eta
%           pparams.kappa
%       nparams     a structure array containing the numerical params:
%           nparams.tend    - End time
%           nparams.xend    - Domain size: 0 < x < xend
%           nparams.dt      - Timestep
%           nparams.nx      - Number points in the spatial discretisation
%       picard      (optional) a structure array containing picard params:
%           picard.iters    - maximum iterations in picard iteration
%           picard.tol      - tolerance in picard iteration
%
%   Outputs:
%       U           displacement matrix
%       W           velocity matrix
%       A           biological mechanism matrix
%       T           time vector
%       F           location of the front (estimate)
%
% Caution:  Outputs [T,F] may behave unexpectedly if code is modified to
%           account for multiple signals.

switch nargout
    case 1
        error("Incorrect number of outputs requested!");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PARAMETERS %%%%%%%%%%%%%%%%

%%%% Physical Parameters

    nu      = pparams.nu;
    eta     = pparams.eta;
    kappa   = pparams.kappa;

%%%% Numerical Parameters

    % Displacement at front
    frontloc = 0;

    % Time domain
    t0      = 0;
    tend    = nparams.tend;
    dt      = nparams.dt;

    % Spacial domain
    x0      = 0;
    xend    = nparams.xend;
    M       = nparams.nx;
    dx      = (xend - x0) / M;

    % Picard Iteration Parameters
    if nargin == 3
        maxiters = varargin{1}.iters;   % Maximum iterations catch
        tol      = varargin{1}.tol;     % Solution tolerance
    else
        maxiters = 100;                 % Maximum iterations catch
        tol     = 1e-6;                 % Solution tolerance
    end
    
    % Initial condition
    Uic = @(x) 1 * (x <= 1) + -1 * (x > 1); 
    Wic = @(x) 0 * x;
    Aic = @(x) -1/eta * ones(size(x));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SOLVE %%%%%%%%%%%%%%%%%%%%%

%%%% Setup Problem

    % Spatial and time vectors
    X = x0:dx:xend;
    T = t0:dt:tend;

    % State variables store
    U = Uic(X);
    W = Wic(X);
    A = Aic(X);

    % Save at each timestep (if necessary)
    if nargout == 3 || nargout == 5

        % Save at each timestep
        Umat = zeros(M+1,(tend - t0)/dt + 1);
        Umat(:,1) = U;

        Wmat = zeros(M+1,(tend - t0)/dt + 1);
        Wmat(:,1) = W;

        Amat = zeros(M+1,(tend - t0)/dt + 1);
        Amat(:,1) = A;
    end

    % Save location of front
    Front = zeros(size(T));
    Front(1) = X(find(U >= frontloc,1,'last'));

%%%% Integrate through time
    for i = 1:(tend - t0)/dt

        % Solution at last time point
        Ulast = U;
        Wlast = W;
        Alast = A;

        diffUlast = diff(Ulast);

        % Solve nonlinear equation Y = G(Y) to get updated state
        b1 = dt / (2 * dx^2);
        b2 = dt / (2);
        b3 = nu * dt / (2);
        b4 = kappa * dt / 2;   

        % Constant contribution from Ylast
        dUlastadd           = Ulast + dt / 2 * Wlast;
        dWlastadd           = Wlast - b2 * Wlast - b3 * Vd(Ulast,Alast,1);
        dWlastadd(1)        = dWlastadd(1) + b1 * diffUlast(1);
        dWlastadd(2:end-1)  = dWlastadd(2:end-1) + b1 * (diffUlast(2:end) - diffUlast(1:end-1));
        dWlastadd(end)      = dWlastadd(end) - b1 * diffUlast(end);
        dAlastadd           = Alast + b4 * (Ulast / eta - Alast);

        % Picard iteration
        converged = false;
        iters = 0;
        while ~converged

            % Current iterations
            Uold = U;
            Wold = W;
            Aold = A;

            diffUold = diff(Uold);

            % Next iteration
            U           = dUlastadd + dt / 2 * Wold;
            W           = dWlastadd - b2 * Wold - b3 * Vd(Uold,Aold,1);
            W(1)        = W(1)      + b1 * diffUold(1);
            W(2:end-1)  = W(2:end-1) + b1 * (diffUold(2:end) - diffUold(1:end-1));
            W(end)      = W(end)    - b1 * diffUold(end);
            A           = dAlastadd + b4 * (Uold / eta - Aold);

            % Check for convergence
            if norm(U - Uold) + norm(W - Wold) + norm(A - Aold) < tol
               converged = true; 
            end

            iters = iters + 1;
            if iters == maxiters
               error('Picard iteration did not converge!'); 
            end

        end

        % If appropriate, store [U,W,A]
        if nargout == 3 || nargout == 5
            Umat(:,i+1) = U;
            Wmat(:,i+1) = W;
            Amat(:,i+1) = A;
        end

        % If appropriate, store [T,F]
        frontlocindex = find(U > frontloc,1,'last');
        if ~isempty(frontlocindex)
            if frontlocindex == M + 1
                Front(i+1) = xend;
            else
                Ua = U(frontlocindex);      Xa = X(frontlocindex);
                Ub = U(frontlocindex + 1);  Xb = X(frontlocindex + 1);
                Front(i+1) = Xa - (Xb - Xa) * Ua / (Ub - Ua);
            end
        else
           Front(i+1) = 0; 
        end

    end % End picard through time
    
    % Handle outputs
    switch nargout
        case 2
            varargout{1} = T;
            varargout{2} = Front;
        case 3
            varargout{1} = Umat;
            varargout{2} = Wmat;
            varargout{3} = Amat;
        case 5
            varargout{1} = T;
            varargout{2} = Front;
            varargout{3} = Umat;
            varargout{4} = Wmat;
            varargout{5} = Amat;
    end

end

function out = Vd(U,A,d)

    out = (U - d) .* (U + d) .* (U - A);

end