function [Umat,Wmat,Amat] = SpringsDiscrete(pparams,nparams,varargin)
%SpringsDiscrete Solve the continuous metamaterial model described by
%Eqs 2.1 to Eq 2.3 in the main document.
%
%   Usage options: 
%       [U,W,A]     = SpringsDiscrete(pparams,nparams)
%       [U,W,A]     = SpringsDiscrete(pparams,nparams,signalback)
%
%   Inputs:
%       pparams     a structure array containing the physical params:
%           pparams.m
%           pparams.k
%           pparams.gamma
%           pparams.Delta
%           pparams.eta
%           pparams.v
%           pparams.N
%           pparams.eps
%       nparams     a structure array containing the numerical params:
%           nparams.tend    - End time
%           nparams.dt      - Timestep
%       signalback  (optional) time to send signal back from RHS 
%
%   Outputs:
%       U           displacement matrix
%       W           velocity matrix
%       A           biological mechanism matrix
%       T           time vector
%       F           location of the front (estimate)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PARAMETERS %%%%%%%%%%%%%%%%

%%%% Physical Parameters

    m       = pparams.m;
    k       = pparams.k;
    g       = pparams.gamma;
    D       = pparams.Delta;
    eta     = pparams.eta;
    v       = pparams.v;
    N       = pparams.N;
    eps     = pparams.eps;


%%%% Numerical Parameters

    t0      = 0;
    dt      = nparams.dt;
    tend    = nparams.tend;
    
%%%% signalback

    if nargin == 3
        signalback = varargin{1};
    else
        signalback = tend + 1;
    end
    sentback = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SOLVE %%%%%%%%%%%%%%%%%%%%%

%%%% Setup Problem

    % State variables store
    U = zeros(N,1);
    W = zeros(N,1);
    A = zeros(N,1);

    % Initial Condition
    U(1)        = +1;
    U(2:end)    = -1;
    A(:)        = -1/eta;

    % Save at each timestep
    Umat = zeros(N,(tend - t0)/dt + 1);
    Wmat = zeros(N,(tend - t0)/dt + 1);
    Amat = zeros(N,(tend - t0)/dt + 1);

    Umat(:,1) = U;
    Wmat(:,1) = W;
    Amat(:,1) = A;

%%%% Solve using forward euler
    for i = 1:(tend - t0)/dt

        diffU = diff(U);

        dU          = W;
        dW          = -1/m * (g * W + v * (U - 1).*(U + 1).*(U - A));
        dW(1)       = dW(1) + k/m * (U(2) - U(1));
        dW(2:end-1) = dW(2:end-1) + k/m * (diffU(2:end) - diffU(1:end-1));
        dW(end)     = dW(end) - k/m * (U(end) - U(end-1));
        dA          = eps * (U / eta - A);

        U = U + dU * dt;
        W = W + dW * dt;
        A = A + dA * dt;   

        Umat(:,i+1) = U;
        Wmat(:,i+1) = W;
        Amat(:,i+1) = A;
        
        % Send back?
        if i * dt > signalback && ~sentback
           
            U(end)      = -1;
            sentback    = true;
            
        end

    end

end