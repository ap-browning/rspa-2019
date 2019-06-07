%% EXAMPLE USAGE
%  Example usage for SpringsDiscrete( .. ) and SpringsContinuous( .. )
%  which solve the mathematical models for the metamaterial presented in
%  Browning et al. (2019).

%  Each function contains a descriptive help text.
%  Run each section below to see:
%       1) A solution to the discrete model, corresponding to 
%          Figure 1e and Figure 1f of the main text
%       2) A solution to the discrete model, corresponding to 
%          Figure 3b and Figure 3f of the main text. Note that the
%          numerical parameters must be refined as per table S2 of the
%          supporting material document to obtain accurate figures and an
%          accurate estimate of the wavespeed.

%% DISCRETE MODEL

% Physical parameters
pparams_d.m         = 1 / 1000;
pparams_d.k         = 1 / 1000;
pparams_d.gamma     = 1 / 1000;
pparams_d.Delta     = 0.002;
pparams_d.eta       = 2;
pparams_d.v         = 1 / 1000;
pparams_d.N         = 101;
pparams_d.eps       = 0.01;

% Numerical parameters
nparams_d.dt        = 0.01;
nparams_d.tend      = 900;

% Time of retransmission
signalback          = 400;

% Solve
[U,~,A]             = SpringsDiscrete(pparams_d,nparams_d,signalback);

%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X               = 1:pparams_d.N;
T               = 0:nparams_d.dt:nparams_d.tend;

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesFontSize',16);

% Plot U
subplot(1,2,1);
imagesc(X,T,U'); set(gca,'YDir','normal'); pbaspect([1 1 1]);
cbar = colorbar; caxis([-1,1]); cbar.Label.Interpreter = 'latex';
xlabel('$i$'); ylabel('$t$'); title('$u_i(t)$');

% Plot A
subplot(1,2,2);
imagesc(X,T,A'); set(gca,'YDir','normal'); pbaspect([1 1 1]);
cbar = colorbar; caxis([-0.5,0.5]); cbar.Label.Interpreter = 'latex';
xlabel('$i$'); ylabel('$t$'); title('$a_i(t)$');


%% CONTINUOUS MODEL

% Physical parameters
pparams_c.nu    = 4;
pparams_c.eta   = 2;
pparams_c.kappa = 0.01;

% Numerical parameters
nparams_c.tend  = 1000;
nparams_c.xend  = 1000;
nparams_c.dt    = 0.01;
nparams_c.nx    = 5000;

% Solve
[T,F,U,~,A]     = SpringsContinuous(pparams_c,nparams_c);

% Estimate wavespeed
cest            = EstimateWavespeed(T,F);

%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X               = linspace(0,nparams_c.xend,nparams_c.nx+1);
Xplot           = X(1:5:end);
Tplot           = T(1:20:end);
Uplot           = U(1:5:end,1:100:end);
Aplot           = A(1:5:end,1:100:end);

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesFontSize',16);

% Plot U
subplot(1,3,1);
imagesc(Tplot,Xplot,Uplot); set(gca,'YDir','normal'); pbaspect([1 1 1]);
cbar = colorbar; caxis([-1,1]); cbar.Label.Interpreter = 'latex';
xlabel('$\hat{t}$'); ylabel('$\hat{x}$'); title('$\hat{u}(\hat{x},\hat{t})$');

% Plot A
subplot(1,3,2);
imagesc(Tplot,Xplot,Aplot); set(gca,'YDir','normal'); pbaspect([1 1 1]);
cbar = colorbar; caxis([-0.5,0.5]); cbar.Label.Interpreter = 'latex';
xlabel('$\hat{t}$'); ylabel('$\hat{x}$'); title('$\hat{a}(\hat{x},\hat{t})$');

% Plot Front Location
subplot(1,3,3);
plot(T,F,'LineWidth',3); ylim([0 1000]); pbaspect([1 1 1]);
xlabel('$\hat{t}$'); ylabel('$x \: : \: \hat{u} = 0$'); title(['Front Location, $c \approx\: $',num2str(cest)]);


