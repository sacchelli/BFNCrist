clear all;
close all;

addpath('./Functions/');

x_min = 0;
x_max = 2;
Nx = 2*10^2;
Nl = Nx;
nb = 10^2;
sol = true;
X = linspace(0, x_max, Nx)';
Delta = [5, 3, 1];
e = 2;

% define figure properties
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'Images/';
opts.width      = 8;
opts.height     = 6;
opts.fontType   = 'Times';
opts.fontSize   = 9;

% create new figure
fig = figure(1);
fig2 = figure(2);
fig3 = figure(3);

for i = 1:3
    
    delta = 10^(-Delta(i));
    
    
    L = linspace(0, 2*max(e, 1)*x_max, Nl)';
    Densite = corde_rep(X, L, nb, e);
    C = d_mat(X) ./ d_mat(L) .* Densite';
    N = exp(-30*(X-1/2).^2)+exp(-30*(X-3/2).^2);
    N = N / (sum(N)*(X(2)-X(1)));
    
    Q_cumul = C * N;
    
    Q = diff(Q_cumul);
    Q = [Q_cumul(1); Q];
    
    % Generation du bruit
    sigma = 0.02; % variance du bruit
    moy = 0;   % moyenne
    bruit = moy + sigma*randn(Nx,1)*max(Q);
    Qb = Q + bruit;
    
    Q_cumul = cumsum(Qb);

    Ninv = pb_inv(Q_cumul, X, L, Densite, delta);
    
    Q = Q * (L(2)-L(1));
    Qb = Qb * (L(2)-L(1));
    Q_cumul = Q_cumul * (L(2)-L(1));
    

    figure(i)
    plot((1+X)/10, 10*N, 'b-')
    hold on
    plot((1+X)/10, 10*Ninv, 'r-')
    xlabel('$r$ (in mm)', 'interpreter', 'latex')
    ylabel('$\bar\psi(r)$ (in mm$^{-1}$)', 'interpreter', 'latex')
    title(['$\delta = 10^{-' num2str(Delta(i)) '}$'], 'interpreter', 'latex')
end

% scaling
fig.Units               = 'centimeters';
fig.Position(3)         = opts.width;
fig.Position(4)         = opts.height;

% set text properties
set(fig.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);

% remove unnecessary white space
figure(1);
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
% 
% scaling
fig2.Units               = 'centimeters';
fig2.Position(3)         = opts.width;
fig2.Position(4)         = opts.height;

% set text properties
set(fig2.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);

% remove unnecessary white space
figure(2);
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))

% 
% scaling
fig3.Units               = 'centimeters';
fig3.Position(3)         = opts.width;
fig3.Position(4)         = opts.height;

% set text properties
set(fig3.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);

% remove unnecessary white space
figure(3);
% legend({'Actual PSD', 'Estimated PSD'}, 'Interpreter', 'latex', 'Location', 'southoutside');
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
