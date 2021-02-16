clear all;
close all;

addpath('./Functions/');

x_min = 0;
x_max = 1;
Nx = 2*10^2;
Nl = Nx;
nb = 10^2;
X = linspace(0, x_max, Nx)';

E = [1/2, 1, 2]; % Eccentricities to plot


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
    

    e = E(i);
    L = linspace(0, 2*max(e, 1)*x_max, Nl)';
    Densite = corde_rep(X, L, nb, e);
    C = d_mat(X) ./ d_mat(L) .* Densite';
    N = zeros(Nx, 1);
    N(end) = 1/(X(2)-X(1));
    N = N / (sum(N)*(X(2)-X(1)));
    
    Q_cumul = C * N;
    
    Q = diff(Q_cumul);
    Q = [Q_cumul(1); Q];
    

    Q = Q * (L(2)-L(1));
    Q_cumul = Q_cumul * (L(2)-L(1));
    
    figure(i)
    plot(L/10, Q_cumul, 'k-')
    xlabel('$\ell$ (in mm)', 'interpreter', 'latex')
    ylabel('$\bar Q(\ell)$ (dimensionless)', 'interpreter', 'latex')
    title(['$\eta =' num2str(e) '$'], 'interpreter', 'latex')

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
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))