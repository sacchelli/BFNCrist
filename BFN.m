clear all;
close all;

addpath('./Functions/');

x_min = 1;
x_max = 2;
t_max = x_max - x_min;
Nx = 10^2;
Nt = Nx;
Nl = Nx;
nb = 10^2;
e1 = 2; % Eccentricity
e2 = 1; % Eccentricity

nb_iter = 100; % Number of iterations

mu = 2; % Speed ratio

affiche = false; % To plot the movie during computations
film = false; % To save the movie in folder "Videos"

cpt_img = 1;

Err = zeros(1, nb_iter+1);

T = linspace(0, t_max, Nt+1)';
dt = T(2) - T(1);


X = linspace(x_min, x_max, Nx+1)';
dx = X(2) - X(1);
X1 = linspace(x_min-t_max, x_max, 2*Nx+1)';
dx1 = X1(2) - X1(1);
X2 = linspace(x_min-2*t_max+dx, x_max-dx, 3/2*Nx)';
dx2 = X2(2) - X2(1);

mu = mu*dx;

N1 = zeros(2*Nx+1, 2*nb_iter*Nt);
N2 = zeros(3/2*Nx, 2*nb_iter*Nt);

N1(:, 1) = exp(-30*(X1-0.5).^2);
N2(:, 1) = exp(-30*(X2+0.5).^2);


L1 = linspace(0, 2*max(e1, 1)*x_max, Nx+1)';
L2 = linspace(0, 2*max(e2, 1)*x_max, Nx/2+1)';

Densite1 = corde_rep(X, L1, nb, e1);
Densite2 = corde_rep(X(1:2:end), L2, nb, e2);

N1hat = zeros(2*Nx+1, 2*nb_iter*Nt);
N2hat = zeros(3/2*Nx, 2*nb_iter*Nt);

Err(1) = sqrt(norm(N1(:, 1) - N1hat(:, 1))^2*dx1 + norm(N2(:, 1) - N2hat(:, 1))^2*dx2);

i = 1;

for iter = 1:nb_iter

   

if iter == 1
    imin = 2;
else
    imin = 1;
end

for k = imin:(Nt+1)
    
    i = i+1;

    N1(2:end, i) = N1(2:end, i-1) - dt/dx1 * diff(N1(:, i-1));
    N1(1, i) = N1(1, i-1) - dt/dx1 * (N1(1, i-1)-N1(end, i-1));
    
    N2(2:end, i) = N2(2:end, i-1) - 2*dt/dx2 * diff(N2(:, i-1));
    N2(1, i) = N2(1, i-1) - 2*dt/dx2 * (N2(1, i-1)-N2(end, i-1));    
    
    
    Q_cumul1 = Densite1'*N1(end-Nx:end, i-1);
    Q_cumul2 = Densite2'*N2(end-Nx/2:end, i-1);
    
    Q1_sum = Q_cumul1 + interp1(L2, Q_cumul2, L1, 'linear', 0);
    Q2_sum = Q_cumul2 + interp1(L1, Q_cumul1, L2, 'linear', 0);
    
    Q_cumul1hat = Densite1'*N1hat(end-Nx:end, i-1);
    Q_cumul2hat = Densite2'*N2hat(end-Nx/2:end, i-1);
    
    Q1hat_sum = Q_cumul1hat + interp1(L2, Q_cumul2hat, L1, 'linear', 0);
    Q2hat_sum = Q_cumul2hat + interp1(L1, Q_cumul1hat, L2, 'linear', 0);


    Ceps1 = Q1hat_sum - Q1_sum;
    Ceps2 = Q2hat_sum - Q2_sum;
    
    N1hat(2:end, i) = N1hat(2:end, i-1) - dt/dx1 * diff(N1hat(:, i-1));
    N1hat(1, i) = N1hat(1, i-1) - dt/dx1 * (N1hat(1, i-1)-N1hat(end, i-1));
    
    N1hat(end-Nx:end, i) = N1hat(end-Nx:end, i) - mu*dt * (Densite1*Ceps1);
    
    N2hat(2:end, i) = N2hat(2:end, i-1) - 2*dt/dx2 * diff(N2hat(:, i-1));
    N2hat(1, i) = N2hat(1, i-1) - 2*dt/dx2 * (N2hat(1, i-1)-N2hat(end, i-1));
    
    N2hat(end-Nx/2:end, i) = N2hat(end-Nx/2:end, i) - mu*dt * (Densite2*Ceps2);
    
    
    if affiche
        
    figure(1)
    plot(X1, N1(:, i))
    axis([x_min-2*t_max, x_max, 0, 1])
    hold on
    plot(X1, N1hat(:, i))
    plot(X2, N2(:, i))
    plot(X2, N2hat(:, i))
    plot([x_min, x_min],[0, 1], '--k');
    legend({'State $\eta = 2$', 'Observer $\eta = 2$', 'State $\eta = 1$', 'Observer $\eta = 1$'}, 'interpreter', 'latex', 'location', 'northwest')
    xlabel('r')
    
    hold off
    drawnow
    
    if film
        print(['Images/Fig_' num2str(cpt_img)], '-dpng')
        cpt_img = cpt_img+1;
    end
    
    end

    
end

if iter == 20
    i20 = i;
end


for k = 1:(Nt+1)
    
    i = i+1;

    N1(1:end-1, i) = N1(1:end-1, i-1) + dt/dx1 * diff(N1(:, i-1));
    N1(end, i) = N1(end, i-1) + dt/dx1 * (N1(1, i-1)-N1(end, i-1));
    
    N2(1:end-1, i) = N2(1:end-1, i-1) + 2*dt/dx2 * diff(N2(:, i-1));
    N2(end, i) = N2(end, i-1) + 2*dt/dx2 * (N2(1, i-1)-N2(end, i-1));    
    
    Q_cumul1 = Densite1'*N1(end-Nx:end, i-1);
    Q_cumul2 = Densite2'*N2(end-Nx/2:end, i-1);
    
    Q1_sum = Q_cumul1 + interp1(L2, Q_cumul2, L1, 'linear', 0);
    Q2_sum = Q_cumul2 + interp1(L1, Q_cumul1, L2, 'linear', 0);
    

    Q_cumul1hat = Densite1'*N1hat(end-Nx:end, i-1);
    Q_cumul2hat = Densite2'*N2hat(end-Nx/2:end, i-1);
    
    Q1hat_sum = Q_cumul1hat + interp1(L2, Q_cumul2hat, L1, 'linear', 0);
    Q2hat_sum = Q_cumul2hat + interp1(L1, Q_cumul1hat, L2, 'linear', 0);


    Ceps1 = Q1hat_sum - Q1_sum;
    Ceps2 = Q2hat_sum - Q2_sum;
    
    N1hat(1:end-1, i) = N1hat(1:end-1, i-1) + dt/dx1 * diff(N1hat(:, i-1));
    N1hat(end, i) = N1hat(end, i-1) + dt/dx1 * (N1hat(1, i-1)-N1hat(end, i-1));
    
    N1hat(end-Nx:end, i) = N1hat(end-Nx:end, i) - mu*dt * (Densite1*Ceps1);
    
    N2hat(1:end-1, i) = N2hat(1:end-1, i-1) + 2*dt/dx2 * diff(N2hat(:, i-1));
    N2hat(end, i) = N2hat(end, i-1) + 2*dt/dx2 * (N2hat(1, i-1)-N2hat(end, i-1));
    
    N2hat(end-Nx/2:end, i) = N2hat(end-Nx/2:end, i) - mu*dt * (Densite2*Ceps2);
    
    
    if affiche
        
    figure(1)
    plot(X1, N1(:, i))
    axis([x_min-2*t_max, x_max, 0, 1])
    hold on
    plot(X1, N1hat(:, i))
    plot(X2, N2(:, i))
    plot(X2, N2hat(:, i))
    plot([x_min, x_min],[0, 1], '--k');
    legend({'State $\eta = 2$', 'Observer $\eta = 2$', 'State $\eta = 1$', 'Observer $\eta = 1$'}, 'interpreter', 'latex', 'location', 'northwest')
    xlabel('r')

    hold off
    drawnow
    
    if film
        print(['Images/Fig_' num2str(cpt_img)], '-dpng')
        cpt_img = cpt_img+1;
    end
    
    end

    
    

    
end

    Err(iter+1) = sqrt(norm(N1(:, i) - N1hat(:, i))^2*dx1 + norm(N2(:, i) - N2hat(:, i))^2*dx2);


end

% define figure properties
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 8;
opts.height     = 6;
opts.fontType   = 'Times';
opts.fontSize   = 9;

% scatter(1:2*Nt*nb_iter, Err);
fig2 = figure(2);
% scatter(0:nb_iter, Err, '*');
% fin = 2*nb_iter+1;
plot(0:nb_iter, Err, '-k');
start = 30;
P1 = polyfit(start:nb_iter,log(Err(start+1:nb_iter+1)),1)
% P2 = polyfit(log(start:nb_iter),log(Err(start+1:nb_iter+1)),1)
set(gca,'yscale','log');
hold on
plot(0:nb_iter, exp(P1(1)*(0:nb_iter) + P1(2)), '--k')
legend({'$\sqrt{\|\psi_1(t)-\hat\psi_1^{2n}(t)\|^2+\|\psi_2(t)-\hat\psi_2^{2n}(t)\|^2}$', ['$y = '  num2str(exp(P1(2))) ' \times ' num2str(exp(P1(1))) '^n$' ]}, 'Interpreter', 'latex', 'Location', 'northoutside');
xlabel('Number of iterations $2n$ of BFN', 'Interpreter', 'latex')
ylabel('$L^2$-error $\psi-\hat\psi^{2n}$', 'Interpreter', 'latex')
% plot(0:nb_iter, exp(P2(2)) * ((0:nb_iter)^P2(1)));
% set(gca,'xscale','log');
% legend('Err', ['y = '  num2str(exp(P2(2))) ' * n^{' num2str(P2(1)) '}'], 'c')

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


i = i - (Nt+1);


fig = figure(3);
figure(3)
subplot(1, 2, 2)
plot(X1/10, N1(:, i), '-b')
hold on
plot(X1/10, N1hat(:, i20), '-magenta')
hold on
plot(X1/10, N1hat(:, i), '-r')
title('$\eta_2 = 2$', 'Interpreter', 'latex')
axis([x_min/10, x_max/10, 0, 1])
legend({'$\psi_2(t, r)$', '$\hat\psi_2^{20}(t, r)$', '$\hat\psi_2^{100}(t, r)$'}, 'Interpreter', 'latex', 'Location', 'northoutside')
xlabel('$r$ (in mm)', 'Interpreter', 'latex')
subplot(1, 2, 1)
plot(X2/10, N2(:, i), '-b')
hold on
plot(X2/10, N2hat(:, i20), '-magenta')
hold on
plot(X2/10, N2hat(:, i), '-r')
axis([x_min/10, x_max/10, 0, 1])
xlabel('$r$ (in mm)', 'Interpreter', 'latex')
ylabel('PSD (in mm$^{-1}$.m$^{-3}$)', 'Interpreter', 'latex')
title('$\eta_1 = 1$', 'Interpreter', 'latex')
legend({'$\psi_1(t, r)$', '$\hat\psi_1^{20}(t, r)$', '$\hat\psi_1^{100}(t, r)$'}, 'Interpreter', 'latex', 'Location', 'northoutside')


% scaling
fig.Units               = 'centimeters';
fig.Position(3)         = opts.width;
fig.Position(4)         = 1.5*opts.height;

% set text properties
set(fig.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);

% remove unnecessary white space
figure(3);
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))

% scaling
fig4.Units               = 'centimeters';
fig4.Position(3)         = opts.width;
fig4.Position(4)         = 1.5*opts.height;

% set text properties
set(fig4.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);

% remove unnecessary white space
figure(4);
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))


% Save the movie :
if film
    cpt_img = cpt_img - 1;

    delta_t = 1/(cpt_img);
    images = cell(cpt_img,1);
    for i = 1:cpt_img
        images{i} = imread(char(['Images/Fig_' num2str(i) '.png']));
    end
    writerObj = VideoWriter(char('Videos/reconstruit.avi'));
    writerObj.FrameRate = 0.05/delta_t; % Dur?e de la vid?o

    open(writerObj);
    for i = 1:cpt_img
        frame = im2frame(images{i});
        writeVideo(writerObj, frame);
    end
    close(writerObj);
end













