clear all;

addpath('./Functions/');

x_min = 1;
x_max = 2;
t_max = x_max - x_min;
Nx = 10^2;
Nt = Nx;
Nx1 = Nx;
Nx2 = Nx;
Nl = Nx;
nb = 10^2;

% -------------------
% Free parameters: 

e1 = 0.5; % Eccentricity
e2 = 2; % Eccentricity

nb_iter = 1000; % Number of iterations


affiche = false; % To plot the short movie during computations
affichetout = false; % To plot the long movie during computations
film = false; % To save the long movie in folder "Videos"

% ----------------


v1 = 2; % Positive growth rate
v2 = -1; % Negative growth rate

mu = 0.001; % observer gain

cpt_img = 1;

Err = zeros(1, nb_iter+1);

T = linspace(0, t_max, Nt+1)';
dt = T(2) - T(1);

X1 = [linspace(x_min-t_max/x_min^2, x_min-dt, Nx1/2)'; linspace(x_min, x_max, Nx1/2)'];
X2 = [linspace(x_min, x_max, Nx2/2)'; linspace(x_max+2*dt, x_max+2*t_max/x_max^2, Nx2/2)'];

X = linspace(x_min, x_max, Nx)';
dx = X(2) - X(1);

X1t = zeros(Nx1, Nt+1);
X2t = zeros(Nx2, Nt+1);

X1t(:, 1) = X1;
X2t(:, 1) = X2;

N1 = exp(-60*(X1-0.5).^2);
N2 = exp(-60*(X2-1.5).^2);
% N1 = exp(-90*(X1-0.5).^2);
% N2 = exp(-120*(X2-2.25).^2);

N1t = zeros(Nx1, 2*nb_iter*Nt);
N2t = zeros(Nx2, 2*nb_iter*Nt);

N1t(:, 1) = interp1(X1, N1, X);
N2t(:, 1) = interp1(X2, N2, X);

L = linspace(0, 2*max([e1, e2, 1])*x_max, Nx)';

Densite1 = corde_rep(X, L, nb, e1);
Densite2 = corde_rep(X, L, nb, e2);

N1hat = zeros(Nx1, 2*nb_iter*Nt);
N2hat = zeros(Nx2, 2*nb_iter*Nt);

Vit1 = zeros(Nx1, Nt);
Vit2 = zeros(Nx2, Nt);

Err(1) = sqrt(norm(N1)^2*dx + norm(N2)^2*dx);

i = 1;

for iter = 1:nb_iter


for k = 1:Nt
    
    i = i+1;
        
    if iter == 1
    
        Vit1(:, i-1) = (X1t(:, i-1)>x_min).*(X1t(:, i-1)<x_max)./X1t(:, i-1).^2;
        Vit2(:, i-1) = (X2t(:, i-1)>x_min).*(X2t(:, i-1)<x_max)./X2t(:, i-1).^2;
        
        Vit1(isnan(Vit1(:, i-1)), i-1) = 0;
        Vit2(isnan(Vit2(:, i-1)), i-1) = 0;

        Vit1(:, i-1) = (X1t(:, i-1)<=x_min)/x_min^2 + (X1t(:, i-1)>=x_max)/x_max^2  + Vit1(:, i-1);
        Vit2(:, i-1) = (X2t(:, i-1)<=x_min)/x_min^2 + (X2t(:, i-1)>=x_max)/x_max^2  + Vit2(:, i-1);
        
    end
    
    
    X1t(:, i) = X1t(:, i-1) + v1*dt*Vit1(:, k);
    X2t(:, i) = X2t(:, i-1) + v2*dt*Vit2(:, k);
    
    N1t(:, i) = interp1(X1t(:, i), N1, X);
    N2t(:, i) = interp1(X2t(:, i), N2, X);
    
    N1t(isnan(N1t(:, i)), i) = 0;
    N2t(isnan(N2t(:, i)), i) = 0;
    
    Q_cumul1 = Densite1'*N1t(:, i-1);
    Q_cumul2 = Densite2'*N2t(:, i-1);
    
    Q_sum = Q_cumul1 + Q_cumul2;
    
    N1hattemp = interp1(X1t(:, i), N1hat(:, i-1), X);
    N2hattemp = interp1(X2t(:, i), N2hat(:, i-1), X);
    
    N1hattemp(isnan(N1hattemp)) = 0;
    N2hattemp(isnan(N2hattemp)) = 0;
    
    Q_cumul1hat = Densite1'*N1hattemp;
    Q_cumul2hat = Densite2'*N2hattemp;
    
    Qhat_sum = Q_cumul1hat + Q_cumul2hat;

    Ceps = Qhat_sum - Q_sum;
    
    CC1 = mu*dt*interp1(X, (Densite1*Ceps), X1t(:, i-1));
    CC2 = mu*dt*interp1(X, (Densite2*Ceps), X2t(:, i-1));
    
    CC1(isnan(CC1)) = 0;
    CC2(isnan(CC2)) = 0;
    
    N1hat(:, i) = N1hat(:, i-1) - CC1;
    N2hat(:, i) = N2hat(:, i-1) - CC2;
    
    if affichetout || (affiche && mod(i, 2*Nt)==0)
        
    figure(1)
    plot(X1t(:, i), N1)
    axis([x_min-t_max/x_min^2, x_max+2*t_max/x_max^2, 0, 1])
    hold on
    plot(X1t(:, i), N1hat(:, i))
    plot(X2t(:, i), N2)
    plot(X2t(:, i), N2hat(:, i))
    plot([x_min, x_min],[0, 1], '--k');
    plot([x_max, x_max],[0, 1], '--k');
    legend({'State $\eta = 2$', 'Observer $\eta = 2$', 'State $\eta = 1$', 'Observer $\eta = 1$'}, 'interpreter', 'latex', 'location', 'northwest')
    xlabel('r')
    hold off
    drawnow
%     pause(3)
    
%     figure(2)
%     plot(L, Q_cumul1)
%     hold on
%     plot(L, Q_cumul2)
% %     plot(L, Q_cumul1hat) 
% %     plot(L, Q_cumul2hat)
%     hold off
%     legend('1', '2')%, 'hat1', 'hat2')
%     drawnow
% %     pause(1)
    
    if film
        print(['Images/Fig_' num2str(cpt_img)], '-dpng')
        cpt_img = cpt_img+1;
    end
    
    end

    
end

% if iter == 20
%     i20 = i;
% end


for k = Nt:(-1):1
    
    i = i+1;
    
    X1t(:, i) = X1t(:, i-1) - v1*dt*Vit1(:, k);
    X2t(:, i) = X2t(:, i-1) - v2*dt*Vit2(:, k);
    
    N1t(:, i) = interp1(X1t(:, i), N1, X);
    N2t(:, i) = interp1(X2t(:, i), N2, X);
    
    N1t(isnan(N1t(:, i)), i) = 0;
    N2t(isnan(N2t(:, i)), i) = 0;
    
    Q_cumul1 = Densite1'*N1t(:, i-1);
    Q_cumul2 = Densite2'*N2t(:, i-1);
    
    Q_sum = Q_cumul1 + Q_cumul2;
    
    N1hattemp = interp1(X1t(:, i), N1hat(:, i-1), X);
    N2hattemp = interp1(X2t(:, i), N2hat(:, i-1), X);
    
    N1hattemp(isnan(N1hattemp)) = 0;
    N2hattemp(isnan(N2hattemp)) = 0;
    
    Q_cumul1hat = Densite1'*N1hattemp;
    Q_cumul2hat = Densite2'*N2hattemp;
    
    Qhat_sum = Q_cumul1hat + Q_cumul2hat;

    Ceps = Qhat_sum - Q_sum;
    
    CC1 = mu*dt*interp1(X, (Densite1*Ceps), X1t(:, i-1));
    CC2 = mu*dt*interp1(X, (Densite2*Ceps), X2t(:, i-1));
    
    CC1(isnan(CC1)) = 0;
    CC2(isnan(CC2)) = 0;
    
    N1hat(:, i) = N1hat(:, i-1) - CC1;
    N2hat(:, i) = N2hat(:, i-1) - CC2;
    
    
    if affichetout || (affiche && mod(i, 2*Nt)==0)
        
    figure(1)
    plot(X1t(:, i), N1)
    axis([x_min-t_max/x_min^2, x_max+2*t_max/x_max^2, 0, 1])
    hold on
    plot(X1t(:, i), N1hat(:, i))
    plot(X2t(:, i), N2)
    plot(X2t(:, i), N2hat(:, i))
    plot([x_min, x_min],[0, 1], '--k');
    plot([x_max, x_max],[0, 1], '--k');
    legend({'State $\eta = 2$', 'Observer $\eta = 2$', 'State $\eta = 1$', 'Observer $\eta = 1$'}, 'interpreter', 'latex', 'location', 'northwest')
    xlabel('r')
    hold off
    drawnow
    
    if film
        print(['Images/Fig_' num2str(cpt_img)], '-dpng')
        cpt_img = cpt_img+1;
    end
    
    end

    
    if iter == 10
        i20 = i;
    end

    
end

    Err(iter+1) = sqrt((N1 - N1hat(:, i))' * diag(d_mat(X1t(:, i))) * (N1 - N1hat(:, i)) + (N2 - N2hat(:, i))'* diag(d_mat(X2t(:, i))) *(N2 - N2hat(:, i)));


end

%% Figures

close all;

% define figure properties
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 8;
opts.height     = 6;
opts.fontType   = 'Times';
opts.fontSize   = 9;

% % scatter(1:2*Nt*nb_iter, Err);
% fig2 = figure(2);
% % scatter(0:nb_iter, Err, '*');
% % fin = 2*nb_iter+1;
% plot(0:nb_iter, Err, '-k');
% % start = 30;
% % P1 = polyfit(start:nb_iter,log(Err(start+1:nb_iter+1)),1)
% % P2 = polyfit(log(start:nb_iter),log(Err(start+1:nb_iter+1)),1)
% set(gca,'yscale','log');
% hold on
% % plot(0:nb_iter, exp(P1(1)*(0:nb_iter) + P1(2)), '--k')
% % legend({'$\sqrt{\|\psi_1(t)-\hat\psi_1^{2n}(t)\|^2+\|\psi_2(t)-\hat\psi_2^{2n}(t)\|^2}$', ['$y = '  num2str(exp(P1(2))) ' \times ' num2str(exp(P1(1))) '^n$' ]}, 'Interpreter', 'latex', 'Location', 'northoutside');
% xlabel('Number of iterations $2n$ of BFN', 'Interpreter', 'latex')
% ylabel('$L^2$-error $\psi-\hat\psi^{2n}$', 'Interpreter', 'latex')
% % plot(0:nb_iter, exp(P2(2)) * ((0:nb_iter)^P2(1)));
% % set(gca,'xscale','log');
% % legend('Err', ['y = '  num2str(exp(P2(2))) ' * n^{' num2str(P2(1)) '}'], 'c')
% 
% % scaling
% fig2.Units               = 'centimeters';
% fig2.Position(3)         = opts.width;
% fig2.Position(4)         = opts.height;
% 
% % set text properties
% set(fig2.Children, ...
%     'FontName',     'Times', ...
%     'FontSize',     9);
% 
% % remove unnecessary white space
% figure(2);
% set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))



fig = figure(3);
figure(3)
% subplot(1, 2, 1)
plot(X1t(:, i)/10, N1, ':k')
hold on
plot(X1t(:, i20)/10, N1hat(:, i20), '--k')
hold on
plot(X1t(:, i)/10, N1hat(:, i), '-k')
title('$\eta_1 = 0.5$', 'Interpreter', 'latex')
axis([x_min-t_max/x_min^2, x_max+2*t_max/x_max^2, 0, 10]/10)
legend({'$\psi_1(t, r)$', '$\hat\psi_1^{10}(t, r)$', '$\hat\psi_1^{1000}(t, r)$'}, 'Interpreter', 'latex', 'Location', 'northeast')
xlabel('$r$', 'Interpreter', 'latex')
ylabel('PSD', 'Interpreter', 'latex')

% subplot(1, 2, 2)
fig4 = figure(4);
plot(X2t(:, i)/10, N2, ':k')
hold on
plot(X2t(:, i20)/10, N2hat(:, i20), '--k')
hold on
plot(X2t(:, i)/10, N2hat(:, i), '-k')
axis([x_min-t_max/x_min^2, x_max+2*t_max/x_max^2, 0, 10]/10)
xlabel('$r$', 'Interpreter', 'latex')
ylabel('PSD', 'Interpreter', 'latex')
title('$\eta_2 = 2$', 'Interpreter', 'latex')
legend({'$\psi_2(t, r)$', '$\hat\psi_2^{10}(t, r)$', '$\hat\psi_2^{1000}(t, r)$'}, 'Interpreter', 'latex', 'Location', 'northwest')


% scaling
fig.Units               = 'centimeters';
fig.Position(3)         = opts.width;
fig.Position(4)         = opts.height;

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
fig4.Position(4)         = opts.height;

% set text properties
set(fig4.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);

% remove unnecessary white space
figure(4);
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))

%% Save the movie :
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













