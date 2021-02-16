function film_inv(sol, N, Ninv, Q, X, T, L, Densite)
    Nt = size(T, 1);
         
    if sol
    figure(1)
    for i = 1:Nt
        plot(X, N(:, i))
        title('Distribution en taille des cristaux')
        hold on
        plot(X, Ninv(:, i))
        axis([X(1), X(end), min(min(N)), max(max(N))])
        xlabel('x')
        ylabel('n(t, x)')
%         legend('n simu', 'n reconstruit')
        drawnow
        hold off
        pause(0.01)
    end
    
    Q = corde_dist(N, Densite, L, X);
    
    figure(2)
    for i = 1:Nt
        plot(L, Q(:, i))
        title('Distribution en taille des cordes')
        hold on
        plot(L, Q(:, i))
        axis([L(1), L(end), min(min(Q)), max(max(Q))])
        xlabel('l')
        ylabel('q(t, l)')
        drawnow
        hold off
        pause(0.01)
    end
    
    end
        
%     Y_approx = H(N_approx, X);
%     
%     figure
%     plot(T, Y_exact)
%     hold on
%     plot(T, Y_approx)
%     xlabel('t')
%     ylabel('y(t)')
%     legend('Y_{exact}', 'Y_{approx}')
%     title('Mesure')
% 
%     Err_conc = abs(Y_approx - Y_exact)./abs(Y_exact);
%     
%     figure
%     subplot(1, 2, 1)
%     plot(T, Err_conc)
%     xlabel('t')
%     ylabel('abs(Y_{approx} - Y_{exact})./abs(Y_{exact})')
%     title('Erreur relative')
%     subplot(1, 2, 2)
%     plot(T, log(Err_conc))
%     xlabel('t')
%     ylabel('log(abs(Y_{approx} - Y_{exact})./abs(Y_{exact}))')
%     title('Erreur relative en ?chelle semi-log')
%         
%     Err_fbrm = zeros(Nt, 1);
%     for i = 1:Nt
%         Err_fbrm(i) = norm(Q_exact(:,i) - Q_approx(:,i))./norm(Q_exact(:,i));
%     end
%     
%     figure
%     subplot(1, 2, 1)
%     plot(T, Err_fbrm)
%     xlabel('t')
%     ylabel('norm(Q_{exact} - Q_{approx})./norm(Q_{exact})')
%     title('Erreur relative')
%     subplot(1, 2, 2)
%     plot(T, log(Err_fbrm))
%     xlabel('t')
%     ylabel('log(norm(Q_{exact} - Q_{approx})./norm(Q_{exact}))')
%     title('Erreur relative en ?chelle semi-log')
%     
%     Err_etat = zeros(Nt, 1);
%     for i = 1:Nt
%         Err_etat(i) = norm(N_exact(:,i) - N_approx(:,i))./norm(N_exact(:,i));
%     end
%     
%     figure
%     subplot(1, 2, 1)
%     plot(T, Err_etat)
%     xlabel('t')
%     ylabel('norm(N_{exact} - N_{approx})./norm(N_{exact})')
%     title('Erreur relative')
%     subplot(1, 2, 2)
%     plot(T, log(Err_etat))
%     xlabel('t')
%     ylabel('log(norm(N_{exact} - N_{approx})./norm(N_{exact}))')
%     title('Erreur relative en ?chelle semi-log')
    
end