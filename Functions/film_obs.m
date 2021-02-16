function film_obs(N, Q, X, T, L, sol)
    Nt = size(T, 1);
    
    if sol
        
        figure
        u = N(:,1);
        Trace = plot(X, u);
        axis([X(1), X(end), min(min(N)), max(max(N))])
        xlabel('x')
        ylabel('n(t, x)')
        title('Approximation de n a t = 0')
        drawnow
        pause(0.01)
        for i = 2:Nt
            set(Trace, 'Ydata', N(:,i))
            title(['Approximation de n a t = ' num2str(T(i))])
            drawnow
            pause(0.01)
        end
        
        figure
        q = Q(:,1);
        Trace = plot(L, q);
        axis([L(1), L(end), min(min(Q)), max(max(Q))])
        xlabel('l')
        ylabel('q(t, l)')
        title('Approximation de q a t = 0')
        drawnow
        pause(0.01)
        for i = 2:Nt
            set(Trace, 'Ydata', Q(:,i))
            title(['Approximation de q a t = ' num2str(T(i))])
            drawnow
            pause(0.01)
        end
        
    end
    
end