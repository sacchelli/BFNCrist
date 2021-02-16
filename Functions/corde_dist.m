function [Q, Q_cumul] = corde_dist(N, Densite, L, X)
    Nt = size(N, 2);
    Nl = size(L, 1);
    
    C = d_mat(X) ./ d_mat(L) .* Densite';
    
    Q_cumul = zeros(Nl, Nt);
    Q = zeros(Nl, Nt);
    for i = 1:Nt
        Q_cumul(:, i) = C * N(:, i);
    end
    
    Q(1, :) = Q_cumul(1, :);
    for k = 2:Nl
        Q(k, :) = Q_cumul(k, :) - Q_cumul(k-1, :);
    end
    
end