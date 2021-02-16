function Ninv = pb_inv(Q_cumul, X, L, Densite, delta)
 
Nx = size(X, 1);
K = d_mat(X) ./ d_mat(L) .* Densite';

options =  optimoptions(@quadprog, 'Display','off');

H = K' * K + delta*eye(Nx);        
f = -K' * Q_cumul;
Ninv = quadprog(H,f,[],[],[],[],zeros(Nx, 1),Inf*ones(Nx, 1), zeros(Nx, 1), options);


end