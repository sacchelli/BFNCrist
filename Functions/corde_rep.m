function Densite = corde_rep(X, L, nb, e)
    Nx = size(X, 1);
    Nl = size(L, 1);
    Densite = zeros(Nx, Nl);
    Phi = linspace(0, 2*pi, nb+1);
    Theta = linspace(0, pi, nb+1);
    
    for i = 1:Nx
        for j = 1:Nl
            
            x = X(i);
            l = L(j);
            
            k = 0;
            
            for it1 = 1:nb
                for it2 = 1:nb
                    
                    phi = Phi(it1);
                    theta = Theta(it2);
                    
                    
                    alpha = cos(phi)^2/(cos(theta)^2+e^2*sin(theta)^2) +sin(phi)^2;

                    temp = 1 - (l/(2*x))^2*alpha;
                    if temp>=0 && not(isnan(temp)) && not(isinf(temp))
                        
                        k = k + sqrt(temp)*sin(theta);
                        
                    end
                    
                end
            end

            k = 1 - k*pi/(2*nb*nb);
            
            
            Densite(i, j) = k;
                        
        end
    end
end