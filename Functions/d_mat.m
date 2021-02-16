function dx = d_mat(X)
    dx = diff(X);
    dx = [dx; dx(end)];
end