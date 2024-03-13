function [phi] = deim(uIn)
%
% The classical DEIM algorithm of 
% Chaturantabut, S., Sorensen, D.: Nonlinear model reduction via discrete
% empirical interpolation. SIAM Journal on Scientific Computing 32(5),
% 2737–2764 (2010)
%
% The output is an index list, not a mask matrix
%
[n,m] = size(uIn);
[~, phi_opt] = max(abs(uIn(:, 1)));
phi = [phi_opt];
U(:, 1) = uIn(:, 1);
for l = 2:m
    uL = uIn(:, l);
    c = U(phi, :)\uL(phi);
    r = uL - U*c;
    [~, ind_l_opt] = max(abs(r));
    U(:, l) = uL;
    phi = [phi, ind_l_opt];
end
phi = sort(phi);
end