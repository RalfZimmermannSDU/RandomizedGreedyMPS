function [Dxx, q0, p0, x_space] = prepare_FOM(N, L, c)
%--------------------------------------------------------------------------
% set up full ordeer model for Schrodinger equation
%
% INPUTS
%  N  : spatial dimension
%  L  : controls length of spatial interval
%  c  : wave speed constant
%OUTPUTS
% Dxx     : fnite difference operator associated with dxx
% q0      : initial q-state
% p0      : initial p-state
% x_space : discrete representation of spatial interval
%--------------------------------------------------------------------------
dx= L/N;        % spatial resolution
% set up FD operator
e = ones(N,1);
Dxx = spdiags([e -2*e e], -1:1, N, N);
Dxx(N,1) = 1;
Dxx(1,N) = 1;
Dxx = Dxx/dx^2;
% spatial discretization
x0 = 0.0;
x_space = linspace(x0-L/2,x0+L/2,N+1)';
% remove last point for use of periodic boundary conditions
x_space(end) = [];
% initial condition
u0 = (sqrt(2.0)./cosh(x_space-x0)).*exp(i*(c/2)*(x_space-x0));
q0 = imag(u0);
p0 = real(u0);
end