function [grad_qHqp] = grad_qH(q, p, A, para)
%--------------------------------------------------------------------------
% evaluation of the nonlinear right-hand side function
% grad_qH(q,p) of nonlinear Schrodinger 
% depending on a scalar parameter
%
% full-order version corresponding to eq. (13) in the associated paper
% "Randomized greedy magic point selection schemes for 
%  nonlinear model reduction"
% Ralf Zimmermann and Kai Cheng, ACOM
%
% INPUTS:
%  q        : q-state vector
%  p        : p-state vector
%  A        : FD operator A=Dxx
%  para     : scalar parameter (here epsilon)
%
% OUTPUTS
%  grad_qHqp : q-gradient from right hand side of eq. (13) in paper
%--------------------------------------------------------------------------
a_nonlin  = (q.*q + p.*p).*q;
grad_qHqp = A*q + para*a_nonlin;

return;
end

