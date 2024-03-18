function [grad_pHqp] = grad_pH(q, p, A, para)
%--------------------------------------------------------------------------
% evaluation of the nonlinear right-hand side function
% grad_pH(q,p) of nonlinear Schrodinger 
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
%  grad_pHqp : p-gradient from right hand side of eq. (13) in paper
%--------------------------------------------------------------------------
b_nonlin  = (q.*q + p.*p).*p;
grad_pHqp = A*p + para*b_nonlin;

return;
end

