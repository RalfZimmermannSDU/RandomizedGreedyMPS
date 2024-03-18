function [grad_pHqp] = grad_pH_POD(qr, pr, A, para, Uq, Up)
%--------------------------------------------------------------------------
% evaluation of the nonlinear right-hand side function
% grad_pH(q,p) of nonlinear Schrodinger 
% depending on a scalar parameter
%
% POD version corresponding to eq. (15) in the associated paper
% "Randomized greedy magic point selection schemes for 
%  nonlinear model reduction"
% Ralf Zimmermann and Kai Cheng, ACOM
%
% INPUTS:
%  qr        : reduced q-state vector
%  pr        : reduced p-state vector
%  A         : projected FD operator A=Dxxr=Uq'*(Dxx*Up)
%  para      : scalar parameter (here epsilon)
%  Uq        : POD basis for q-states
%  Up        : POD basis for p-states
%
% OUTPUTS
%  grad_pHqp : gradient from right hand side of eq. (15) in paper
%--------------------------------------------------------------------------
q = Uq*qr;
p = Up*pr; % lift reduced states to FOM dimension
grad_pHqp = A*pr + para*Uq'*((q.*q + p.*p).*p);

return;
end

