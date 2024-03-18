function [grad_qHqp] = grad_qH_POD(qr, pr, A, para, Uq, Up)
%--------------------------------------------------------------------------
% evaluation of the nonlinear right-hand side function
% grad_qH(q,p) of nonlinear Schrodinger 
% depending on a scalar parameter
%
% POD version corresponding to eq. (16) in the associated paper
% "Randomized greedy magic point selection schemes for 
%  nonlinear model reduction"
% Ralf Zimmermann and Kai Cheng, ACOM
%
% INPUTS:
%  qr        : reduced q-state vector
%  pr        : reduced p-state vector
%  A         : projected FD operator A=Dxxr'=Up'*(Dxx*Uq)
%  para      : scalar parameter (here epsilon)
%  Uq        : POD basis for q-states
%  Up        : POD basis for p-states
%
% OUTPUTS
%  grad_qHqp : gradient from right hand side of eq. (16) in paper
%--------------------------------------------------------------------------
q = Uq*qr;
p = Up*pr; % lift reduced POD states to FOM dimension
grad_qHqp = A*qr + para*Up'*((q.*q + p.*p).*q);

return;
end

