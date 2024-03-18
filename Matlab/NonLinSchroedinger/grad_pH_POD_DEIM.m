function [grad_pHqp] = grad_pH_POD_DEIM(qr, pr, A, para, Uq, Up, b_DEIM_op, Pb)
%--------------------------------------------------------------------------
% evaluation of the nonlinear right-hand side function
% grad_pH(q,p) of nonlinear Schrodinger 
% depending on a scalar parameter
%
% POD-DEIM version corresponding to eq. (19) in the associated paper
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
%  b_deim_op : DEIM operator for nonlinear b-terms, see eq. (17) in paper
%  Pb        : list of magic points collected for b-terms
%
% OUTPUTS
%  grad_pHqp : gradient from right hand side of eq. (19) in paper
%--------------------------------------------------------------------------
q = Uq(Pb,:)*qr;  % this realizes the operation Pb'*Uq*qr
p = Up(Pb,:)*pr;  % this realizes the operation Pb'*Up*pr

b_nonlin  = (q.*q + p.*p).*p;
% apply DEIM operator
b_DEIM = b_DEIM_op*b_nonlin;

grad_pHqp = A*pr + para*b_DEIM;  % see eq. (19)

return;
end

