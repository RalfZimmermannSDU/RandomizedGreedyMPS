function [grad_qHqp] = grad_qH_POD_DEIM(qr, pr, A, para, Uq, Up, a_DEIM_op, Pa)
%--------------------------------------------------------------------------
% evaluation of the nonlinear right-hand side function
% grad_qH(q,p) of nonlinear Schrodinger 
% depending on a scalar parameter
%
% POD-DEIM version corresponding to eq. (20) in the associated paper
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
%  a_deim_op : DEIM operator for nonlinear a-terms, see eq. (18) in paper
%  Pa        : list of magic points collected for a-terms
%
% OUTPUTS
%  grad_qHqp : gradient from right hand side of eq. (20) in paper
%--------------------------------------------------------------------------

q = Uq(Pa,:)*qr;  % this realizes the operation Pa'*Uq*qr
p = Up(Pa,:)*pr;  % this realizes the operation Pa'*Up*pr

a_nonlin  = (q.*q + p.*p).*q;
% apply DEIM operator
a_DEIM = a_DEIM_op*a_nonlin;

grad_qHqp = A*qr + para*a_DEIM;  % see eq. (20)

return;
end

