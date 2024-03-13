function [grad_pHqp] = grad_pH_POD_DEIM(qr,pr,A, para, Uq, Up, b_DEIM_op, Pb)
% evaluation of the nonlinear right-hand side function
% grad_pH(q,p), depending on a scalar parameter
%
% A is the finite differences operator projectes as A=Dxxr=Uq'*(Dxx*Up)

q = Uq(Pb,:)*qr;  % = Pb'*Uq*qr
p = Up(Pb,:)*pr;  % = Pb'*Uq*pr

b_nonlin  = (q.*q + p.*p).*p;
% apply DEIM operator
b_DEIM = b_DEIM_op*b_nonlin;

grad_pHqp = A*pr + para*b_DEIM;

return;
end

