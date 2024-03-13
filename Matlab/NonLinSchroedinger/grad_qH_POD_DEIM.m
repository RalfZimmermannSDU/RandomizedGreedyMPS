function [grad_qHqp] = grad_qH_POD_DEIM(qr,pr,A, para, Uq, Up, a_DEIM_op, Pa)
% evaluation of the nonlinear right-hand side function
% grad_qH(q,p), depending on a scalar parameter
%
% A is the finite differences operator projectes as A=Dxxr=Up'*(Dxx*Uq)

q = Uq(Pa,:)*qr;  % = Pa'*Uq*qr
p = Up(Pa,:)*pr;  % = Pa'*Uq*pr

a_nonlin  = (q.*q + p.*p).*q;
% apply DEIM operator
a_DEIM = a_DEIM_op*a_nonlin;

grad_qHqp = A*qr + para*a_DEIM;

return;
end

