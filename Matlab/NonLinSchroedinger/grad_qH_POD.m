function [grad_qHqp] = grad_qH_POD(qr,pr,A, para, Uq, Up)
% evaluation of the nonlinear right-hand side function
% grad_qH(q,p), depending on a scalar parameter
%
% POD-DEIM-ROM version
% A is the finite differences operator projectes as 

q = Uq*qr;
p = Up*pr;
grad_qHqp = A*qr + para*Up'*((q.*q + p.*p).*q);

return;
end

