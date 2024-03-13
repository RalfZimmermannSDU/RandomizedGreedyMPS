function [grad_pHqp] = grad_pH_POD(qr,pr,A, para, Uq, Up)
% evaluation of the nonlinear right-hand side function
% grad_pH(q,p), depending on a scalar parameter
%
% POD-ROM version
% A is the finite differences operator

q = Uq*qr;
p = Up*pr;
grad_pHqp = A*pr + para*Uq'*((q.*q + p.*p).*p);

return;
end

