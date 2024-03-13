function [grad_pHqp] = grad_pH(q,p,A, para)
% evaluation of the nonlinear right-hand side function
% grad_pH(q,p), depending on a scalar parameter
%
% A is the finite differences operator
b_nonlin  = (q.*q + p.*p).*p;
grad_pHqp = A*p + para*b_nonlin;

return;
end

