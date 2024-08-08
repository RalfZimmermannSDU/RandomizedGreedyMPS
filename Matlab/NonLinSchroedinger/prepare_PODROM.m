function [ROM] = prepare_PODROM(FOM, r)
%
% INPUT  : FOM struct
% OUTPUT : ROM struct
% side effects: opens a figure to show singular value decay
%
% FOM.snapshots contains the q- and p-snapshots stacked as
%     |q1,q2,...,qm|
% Y = |p1,p2,...,pm| 
%
N          = FOM.dimN;
ROM.dimPOD = r;
ROM.Yq     = FOM.snapshots(1:N,:);
ROM.Yp     = FOM.snapshots(N+1:2*N,:);

% create separate POD bases for Yq and Yp
% set up POD states: Uq*qr(t) ~ q(t),  Up*pr(t) ~ p(t)
[Uq,Sq,~] = svd(ROM.Yq);
ROM.Uq     = Uq(:,1:r);
[Up,Sp,~] = svd(ROM.Yp);
ROM.Up     = Up(:,1:r);
dSq = diag(Sq);
dSp = diag(Sp);
figure;
semilogy(1:N, diag(Sq), 1:N, diag(Sq), 'k-', 1:10:N, dSq(1:10:N), 'k*',...
         1:N, diag(Sp), 1:N, diag(Sp), 'b--', 5:10:N, dSp(5:10:N), 'bs');
legend('Singular values of Yq','Singular values of Yp')

% the projected ODE is
% qr(t) =  (Uq'*Dxx*Up)*pr + eps*Uq'*b(Uq*qr(t),Up*pr(t)) = \nabla_pH
% pr(t) = -(Up'*Dxx*Uq)*qr + eps*Up'*a(Uq*qr(t),Up*pr(t)) =-\nabla_qH
% see eqs. (15),(16) in corrsponding paper

% The linear operator is projected as
ROM.lin_op = ROM.Up'*(FOM.lin_op*ROM.Uq); % this is Dxxr
% Dxx is sym. => Dxxr' = (Up'*Dxx*Uq)
% Dxxr  features in grad_qH, Dxxr' features in grad_pH

end