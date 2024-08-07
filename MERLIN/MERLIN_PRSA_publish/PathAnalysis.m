function [Uhis,load_his,Data] = PathAnalysis(truss,angles,F,b_lambda,MaxIcr)
tol = 1e-6; MaxIter = 50; 
Node = truss.Node;
AllDofs = [1:3*size(Node,1)];
U = truss.U0;
Uhis = zeros(3*size(Node,1),MaxIcr);
Data.Exbar = zeros(size(truss.Bars,1),MaxIcr); 
Data.FdAngle = zeros(size(angles.fold,1),MaxIcr); Data.LFdAngle = Data.FdAngle;
Data.BdAngle = zeros(size(angles.bend,1),MaxIcr); Data.LBdAngle = Data.BdAngle;

FreeDofs = setdiff(AllDofs,truss.FixedDofs);
lmd = 0; icrm = 0; MUL = [U,U];
load_his = zeros(MaxIcr,1);
F = F(:,1);
while icrm<MaxIcr
    icrm = icrm+1;
    iter = 0; err = 1;
    fprintf('icrm = %d, lambda = %6.4f\n',icrm,lmd);
    while err>tol && iter<MaxIter
        iter = iter+1; 
        % There are two versions of the function for assembling global
        % stiffness matrix:
        % 'GlobalK_fast_ver': speed-optimized version;
        % 'GlobalK_edu_ver': easy-to-read educational version.
        [IF,K] = GlobalK_fast_ver(U,Node,truss,angles);
        % [IF,K] = GlobalK_edu_ver(U,Node,truss,angles);
        R = lmd*F-IF;   MRS = [F,R];
        MUL(FreeDofs,:) = K(FreeDofs,FreeDofs)\MRS(FreeDofs,:);
        dUp = MUL(:,1); dUr = MUL(:,2);
         if iter==1, dUr = 0*dUr; end
        dlmd=nlsmgd(icrm,iter,dUp,dUr,b_lambda);
        dUt = dlmd*dUp+dUr;
        U = U+dUt;
        err = norm(dUt(FreeDofs));
        lmd = lmd+dlmd;
        fprintf('    iter = %d, err = %6.4f, dlambda = %6.4f\n',iter,err,dlmd);
        if err > 1e8, disp('Divergence!'); break; end
    end

    if iter>15
        b_lambda = b_lambda/2;
        disp('Reduce constraint radius!')
        icrm = icrm-1;
        U = Uhis(:,max(icrm,1));
        lmd = load_his(max(icrm,1));
    elseif iter<3
        disp('Increase constraint radius!')
        b_lambda = b_lambda*1.5;
        Uhis(:,icrm) = U;
        load_his(icrm) = lmd; 
        [Exbari,FdAnglei,BdAnglei,LFdAnglei,LBdAnglei] = GetData(U,Node,truss,angles);
        Data.Exbar(:,icrm) = Exbari; 
        Data.FdAngle(:,icrm) = FdAnglei; Data.LFdAngle(:,icrm) = LFdAnglei;
        Data.BdAngle(:,icrm) = BdAnglei; Data.LBdAngle(:,icrm) = LBdAnglei;
    else
        Uhis(:,icrm) = U;
        load_his(icrm) = lmd; 
        [Exbari,FdAnglei,BdAnglei,LFdAnglei,LBdAnglei] = GetData(U,Node,truss,angles);
        Data.Exbar(:,icrm) = Exbari; 
        Data.FdAngle(:,icrm) = FdAnglei; Data.LFdAngle(:,icrm) = LFdAnglei;
        Data.BdAngle(:,icrm) = BdAnglei; Data.LBdAngle(:,icrm) = LBdAnglei;
    end
end

icrm = icrm+1;
Uhis(:,icrm:end) = [];
load_his(icrm:end,:) = [];
Data.Exbar(:,icrm:end) = []; 
Data.FdAngle(:,icrm:end) = []; Data.LFdAngle(:,icrm:end) = [];
Data.BdAngle(:,icrm:end) = []; Data.LBdAngle(:,icrm:end) = [];