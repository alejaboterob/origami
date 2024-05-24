function dl=nlsmgd(step,ite,dup,dur,cmp)
% Modified generalized displacement control method.

% Define the incremental load factor signal
global dupp1 sinal

if step==1 && ite==1
    sinal=1;
    dupp1=dup;
end

if ite==1
    sinal=sinal*sign(dot(dupp1,dup));
    dupp1=dup;
end

% Calculate the incremental load factor
global dupc1 numgsp

if ite==1
    if step==1
        dl=cmp;
        numgsp=dot(dup,dup);
        dupc1=dup;
    else
        gsp=numgsp/dot(dup,dup);
        dl=sinal*cmp*sqrt(gsp);
        dupc1=dup;
    end
else
    dl=-dot(dupc1,dur)/dot(dupc1,dup);
end