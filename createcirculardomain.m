function [newdomain]=createcirculardomain(diam,sizesimspace)
% function [newdomain]=createcirculardomain(diam,sizesimspace)
% simulation space: sizesimspace x sizesimspace pixels de 1 nm
% 1 cercle with diameter diam
%
% Marianne Renner 04/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newdomain=zeros(sizesimspace,sizesimspace);
poscentro=round(sizesimspace/2);

centro=round(diam/2);
radio=diam/2;

synapse=zeros(diam+1,diam+1); %total image
for i=1:diam+1
    for j=1:diam+1
        distcentro=sqrt((centro-i)^2+(centro-j)^2);
        if distcentro<diam/2
            synapse(i,j)=1;
        end
    end
end
t=1;
for i=poscentro-radio : poscentro+radio
    k=1;
    for j=poscentro-radio:poscentro+radio
        newdomain(i,j)=synapse(t,k);
        k=k+1;
    end
    t=t+1;
end

%figure
%imshow(newdomain,'InitialMagnification','fit')

clear synapse synapse2 numberpergroup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
