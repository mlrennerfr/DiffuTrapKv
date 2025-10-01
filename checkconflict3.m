function forbidden=checkconflict3(matrixposmol, nro, firstpos,diammol,group)
%function forbidden=checkconflict3(matrixposmol, nro, firstpos,diammol,group);
%
% called by DoSimulate.m
% compares the next position of the molecule nro with the actual positions
% of the molecules nearby (diammol= molecules diameter).
% forbidden=1 if the next position is not possible (excluded volume due to
% another molecule)
%
% Marianne Renner nov 12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<5
    group=1;
end

forbidden=0;
if size(diammol,2)<2
    diam=diammol;
else
    diam=(diammol(1)+diammol(2))/2; %in case of two sizes, it chooses the larger one
end

if firstpos==1
    dist=diam*1.5; %first position, not in touch with others
else
    dist=diam; % min distance: 10 nm
end

posx=matrixposmol(nro,2);
posy=matrixposmol(nro,3);

indexx=find(abs(matrixposmol(:,2)-posx)<dist); % first approximate
if isempty(indexx)==0
    indexy=find(abs(matrixposmol(indexx,3)-posy)<dist); % first approximate
    if isempty(indexy)==0
        indexmol=find(matrixposmol(indexx(indexy),1)~= nro) ;
        if isempty(indexmol)==0 % there is(are) another molecules(s) too close
            % check real distance
            for m=1:size(indexmol,1)
                
                difx=abs(matrixposmol(indexx(indexy(indexmol(m))),2)-posx);
                dify=abs(matrixposmol(indexx(indexy(indexmol(m))),3)-posy);

                if isempty(group)==0
                    if nro>size(group,2)
                      %  mindist=min(diammol)
                       mindist=dist;
                    else
                        nromolcheck=matrixposmol(indexx(indexy(indexmol(m))),1);
                       % mindist=diammol(group(nro))/2 +diammol(group(nromolcheck))/2 ; %sum of two radius
                        mindist=diam;
                    end
                end
                if sqrt(difx^2+dify^2)<mindist 
                    forbidden=1;
                end
            end
        end
    end
end

%%%%%%eof %%%%%%%
