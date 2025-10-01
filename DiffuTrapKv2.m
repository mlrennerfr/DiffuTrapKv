function DiffuTrapKv2(tours)
%function DiffuTrapKv2(tours)
% it creates un space with 1 circular domain where molecules can undergo
% immobilizations due to interactions
% each molecule may have 0 to 4 interaction sites (Cter)
% it creates nrotraj random trajectories of length lengthtrc, with given D
% steps of 1 ms initially, converted to the "real" till given
%
% Marianne Renner 04/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<1
    tours=1;
end

current_dir=cd;

%%%%%%%%%%%%%
till=1; %1 ms
%%%%%%%%%%%

% voir
Dslot=0;
Pextra=0;
Psyn=0;
distconfextra=50;
distconfsyn=10;

prompt = {'Number of Cter - number of molecules:','Diffusion coefficient:',...
    'Trajectory length (frames):','Cercle diameter (nm):','Simulation space (pixels):',...
    'kon (s-1)','koff (s-1):','Time between images (ms):','Pixel size (nm):','Localization accuracy (nm):','Molecule diameter (nm):'};
num_lines= 1;
dlg_title = 'Simulation parameters 1';
def = {'1 100','0.1','50','700','8','0.7 0.7','0.15','75','167','15','10'}; % default values

answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0
       return; 
end
ctergroups=str2num(answer{1});     % several values possible

Dvalue=str2num(answer{2})  ;      %
lengthtrcframes=str2num(answer{3}); % in simulation frames
cerclediam=str2num(answer{4});
sizesimpx=str2num(answer{5});

%valkon=str2num(answer{6}); 
allkon=str2num(answer{6}); 
valkon=allkon(1);
valkonout=allkon(2);

valkoff=str2num(answer{7});

realtill=str2num(answer{8});
szpx=str2num(answer{9});
paccu=str2num(answer{10});
diammol=str2num(answer{11}); 

lengthtrc=lengthtrcframes* realtill/till; % in simulation frames
sizesimspace=sizesimpx*szpx;

probinterac=valkon/1000 ; %prob pour till = 1ms
probinteracout=valkonout/1000 ; %prob pour till = 1ms

probfree=valkoff/1000;

  %  ctergroups(1)=3 %4;
  %  ctergroups(2)=val1(seriescounter);
  %  ctergroups(3)= 1 % 2% 1  % 
  %  ctergroups(4)=val2(seriescounter);
    
numbergroups=size(ctergroups,2)/2;
nrotraj=0;
listecter=[];
reportcter=[];
poscter=1;

for i=1:numbergroups
    poscter=i*2-1;
    for j=1:ctergroups(poscter+1) % for all the molecules of this group
        nrotraj=nrotraj+1;
        listecter=[listecter; nrotraj ctergroups(poscter) 0 0 0 0 ]; % mol, C ter , colonnes posible Cter for interaction status
    end
    reportcter=[reportcter; ctergroups(poscter) ctergroups(poscter+1)];
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

waitbarhandle=waitbar( 0,'Please wait...','Name','Simulating trajectories') ;
 
%--------------------------------------------------------------------------
% create domain for interactions
domain=createcirculardomain(cerclediam,sizesimspace);
%figure
%imshow(domain,'InitialMagnification','fit')

%------------------------------------------------------------------------
% initialize
pos{nrotraj,2}=[]; % positions all trajectories, all simulation frames 
newtrcdata=[]; %preallocate!!!! % trc files, in wanted frame rate
matrixposmol=zeros(nrotraj,5);  % positions mol at each frame and interaction status

%-------------------------------------------------------------------------

% first point 
nro=1;
firstpos=1;
prob=zeros(nro,1);
D=zeros(nro,1);
sizeconf=zeros(nro,1);
first=zeros(nro,1);
probconf=zeros(nro,1);

while nro<nrotraj+1
    if exist('waitbarhandle') %#ok<EXIST>
        waitbar(1/lengthtrcframes,waitbarhandle,'Initializing...');
    end
    
    % chooses an initial position free of molecules -----------------------
    forbidden=1;
    while forbidden
           pos{nro}(1,1)=((rand-0.5)*sizesimspace)+sizesimspace/2; %initial pos in x at random
           pos{nro}(1,2)=((rand-0.5)*sizesimspace)+sizesimspace/2; %initial pos in y at random 
           matrixposmol(nro,1)=nro;
           matrixposmol(nro,2)=pos{nro}(1,1);
           matrixposmol(nro,3)=pos{nro}(1,2);
           matrixposmol(nro,4)=0; %not interacting for now
           matrixposmol(nro,5)=listecter(nro,2); % Cter  VOIR not needed?
           
           if nro>1
               forbidden=checkconflict3(matrixposmol, nro, firstpos,diammol);
           else
               forbidden=0;
           end
    end
    
    % inside domain at the beginning? -----------------------------------------------------------
                posy=floor(pos{nro}(1,2));
            if posy==0; posy=1; end
            posx=floor(pos{nro}(1,1));
            if posx==0; posx=1; end

    if domain(posy,posx)==0 %transposed!!!!!!!!
        %extra
        D(nro)=Dvalue;        
        prob(nro)=Pextra; sizeconf(nro)=distconfextra;
    else
       % if domain(round(pos{nro}(1,2)),round(pos{nro}(1,1)))==1 % inside domain
           % disp(nro)
            D(nro)=Dvalue;
            prob(nro)=Psyn; sizeconf(nro)=distconfsyn;
            
            % interactions : calculation of probabilites
            valCter=matrixposmol(nro,5); % how many interactions sites
            if valCter>0
                for ct=1:valCter
                    % calculate binding probability for each Cter
                    binding=rand; 
                    if probinterac>binding %has to interact
                        listecter(nro,ct+2)=1;
                    end
                end
            end
            matrixposmol(nro,4)=max(listecter(nro,3:6)); % if only one site interacting, the status in interacting
            if matrixposmol(nro,4)==1
                D(nro)=Dslot;
                prob(nro)=1;
            end
    
       % end %inside domain
    end   %check localization
    
    first(nro)=1;
    probconf(nro)=prob(nro);
    if pos{nro}(1,1)<1;         pos{nro}(1,1)=1; end
    if pos{nro}(1,2)<1;         pos{nro}(1,2)=1; end
   % newtrcdata=[newtrcdata; nro 1 pos{nro}(1,1) pos{nro}(1,2) D(nro) domain(round(pos{nro}(1,2)),round(pos{nro}(1,1))) matrixposmol(nro,4) matrixposmol(nro,5)];
    nro=nro+1;
    
end % loop trajectories first point
%-----------------------------------------------------------------------
firstpos=0;

% loop points 

if exist('waitbarhandle') %#ok<EXIST>
    waitbar(2/lengthtrcframes,waitbarhandle,'Frame # 1');
end

nexttill=realtill/till+1; % first point
counttill=2;

for i=2: lengthtrc % loop all frames from frame=2
    
    %  all trajectories at frame=i
    for nro=1:nrotraj              
        
        % ver
        if probconf(nro)==prob(nro) % same conf than before
        else
            if prob(nro)==1
                first(nro)=1; %start of strong confinement
            else
                first(nro)=0; % not strong confinement
            end
        end
        probconf(nro)=prob(nro);
        
        % displacements - pixel simulation=1 nm
        ang=rand*360/57.29; %radians
        Displ=abs(normrnd(0,sqrt(2*D(nro)*(till/1000)/((1/1000)^2)))); % pixel/frame   
        
        % confinement or obstacles
        conf=rand; %random probability  
        
        if conf<probconf(nro)  
            if probconf(nro)==1 % immo or ++ conf: stop periods
                if first(nro)==1
                    center{nro}=pos{nro}(1,:); %#ok<*AGROW> % center of the confinement area
                    first(nro)=0;
                end
            else % confinement soft: obstacles
                center{nro}=[];
            end
            if Displ>sizeconf(nro) % displacement in x outside of the confinement area/obstacle
                veces=floor(Displ/sizeconf(nro));
                Displ=sizeconf(nro)-(Displ-veces*sizeconf(nro));
                controldepas=1;
            else
                controldepas=0;
            end
        else
            first(nro)=1;        
            center{nro}=[];
        end %conf   
        
        % new positions
        if isempty(center{nro})==0 && controldepas==1  % confined
            newposx=center{nro}(1)+Displ*sin(ang);
            newposy=center{nro}(2)+Displ*cos(ang);
        else  % not confined
            newposx=pos{nro}(1,1)+Displ*sin(ang);
            newposy=pos{nro}(1,2)+Displ*cos(ang);
        end
        matrixposmol(nro,2)=newposx;
        matrixposmol(nro,3)=newposy;
        
        %check next point
       % [possynx,possyny,newposx, newposy]=defineposmol2(newposx, newposy, sizesimspace,ang);
        [possynx,possyny,newposx, newposy]=defineposmol2(newposx, newposy, sizesimspace);
        
        %----------------------------------------------------------------
        % check forbidden: positions around    
        forbidden=checkconflict3(matrixposmol, nro, firstpos,diammol);
        
        if forbidden==1
             % att the same mol!!!
                % new position
                if D(nro)==Dslot % does not move
                    newposx=center{nro}(1);
                    newposy=center{nro}(2);
                else
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    newdispl=Displ-max(Displ/2,15); %recul min 15 nm
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if isempty(center{nro})==0 && controldepas==1  % confined
                        newposx=center{nro}(1)+newdispl*sin(ang);
                        newposy=center{nro}(2)+newdispl*cos(ang);
                    else  % not confined
                        newposx=pos{nro}(1,1)+newdispl*sin(ang);
                        newposy=pos{nro}(1,2)+newdispl*cos(ang);
                    end
                end
               % [possynx,possyny,newposx, newposy]=defineposmol2(newposx, newposy, sizesimspace,ang);
                [possynx,possyny,newposx, newposy]=defineposmol2(newposx, newposy, sizesimspace);
                
                matrixposmol(nro,2)=newposx;
                matrixposmol(nro,3)=newposy;
                
                 % control again
                forbidden=checkconflict3(matrixposmol, nro, firstpos,diammol);
                if forbidden==1                    
                        % cannot move
                        newposx=floor(pos{nro}(1,1));
                        newposy=floor(pos{nro}(1,2)); 
                        %[possynx,possyny,newposx, newposy]=defineposmol2(newposx, newposy, sizesimspace,ang);
                        [possynx,possyny,newposx, newposy]=defineposmol2(newposx, newposy, sizesimspace);
                        matrixposmol(nro,2)=newposx;
                        matrixposmol(nro,3)=newposy;
                end % forbidden 2
         end % forbidden 1
        %--------------------------------------------------------------
        
        % check ongoing interactions
        
         %   contabilize interacting time
        valCter=matrixposmol(nro,5);
        if matrixposmol(nro,4)==1
            for ct=1:valCter
                %calculates all the sites independently
                if listecter(nro,ct+2)>0 %interacting
                    % calculate Pfree
                    free=rand; 
                    if probfree>free %getps free
                        listecter(nro,ct+2)=0;
                    end
                end
            end
            matrixposmol(nro,4)=max(listecter(nro,3:6)); % if at least one site interacting, the status in interacting
            
        else % not interacting
            posy=floor(pos{nro}(1,2));
            if posy==0; posy=1; end
            posx=floor(pos{nro}(1,1));
            if posx==0; posx=1; end
            if domain(posy,posx)==1 % inside domain
                %calculate Pbind
                for ct=1:valCter
                    %calculates all the sites independently
                    binding=rand; 
                    if probinterac>binding % interacts
                        listecter(nro,ct+2)=1;
                    end
                end
                matrixposmol(nro,4)=max(listecter(nro,3:6)); % if at least one site interacting, the status in interacting
            else %out domain
                %calculate Pbind outside domain
                for ct=1:valCter
                    %calculates all the sites independently
                    binding=rand; 
                    if probinteracout>binding % interacts
                        listecter(nro,ct+2)=1;
                    end
                end
                matrixposmol(nro,4)=max(listecter(nro,3:6)); % if at least one site interacting, the status in interacting
            end
        end
        
        if matrixposmol(nro,4)==0 %free
            D(nro)=Dvalue;
            prob(nro)=Psyn; sizeconf(nro)=distconfsyn;
        else
            
            D(nro)=Dslot;
            prob(nro)=1;
        end
  
        
        pos{nro}(1,1)=newposx;
        pos{nro}(1,2)=newposy;
        
         matrixposmol(nro,2)=newposx;
         matrixposmol(nro,3)=newposy;

        if i==nexttill % VER
        %if i>75000 && i==nexttill % VER
           newtrcdata=[newtrcdata; nro counttill newposx newposy D(nro) domain(possynx,possyny) matrixposmol(nro,4) matrixposmol(nro,5)];
        end
          
        
    end % trajectories

   
    %-------------------------------------------------------------------
    if i==nexttill % VER

        counttill=counttill+1;
        if exist('waitbarhandle') %#ok<EXIST>
            waitbar(counttill/lengthtrcframes,waitbarhandle,['Frame # ',num2str(counttill)]);
        end
        nexttill=counttill*realtill/till; % VER
    end
    %---------------------------------------------------------------------
    
end% points trajectories

if exist('waitbarhandle') %#ok<EXIST>
   close(waitbarhandle);
end


%--------------------------------------------------------------------------
% prepare and save data trajectories


% save data
if isfolder('trc'); else ; mkdir('trc'); end
currentdir=cd;

cd('trc')
d=dir('*-all.con.trc*');
st = {d.name};
nroorder=size(st,2)+1;
savename = 'sim-' ;

% noise
for t=1:size(newtrcdata,1)
    noisex=normrnd(0,paccu/2);
    noisey=normrnd(0,paccu/2);   
    newtrcdata(t,3)=newtrcdata(t,3)+noisex ;
    newtrcdata(t,4)=newtrcdata(t,4)+noisey ;
end



% sort by group
for gg=0:max(newtrcdata(:,8))
    index=find(newtrcdata(:,8)==gg);
    if isempty(index)==0
        aux=newtrcdata(index,:);
        save(['cter',num2str(gg),'-',num2str(nroorder),'.raw.trc'],'aux','-ascii'); % only syn loc
    end
    clear aux
end
% save all
%save([savename,'.raw.trc'],'newtrcdata','-ascii'); % in nm, for plotting

% convertion to desired pixel size 
for t=1:size(newtrcdata,1)
    newtrcdata(t,3)=newtrcdata(t,3)/szpx ;
    newtrcdata(t,4)=newtrcdata(t,4)/szpx ;
end

%no localization
noloctrc=newtrcdata;
noloctrc(:,6)=zeros(size(noloctrc,1),1);

% sort by group
for gg=0:max(newtrcdata(:,8))
    index=find(newtrcdata(:,8)==gg);
    if isempty(index)==0
        aux=newtrcdata(index,:);
        auxnoloc=noloctrc(index,:);
        counter=1;
        for nn=1:max(aux(:,1))
            indexmol=find(aux(:,1)==nn);
            if isempty(indexmol)==0
                aux(indexmol,1)=counter;
                auxnoloc(indexmol,1)=counter;
                counter=counter+1;
            end
        end
       % save(['cter',num2str(gg),'-',num2str(nroorder),'.loc.trc'],'aux','-ascii'); %  with loc
        save(['cter',num2str(gg),'-',num2str(nroorder),'.con.trc'],'auxnoloc','-ascii'); %  without loc
    end
end

% save trc with all groups
%save([savename,'-all.loc.trc'],'newtrcdata','-ascii'); % with loc
%save([savename,'-all.con.trc'],'noloctrc','-ascii'); % without loc

cd(currentdir)


%-------------------------------------------------------------------------
%report simulation

name=['report',savename,'.txt'];
fi = fopen(name,'wt');
if fi<3
    error('File not found or readerror.');
end
report{1}= ['Simulation number ',num2str(nroorder),' - ',num2str(nrotraj),'  trajectories'];
report{2}= '  ';
report{3}= ['D: ',num2str(Dvalue),'(µm2/s)'];
report{4}= ['Trajectory length (frames): ',num2str(lengthtrcframes),'. Simulation space : ',num2str(sizesimpx),' pixels.'];
report{5}= ['Domain diameter: ',num2str(cerclediam),' nm'];
report{6}= ['koff: ',num2str(valkoff),' s-1. kon :', num2str(valkon)];
report{7}= ['Real till: ',num2str(realtill),' nm    Real pixel size: ',num2str(szpx),' nm   Loc accuracy : ', num2str(paccu),'  nm'];
report{8}= [num2str(numbergroups),' types of molecules'];
nroorder=8;
for ct=1:size(reportcter,1)
    nroorder=nroorder+1;
    report{nroorder}=[num2str(reportcter(ct,2)),' molecules with ',num2str(reportcter(ct,1)),' Cter'];
    ct=ct+1;
end
for celda=1:nroorder
    fseek(fi,10,0);
    fprintf(fi,'%10s\n',report{celda});
end
fclose(fi);

clear trcdata newtrc newtrcdata synapse newsynapse newsynapse2 report occupation
cd(current_dir);
