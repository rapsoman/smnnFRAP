function [trec,xrec,qrec,brec]=simulation_code(d,res_time,Fimm)
%code to simulate FRAP experiments

%INPUTS:
%d=diffusion sigma, res_time = residence time, Fimm=bound fraction

%OUTPUTS:
%trec: 1 x nobs vector of simulation times (for now set equal to tobs)
%qrec: N x nobs matrix of binding status (0=unbound,1=bound) of the N particles at times trec
%xrec: 3*N x nobs matrix of the 3D position of the N particles at times trec
%brec: N x nobs matrix of bleaching status of the N particles at times trec


plotmode=1; 
%-----------------time recording parameters--------------------------------
Dtobs=0.0680; 
tobs=(0:Dtobs:(Dtobs*301)); 
nobs=length(tobs);
%bleaching interval
bleachint=[49.5*Dtobs,51.5*Dtobs];

%----- ellipsoidal domain G: (inv(UG)x-cG)'inv(DG^2)(inv(UG)x-cG)=1--------
DG=diag([5,4,4]); %radii
UG=eye(3); %rotation
cG=[0;0;0]; %center
DG2inv=inv(DG^2);
diagDG2inv=diag(DG2inv)';
UGinv=inv(UG);
[XG1,YG1]=Create2DEllipsoid(cG(1:2),diag(DG(1:2,1:2)),eye(2)); %x,y axis - apo panw
[YG2,ZG2]=Create2DEllipsoid(cG(2:3),diag(DG(2:3,2:3)),eye(2)); %x,z axis - apo plai

%----------------------domain of FRAP experiment---------------------------
%bleaching shape
DF=diag([2,2,2]);
UF=eye(3);
cF=[0;0;0];
DF2inv=inv(DF^2);
UFinv=inv(UF);
[XF1,YF1]=Create2DEllipsoid(cF(1:2),diag(DF(1:2,1:2)),eye(2));
[YF2,ZF2]=Create2DEllipsoid(cF(2:3),diag(DF(2:3,2:3)),eye(2));

%----------------------------Approximation parameters----------------------
h1=0.1; %approximation parameter
n=3; %dimension of the continuous-state space 
N=30000; % number of particles (dimensionless)
rho=(1/(n*d^2));
pmax=0.9;
alpha=0.2;
%lamdas
lambdarelease=1/res_time; %propensity of release of a bound molecule
lambdabind=lambdarelease*Fimm/(1-Fimm); %binding propensity of an unbound molecule
%maximum binding/unbinding propensity: needed in the simulation to
%determine the rate at which events happen
lambdamax=max([lambdarelease,lambdabind]); %(>= max_q max_x SwitchingPropensity(x))
h2=sqrt(-(log(1-pmax))/(lambdamax*rho));
h3=sqrt((alpha*Dtobs)/rho);
h=min([h1 h2 h3]);
%average time interval
EDeltat=rho*h^2;
%Efficiency of bleaching
bleacheff=100000; 
%j_k propensity (for a given particle, maximum probability that a binding/unbinding event happens in a time interval EDeltat)
pevent=1-exp(-lambdamax*EDeltat); 

%--------probabilities of moves towards immediate neighbors (approximation of diffusion on a grid)
% (x,x+-[1;0;0],x+-[0;1;0],x+-[0;0;1]) 
moves=[0,0,0;1,0,0;-1,0,0;0,1,0;0,-1,0;0,0,1;0,0,-1]';
xsicoeff=zeros(1,n+1);
xsicoeff(1)=2/(rho*d^2)-2*n;
xsicoeff(2:(n+1))=0/d^2;    
ccoeff=1/(2*sum(cosh(h*xsicoeff(2:(n+1))))+xsicoeff(1));
pmoves=zeros(1,(2*n+1));
pmoves(1)=xsicoeff(1);
pmoves(2)=exp(h*xsicoeff(2));
pmoves(3)=exp(-h*xsicoeff(2));
pmoves(4)=exp(h*xsicoeff(3));
pmoves(5)=exp(-h*xsicoeff(3));
pmoves(6)=exp(h*xsicoeff(4));
pmoves(7)=exp(-h*xsicoeff(4));

pmoves=pmoves*ccoeff;

%canonical reflection cones (defined here for speed but used in subroutine ReflectionScheme below)
global ref_A1 ref_A2 ref_A3 ref_A4 ref_A5 ref_A6 ref_A7 ref_A8
ref_A1=[1,1,-1; 1,0,-1;0,1,-1]';
ref_A2=[0,0,-1; 1,0,-1;0,1,-1]';
ref_A3=[1,-1,-1; 1,0,-1;0,-1,-1]';
ref_A4=[0,0,-1; 1,0,-1;0,-1,-1]';
ref_A5=[-1,-1,-1; 0,-1,-1;-1,0,-1]';
ref_A6=[0,0,-1; 0,-1,-1;-1,0,-1]';
ref_A7=[-1,1,-1;-1,0,-1;0,1,-1]';
ref_A8=[0,0,-1;-1,0,-1;0,1,-1]';

%-----------------------------Simulation code------------------------------
trec=zeros(1,nobs);
qrec=zeros(N,nobs);
xrec=zeros(N*n,nobs);
brec=zeros(N,nobs);

%initialization - first recording point
t=0;
x=zeros(N*n,1);
q=zeros(N,1);
isbleached=zeros(N,1);
xnew=x; %for allocation only
qnew=q;
ntried=0;
nplaced=0;

%-------------- Random particle positioning within boundaries--------------
disp('Random initial positioning...');
while nplaced<N,
    ntried=ntried+1;
    xdraw=UG*DG*(2*rand(3,1)-1)+cG; %random position of a particle
    if InG(xdraw,cG,DG2inv,UGinv), %is it in the nucleus?
        x(nplaced*n+[1,2,3])=xdraw;
        q(nplaced+1)=randsrc(1,1,[0,1;[lambdarelease,lambdabind]/(lambdabind+lambdarelease)]);
        nplaced=nplaced+1;
    end
end
disp('Start simulation...');
%saving initial state if part of observation times
if tobs(1)==0,   
    iobs=1;    
    trec(1)=t;
    xrec(:,1)=x;
    qrec(:,1)=q;
    brec(:,1)=isbleached;       
else
    iobs=0;
end
nevents=0;
neventstot=0;
refcnt=0;
refcnttot=0;
nswcycles=0;
tic

%repeat discrete-time simulation until latest observation time is reached
while iobs<nobs    
%random time step
Deltat_k=rande(1/EDeltat,1);     
%--------simulate FRAP-----------------------
    if (t>=bleachint(1)) && (t<bleachint(2))
        randbleach=rand(1,N);
        for np=1:N
            if InG(x((np-1)*n+[1,2,3]),cF,DF2inv,UFinv) %if in bleaching region
                isbleached(np)=isbleached(np) || (randbleach(np)<Deltat_k*bleacheff/(1+Deltat_k*bleacheff));
            end
        end
    end
%--------------------------------------------------
    % move time ahead
    tnew=t+Deltat_k;   
    % random numbers for simulation of binding
    j_k_all=randsrc(N,1,[0,1;(1-pevent),pevent]);
    %indices of particles that may undergo a switching (binding/unbinding) event
    swidx=find(j_k_all==1);       
    mvidx=find(q==0); %find indices of particles that are not immobilized
    mvnum=length(mvidx);
    move_next_all=randsrc(1,mvnum,[1:size(moves,2);pmoves]);
    ref_all=rand(1,mvnum);
    qnew=q;
    xnew=x; %unbound particles will be moved next
    
        for i=1:mvnum %for all mobile molecules
            npidx=n*(mvidx(i)-1)+[1:3];          
            move_next=move_next_all(i);
            xnew(npidx)=x(npidx)+moves(:,move_next)*h;          
            %check for reflection       
            z=UGinv*xnew(npidx)-cG;         
            %If the particle landed out of boundaries, send it back in by a
            %randomized technique
            if diagDG2inv*z.^2>1%z'*DG2inv*z>1,%~InG(xnew(npidx),cG,DG2inv,UGinv),
            [Y,refprobs]=ReflectionScheme(xnew(npidx),cG,diagDG2inv,UGinv,h);
            refmove=find(ref_all(i)<cumsum(refprobs));% randsrc(1,1,[1:size(Y,2);refprobs']);
            y=Y(:,refmove(1));
            xnew(npidx)=xnew(npidx)+y;
            refcnt=refcnt+1;
            refcnttot=refcnttot+1;            
            end
        end

    %if a switch may have occurred, then simulate exactly what happened
    if ~isempty(swidx)       
        nswcycles=nswcycles+1;                    
        swnum=length(swidx);
        qk_all=rand(1,swnum);       
        for i=1:swnum           
            np=swidx(i);
            npidx=n*(np-1)+[1:3];
            Lambda=SwitchingPropensity(x(npidx),q(np),cG,diagDG2inv,UGinv,lambdabind,lambdarelease)/lambdamax;                    
            Iq=zeros(1,2);
            Iq(q(np)+1)=1;          
            %simulate random events using random numbers generated above
            qnew(np)=(qk_all(i)<(Iq(2)+Lambda(2)));            
            %for debugging/checking purposes, count number of events
            if qnew(np)~=q(np), nevents=nevents+1; neventstot=neventstot+1; end          
        end         
    end
     
    %update state for next iteration
    x=xnew;
    q=qnew;    
    told=t; %used for visualization only (below)
    t=tnew;
    
    %if an observation time point is reached, record the system state
    if t>tobs(iobs+1)       
        toc
        iobs=iobs+1;       
        iobs       
        trec(iobs)=t;
        xrec(:,iobs)=x;
        qrec(:,iobs)=q;
        brec(:,iobs)=isbleached;
        nevents=0;
        refcnt=0;
        nswcycles=0;     
        %display the particles state in graphical form, if requested
        if (plotmode==1)
            showcell(YG1,XG1,YG2,ZG2,YF1,XF1,YF2,ZF2,q,x,isbleached,(told>=bleachint(1)) && (told<bleachint(2)));
        end     
        tic      
    end       
end
toc

%------------------------ end of simulation-------------------------
%Additional functions needed

%random generations
function [sample]=rande(lambda,niter)
sample=-log(rand(1,niter))/lambda;

%functions for binding domain and probability
function [Lambda]=SwitchingPropensity(x,q,cU,diagDU2inv,UUinv,lambdabind,lambdarelease)
%return the appropriate switching propensity for the given particle
%position and the specified binding domain
%(ellipsoidal binding subdomain U)
z=UUinv*x-cU;
if diagDU2inv*z.^2<=1
    Lambda=[-lambdabind,lambdabind;lambdarelease,-lambdarelease];
else
    Lambda=zeros(2);
end
Lambda=Lambda(q+1,:);


function [isin]=InG(x,c,D2inv,Uinv)
%evaluates if x belongs to G
z=Uinv*x-c;
isin=z'*D2inv*z<=1;

%returns candidate reflection moves and corresponding probabilities for a
%particle out of boudary. See comments next to its call above.
function [Y,refprobs]=ReflectionScheme(x,c,diagD2inv,Uinv,h)
refdir=-((Uinv*x-c)'.*diagD2inv*Uinv)';
rx=refdir(1);
ry=refdir(2);
rz=refdir(3);
absrx=abs(rx);
absry=abs(ry);
absrz=abs(rz);
if (rz<0) && (absrx<absrz) && (absry<absrz) %z=-h quadrant
    ref_A=RefConeSel(rx,ry,rz);
elseif (rz>0) && (absrx<absrz) && (absry<absrz) %z=h quadrant
    ref_A=RefConeSel(rx,-ry,-rz);
    ref_A(2:3,:)=-ref_A(2:3,:);
elseif (ry<0) && (absrx<absry) && (absrz<absry) %y=-h quadrant
    ref_A=RefConeSel(rx,-rz,ry);
    ref_A(2:3,:)=[ref_A(3,:);-ref_A(2,:)];
elseif (ry>0) && (absrx<absry) && (absrz<absry) %y=h quadrant
    ref_A=RefConeSel(rx,rz,-ry);
    ref_A(2:3,:)=[-ref_A(3,:);ref_A(2,:)];
elseif (rx<0) && (absry<absrx) && (absrz<absrx) %x=-h quadrant
    ref_A=RefConeSel(-rz,ry,rx);
    ref_A([1,3],:)=[ref_A(3,:);-ref_A(1,:)];
else %(rx>0) && (absry<absrx) && (absrz<absrx), %x=h quadrant
    ref_A=RefConeSel(rz,ry,-rx);
    ref_A([1,3],:)=[-ref_A(3,:);ref_A(1,:)];
end
Y=ref_A*h;
tmp1=(Uinv*(x+Y(:,1))-c).^2;
tmp2=(Uinv*(x+Y(:,2))-c).^2;
tmp3=(Uinv*(x+Y(:,3))-c).^2;

%if short moves are not enough for an appropriate simulation of reflection
%(direction orthogonal to the boundary is not in the convex hull of the canonical moves toward the interior, 
%because at least one of these is out of boundary) then consider larger
%moves
if (diagD2inv*tmp1>=1) || (diagD2inv*tmp2>=1) || (diagD2inv*tmp3>=1)
    Y=2*Y;
end
refprobs=Y\refdir;
refprobs=refprobs/sum(refprobs);

%If larger moves are not enough, give error. This never happened, but I
%never proved mathematically that Cube 2 is enough...
if ~isempty(find((refprobs<0) | (refprobs==NaN) | (abs(refprobs)==Inf) ))
    disp('Error!')
    pause;
end

function [ref_A]=RefConeSel(rx,ry,rz)
%choses the right reflection moves depending on particle position relative
%to a "centered and normalized" diffusion domain
global ref_A1 ref_A2 ref_A3 ref_A4 ref_A5 ref_A6 ref_A7 ref_A8
xP=-rx/rz; yP=-ry/rz; xPyP=xP+yP; yP_xP=yP-xP;
if xPyP>=1 %region 1
    ref_A=ref_A1;
elseif (xPyP<=1) && (xP>=0) && (yP>=0) %region 2
    ref_A=ref_A2;
elseif yP-xP<=-1 %region 3
    ref_A=ref_A3;
elseif (yP-xP>=-1) && (xP>=0) && (yP<=0) %region 4
    ref_A=ref_A4;
elseif xPyP<=-1 %region 5
    ref_A=ref_A5;
elseif (xPyP>=-1) && (xP<0) && (yP<0) %region 6
    ref_A=ref_A6;
elseif yP_xP>=1 %region 7
    ref_A=ref_A7;
else %(yP_xP<=1) && (xP<0) && (yP>0), %region 8
    ref_A=ref_A8;
end



