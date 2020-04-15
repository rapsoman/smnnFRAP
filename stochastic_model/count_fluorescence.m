
function [RFI,Ifocus,Ifocus0,Itot,Itot0]=count_fluorescence(trec,xrec,qrec,brec)
%function that produces a FRAP curve from the simulation output

%ellipsoidal domain G: (inv(UG)x-cG)'inv(DG^2)(inv(UG)x-cG)=1
plotmode=1;
DG=diag([5,4,4]); %radii
UG=eye(3); %rotation
cG=[0;0;0]; %center
DG2inv=inv(DG^2);
diagDG2inv=diag(DG2inv)';
UGinv=inv(UG);
[XG1,YG1]=Create2DEllipsoid(cG(1:2),diag(DG(1:2,1:2)),eye(2));
[YG2,ZG2]=Create2DEllipsoid(cG(2:3),diag(DG(2:3,2:3)),eye(2));


%photobleaching geometry (spherical spot for now)
%MAKE SURE THAT THE GEOMETRIES (WHERE YOU SIMULATE FRAP AND WHERE YOU COUNT FLUORESCENCE) MATCH!

%bleaching shape
DF=diag([1,1,1]);
UF=eye(3);
cF=[4;0;0];
DF2inv=inv(DF^2);
UFinv=inv(UF);
[XF1,YF1]=Create2DEllipsoid(cF(1:2),diag(DF(1:2,1:2)),eye(2));
[YF2,ZF2]=Create2DEllipsoid(cF(2:3),diag(DF(2:3,2:3)),eye(2));

DN=diag([1,1,1]); %radii
UN=eye(3); %rotation
cN=[4;0;0]; %center
DN2inv=inv(DN^2);
diagDN2inv=diag(DN2inv)';
UNinv=inv(UN);
[XN1,YN1]=Create2DEllipsoid(cN(1:2),diag(DN(1:2,1:2)),eye(2));
[YN2,ZN2]=Create2DEllipsoid(cN(2:3),diag(DN(2:3,2:3)),eye(2));


[N,nobs]=size(qrec);


%identification of prebleach images
n0=find(sum(brec)==0);
n0=n0(end);
Itot=zeros(1,nobs);
Ifocus=zeros(1,nobs);

%fluorescence intensity (number of fluorescent particles)
Itot=sum(1-brec);
%fluorescence intensity (count of fluorescent particles) in F (bleaching region)
for iobs=1:nobs
    for np=1:N
        xnp=xrec(3*(np-1)+[1,2,3],iobs);
        Ifocus(iobs)=Ifocus(iobs)+InG(xnp,cF,DF2inv,UFinv)*(1-brec(np,iobs));
    end
end

%average (over time) of prebleach global fluorescence intensity 
Itot0=mean(Itot(1:n0));
%average (over time) of prebleach fluorescence intensity inside
%the bleaching region
Ifocus0=mean(Ifocus(1:n0));
%Relative fluorescence intensity in the bleaching region
RFI=(Ifocus./Itot)/(Ifocus0/Itot0);

if plotmode==1
    for iobs=1:nobs
        iobs
        disp(num2str(trec(iobs)));
        %showcellposter(YG1,XG1,YG2,ZG2,YF1,XF1,YF2,ZF2,qrec(:,iobs),xrec(:,iobs),brec(:,iobs),0);
        showcellposter(YG1,XG1,YG2,ZG2,YN1,XN1,YN2,ZN2,YF1,XF1,YF2,ZF2,qrec(:,iobs),xrec(:,iobs),brec(:,iobs),0);

        Fmov(iobs)=getframe(gcf);
    end
end


function [isin]=InG(x,c,D2inv,Uinv)
%evaluates if x belongs to G

z=Uinv*x-c;
isin=z'*D2inv*z<=1;




