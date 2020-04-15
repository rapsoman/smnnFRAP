function [X,Y]=Create2DEllipsoid(c,D,U)

Npts=100;
phi=[-Npts:(Npts)]/Npts*pi;

X=zeros(1,2*Npts);
Y=zeros(1,2*Npts);


for ph=1:2*Npts

        pt=zeros(2,1);
        pt(1)=D(1)*cos(phi(ph));
        pt(2)=D(2)*sin(phi(ph));

        pt=U*pt+c;
        
        X(ph)=pt(1);
        Y(ph)=pt(2);

end

X=X(:);
Y=Y(:);
