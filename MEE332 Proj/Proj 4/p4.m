%parameters

W=50; %lb
r=6; %in
b=.75; % also height of trapezoid
l=.312; % also outer width of trapezoid 
d=.312;
F=W./2;
a=1.62; %trapezoid inner width
N=5.*10^5;

%% rectangular stress
ro=r+l;
rc=r+l./2;
A=l.*b;
M=F.*rc;
rn=(ro-r)./(log(ro./r));
e=rc-rn;
ci=rn-r;
co=ro-rn;
sir=(M./(e.*A)).*(ci./r)+(F./A); %psi
sor=(-M./(e.*A)).*(co./ro)+(F./A); %psi

%% circular stress
ro=r+d;
rc=r+(d./2);
A=(pi./4).*(d.^2);
M=F.*rc;
rn=((d./2).^2)./(2.*(rc-sqrt(rc.^2-(d./2).^2)));
e=rc-rn;
ci=rn-r;
co=ro-rn;
sic=(M./(e.*A)).*(ci./r)+(F./A); %psi
soc=(-M./(e.*A)).*(co./ro)+(F./A); %psi

%% trapezoid stress
h=b;
ro=r+h;
rc=r+((h./3).*((2.*l+a)./(a+l)));
A=(h.*(a+l))./2;
dAr1=(a+r.*((a-l)./(ro-r))).*log(ro)-(((a-l)./(ro-r)).*ro);
dAr2=(a+r.*((a-l)./(ro-r))).*log(r)-(((a-l)./(ro-r)).*r);
dAr=dAr1-dAr2;
rn=A./dAr;
e=rc-rn;
ci=rn-r;
co=ro-rn;
sit=(M./(e.*A)).*(ci./r)+(F./A); %psi
sot=(-M./(e.*A)).*(co./ro)+(F./A); %psi

%% sort stresses into array

sr=max(abs(sir),abs(sor));
st=max(abs(sit),abs(sot));
sc=max(abs(sic),abs(soc)); 
s=[sr sc st];
      
%% 1020 Stl hot rolled
Sy=30.*(10.^3);
Sut=55.*(10.^3);
Sep=.5.*Sut;

Cload=1; %bending

A95r=.05.*b.*l; %rectangular
A95c=.010462.*(d.^2); %circular
A95t=.05.*((h.*(a+l))./2); %trapezoidal
A95=[A95r A95c A95t];

for i=1:3
    
    smin=0;
    sa=(s(i)-smin)./2; %alternating
    sm=(s(i)+smin)./2; %mean
    
    deq=sqrt(A95(i)./.0766);

    if deq <= .3
        Csize=1;
    elseif deq <= 10
        Csize=.869.*(deq.^(-.097));
    else
        Csize=.6;
    end


    A1=38545; %forged
    bb=-.995;
    Csurf=A1.*(Sut).^bb;

    if Csurf >= 1
        Csurf=1;
    end

    Ctemp=1;
    Creliab=.702;
    Se=Sep.*Cload.*Csize.*Csurf.*Ctemp.*Creliab;

    %cycle
    Sm=.9.*Sut;
    z=log10(1000)-log10(N);
    b1=(1./z).*log10(Sm./Se);
    a1=10.^(log10(Sm)-(3.*b1));
    SN=a1.*(N.^b1);

    %Fatigue Failure
    NfFs(i)=(SN.*Sut)./(sa.*Sut+sm.*SN); %sa/sm constant case 3

    %Static Failure
    %DE
    NfSs(i)=Sy./s(i);
end

%% Class 40 Cast Iron
Suc=140.*(10.^3);
Sut=42.*(10.^3);
Sfp=.4.*Sut;

Cload=1; %bending

A95r=.05.*b.*l; %rectangular
A95c=.010462.*(d.^2); %circular
A95t=.05.*((h.*(a+l))./2); %trapezoidal
A95=[A95r A95c A95t];

for i=1:3
    
    smin=0;
    sa=(s(i)-smin)./2; %alternating
    sm=(s(i)+smin)./2; %mean
    
    deq=sqrt(A95(i)./.0766);

    if deq <= .3
        Csize=1;
    elseif deq <= 10
        Csize=.869.*(deq.^(-.097));
    else
        Csize=.6;
    end


    Csurf=1; %cast iron

    Ctemp=1;
    Creliab=.702;
    Sf=Sfp.*Cload.*Csize.*Csurf.*Ctemp.*Creliab;

    %cycle
    Sm=.9.*Sut;
    z=log10(1000)-log10(N);
    b1=(1./z).*log10(Sm./Sf);
    a1=10.^(log10(Sm)-(3.*b1));
    SN=a1.*(N.^b1);

    %Fatigue Failure
    NfFc(i)=(SN.*Sut)./(sa.*Sut+sm.*SN); %sa/sm constant case 3

    %Static Failure
    %principal stresses
    s1=(s(i)./2)+sqrt((s(i)./2).^2 + 0); %shear = 0
    s3=(s(i)./2)-sqrt((s(i)./2).^2 + 0);
    %both s1, s3 positive so zone 1 of Modified Mohr
    NfSc(i)=Sut./s1;
end

%% 2024 Al
Sy=11.*(10.^3);
Sut=26.*(10.^3);
Sfp=.4.*Sut;

Cload=1; %bending

A95r=.05.*b.*l; %rectangular
A95c=.010462.*(d.^2); %circular
A95t=.05.*((h.*(a+l))./2); %trapezoidal
A95=[A95r A95c A95t];

for i=1:3
    
    smin=0;
    sa=(s(i)-smin)./2; %alternating
    sm=(s(i)+smin)./2; %mean
    
    deq=sqrt(A95(i)./.0766);

    if deq <= .3
        Csize=1;
    elseif deq <= 10
        Csize=.869.*(deq.^(-.097));
    else
        Csize=.6;
    end


    A1=38545; %forged
    bb=-.995;
    Csurf=A1.*(Sut).^bb;

    if Csurf >= 1
        Csurf=1;
    end

    Ctemp=1;
    Creliab=.702;
    Sf=Sfp.*Cload.*Csize.*Csurf.*Ctemp.*Creliab;

    %cycle
    Sm=.9.*Sut;
    z=log10(1000)-log10(N);
    b1=(1./z).*log10(Sm./Sf);
    a1=10.^(log10(Sm)-(3.*b1));
    SN=a1.*(N.^b1);

    %Fatigue Failure
    NfFa(i)=(SN.*Sut)./(sa.*Sut+sm.*SN); %sa/sm constant case 3

    %Static Failure 
    %DE
    NfSa(i)=Sy./s(i); 
end


%% print result
Name={'rectangular';'circular';'trapezoidal'};
FatigueSF=NfFs.';
StaticSF=NfSs.';
Steel=table(Name,FatigueSF,StaticSF)
FatigueSF=NfFc.';
StaticSF=NfSc.';
Castiron=table(Name,FatigueSF,StaticSF)
FatigueSF=NfFa.';
StaticSF=NfSa.';
Aluminium=table(Name,FatigueSF,StaticSF)

