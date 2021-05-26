    
    % This code produces plots of distance from ceiling vs stress in cables
    % AB, AC, and AD

    %variable declaring
    sAB=1;
    sAC=1;
    sAD=1;
    
    % for loop counter variable
    a=1;

    % Gravity (m/s^2)
    g= 9.81;
    
    % Masses (kg)
    m=[1000 2000];
    
    % Known Dimensions (m)
    x1= 3.2;
    x2= 2.7;
    y1= 4.0;
    y2= 3.6;
    z= linspace(0.1,6,60); % +.1 increment
   
    % Diameter (mm)
    d= linspace(10,30,3); %10,20,30 
                            
    % Saftey Factor (=UTS/Tensile Stress)
    SF= 3;
    
    % Ultimate Tensile Strength (Al 2024) (Pa)
    UTS= 469000000;    
    
    % Max Allowable Stress (sAB, sAC, sAD need to be less than this)
    s= UTS/SF;
    
    for n=1:length(m) %loop for each Masses
        
        m1=m(n);
        
        for k=1:length(d) %loop for diameters 10,20,30

            d1= d(k);

            %cross-sectional area
            Ac= (pi.*(d1.^2))./4;
            Acm= Ac./(1000^2);       % (m^2)

            for i=1:length(z) %loop for plotting z v. d
                %cable lengths
                LAC= sqrt((x1^2)+(z(i).^2));
                LAB= sqrt((x2^2)+(y2^2)+(z(i).^2));
                LAD= sqrt((x2^2)+(y1^2)+(z(i).^2));

                %tensions
                A= [(-x2)./LAB x1./LAC (-x2)./LAD; (-y2)./LAB 0 y1./LAD; (-z(i))./LAB (-z(i))./LAC (-z(i))./LAD];
                B= [0; 0; -m1*g];
                T= A\B;
                TAB= T(1,1);
                TAC= T(2,1);
                TAD= T(3,1);

                %stresses
                sAB(i)=TAB./Acm;
                sAC(i)=TAC./Acm;
                sAD(i)=TAD./Acm;
                
                %safety factor
                SF1(i)= UTS./sAC(i);
                
            end
            
            % plot z v. stress
            figure(a);
            grid on
            title("Mass = " + m1 + "kg" + newline + "d = " + d1 + "mm")
            hold on
            yline(s,'-.b');
            pAB=plot(z,sAB);
            pAC=plot(z,sAC);
            pAD=plot(z,sAD);
            xlabel('z distance from ceiling(m)') ;
            ylabel('Tensile Stress (Pa)') ;
            legend('Max Allowable Stress for SF=3','Cable AB','Cable AC','Cable AD');
            hold off
            
            % plot z v. SF
            figure(7);
            grid on
            hold on
            pSF=plot(z,SF1);
            title("z vs SF")
            xlabel('z distance from ceiling(m)') ;
            ylabel('Safety Factor') ;
            hold off
            legend('m=1000kg d=10mm','m=1000kg d=20mm','m=1000kg d=30mm','m=2000kg d=10mm','m=2000kg d=20mm','m=2000kg d=30mm');
            a=a+1;
        end

    end
    