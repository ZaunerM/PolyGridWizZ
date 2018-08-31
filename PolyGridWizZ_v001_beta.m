%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  PolyGridWizZ (Polynomial Gridgeneration Wizzard)  %%%
%%%               (Beta-version 0.0.1               )  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Developed by Markus Zauner and Neil Sandham     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           UNIVERSITY OF SOUTHAMPTON                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Corresponding author:      Markus Zauner        %%%
%%%                           m.zauner@soton.ac.uk     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  developed for Matlab 8.6.0.267246 (R2015b)        %%%
%%%     tested for Matlab 9.4.0.813654 (R2018a)        %%%
%% INTRODUCTION
%
% Supporting comments by Markus Zauner 
% ("Please excuse my English and the messy coding - I am not a software engineer ;)")
%
% This code is a beta version for airfoils with a blunt and sharp trailing edge
% for Matlab 8.6.0.267246 (R2015b)
%
% The version 1.0.0 is supposed to contain following features:
% - parallelised 2D grid generation
% - optimised modular structure
% - open-source version written in Python (www.Jupyter.org)
% - gridline-intersection check
% _________________________________________________________________________
% Please cite this code referring to: 
% Proceedings of the 10th International Conference on Computational Fluid Dynamics (ICCFD10)
% "Multiblock structured grids for direct numerical simulations of transonic wing sections"
% Zauner, M. and Sandham, N. D., 2018
% Paper ID: ICCFD10-2018-129
% _________________________________________________________________________
%% Guide:
%-1. It might be useful to use version control (e.g. Bitbucket) to track changes of your grid
% 0. First at all, we recommend to fold all sections (CTRL+=)
% 1. Set parameters in P1, P2 & P3 
%    In the case of NACA profiles, the four digits are specified in P1
% 2. Check inputs in I1 and correct them by commenting A(i) out
%    After the final input file was safed before generating the full grid, 
%    A(i) can uncommented again.
% 3. In case you add new parameters you can either add them to I2 or extend
% the input file in I1 and F1 -> see sec I1 for more details
%    I2 also allows to overrule parameters in case you want to make sure
%    they won't change. In that example the eta-spacing at the wall and outer 
%    bound are hardcoded in I2
% 4. Run now every single section separately and check results carefully.
%    It might be necessary to do several iterations to optimise the grid
%    --> That's the price you pay for a fully flexible grid generation.
%        Therefore you should think carefully what you want your grid to
%        look like and what resolution you aim for.
%
%    If you change the resolution or other parameters that are coupled to
%    sections that were executed before, you have to rerun the code. In
%    that case, (1) save the input file by executing the corresponding section 
%    'S1' and (2) put a 'stop' command after the recent section and (3) 
%    run the whole code (F5). When you reach the 'stop', you can continue.
% 5. Contour 1-4 are the segments on the airfoil surface.
%    Functions Poly6 and Poly6_end are applied to calculate spacing.
%      The third derivative is a control parameter and applied at the first point
%      or the last point in Poly6_end and Poly6, respectively.
%    Initially, you set draft='T' and find you best setting by changing
%      parameters within the corresponding if-condition
%    Later, when you have written out the final input-file, you can change
%      it back to draft='F'. 
%    Please mind that each segment (C1-C4) is coupled with others. It might
%      take you some iterations to find the best settings for C1-C4
% 6. The eta-gridlines are designed at representative locations either by
%    one polynomial (Style 1&2) ord by more polynomials (Style 3). Please
%    find an example for the design process in the corresponding paper.
% 7. Before generating the full grid, test the C-block in section T1
% 8. Generate Block 2 by executing E1
% 9. Then execute 
%      - E2a for a blunt trailing edge
%      - E2b for a sharp trailing edge
%10. For checking, copy all blocks from the created directory to the
%      directory of the matlab script and delete the first line of the .dat-files
%    Then run the corresponding Testing section at the end of the script
%% Reset everything & calc Carpenter coefficients
clearvars
close all
global bc4
bc4 = ones(4,6);
CarpCoeff
%% P1: Read-in settings          (User input required)

readfiles=false       % true  -> Read in ascii-files of airfoil
                      % false -> Airfoil is defined by function
                        
%Don't forget to set sharp=true in case of a sharp TE profile
sharp=true            % true  -> sharp trailing edge
                      % false -> blunt trailing edge
                        
% READ-IN SETTINGS
    Contour1='Contour1.dat' % Upper surface rear part
    Contour2='Contour2.dat' % Lower surface rear part
    Contour3='Contour3.dat' % Upper surface front part
    Contour4='Contour4.dat' % Lower surface front part
% FUNCTION SETTINGS (example cut-off NACA airfoils)
    TE_thickness = 0.005 % Trailing Edge Thickness (Profile gets cut-off)
    if sharp==true
        TE_thickness = 0
    end
    m = 4/100;           % First digit
    p = 4/10;            % Second digit
    t=12/100;            % Third and fourth digits
    cx = 1.0;             % Chord lenght
    ncheck = 100000;     % Resolution of raw profile -----> Check later variable 'space'<<smallest cell
    sang = 0.0;          % Sweep angle (degrees) <-- NOT TESTED
    sang = sang/180*pi;
%   Profile function gets defined between line 32-96
%
% Run Calculation 
if readfiles
    C1 = load(Contour1); %upper TE
    C2 = load(Contour2); %lower TE
    C3 = load(Contour3); %upper LE
    C4 = load(Contour4); %lower LE
    'Contour read in'
end
if sharp==false && readfiles==false
    xck = linspace(0,1,ncheck);
    xuck = linspace(0,1,ncheck);
    yuck = linspace(0,1,ncheck);
    suck = linspace(0,1,ncheck);
    xlck = linspace(0,1,ncheck);
    ylck = linspace(0,1,ncheck);
    slck = linspace(0,1,ncheck);
    
    ytck = t/0.2*( 0.2969*sqrt(xck)+(-0.1260*xck)+(-0.3516*(xck).^2)+0.2843*(xck).^3+(-0.1015)*(xck).^4 );

    for i=2:ncheck
        if(xck(i) <= p)
            yc = m*xck(i)/(p^2)*(2*p-xck(i));
            dycdx = 2*m/(p^2)*(p-xck(i));
        else
            yc = m/((1-p)^2)*(1-2*p+2*p*xck(i)-xck(i).^2); %(c-xck(i))*(1+xck(i)/c-2*p);
            dycdx = 2*m/((1-p)^2)*(p-xck(i));
        end

        yuck(i) = yc+ytck(i)*cos(atan(dycdx));
        xuck(i) = xck(i)*cos(sang)-ytck(i)*sin(atan(dycdx));
        suck(i) = suck(i-1)+sqrt((yuck(i)-yuck(i-1))^2+(xuck(i)-xuck(i-1))^2);
        %
        ylck(i) = yc-ytck(i)*cos(atan(dycdx));
        xlck(i) = xck(i)*cos(sang)+ytck(i)*sin(atan(dycdx));
        slck(i) = slck(i-1)+sqrt((ylck(i)-ylck(i-1))^2+(xlck(i)-xlck(i-1))^2);
    end
    
    D=ytck-ylck;
    i=round(ncheck/2);
    while D(i)>=TE_thickness && i<=ncheck
        n_max=i;
        i=i+1;
    end
            
    [minxu1,minxu2]=min(xuck)
    [minxl1,minxl2]=min(xlck)
    if minxu1 < minxl1
        xu_fin=xuck(minxu2:end);
        yu_fin=yuck(minxu2:end);
        xl_fin=xuck(minxu2:-1:1);
        yl_fin=yuck(minxu2:-1:1);
        xl_fin(minxu2:minxu2+size(xlck,2)-1)=xlck(1:end);
        yl_fin(minxu2:minxu2+size(ylck,2)-1)=ylck(1:end);
    else
        'allocation of arrays still needs to be done accordingly'
        stop
    end

    i=1
    while xu_fin(i)<=0.5
        nu_break=i;
        i=i+1;
    end
    while xl_fin(i)<=0.5
        nl_break=i;
        i=i+1;
    end

    C1(:,1)=xu_fin(nu_break:end).*cx;
    C1(:,3)=yu_fin(nu_break:end).*cx;
    C2(:,1)=xl_fin(nl_break:end).*cx;
    C2(:,3)=yl_fin(nl_break:end).*cx;
    C3(:,1)=xu_fin(1:nu_break).*cx;
    C3(:,3)=yu_fin(1:nu_break).*cx;
    C4(:,1)=xl_fin(1:nl_break).*cx;
    C4(:,3)=yl_fin(1:nl_break).*cx;
    
    figure
    plot(C1(:,1),C1(:,3),'r')
    hold on
    plot(C2(:,1),C2(:,3),'g')
    plot(C3(:,1),C3(:,3),'b')
    plot(C4(:,1),C4(:,3),'k')
    daspect([1 1 1])
end
if sharp==true && readfiles==false
    xck = linspace(0,1,ncheck);
    xuck = linspace(0,1,ncheck);
    yuck = linspace(0,1,ncheck);
    suck = linspace(0,1,ncheck);
    xlck = linspace(0,1,ncheck);
    ylck = linspace(0,1,ncheck);
    slck = linspace(0,1,ncheck);

    ytck = t/0.2*( 0.2969*sqrt(xck)+(-0.1260*xck)+(-0.3516*(xck).^2)+0.2843*(xck).^3+(-0.1036)*(xck).^4 );
    
    for i=2:ncheck
        if(xck(i) <= p)
            yc = m*xck(i)/(p^2)*(2*p-xck(i));
            dycdx = 2*m/(p^2)*(p-xck(i));
        else
            yc = m/((1-p)^2)*(1-2*p+2*p*xck(i)-xck(i).^2); %(c-xck(i))*(1+xck(i)/c-2*p);
            dycdx = 2*m/((1-p)^2)*(p-xck(i));
        end

        yuck(i) = yc+ytck(i)*cos(atan(dycdx));
        xuck(i) = xck(i)*cos(sang)-ytck(i)*sin(atan(dycdx));
        suck(i) = suck(i-1)+sqrt((yuck(i)-yuck(i-1))^2+(xuck(i)-xuck(i-1))^2);
        %
        ylck(i) = yc-ytck(i)*cos(atan(dycdx));
        xlck(i) = xck(i)*cos(sang)+ytck(i)*sin(atan(dycdx));
        slck(i) = slck(i-1)+sqrt((ylck(i)-ylck(i-1))^2+(xlck(i)-xlck(i-1))^2);
    end
    
    D=ytck-ylck;
    i=round(ncheck/2);
    while D(i)>=TE_thickness && i<=ncheck
        n_max=i;
        i=i+1;
    end
            
    [minxu1,minxu2]=min(xuck)
    [minxl1,minxl2]=min(xlck)
    if minxu1 < minxl1
        xu_fin=xuck(minxu2:end);
        yu_fin=yuck(minxu2:end);
        xl_fin=xuck(minxu2:-1:1);
        yl_fin=yuck(minxu2:-1:1);
        xl_fin(minxu2:minxu2+size(xlck,2)-1)=xlck(1:end);
        yl_fin(minxu2:minxu2+size(ylck,2)-1)=ylck(1:end);
    else
        'allocation of arrays still needs to be done accordingly'
        stop
    end

    i=1
    while xu_fin(i)<=0.5
        nu_break=i;
        i=i+1;
    end
    while xl_fin(i)<=0.5
        nl_break=i;
        i=i+1;
    end

    C1(:,1)=xu_fin(nu_break:end).*cx;
    C1(:,3)=yu_fin(nu_break:end).*cx;
    C2(:,1)=xl_fin(nl_break:end).*cx;
    C2(:,3)=yl_fin(nl_break:end).*cx;
    C3(:,1)=xu_fin(1:nu_break).*cx;
    C3(:,3)=yu_fin(1:nu_break).*cx;
    C4(:,1)=xl_fin(1:nl_break).*cx;
    C4(:,3)=yl_fin(1:nl_break).*cx;
    
    figure
    plot(C1(:,1),C1(:,3),'r')
    hold on
    plot(C2(:,1),C2(:,3),'g')
    plot(C3(:,1),C3(:,3),'b')
    plot(C4(:,1),C4(:,3),'k')
    daspect([1 1 1])
end
'Raw contour generated'
a_rad_offset=atan((C1(end,1)-C2(end,1))/(C1(end,3)-C2(end,3)))
%
% Redistribute points along surface
CC=flipud(C2);
CC(size(C2,1):size(C2,1)+size(C4,1)-1,:)=flipud(C4);
CC(size(C2,1)+size(C4,1)-1:size(C2,1)+size(C4,1)-1+size(C3,1)-1,:)=C3;
CC(size(C2,1)+size(C4,1)-1+size(C3,1)-1:size(C2,1)+size(C4,1)-1+size(C3,1)-1+size(C1,1)-1,:)=C1;
%
sCC   = calc_s(CC(:,1),CC(:,3));
space = sCC(1,end)/size(CC,1)
sCCnew= linspace(0,sCC(end),size(CC,1));
CCnew = Interp_onto(size(CC,1),CC,sCC',sCCnew,'t');
%
C2=flipud(CC(1:size(C2,1),:));
C4=flipud(CC(size(C2,1):size(C2,1)+size(C4,1)-1,:));
C3=CC(size(C2,1)+size(C4,1)-1:size(C2,1)+size(C4,1)-1+size(C3,1)-1,:);
C1=CC(size(C2,1)+size(C4,1)-1+size(C3,1)-1:size(C2,1)+size(C4,1)-1+size(C3,1)-1+size(C1,1)-1,:);
clearvars CC CCnew sCCnew sCC

'Equidistant resolution --->',space
%% P2: Setting angle of attack   (User settings required)

% Rotate Profile
a_deg = 4.0;                  % Angle of attack (degrees)
xor   = 0.5;                  % X-Position of rotation axis
yor   = 0.0;                  % Y-Position of rotation axis
% Perform calculation
if a_deg ~= 0.0
    % Rotate Profile
    a_deg=a_deg/360*2*pi;
    C1_(:,1)=C1(:,1)-xor;
    C1_(:,2)=C1(:,3)-yor;
    C2_(:,1)=C2(:,1)-xor;
    C2_(:,2)=C2(:,3)-yor;
    C3_(:,1)=C3(:,1)-xor;
    C3_(:,2)=C3(:,3)-yor;
    C4_(:,1)=C4(:,1)-xor;
    C4_(:,2)=C4(:,3)-yor;
    M = [cos(a_deg), sin(a_deg); -sin(a_deg), cos(a_deg)];
    
    for i=1:size(C1_,1)
        C1_(i,:)=(M*C1_(i,:)')'+[xor, yor];
    end
    for i=1:size(C2_,1)
        C2_(i,:)=(M*C2_(i,:)')'+[xor, yor];
    end
    for i=1:size(C3_,1)
        C3_(i,:)=(M*C3_(i,:)')'+[xor, yor];
    end
    for i=1:size(C4_,1)
        C4_(i,:)=(M*C4_(i,:)')'+[xor, yor];
    end
    
    figure
    plot(C1_(:,1),C1_(:,2),'-r')
    daspect([1 1 1])
    hold on
    plot(C2_(:,1),C2_(:,2),'-r')
    plot(C3_(:,1),C3_(:,2),'-r')
    plot(C4_(:,1),C4_(:,2),'-r')
    
    plot(C1(:,1),C1(:,3),'-.b')
    daspect([1 1 1])
    hold on
    plot(C2(:,1),C2(:,3),'-.b')
    plot(C3(:,1),C3(:,3),'-.b')
    plot(C4(:,1),C4(:,3),'-.b')
    
    % Overwrite Values for C
    C1(:,1)=C1_(:,1);
    C1(:,3)=C1_(:,2);
    C2(:,1)=C2_(:,1);
    C2(:,3)=C2_(:,2);
    C3(:,1)=C3_(:,1);
    C3(:,3)=C3_(:,2);
    C4(:,1)=C4_(:,1);
    C4(:,3)=C4_(:,2);
    
    % Find new leading edge
    C_LE=flipud(C3);
    C_LE(size(C3,1):size(C3,1)+size(C4,1)-1,:)=(C4);
    dx_LE=deriv(C_LE(:,1),1);
    i_LE=0;
    i=1;
    while dx_LE(i)<=0
        i_LE=i;
        i=i+1;
    end    
    %i_LE=i_LE+1; %Correction
    
    C_LE3(1:i_LE,:)=C_LE(1:i_LE,:);
    C_LE4=C_LE(i_LE:size(C_LE,1),:);
    test=size(C_LE3,1)+size(C_LE4,1);
    clearvars C3 C4
    C3=flipud(C_LE3);
    C4=C_LE4;

    clearvars C1_ C2_ C3_ C4_ C_LE dx_LE C_LE3 C_LE4 y_LE
    
    % Find new Center
    C_up=C3;
    C_up(size(C3,1):size(C3,1)+size(C1,1)-1,:)=C1;
    i_cu=0;
    i=1;
    while C_up(i,1)<=0.5
        i_cu=i;
        i=i+1;
    end
    clearvars C3 C1
    C3(1:i_cu,:)=C_up(1:i_cu,:);
    C1=C_up(i_cu:size(C_up,1),:);
    
    yq=interp1(C_up(:,1),C_up(:,3),0.5);
    C3(end,1)=0.5;
    C3(end,3)=yq;
    C1(1,1)=0.5;
    C1(1,3)=yq;
    clearvars yq
    
    C_lo=C4;
    C_lo(size(C4,1):size(C4,1)+size(C2,1)-1,:)=(C2);
    i_cl=0;
    i=1;
    while C_lo(i,1)<=0.5
        i_cl=i;
        i=i+1;
    end
    clearvars C4 C2
    C4(1:i_cl,:)=C_lo(1:i_cl,:);
    C2=C_lo(i_cl:size(C_lo,1),:);
    
    yq=interp1(C_lo(:,1),C_lo(:,3),0.5);
    C4(end,1)=0.5;
    C4(end,3)=yq;
    C2(1,1)=0.5;
    C2(1,3)=yq;
    
    test=size(C1,1)+size(C2,1)+size(C3,1)+size(C4,1)

    clearvars C1_ C2_ C3_ C4_ C_LE dx_LE C_LE3 C_LE4 C_lo C_up
    
    % Put x-axis trough LE-Point
    y_LE=C3(1,3)
    C1(:,3)=C1(:,3)-y_LE;
    C2(:,3)=C2(:,3)-y_LE;
    C3(:,3)=C3(:,3)-y_LE;
    C4(:,3)=C4(:,3)-y_LE;
    xq=C3(1,1)
    
    plot(C1(:,1),C1(:,3),'-g')
    daspect([1 1 1])
    hold on
    plot(C2(:,1),C2(:,3),'-g')
    plot(C3(:,1),C3(:,3),'-g')
    plot(C4(:,1),C4(:,3),'-g')
end
y_LE=C3(1,3);
xq=C3(1,1);
'Contour rotated to adjust angle of attack'
%% P3: Basic domain settings
%
rad    = 7.3;           % Radius of half circle
x_circ = 0.5;           % X-Coordinate of center of half circle
y_circ = 0.0;           % Y-Coordinate of center of half circle
%
xOut     =5.5;          % X-position of outlet
%
xTEu     = C1(end,1);   % X-Position of upper TE
xTEt     = C1(end,1);   % X-Position of last point of Block 2 in J 
xTEl     = C2(end,1);   % X-Position of lower TE
xTElt    = C2(end,1);   % X-Position of first point of Block 2 in J
yTE1 =  C1(end,3);      % Y-Position of upper TE 
yTE2 =  C2(end,3);      % Y-Position of lower TE
yTE3 =  rad;            % Y-Position of last point of Block 2 in J 
yTE4 = -rad;            % Y-Position of first point of Block 2 in J
%
xcu    = C1(1,1);       % X-Position of upper center
ycu    = C1(1,3);       % Y-Position of upper center
xcl    = C2(1,1);       % X-Position of lower center
ycl    = C2(1,3);       % Y-Position of lower center
%
N_discr = round(5000*rad/7.5);
%% I1: Read input (Very important! -> Main input)
%
%In this section, all the parameters are listed. 
%By default, it reads them in from the Input.txt
%Before creating the full 2D grid, the grid generator updates this file and
%saves a copy in a folder named Grid@HH_MM_dd_mm_jj.
%
%It is not recommended to overrule parameters here....use the next section!
%
%The input_readable.txt file is just to map numbers to the variable-name,
%but effectively only the input.txt containing only values is read in.
%
fileID = fopen('./Input.txt','r');
formatSpec = '%e';
sizeA = [1 Inf];
A = fscanf(fileID,formatSpec,sizeA)
fclose(fileID);
%
dummy			=A(1) %TE
NTEw			=A(2)   % Resolution of trailing edge (TE) - in case of blunt trailing edge
dummy			=A(3) %C1 
Nc1             =A(4)   % Resolution Section 1 (up)
dscu			=A(5)   % Spacing in Xi-direction @ half-chord at upper side
dsTEu			=A(6)   % Spacing in Xi-direction @ TE at upper side
ddsTEu			=A(7)   % 1st deriv of spacing @ TE at upper side
ddscu			=A(8)   % 1st deriv of spacing @ half-chord at upper side
dddsTEu			=A(9)   % 2nd deriv of spacing @ TE at upper side
dddscu			=A(10)  % 2nd deriv of spacing @ half-chord at upper side
phi_center1		=A(11)  % Wall-angle at upper half-chord
dummy			=A(12) %C2
Nc2         	=A(13)   % Resolution Section 2 (low)
dscl			=A(14)   % Spacing in Xi-direction @ half-chord at lower side
dsTEl			=A(15)   % Spacing in Xi-direction @ TE at lower side
ddsTEl			=A(16)   % 1st deriv of spacing @ TE
ddscl			=A(17)   % 1st deriv of spacing @ half-chord
dddsTEl			=A(18)   % 2nd deriv of spacing @ TE
dddscl			=A(19)   % 2nd deriv of spacing @ half-chord
phi_center2		=A(20)   % Wall-angle at lower half-chord
dummy			=A(21) %C3
Nc3         	=A(22)   % Resolution Section 3 (up)
dsLE			=A(23)   % Spacing in Xi-direction @ LE
ddsLE			=A(24)   % 1st deriv of spacing @ LE
dddsLE			=A(25)   % 2nd deriv of spacing @ LE
dummy			=A(26) %C4
Nc4         	=A(27)   % Resolution Section 4 (low)
%dsLE			=A(28)   % Spacing in Xi-direction @ LE
%ddsLE			=A(29)   % 1st deriv of spacing @ LE
dddsLE			=A(30)   % 2nd deriv of spacing @ LE
dummy			=A(31) %Airfoil extension
Nbuffer			=A(32)   % Points in Xi used for transitioning from TE-continuation to horizontal grid line
Nw              =A(33)   % Total number of grid points in Xi in B1&B3
xTEc_save		=A(34)   % Approx. distance where transition to horizontal xi-gridlines should be finished (initial setting gets overwritten later)
StretchW		=A(35)   % Stretching factor that causes the upper and lower contiuation to diverge in order to get a more uniform resolution in the wake-region away from the airfoil.
                         % For sharp TE this factor defines the stretching of grid cells along B1/B3 interface in eta-dirextion
Nk              =A(36)   % xTE/Nk -> Fraction of contour-blending (to steer converging wall-gridlines in the wake)
dsTEc			=A(37)   % Xi-Spacing at the end of transient
ddsTEc			=A(38)   % 1st derivative of Xi-Spacing at the end of transient
dddsTEc			=A(39)   % 2nd derivative of Xi-Spacing at the end of transient
dummy			=A(40) %Outlet
dsOut			=A(41)   % Xi-Spacing at outlet
ddsOut			=A(42)   % 1st derivative of Xi-Spacing
dddsOut			=A(43)   % 2nd derivative of Xi-Spacing 
dummy			=A(44) %Upper top boundaries
dsTEut			=A(45)   % Xi-Spacing at TE
ddsTEut			=A(46)   % 1st derivative of Xi-Spacing
dddsTEut		=A(47)   % 2nd derivative of Xi-Spacing 
dscut			=A(48)   % Xi-Spacing at half-chord
ddscut			=A(49)   % 1st derivative of Xi-Spacing
dddscut			=A(50)   % 2nd derivative of Xi-Spacing
dummy			=A(51) %Lower lower boundaries
dsTElt			=A(52)   % Xi-Spacing at TE
ddsTElt			=A(53)   % 1st derivative of Xi-Spacing
dddsTElt		=A(54)   % 2nd derivative of Xi-Spacing
dsclt			=A(55)   % Xi-Spacing at half-chord
ddsclt			=A(56)   % 1st derivative of Xi-Spacing
dddsclt			=A(57)   % 2nd derivative of Xi-Spacing
dummy			=A(55) %Angles at the outer surface are approximated by polynomials (see TOP of Contour 3)
N3tphic			=A(59)   % Control Point for approximating angles
ddphi3t_LEc		=A(60)   % Second derivative of the angle at N3tphic
dLE_tan			=A(61)   % Spacing at leading edge top in xi-direction (tangential) 
ddLE_tan		=A(62)   % 1st derivative of spacing at leading edge top in xi-direction (tangential) 
dddLE_tan		=A(63)   % 2nd derivative of spacing at leading edge top in xi-direction (tangential) 
dummy			=A(64) %Angles at the outer surface are approximated by polynomials (see TOP of Contour 4)
N4tphic			=A(65)   % Control Point for approximating angles
ddphi4t_LEc		=A(66)   % Second derivative of the angle at N4tphic
dummy           	=A(67) %Eta gridlines
dummy               =A(68) %Define @ TE, half-chord and LE
Ny1                 =A(69)   % Points near wall
Ny2                 =A(70)   % Points far from wall
NyTOT               =A(71)   % Total # of points in eta (Ny1+Ny2-1)
control_dist1		=A(72)   % Approximate distance of the border of the near-wall section @ TE
control_dist2		=A(73)   % Approximate distance of the border of the near-wall section @ Half-chord
control_dist3		=A(74)   % Approximate distance of the border of the near-wall section @ LE
prof_deriv1_top1	=A(75)  % Eta-spacing at outer boundary @ TE
prof_deriv2_top1	=A(76)   % 1st derivative of eta-spacing at outer boundary
prof_deriv3_top1	=A(77)   % 2nd derivative of eta-spacing at outer boundary
prof_deriv1_wall1	=A(78)   % Eta-spacing at airfoil surface
prof_deriv2_wall1	=A(79)   % 1st derivative of eta-spacing at airfoil surface
prof_deriv3_wall1	=A(80)   % 2nd derivative of eta-spacing at airfoil surface
prof_deriv1_cont1	=A(81)   % Eta-spacing at interface
prof_deriv2_cont1	=A(82)   % 1st derivative of eta-spacing at interface
dummy           	=A(83)   %
prof_deriv1_top2	=A(84)  % Eta-spacing at outer boundary @ half chord
prof_deriv2_top2	=A(85)   % 1st derivative of eta-spacing at outer boundary
prof_deriv3_top2	=A(86)   % 2nd derivative of eta-spacing at outer boundary
prof_deriv1_wall2	=A(87)   % Eta-spacing at airfoil surface
prof_deriv2_wall2	=A(88)   % 1st derivative of eta-spacing at airfoil surface
prof_deriv3_wall2	=A(89)   % 2nd derivative of eta-spacing at airfoil surface
prof_deriv1_cont2	=A(90)   % Eta-spacing at interface
prof_deriv2_cont2	=A(91)   % 1st derivative of eta-spacing at interface
dummy           	=A(92)   %
prof_deriv1_top3	=A(93)  % Eta-spacing at outer boundary @ LE
prof_deriv2_top3	=A(94)   % 1st derivative of eta-spacing at outer boundary
prof_deriv3_top3	=A(95)   % 2nd derivative of eta-spacing at outer boundary
prof_deriv1_wall3	=A(96)   % Eta-spacing at airfoil surface
prof_deriv2_wall3	=A(97)   % 1st derivative of eta-spacing at airfoil surface
prof_deriv3_wall3	=A(98)   % 2nd derivative of eta-spacing at airfoil surface
prof_deriv1_cont3	=A(99)   % Eta-spacing at interface
prof_deriv2_cont3	=A(100)  % 1st derivative of eta-spacing at interface
dummy           	=A(101) %Additional Parameters
const_line      	=A(102)  % Force number of points being on wall-optimised line (Eta-gridlines)
const_line2     	=A(103)  % Force number of points being on top-optimised line (Eta-gridlines)
blend_C         	=A(104)  % Number of points left & right of LE which are uesed for a blending over horizontal gridlines (Eta-gridlines)
%% S1: Save input-file
safe_input='on'
if safe_input=='on'
    vars=[
        dummy %TE
        NTEw
        dummy %C1
        Nc1
        dscu
        dsTEu
        ddsTEu
        ddscu
        dddsTEu
        dddscu
        phi_center1
        dummy %C2
        Nc2
        dscl
        dsTEl
        ddsTEl
        ddscl
        dddsTEl
        dddscl
        phi_center2
        dummy %C3
        Nc3
        dsLE
        ddsLE
        dddsLE
        dummy %C4
        Nc4
        dsLE
        ddsLE
        dddsLE
        dummy %Aifoil extension
        Nbuffer
        Nw
        xTEc_save
        StretchW
        Nk
        dsTEc
        ddsTEc
        dddsTEc
        dummy %Outlet
        dsOut
        ddsOut
        dddsOut
        dummy %TC1
        dsTEut
        ddsTEut
        dddsTEut
        dscut
        ddscut
        dddscut
        dummy %TC2
        dsTElt
        ddsTElt
        dddsTElt
        dsclt
        ddsclt
        dddsclt
        dummy %TC3
        N3tphic
        ddphi3t_LEc
        dLE_tan
        ddLE_tan
        dddLE_tan
        dummy %TC4
        N4tphic
        ddphi4t_LEc
        dummy %Eta gridlines
        dummy %TE
        Ny1
        Ny2
        NyTOT
        control_dist1
        control_dist2
        control_dist3
        prof_deriv1_top1
        prof_deriv2_top1
        prof_deriv3_top1
        prof_deriv1_wall1
        prof_deriv2_wall1
        prof_deriv3_wall1
        prof_deriv1_cont1
        prof_deriv2_cont1
        dummy %center
        prof_deriv1_top2
        prof_deriv2_top2
        prof_deriv3_top2
        prof_deriv1_wall2
        prof_deriv2_wall2
        prof_deriv3_wall2
        prof_deriv1_cont2
        prof_deriv2_cont2
        dummy %LE
        prof_deriv1_top3
        prof_deriv2_top3
        prof_deriv3_top3
        prof_deriv1_wall3
        prof_deriv2_wall3
        prof_deriv3_wall3
        prof_deriv1_cont3
        prof_deriv2_cont3
        dummy %Final preparation and additiona parameters
        const_line
        const_line2
        blend_C
        ]
    fid = fopen('./Input.txt','w');
    for i=1:size(vars,1)
        fprintf( fid, '%e \n', vars(i,1));
    end
    fclose(fid);
    'input saved'
end
%% I2: Additional user-input (Always check)
% List here new parameters or add them to the input-file.
% If you add them to the input file, make sure you add them:
% - at the end of section I1
% - in section F1 ('Final preparation for creating Blocks') in the 'safe_input' if-condition
%
% Here you can also define some relations: e.g. dsTEu=dsTEut (the xi-spacing at the
% corner of the top boundary should equal the spacing at the upper trailing edge corner)
%--------------------------------------------------------------------------
dsTEu=dsTEut
%--------------------------------------------------------------------------
% Here some default settings
phi_center2=-a_deg
%
N3tphic=round(Nc3*0.75);       % Control Point for approximating angles at outer boundaries
ddphi3t_LEc=0;                 % Second derivative of the boundary angle at N3tphic
N4tphic=round(Nc4*0.75);       % Control Point for approximating angles
ddphi4t_LEc=0;                 % Second derivative of the boundary angle at N4tphic
%
const_line=1;                       % Force number of points being on wall-optimised line
const_line2=round(N_discr*0.1);     % Force number of points being on top-optimised line
blend_C=round(Nc2*0.25)             % Number of points left & right of LE which are uesed for a blending
%--------------------------------------------------------------------------
%Overrule parameters here
prof_deriv1_wall1=3e-4
prof_deriv1_wall2=5e-4
prof_deriv1_wall3=1e-3
prof_deriv1_top1=0.02
prof_deriv1_top2=0.02
prof_deriv1_top3=0.02

%Those parameters are just for the default blunt trailing edge
% Nk=20 %recommended as default for sharp trailing edge: Nk=10
% StretchW_=5 % default as 1 means no stretching (used for sharp TE)
% StretchW=5  % default as 1 means no stretching (used for sharp TE)
% Nk_low=Nk*5 % For blunt trailing edge the control factor can be adjusted separately

%Those parameters are just for the default sharp trailing edge
Nk=10
StretchW_=1
StretchW=1  
Nk_low=Nk
%--------------------------------------------------------------------------
%% Start calculation of airfoil surface
%%   Discretisation of trailing edge (TE) 
if sharp==false
    a_deg =  a_deg+a_rad_offset;
    y_TE  =  linspace(yTE2,yTE1,NTEw);
    x_TE  =  xTEl+(xTEu-xTEl)/(yTE1-yTE2)*(y_TE-yTE2);
    s_TE  =  sqrt((xTEu-xTEl)^2+(yTE1-yTE2)^2);
    
    deriv_yTE=deriv(y_TE',1);
    deriv_xTE=deriv(x_TE',1);
    phi_TE=-atan((xTEu-xTEl)/(yTE1-yTE2))
    'spacing in J at TE set to'
    djTE=s_TE/(NTEw-1)
    sTE=calc_s(x_TE',y_TE')
    dsTE=deriv(sTE',1)
    'Trailing edge generated'
else
    'sharp TE'
    phi_TE=linspace(1,2,NTEw)*0-a_deg;
    djTE=5*dsTEu;
    dsTE=djTE;
    phi_TE=phi_TE.*0+(-a_deg);
end
%%   Contour 1
%_______________________________
% Calculate Spacing
sC1  = calc_s(C1(:,1),C1(:,3));

% Drafting
draft='F'
if draft=='T' % Set draft flag to 'T' to adjust distribution of points
              % Execute section by CTRL+ENTER. When finished, change 
              % parameters in I3&I4 and set draft flag to 'F'
    % Resolution
    Nc1=440 %Nc1
    %
    dscu=0.0015 %dscu
    dsTEu=dsTEu
    %
    ddsTEu=ddsTEu
    ddscu=ddscu
    %     
    dddsTEu=dddsTEu % used for Poly6_end (lines below)
    dddscu=dddscu % used for Poly6        (lines below)
    %
    phi_center1=0;
end

%s_C1=Poly6(Nc1,0,sC1(end),dscu,dsTEu,ddscu,ddsTEu,dddscu,'t');       % <--Choose
s_C1=Poly6_end(Nc1,0,sC1(end),dscu,dsTEu,ddscu,ddsTEu,dddsTEu,'t'); % <--Choose

% Interpolate spacing onto Contour
Ic1  = Interp_onto(Nc1,C1,sC1,s_C1,'f');
s_Ic1= calc_s(Ic1(:,1),Ic1(:,2));

%_______________________________
% Defining angle for J_Gridlines

% wall-normal
phi1_tan=atan(deriv(Ic1(:,2),1)./deriv(Ic1(:,1),1));

% blending kernel
blendF=Poly6(size(phi1_tan,1),0,1,0,0,0,0,0,'f').^1;
blendF=1-fliplr(blendF);
figure
plot(blendF)

% blend from phi_w at center to AoA
phi1_tan_s=phi_center1+blendF'.*(-a_deg-phi_center1);

% check function
dphi1_tan_s=deriv(phi1_tan_s,1);
ddphi1_tan_s=deriv(dphi1_tan_s,1);
if draft=='T'
    figure
    plot(Ic1(:,1),deriv(s_C1',1))
    
    figure
    plot(Ic1(:,1),phi1_tan_s)
    hold on
    plot(Ic1(:,1),phi1_tan)
    plot(Ic1(:,1),deriv(deriv(phi1_tan,1),1))
    stop
end
'Contour 1 generated'
%clearvars 'C1' 'sC1'
%%   Contour 2
close all
% _______________________________
% Calculate Spacing
sC2 = calc_s(C2(:,1),C2(:,3));

% Drafting
draft='F'
if draft=='T' % Set draft flag to 'T' to adjust distribution of points
              % Execute section by CTRL+ENTER. When finished, change 
              % parameters in I3&I4 and set draft flag to 'F'
    % Resolution
    Nc2=Nc2
    %
    dscl=0.0013 %dscu
    dsTEl=dsTEu
    %
    ddsTEl=ddsTEu
    ddscl=ddscu
    %     
    dddsTEl=dddsTEu % used for Poly6_end (lines below)
    dddscl=dddscu   % used for Poly6     (lines below)
    %
end

%s_C2=Poly6(Nc2,0,sC2(end),dscl,dsTEl,ddscl,ddsTEl,dddscl,'f');
s_C2=Poly6_end(Nc2,0,sC2(end),dscl,dsTEl,ddscl,ddsTEl,dddsTEl,'t');

% Interpolate spacing onto Contour
Ic2 = Interp_onto(Nc2,C2,sC2,s_C2,'f');
s_Ic2=calc_s(Ic2(:,1),Ic2(:,2));

%_______________________________
%Defining angle for J_Gridlines

%calc wall normal
phi2_tan=atan(deriv(Ic2(:,2),1)./deriv(Ic2(:,1),1));

%Define function
phi2_tan_s=Poly5_phi(size(phi2_tan,1),C2(1,1),C2(end,1),phi_center2,-a_deg,0,0,0,0,'f');
phi2_tan_s=phi2_tan_s';

%Define function
dphi2_tan_s=deriv(phi2_tan_s,1);
ddphi2_tan_s=deriv(dphi2_tan_s,1);

if draft=='T'
    figure
    plot(Ic2(:,1),deriv(s_C2',1))
    
    figure
    plot(Ic2(:,1),phi2_tan_s)
    hold on
    plot(Ic2(:,1),phi2_tan)
    stop
end

'Contour 2 generated'
%clearvars 'C2' 'sC2'
%%   Contour 3
close all
% Calculate Spacing
sC3 = calc_s(C3(:,1),C3(:,3));
draft='F'
if draft=='T' % Set draft flag to 'T' to adjust distribution of points
              % Execute section by CTRL+ENTER. When finished, change 
              % parameters in I3&I4 and set draft flag to 'F'
    % Resolution
    Nc3=290
    dsLE=0.002 %dsLE
    ddsLE=ddsLE  
    dddsLE=dddsLE % used for Poly6_end (lines below)
end
%s_C3=Poly6_end(Nc3,0,sC3(end),dsLE,dscu,ddsLE,ddscu,dddscu,'t');
s_C3=Poly6(Nc3,0,sC3(end),dsLE,dscu,ddsLE,ddscu,dddsLE,'t');

% Interpolate spacing onto Contour
Ic3 = Interp_onto(Nc3,C3,sC3,s_C3,'f');
s_Ic3=calc_s(Ic3(:,1),Ic3(:,2));

if draft=='T'
    figure
    plot(Ic3(:,1),deriv(s_C3',1))
    hold on
    stop
end

'Contour 3 generated'

%clearvars 'C3' 'sC3'
%%   Contour 4
close all
% Calculate Spacing
sC4 = calc_s(C4(:,1),C4(:,3));
draft='F'
if draft=='T' % Set draft flag to 'T' to adjust distribution of points
              % Execute section by CTRL+ENTER. When finished, change 
              % parameters in I3&I4 and set draft flag to 'F'
    % Resolution
    Nc4=310
    dsLE=dsLE
    ddsLE=ddsLE  
    dddsLE=dddsLE % used for Poly6_end (lines below)
end
%s_C4=Poly6_end(Nc4,0,sC4(end),dsLE,dscl,ddsLE,ddscl,dddscl,'t');
s_C4=Poly6(Nc4,0,sC4(end),dsLE,dscl,ddsLE,ddscl,dddsLE,'t');

% Interpolate spacing onto Contour
Ic4 = Interp_onto(Nc4,C4,sC4,s_C4,'f');
s_Ic4=calc_s(Ic4(:,1),Ic4(:,2));

if draft=='T'
    figure
    plot(Ic4(:,1),deriv(s_C4',1))
    hold on
    stop
end
'Contour 4 generated'

%clearvars 'C3' 'sC3'
%%   Define wall-angels of eta-gridlines <- change blending-speed by modifying 'fac'
close all
%calculate wall normals and correct sign
Ic34=vertcat(flipud(Ic4),Ic3);
phi34=atan(deriv(Ic34(:,2),1)./deriv(Ic34(:,1),1));
figure
plot(phi34(:,1))
turn=1
for i=2:size(phi34(:,1),1)
    if abs(phi34(i)-phi34(i-1))>2        
        abs(phi34(i)-phi34(i-1))
        turn=-1
    end
    phi34(i)=phi34(i)*turn;
end
hold on
plot(phi34(:,1))
phi34(:,1)=-phi34(:,1)
phi4_tan_=phi34(1:size(Ic4(:,1),1))
phi4_tan_=-flipud(phi4_tan_)
phi3_tan_=phi34(size(Ic4(:,1),1)+1:end)
plot(phi4_tan_)
%Section C3 
fac=1
blendF=Poly6(size(phi3_tan_,1),0,1,0,0,0,0,0,'f').^fac;
%blendF=1-fliplr(blendF);
phi3_tan_s=phi3_tan_+blendF'.*(phi_center1-phi3_tan_);
plot(blendF)
plot(s_Ic3,blendF)
figure
plot(s_Ic3,phi3_tan_s)
hold on
plot(s_Ic3,phi3_tan_)

%Check carefully.....harsh changes in the wall-angels can lead to gridline intersections
%
% figure
% plot(s_Ic3,deriv(deriv(phi3_tan_s,1),1))

%Section C4
fac=10
blendF=Poly6(size(phi4_tan_,1),0,1,0,0,0,0,0,'f').^fac;
phi4_tan_s=phi4_tan_+blendF'.*(phi_center2-phi4_tan_);
plot(blendF)
plot(s_Ic4,blendF)
figure
plot(s_Ic4,phi4_tan_s)
hold on
plot(s_Ic4,phi4_tan_)

%check overall

% figure
% plot(Ic3(:,1),deriv(deriv(phi3_tan_,1),1),'--r')
% hold on
% plot(Ic4(:,1),deriv(deriv(phi4_tan_,1),1),'--b')

% plot(Ic4(:,1),phi4_tan_s,'b')
% hold on
% plot(Ic4(:,1),phi4_tan_,'--r')

'Angles for contour 3&4 generated'

%clearvars 'C3' 'sC3'
%%   Merge sections and check contour (no user input needed)
close all

C_tot=flipud(Ic2);
C_tot(size(Ic2,1):size(Ic2,1)+size(Ic4,1)-1,:)=flipud(Ic4);
C_tot(size(Ic2,1)+size(Ic4,1)-1:size(Ic2,1)+size(Ic4,1)-1+size(Ic3,1)-1,:)=Ic3;
C_tot(size(Ic2,1)+size(Ic4,1)-1+size(Ic3,1)-1:size(Ic2,1)+size(Ic4,1)-1+size(Ic3,1)-1+size(Ic1,1)-1,:)=Ic1;
sC_tot = calc_s(C_tot(:,1),C_tot(:,2));
sC_tot = sC_tot';

figure
plot(C_tot(:,1),C_tot(:,2))

figure
plot(sC_tot)
plot(C_tot(:,1),deriv(sC_tot,1))
%plot(deriv(deriv(sC_tot,1),1))
 
%% Calculate transition region downstream of the airfoil
%%   Create shape of airfoil extension (Needs a bit of playing and depends strongly on the airfoil)
%  --> often trade-off between smooth continuation of the airfoil, but
%  avoiding to squeeze the region between upper and lower curve too much  
%  --> the best would be to rewrite this whole section in respect to
%  investigated airfoil.
close all
clearvars C_test C_new C_test2 C_new2
dxcont=C1(end,1)-C1(end-1,1);
xTEc=xTEc_save;
N_discrW=round(xTEc/dxcont,3,'significant');
StretchW_=StretchW
if sharp==true
    StretchW=0;
    s_TE=0;
    djTE=0;
end
s_TE=sqrt((xTEu-xTEl)^2+(yTE1-yTE2)^2);
yTEuC=yTE1+StretchW*s_TE/2 %djTE/2*NTEw-djTE/2*NTEw
yTElC=yTE2-StretchW*s_TE/2 %djTE/2*NTEw-djTE/2*NTEw

% It is tried to adapt the raw-resolution of the airfoil extension to the
%   raw resolution of the airfoil
x_c1  = C1(:,1);
dxC1  = deriv(x_c1,1);
'Overwriting xTEc'
xTEc  = C1(end,1)+(N_discrW-1)*dxC1(end); %xTEc assuming linear spacing corresponding to the last spacing of the raw airfoil profile
xTEc_ = C1(end,1)+(N_discrW/Nk-1)*dxC1(end);
y_c1  = C1(:,3);
dyC1  = deriv(y_c1,1);
ddyC1 = deriv(dyC1,1);
dydxC1=dyC1./dxC1;
ddydxC1=dyC1./dxC1.*deriv((1./dxC1),1)+ddyC1./((dxC1).^2);

yw1=yTE1-(xTEc_-xTEu)*atan(a_deg-a_rad_offset);
yw2=yTE2-(xTEc_-xTEl)*atan(a_deg-a_rad_offset);

%Generate polynomial extension
w_xu=linspace(xTEu,xTEc,N_discrW);
dw_xu=deriv(w_xu',1);
w_yu1=Poly6_s(N_discrW/Nk,xTEu,xTEc_,yTE1,yw1,dydxC1(end),0,ddydxC1(end),0,0,'f');
w_yu1(N_discrW/Nk:size(w_xu,2))=w_yu1(N_discrW/Nk);

%Blend y-coordinate of polynomial with a horizontal line going through the
%  corner of the trailing edge
blendW=Poly6_end(N_discrW/Nk,0,1,0,0,0,0,0,'f');
bl(1:N_discrW/Nk)=blendW;
bl(N_discrW/Nk+1:N_discrW)=bl(end);
blendW=bl;% Assign blending function
clearvars lam bl bl_ a b c d blendWL
w_yu = w_yu1 + blendW.*(yTE1-w_yu1);

figure
plot(C1(:,1),C1(:,3),'-k')
hold on
plot(C2(:,1),C2(:,3),'-k')
plot(w_xu,w_yu1,'--')
plot(w_xu,w_yu1.*0+yTE1,'--')
plot(w_xu,w_yu,'--')

%Blend that curve again to generate diverging end
bl=Poly6_end(N_discrW,0,1,0,0,0,0,0,'f');
blendWL=bl;               %% Assign blending function
w_yu= w_yu + blendWL.*(yTEuC-w_yu);
clearvars lam bl bl_ a b c d blendWL
plot(w_xu,w_yu,'-r')
daspect([1 1 1])
%Generate equidistant curve for later
w_yu2=w_yu-s_TE;

% Create lower line
% Generate the x-coordinates of the lower line in the same way as for the
%     upper one
x_c2  = C2(:,1);
dxC2  = deriv(x_c2,1);
xTEcL = C2(end,1)+(N_discrW-1)*dxC2(end);
w_xl_ = linspace(xTEl,xTEcL,N_discrW);
w_xl_2= linspace(xTEl,xTEc,N_discrW);
y_c2  = C2(:,3);
dyC2  = deriv(y_c2,1);
ddyC2 = deriv(dyC2,1);
dydxC2=dyC2./dxC2;
ddydxC2=dyC2./dxC2.*deriv((1./dxC2),1)+ddyC2./((dxC2).^2);
% Blend the generated x-coordinates with the x-coordinates of the upper
%   extension curve, so that the upper and lower extension end with the
%   same x-coordinate in the wake. That ensures that the eta-gridlines
%   after this trasition region will be vertical
bl=Poly6_end(N_discrW,0,1,0,0,0,0,0,'f');
blendWL=bl;% Assign blending function
w_xl= w_xl_ + blendWL.*(w_xl_2-w_xl_);
clearvars lam bl bl_ a b c d blendWL

%Calculate the y-coordinate of the lower extension again by a polinomial
w_yl_=Poly6_sL(N_discrW/Nk,w_xl(1:N_discrW/Nk),yTE2,yw2,dydxC2(end),0,ddydxC2(end),0,0,'f');
w_yl_(N_discrW/Nk:size(w_xl,2))=w_yl_(N_discrW/Nk);
plot(w_xl,w_yl_)
daspect([1 1 1])
%Blend with equidistant curve to upper extension
blendWL=Poly6_end(N_discrW/Nk_low,0,1,0,0,0,0,0,'f');
blendWL(N_discrW/Nk_low:size(w_xl,2))=blendWL(end);
w_yl=w_yl_ + blendWL.*(w_yu2-w_yl_);
plot(w_xl,w_yu2)
plot(w_xl,w_yl)
clearvars lam bl bl_ a b c d blendWL
%Blend that curve again to generate diverging end
bl=Poly6_end(N_discrW,0,1,0,0,0,0,0,'f');
blendWL=bl;%% Assign blending function
w_yl= w_yl + blendWL.*(yTElC-w_yl);
clearvars lam bl bl_ a b c d blendWL
plot(w_xl,w_yl,'-r')

%......and quickly check it in between:
clearvars C_test1
C_test1(:,1)=C1(:,1);
C_test1(:,2)=C1(:,3);
C_test1(size(C1,1)+1:size(C1,1)+size(w_xu,2)-1,1)=w_xu(1,2:end);
C_test1(size(C1,1)+1:size(C1,1)+size(w_xu,2)-1,2)=w_yu(1,2:end);
C_new1(:,1)=w_xu(1,1:end);
C_new1(:,2)=w_yu(1,1:end);
figure
plot(C_test1(:,1),deriv(deriv(C_test1(:,2),1),1)./(deriv(C_test1(:,1),1).^2)+deriv(C_test1(:,2),1)./deriv(C_test1(:,1),1).*deriv(1./(deriv(C_test1(:,1),1)),1),'b')

C_test2(:,1)=C2(:,1);
C_test2(:,2)=C2(:,3);
C_test2(size(C2,1)+1:size(C2,1)+size(w_xl,2)-1,1)=w_xl(1,2:end);
C_test2(size(C2,1)+1:size(C2,1)+size(w_xl,2)-1,2)=w_yl_(1,2:end);
C_new2(:,1)=w_xl(1,1:end);
C_new2(:,2)=w_yl_(1,1:end);
figure
plot(deriv(deriv(C_test2(:,2),1),1)./(deriv(C_test2(:,1),1).^2)+deriv(C_test2(:,2),1)./deriv(C_test2(:,1),1).*deriv(1./(deriv(C_test2(:,1),1)),1),'b')
'continuous airfoil extension generated'
%%   Discretisation of airfoil extension and corresponding outer domain bound
close all
if sharp == false
    
    sCwu = calc_s(w_xu',w_yu');
    sCwl = calc_s(w_xl',w_yl');
    draft='F'
    if draft=='T'
        Nbuffer=180 %Nbuffer
        Nw=280 %Nw
        dsTEc=dsTEc
        ddsTEc=ddsTEc
        dddsTEc=dddsTEc
        dsOut=0.08
        ddsOut=0
        dddsOut=0
    end   
    % Calculate spacing function for airfoil extension
    s_Cwu_end='f'  % <--choose Poly6
    %s_Cwu_end='t' % <--choose Poly6_end
    if s_Cwu_end=='t'
        s_Cwu=Poly6_end(Nbuffer,0,sCwu(end),dsTEu,dsTEc,ddsTEu,ddsTEc,dddsTEc,'t');
        s_Cwl=Poly6_end(Nbuffer,0,sCwl(end),dsTEl,dsTEc,ddsTEl,ddsTEc,dddsTEc,'t');
    else
        s_Cwu=Poly6(Nbuffer,0,sCwu(end),dsTEu,dsTEc,ddsTEu,ddsTEc,dddsTEu,'t');
        s_Cwl=Poly6(Nbuffer,0,sCwl(end),dsTEl,dsTEc,ddsTEl,ddsTEc,dddsTEl,'t');
    end
    Nw-Nbuffer
    %
    % Calculate spacing function for remaining section unil outlet
    s_Cwu2_end='f'  % <--choose Poly6
    %s_Cwu2_end='t' % <--choose Poly6_end
    if s_Cwu2_end=='t'
        s_Cwu2=Poly6_end(Nw-Nbuffer,sCwu(end),xOut,dsTEc,dsOut,ddsTEc,ddsOut,0,'t');
    else
        s_Cwu2=Poly6(Nw-Nbuffer,sCwu(end),xOut,dsTEc,dsOut,ddsTEc,ddsOut,0,'t');
    end
    %Carry out interpolation for upper side
    Cwu(:,1)=w_xu;
    Cwu(:,3)=w_yu;
    Cwl(:,1)=w_xl;
    Cwl(:,3)=w_yl;
    
    Icwu = Interp_onto(Nbuffer,Cwu,sCwu,s_Cwu,'f');
    s_Icwu=calc_s(Icwu(:,1),Icwu(:,2));
    
    % test
    if draft=='F'
        figure
        plot(Icwu(:,1),deriv(s_Icwu',1))
        hold on
        s_CwuTEST=Poly6(Nw-Nbuffer,Icwu(end,1),xOut,dsTEc,dsOut,ddsTEc,ddsOut,0,'f');
        s_CwuTEST=Poly6_end(Nw-Nbuffer,Icwu(end,1),xOut,dsTEc,dsOut,ddsTEc,ddsOut,0,'f');
        % Check remaining piece until outlet
        s_Cwu2=Poly6(Nw-Nbuffer,sCwu(end),xOut,dsTEc,dsOut,ddsTEc,ddsOut,0,'t');
        
        %plot(s_CwuTEST,deriv(s_CwuTEST',1))
        clear vars s_CwuTEST
        stop
    end
        
    %Carry out interpolation for upper side
    Icwl_1 = Interp_onto(Nbuffer,Cwl,sCwl,s_Cwl,'f');
    Icwl_2(:,1)=Icwu(:,1)-Icwu(1,1)+w_xl(1);
    Icwl_2(:,2)=spline(w_xl,w_yl,Icwl_1(:,1));
    
    lam=[1:round(Nbuffer*0.3)]; %Delaying everything a bit
    a=-2/( 2*((1-size(lam,2)^3)-3*(1-size(lam,2))) - ((1-size(lam,2)^2)-2*(1-size(lam,2)))*3*(1+size(lam,2)) );
    b=-a*(3*(1+size(lam,2)))/2;
    c=-3*a-2*b;
    d=-a-b-c;
    bl_=a*lam.^3+b*lam.^2+c*lam+d;
    bl(1)=0;
    bl(2:size(lam,2)+1)=bl_;
    blendWL=bl;               %% Assign blending function
    blendWL(size(bl,2)+1:size(Icwl_1,1))=1;
    Icwl(:,1)=Icwl_1(:,1);%+blendWL'.*(Icwl_2(:,1)-Icwl_1(:,1));
    Icwl(:,2)=Icwl_1(:,2);%+blendWL'.*(Icwl_2(:,2)-Icwl_1(:,2));
    clearvars lam bl bl_ a b c d blendWL
    
    C_tot2=flipud(Icwl);
    C_tot2(size(Icwl,1)+1:size(Icwl,1)+size(C_tot,1)-1,:)=C_tot(2:end,:);
    C_tot2(size(Icwl,1)+size(C_tot,1)-1:size(Icwl,1)+size(C_tot,1)+size(Icwu,1)-1-1,:)=Icwu;
    
    sC_tot2 = calc_s(C_tot2(:,1),C_tot2(:,2));
    sC_tot2 = sC_tot2';
    
    figure 
    plot(Icwl(:,1),Icwl(:,2))
    hold on
    plot(Icwu(:,1),Icwu(:,2))
    
    figure
    plot(deriv(s_Cwl',1))
    plot(Icwl(:,1),deriv(s_Cwl',1))
    hold on
    plot(Icwu(:,1),deriv(s_Cwu',1))
    hold on
    
    figure
    plot(deriv(s_Cwu2',1))

end % Blunt TE
%
if sharp == true
    draft='F'
    if draft=='T'
        Nbuffer=200 %Nbuffer
        Nw=354 %Nw
        dsTEc=0.05
        ddsTEc=(dsTEc-dsTEu)/2*0.1*1;
        dddsTEc=dddsTEc
        ddsOut=0
    end
    w_xu = (w_xu+w_xl)/2;
    w_yu = (w_yu+w_yl)/2;
    sCwu = calc_s(w_xu',w_yu');
    %
    % Calculate spacing function for airfoil extension
    s_Cwu_end='f'  % <--choose Poly6
    %s_Cwu_end='t' % <--choose Poly6_end
    if s_Cwu_end=='t'
        s_Cwu=Poly6_end(Nbuffer,0,sCwu(end),dsTEu,dsTEc,ddsTEu,ddsTEc,dddsTEc,'t');
    else
        s_Cwu=Poly6(Nbuffer,0,sCwu(end),dsTEu,dsTEc,ddsTEu,ddsTEc,dddsTEu,'t');
    end
    Nw-Nbuffer
    %
    % Calculate spacing function for remaining section unil outlet
    s_Cwu2_end='f'  % <--choose Poly6
    %s_Cwu2_end='t' % <--choose Poly6_end
    if s_Cwu2_end=='t'
        s_Cwu2=Poly6_end(Nw-Nbuffer,sCwu(end),xOut,dsTEc,dsOut,ddsTEc,ddsOut,0,'t');
    else
        s_Cwu2=Poly6(Nw-Nbuffer,sCwu(end),xOut,dsTEc,dsOut,ddsTEc,ddsOut,0,'t');
    end
    %
    Cwu(:,1)=w_xu;
    Cwu(:,3)=w_yu;
    
    Icwu = Interp_onto(Nbuffer,Cwu,sCwu,s_Cwu,'f');
    Icwl=Icwu;
    
    if draft=='T'
        figure
        plot(Icwu(:,1),deriv(s_Icwu',1))
        hold on
        s_CwuTEST=Poly6(Nw-Nbuffer,Icwu(end,1),xOut,dsTEc,dsOut,ddsTEc,ddsOut,0,'f');
        s_CwuTEST=Poly6_end(Nw-Nbuffer,Icwu(end,1),xOut,dsTEc,dsOut,ddsTEc,ddsOut,0,'f');
        s_Cwu2=Poly6(Nw-Nbuffer,sCwu(end),xOut,dsTEc,dsOut,ddsTEc,ddsOut,0,'t');
        s_Cwu=Poly6(Nbuffer,0,sCwu(end),dsTEu,dsTEc,ddsTEu,ddsTEc,dddsTEu,'t');
        plot(s_CwuTEST,deriv(s_CwuTEST',1))
        clear vars s_CwuTEST
        stop
    end
    s_Icwu=calc_s(Icwu(:,1),Icwu(:,2));
end  % Sharp TE

'Extension discretised'
%% Define angels so that the eta-gridlines transition to vertical downstream of the airfoil
close all
phi_w_orig=-atan((Icwu(:,1)-Icwl(:,1))./(Icwu(:,2)-Icwl(:,2)));
phi_wu_=Poly5_phi(size(phi_w_orig,1),Icwu(1,1),Icwu(end,1),phi1_tan_s(end),0,0,0,0,0,'t'); 
phi_wl_=Poly5_phi(size(phi_w_orig,1),Icwl(1,1),Icwl(end,1),phi2_tan_s(end),0,0,0,0,0,'t');
phi_wu=phi_wu_;
phi_wu(size(phi_wu_,2)+1:size(phi_w_orig,1))=0;
phi_wl=phi_wl_;
phi_wl(size(phi_wl_,2)+1:size(phi_w_orig,1))=0;

phi_tot=flipud(phi2_tan_s);
phi_tot(size(Ic2,1):size(Ic2,1)+size(Ic4,1)-1,:)=flipud(phi4_tan_s);
phi_tot(size(Ic2,1)+size(Ic4,1)-1:size(Ic2,1)+size(Ic4,1)-1+size(Ic3,1)-1,:)=phi3_tan_s;
phi_tot(size(Ic2,1)+size(Ic4,1)-1+size(Ic3,1)-1:size(Ic2,1)+size(Ic4,1)-1+size(Ic3,1)-1+size(Ic1,1)-1,:)=phi1_tan_s;

phi_tot2=fliplr(phi_wl);
phi_tot2(size(phi_wl,2)+1:size(phi_wl,2)+size(phi_tot,1)-1)=phi_tot(2:end);
phi_tot2(size(phi_wl,2)+size(phi_tot,1)-1:size(phi_wl,2)+size(phi_tot,1)-1+size(phi_wu,2)-1)=phi_wu;

figure
plot(phi_tot2)
plot(C_tot(:,1),deriv(deriv((phi_tot),1),1))
'Angels defined for extension'
%% Calculate gridlines at the outer boundary of the C-block
%%   Top of Contour 1
close all
clearvars It1
draft='F'
if draft=='T'
    dsTEut=dsTEu
    dscut=0.00215
    ddscut=-2.1e-5  %ddscu
    dddscut=3e-7
end

It1(:,1)=Poly6(Nc1,xcu,xTEu,dscut,dsTEut,ddscut,ddsTEut,dddscut,'t');
%It1(:,1)=Poly6_end(Nc1,xcu,xTEu,dscut,dsTEu,ddscut,ddsTEu,dddsTEu,'t');
It1(:,2)=rad;
st1=calc_s(It1(:,1),It1(:,2));
dst1=deriv(st1',1);
ddst1=deriv(dst1,1);
dddst1=deriv(ddst1,1);

phi1_top=dst1(:,1)*0.0;

if draft=='T'
    figure
    plot(It1(:,1),deriv(It1(:,1),1))
    stop
end

'Top of Contour 1 generated'
%%   TOP of Contour 2
close all
clearvars It2
%
%default <-- Using settings of C1t
dsclt=dscut;
ddsclt=ddscut;
dddsclt=dddscut;
dsTElt=dsTEut;
ddsTElt=ddsTEut;
dddsTElt=dddsTEut;
%
draft='F'
if draft=='T'
    dsclt=dscut
    ddsclt=ddscut
    dddsclt=dddscut
end
It2(:,1)=Poly6(Nc2,xcl,xTEl,dsclt,dsTElt,ddsclt,ddsTElt,dddsclt,'t');
%It2(:,1)=Poly6_end(Nc2,xcl,xTEl,dsclt,dsTElt,ddsclt,ddsTElt,dddsTElt,'f');

It2(:,2)=-rad;
st2=calc_s(It2(:,1),It2(:,2));
dst2=deriv(st2',1);
ddst2=deriv(dst2,1);
dddst2=deriv(ddst2,1);
phi2_top=dst2(:,1)*0.0;

if draft=='T'
    figure
    plot(It2(:,1),deriv(It2(:,1),1))
    stop
end

'Top of Contour 2 generated'
%%   TOP of Contour 3
close all
% Angles at the outer surface are approximated by polynomials to steer
% angle close to symmetry-line to avoid intersections of gridlines
  %N3tphic=round(Nc3*0.75);       % Control Point for approximating angles
  %ddphi3t_LEc=0;                 % Second derivative of the boundary angle at N3tphic
close all
draft='F'
if draft=='T'
    dLE_tan=0.08 %dLE_tan
    ddLE_tan=ddLE_tan
    dddLE_tan=dddLE_tan
end

s_t3 = Poly6_end(size(s_Ic3,2),0,rad*pi/2,dLE_tan,dscut,ddLE_tan,ddscut,dddscut,'t');
%s_t3 = Poly6(size(s_Ic3,2),0,rad*pi/2,dLE_tan,dscut,ddLE_tan,ddscut,dddLE_tan,'t');

phi_poly=linspace(0,pi/2,10000000);
x_poly=x_circ-rad.*cos(phi_poly);
y_poly=rad.*sin(phi_poly);
C3t(:,1)=x_poly;
C3t(:,3)=y_poly;
s_poly(1)=0;
dx_poly=deriv(x_poly',1);
dy_poly=deriv(y_poly',1);
for i=2:size(x_poly,2)
    s_poly(i)=s_poly(i-1)+(dx_poly(i)^2+dy_poly(i)^2)^0.5;
end

figure
plot(deriv(s_t3',1))

It3 = Interp_onto(size(s_Ic3,2),C3t,s_poly,s_t3,'f');

phi3t_tan=atan(deriv(It3(:,2),1)./deriv(It3(:,1),1));
phi3t_tan(1)=pi/2;

dphi3t_tan=deriv(phi3t_tan,1);
ddphi3t_tan=deriv(dphi3t_tan,1);

phi3t_tan_s=smooth(phi3t_tan,0.10,'sgolay',4);
phi3t_tan_s(1:N3tphic) = phi3t_tan_s(1:N3tphic);
phi3t_tan_s(N3tphic:end) = Poly5_phi(size(phi3t_tan(N3tphic:end),1),C3t(N3tphic,3),C3t(end,3),phi3t_tan(N3tphic,1),0,dphi3t_tan(N3tphic,1),0,ddphi3t_tan(N3tphic,1),0,'t');
dphi3t_tan_s=deriv(phi3t_tan_s,1);
ddphi3t_tan_s=deriv(dphi3t_tan_s,1);
phi3t_tan_s(1:N3tphic) = Poly5_phi(size(phi3t_tan(1:N3tphic),1),C3t(1,3),C3t(N3tphic,3),pi/2,phi3t_tan_s(N3tphic),0,dphi3t_tan_s(N3tphic),ddphi3t_LEc,ddphi3t_tan_s(N3tphic),'t');

dphi3t_tan_s=deriv(phi3t_tan_s,1);
ddphi3t_tan_s=deriv(dphi3t_tan_s,1);

if draft=='T'
    figure
    %plot(s_t3,deriv(s_t3',1))
    plot(It3(:,1),deriv(s_t3',1))
    
    figure
    plot(phi3t_tan_s)
    hold on
    plot(phi3t_tan)
    stop
end

'Top of Contour 3 generated'
%%   TOP of Contour 4
close all
% Angles at the outer surface are approximated by polynomials to steer
% angle close to symmetry-line to avoid intersections of gridlines
  N4tphic=round(Nc4*0.75);       % Control Point for approximating angles
  ddphi4t_LEc=0;                 % Second derivative of the boundary angle at N3tphic

s_t4 = Poly6_end(size(s_Ic4,2),0,rad*pi/2,dLE_tan,dsclt,ddLE_tan,ddsclt,dddsclt,'t');
%s_t4 = Poly6(size(s_Ic4,2),0,rad*pi/2,dLE_tan,dsclt,ddLE_tan,ddsclt,dddLE_tan,'t');

clearvars x_poly y_poly dx_poly dy_poly s_poly 
phi_poly=linspace(0,pi/2,10000000);
x_poly=x_circ-rad.*cos(phi_poly);
y_poly=-rad.*sin(phi_poly);
C4t(:,1)=x_poly;
C4t(:,3)=y_poly;
s_poly(1)=0;
dx_poly=deriv(x_poly',1);
dy_poly=deriv(y_poly',1);
for i=2:size(x_poly,2)
    s_poly(i)=s_poly(i-1)+(dx_poly(i)^2+dy_poly(i)^2)^0.5;
end

It4 = Interp_onto(size(s_Ic4,2),C4t,s_poly,s_t4,'f');

phi4t_tan=atan(deriv(It4(:,2),1)./deriv(It4(:,1),1));
phi4t_tan(1)=-pi/2;
dphi4t_tan=deriv(phi4t_tan,1);
ddphi4t_tan=deriv(dphi4t_tan,1);

phi4t_tan_s=smooth(phi4t_tan,0.10,'sgolay',4);
phi4t_tan_s(1:N4tphic) = phi4t_tan_s(1:N4tphic);
phi4t_tan_s(N4tphic:end) = Poly5_phi(size(phi4t_tan(N4tphic:end),1),C4t(N4tphic,3),C4t(end,3),phi4t_tan(N4tphic,1),0,dphi4t_tan(N4tphic,1),0,ddphi4t_tan(N4tphic,1),0,'f');
dphi4t_tan_s=deriv(phi4t_tan_s,1);
ddphi4t_tan_s=deriv(dphi4t_tan_s,1);
phi4t_tan_s(1:N4tphic) = Poly5_phi(size(phi4t_tan(1:N4tphic),1),C4t(1,3),C4t(N4tphic,3),-pi/2,phi4t_tan_s(N4tphic),0,dphi4t_tan_s(N4tphic),ddphi4t_LEc,ddphi4t_tan_s(N4tphic),'f');
dphi4t_tan_s=deriv(phi4t_tan_s,1);
ddphi4t_tan_s=deriv(dphi4t_tan_s,1);


figure
%plot(s_t4,deriv(s_t4',1))
plot(It4(:,1),deriv(s_t4',1))

figure
plot(phi4t_tan_s)
hold on
plot(phi4t_tan)


'Top of Contour 4 generated'

%% Merging zones for generating Block 2 (No input needed)
close all
Ct_tot=flipud(It2);
Ct_tot(size(It2,1):size(It2,1)+size(It4,1)-1,:)=flipud(It4);
Ct_tot(size(It2,1)+size(It4,1)-1:size(It2,1)+size(It4,1)-1+size(It3,1)-1,:)=It3;
Ct_tot(size(It2,1)+size(It4,1)-1+size(It3,1)-1:size(It2,1)+size(It4,1)-1+size(It3,1)-1+size(It1,1)-1,:)=It1;
sCt_tot = calc_s(Ct_tot(:,1),Ct_tot(:,2));
sCt_tot = sCt_tot';
% 
phit_tot=flipud(phi2_top);
phit_tot(size(It2,1):size(It2,1)+size(It4,1)-1,:)=flipud(phi4t_tan_s);
phit_tot(size(It2,1)+size(It4,1)-1:size(It2,1)+size(It4,1)-1+size(It3,1)-1,:)=phi3t_tan_s;
phit_tot(size(It2,1)+size(It4,1)-1+size(It3,1)-1:size(It2,1)+size(It4,1)-1+size(It3,1)-1+size(It1,1)-1,:)=phi1_top;
% 
figure
plot(Ct_tot(:,1),deriv(deriv(sCt_tot,1),1))

'All Contours merged'

%% Define polynomial describing boundary conditions for eta-gridlines
%%   Define and Test Eta-Gridline at TE
close all
%Load input parameters
control_dist=control_dist1
prof_deriv1_wall=prof_deriv1_wall1;
prof_deriv2_wall=prof_deriv2_wall1;
prof_deriv3_wall=prof_deriv3_wall1;
prof_deriv1_cont=prof_deriv1_cont1;
prof_deriv2_cont=prof_deriv2_cont1;
prof_deriv1_top=prof_deriv1_top1;
prof_deriv2_top=prof_deriv2_top1;
prof_deriv3_top=prof_deriv3_top1;
%
i=1 % <---- Specify test-location (1 -> Lower TE) 
%
if 1==1 %(No changes here)
    sig_phi=(Ct_tot(i,2)+1e-20)/(abs(Ct_tot(i,2))+1e-20);
    smax=((C_tot(i,1)-Ct_tot(i,1))^2+(C_tot(i,2)-Ct_tot(i,2))^2)^0.5;
    s=linspace(0,smax,N_discr);
    f_pol=Poly6_sn(N_discr,0,N_discr,0,1,0,0,0,0,0,'f');
    f_blend=f_pol;
    f_blend2=1-fliplr(f_pol);
    clearvars f_pol
    xg1=C_tot(i,1)-sig_phi*s.*sin(phi_tot(i));
    yg1=C_tot(i,2)+sig_phi*s.*cos(phi_tot(i));
    xc1=linspace(C_tot(i,1),Ct_tot(i,1),N_discr);
    yc1=linspace(C_tot(i,2),Ct_tot(i,2),N_discr);
    xg=xg1+f_blend2.*(xc1-xg1);
    yg=yg1+f_blend2.*(yc1-yg1);
    xh2=Ct_tot(i,1)-sig_phi*(-s(end)+s).*sin(phit_tot(i));
    yh2=Ct_tot(i,2)+sig_phi*(-s(end)+s).*cos(phit_tot(i));
    xc2=linspace(C_tot(i,1),Ct_tot(i,1),N_discr);
    yc2=linspace(C_tot(i,2),Ct_tot(i,2),N_discr);
    xh=xc2+f_blend.*(xh2-xc2);
    yh=yc2+f_blend.*(yh2-yc2);
    f_blend=f_blend2;
    x_fin=xg+(xh-xg).*f_blend;
    y_fin=yg+(yh-yg).*f_blend;
    GL(:,1)=x_fin;
    GL(:,3)=y_fin;
    dGLdy=deriv(GL(:,1),1)./deriv(GL(:,2),1);
    ddGLdy=deriv(deriv(GL(:,1),1),1)./((deriv(GL(:,2),1)).^2);
    s_GL=calc_s(GL(:,1),GL(:,3));
end % Generate Shape (No input needed)
%
% Specify distribution
Style=1 % <-- Choose Version 1/2/3
%
%Version 1: Simplest version
if Style==1
    draft='f'
    if draft=='t'
        Ny1=170
        Ny2=500
        NyTOT=Ny1+Ny2-1
        prof_deriv1_wall=0.0003 %gets overwritten for blunt TE
        prof_deriv1_top=0.03
        prof_deriv2_top=0
        prof_deriv3_top=0e-4
        prof_deriv2_wall=0;
        prof_deriv3_wall=0;
        prof_deriv1_cont=0;
        prof_deriv2_cont=0;
    end
    if sharp == false
        prof_deriv1_wall=djTE
    end
    s_JJ=Poly6_J(Ny1+Ny2-1,0,s_GL(end),prof_deriv1_wall,prof_deriv1_top,prof_deriv2_wall,prof_deriv2_top,prof_deriv3_wall,'t');
end % <-- Set Spacing
%Version 2: Simple simple & no stretching at both ends (Cheaper)
% --> a good starting point for Version 3
if Style==2 || Style==3
    draft='f'
    if draft=='t'
        NyTOT=100
        prof_deriv1_top=1
        prof_deriv2_top=0e-2
        prof_deriv2_top=0e-4
        prof_deriv1_wall=2e-6
        prof_deriv2_wall=0
        prof_deriv3_wall=0
        prof_deriv1_cont=0;
        prof_deriv2_cont=0;
    end
    if sharp == false
        prof_deriv1_wall=djTE
    end
    s_JJ=Poly6_Jend(NyTOT,0,s_GL(end),prof_deriv1_wall,prof_deriv1_top,prof_deriv2_wall,prof_deriv2_top,prof_deriv3_top,'t');
    figure
    plot(s_JJ,deriv(s_JJ',1))
end % <-- Set Spacing
% First draft for Version 3
if Style==3
    %Calculate starting point based on Style 2
    suggest='t'
    if suggest=='t'
        control_dist=0.05 %control_dist3
        %index=find(abs(s_JJ-control_dist)<0.001) % <-- Default for Ny1
        index=125 % <-- if Ny1 is already fix defined
        Ny1=index(end)
        Ny2=NyTOT-index(end)+1
        control_dist=s_JJ(index(end))
        ds_JJ=deriv(s_JJ',1);
        prof_deriv1_cont=ds_JJ(index(end));
        dds_JJ=deriv(ds_JJ,1);
        prof_deriv2_cont=dds_JJ(index(end));
    end
    % Fine-tuning for Version 3: Takes a bit of playing around
    draft='f'
    if draft=='t'
        index=125
        Ny1=index
        Ny2=NyTOT-Ny1+1
        prof_deriv1_cont=prof_deriv1_cont*0.9
        prof_deriv2_cont=prof_deriv2_cont*1.5
    end
    s_JJ=Poly_J(Ny1,Ny2,0,control_dist,s_GL(end),prof_deriv1_wall,prof_deriv1_cont,prof_deriv1_top,prof_deriv2_wall,prof_deriv2_cont,prof_deriv2_top,'t');
end
figure
plot(s_JJ,deriv(s_JJ',1))
%%   Saving Parameters 
control_dist1=control_dist;
prof_deriv1_wall1=prof_deriv1_wall;
prof_deriv2_wall1=prof_deriv2_wall;
prof_deriv3_wall1=prof_deriv3_wall;
prof_deriv1_cont1=prof_deriv1_cont;
prof_deriv2_cont1=prof_deriv2_cont;
prof_deriv1_top1=prof_deriv1_top;
prof_deriv2_top1=prof_deriv2_top;
prof_deriv3_top1=prof_deriv3_top;
'TE defined'
%%   Define and Test Eta-Gridline at Center
close all
%Load input parameters
control_dist=control_dist2
prof_deriv1_wall=prof_deriv1_wall2;
prof_deriv2_wall=prof_deriv2_wall2;
prof_deriv3_wall=prof_deriv3_wall2;
prof_deriv1_cont=prof_deriv1_cont2;
prof_deriv2_cont=prof_deriv2_cont2;
prof_deriv1_top=prof_deriv1_top2;
prof_deriv2_top=prof_deriv2_top2;
prof_deriv3_top=prof_deriv3_top2;
%
i=Nc2 % <---- Specify test-location (1 -> Lower TE)
%
if 1==1 %(No changes here)
    sig_phi=(Ct_tot(i,2)+1e-20)/(abs(Ct_tot(i,2))+1e-20);
    smax=((C_tot(i,1)-Ct_tot(i,1))^2+(C_tot(i,2)-Ct_tot(i,2))^2)^0.5;
    s=linspace(0,smax,N_discr);
    f_pol=Poly6_sn(N_discr,0,N_discr,0,1,0,0,0,0,0,'f');
    f_blend=f_pol;
    f_blend2=1-fliplr(f_pol);
    clearvars f_pol
    xg1=C_tot(i,1)-sig_phi*s.*sin(phi_tot(i));
    yg1=C_tot(i,2)+sig_phi*s.*cos(phi_tot(i));
    xc1=linspace(C_tot(i,1),Ct_tot(i,1),N_discr);
    yc1=linspace(C_tot(i,2),Ct_tot(i,2),N_discr);
    xg=xg1+f_blend2.*(xc1-xg1);
    yg=yg1+f_blend2.*(yc1-yg1);
    xh2=Ct_tot(i,1)-sig_phi*(-s(end)+s).*sin(phit_tot(i));
    yh2=Ct_tot(i,2)+sig_phi*(-s(end)+s).*cos(phit_tot(i));
    xc2=linspace(C_tot(i,1),Ct_tot(i,1),N_discr);
    yc2=linspace(C_tot(i,2),Ct_tot(i,2),N_discr);
    xh=xc2+f_blend.*(xh2-xc2);
    yh=yc2+f_blend.*(yh2-yc2);
    f_blend=f_blend2;
    x_fin=xg+(xh-xg).*f_blend;
    y_fin=yg+(yh-yg).*f_blend;
    GL(:,1)=x_fin;
    GL(:,3)=y_fin;
    dGLdy=deriv(GL(:,1),1)./deriv(GL(:,2),1);
    ddGLdy=deriv(deriv(GL(:,1),1),1)./((deriv(GL(:,2),1)).^2);
    s_GL=calc_s(GL(:,1),GL(:,3));
end % Generate Shape (No input needed)
%
%Version 1: Simplest version
if Style==1
    draft='f'
    if draft=='t'
        control_dist=control_dist1
        prof_deriv1_wall=prof_deriv1_wall1;
        prof_deriv2_wall=prof_deriv2_wall1;
        prof_deriv3_wall=prof_deriv3_wall1;
        prof_deriv1_cont=prof_deriv1_cont1;
        prof_deriv2_cont=prof_deriv2_cont1;
        prof_deriv1_top=prof_deriv1_top1;
        prof_deriv2_top=prof_deriv2_top1;
        prof_deriv3_top=prof_deriv3_top1;
    end
    s_JJ=Poly6_J(Ny1+Ny2-1,0,s_GL(end),prof_deriv1_wall,prof_deriv1_top,prof_deriv2_wall,prof_deriv2_top,prof_deriv3_wall,'t');
end % <-- Set Spacing
%Version 2: Simple simple & no stretching at both ends (Cheaper)
% --> a good starting point for Version 3
if Style==2 || Style==3
    draft='f'
    if draft=='t'
        prof_deriv1_top=1
        prof_deriv2_top=0e-2
        prof_deriv3_top=0e-4
        prof_deriv1_wall=2e-6
        prof_deriv2_wall=0
        prof_deriv3_wall=0
    end
    s_JJ=Poly6_Jend(NyTOT,0,s_GL(end),prof_deriv1_wall,prof_deriv1_top,prof_deriv2_wall,prof_deriv2_top,prof_deriv3_top,'t');
end % <-- Set Spacing
% First draft for Version 3
if Style==3
    % Fine-tuning for Version 3: Takes a bit of playing around
    suggest='t'
    if suggest=='t'
        control_dist=control_dist
        index=find(abs(s_JJ-control_dist)<0.001) % <-- Default for Ny1
        index=index(end)
        index=125 % <-- Index must agree with other positions
        Ny1=index
        Ny2=NyTOT-index+1
        control_dist=s_JJ(index)
        ds_JJ=deriv(s_JJ',1);
        prof_deriv1_cont=ds_JJ(index);
        dds_JJ=deriv(ds_JJ,1);
        prof_deriv2_cont=dds_JJ(index);
    end
    draft='f'
    if draft=='t'
        Ny1=125
        Ny2=NyTOT-Ny1+1
        prof_deriv1_cont=prof_deriv1_cont
        prof_deriv2_cont=prof_deriv2_cont
    end
    s_JJ=Poly_J(Ny1,Ny2,0,control_dist,s_GL(end),prof_deriv1_wall,prof_deriv1_cont,prof_deriv1_top,prof_deriv2_wall,prof_deriv2_cont,prof_deriv2_top,'t');
end
figure
plot(s_JJ,deriv(s_JJ',1))
%%   Saving Parameters 
control_dist2=control_dist;
prof_deriv1_wall2=prof_deriv1_wall;
prof_deriv2_wall2=prof_deriv2_wall;
prof_deriv3_wall2=prof_deriv3_wall;
prof_deriv1_cont2=prof_deriv1_cont;
prof_deriv2_cont2=prof_deriv2_cont;
prof_deriv1_top2=prof_deriv1_top;
prof_deriv2_top2=prof_deriv2_top;
prof_deriv3_top2=prof_deriv3_top;
'Center defined'
%%   Define and Test Eta-Gridline at LE
close all
%
%Load input parameters
control_dist=control_dist3
prof_deriv1_wall=prof_deriv1_wall3;
prof_deriv2_wall=prof_deriv2_wall3;
prof_deriv3_wall=prof_deriv3_wall3;
prof_deriv1_cont=prof_deriv1_cont3;
prof_deriv2_cont=prof_deriv2_cont3;
prof_deriv1_top=prof_deriv1_top3;
prof_deriv2_top=prof_deriv2_top3;
prof_deriv3_top=prof_deriv3_top3;
%
i=Nc2+Nc4 % <---- Specify test-location (1 -> Lower TE)
%
if 1==1 %(No changes here)
    sig_phi=(Ct_tot(i,2)+1e-20)/(abs(Ct_tot(i,2))+1e-20);
    smax=((C_tot(i,1)-Ct_tot(i,1))^2+(C_tot(i,2)-Ct_tot(i,2))^2)^0.5;
    s=linspace(0,smax,N_discr);
    f_pol=Poly6_sn(N_discr,0,N_discr,0,1,0,0,0,0,0,'f');
    f_blend=f_pol;
    f_blend2=1-fliplr(f_pol);
    clearvars f_pol
    xg1=C_tot(i,1)-sig_phi*s.*sin(phi_tot(i));
    yg1=C_tot(i,2)+sig_phi*s.*cos(phi_tot(i));
    xc1=linspace(C_tot(i,1),Ct_tot(i,1),N_discr);
    yc1=linspace(C_tot(i,2),Ct_tot(i,2),N_discr);
    xg=xg1+f_blend2.*(xc1-xg1);
    yg=yg1+f_blend2.*(yc1-yg1);
    xh2=Ct_tot(i,1)-sig_phi*(-s(end)+s).*sin(phit_tot(i));
    yh2=Ct_tot(i,2)+sig_phi*(-s(end)+s).*cos(phit_tot(i));
    xc2=linspace(C_tot(i,1),Ct_tot(i,1),N_discr);
    yc2=linspace(C_tot(i,2),Ct_tot(i,2),N_discr);
    xh=xc2+f_blend.*(xh2-xc2);
    yh=yc2+f_blend.*(yh2-yc2);
    f_blend=f_blend2;
    x_fin=xg+(xh-xg).*f_blend;
    y_fin=yg+(yh-yg).*f_blend;
    GL(:,1)=x_fin;
    GL(:,3)=y_fin;
    dGLdy=deriv(GL(:,1),1)./deriv(GL(:,2),1);
    ddGLdy=deriv(deriv(GL(:,1),1),1)./((deriv(GL(:,2),1)).^2);
    s_GL=calc_s(GL(:,1),GL(:,3));
end % Generate Shape (No input needed)
%
%Version 1: Simplest version
if Style==1
    draft='f'
    if draft=='t'
        control_dist=control_dist1
        prof_deriv1_wall=prof_deriv1_wall1;
        prof_deriv2_wall=prof_deriv2_wall1;
        prof_deriv3_wall=prof_deriv3_wall1;
        prof_deriv1_cont=prof_deriv1_cont1;
        prof_deriv2_cont=prof_deriv2_cont1;
        prof_deriv1_top=prof_deriv1_top1;
        prof_deriv2_top=prof_deriv2_top1;
        prof_deriv3_top=prof_deriv3_top1;
    end
    s_JJ=Poly6_J(Ny1+Ny2-1,0,s_GL(end),prof_deriv1_wall,prof_deriv1_top,prof_deriv2_wall,prof_deriv2_top,prof_deriv3_wall,'t');
end % <-- Set Spacing
%Version 2: Simple simple & no stretching at both ends (Cheaper)
% --> a good starting point for Version 3
if Style==2 || Style==3
    draft='f'
    if draft=='t'
        prof_deriv1_top=0.01 %prof_deriv1_top
        prof_deriv2_top=prof_deriv2_top
        prof_deriv3_top=prof_deriv3_top
        prof_deriv1_wall=prof_deriv1_wall
        prof_deriv2_wall=prof_deriv2_wall
        prof_deriv3_wall=prof_deriv3_wall
    end
    s_JJ=Poly6_Jend(Ny1+Ny2-1,0,s_GL(end),prof_deriv1_wall,prof_deriv1_top,prof_deriv2_wall,prof_deriv2_top,prof_deriv3_top,'t');
end % <-- Set Spacing
% First draft for Version 3
if Style==3
    suggest='t'
    if suggest=='t'
        %control_dist=0.05
        %index=find(abs(s_JJ-control_dist)<0.01) % <-- Index must agree with other positions
        index=125 %index(end)
        Ny1=index
        Ny2=NyTOT-index+1
        control_dist=s_JJ(index)
        ds_JJ=deriv(s_JJ',1);
        prof_deriv1_cont=ds_JJ(index);
        prof_deriv1_cont=prof_deriv1_cont;
        dds_JJ=deriv(ds_JJ,1);
        prof_deriv2_cont=dds_JJ(index);
    end
    % Fine-tuning for Version 3: Takes a bit of playing around
    draft='t'
    if draft=='t'
          prof_deriv1_wall=prof_deriv1_wall3;
          prof_deriv1_top=prof_deriv1_top3
    end
    s_JJ=Poly_J(Ny1(1),Ny2(1),0,control_dist(1),s_GL(end),prof_deriv1_wall,prof_deriv1_cont(1),prof_deriv1_top,prof_deriv2_wall,prof_deriv2_cont(1),prof_deriv2_top,'t');
end
figure
plot(s_JJ,deriv(s_JJ',1))

figure
plot(s_JJ)
%%   Saving Parameters 
control_dist3=control_dist;
prof_deriv1_wall3=prof_deriv1_wall;
prof_deriv2_wall3=prof_deriv2_wall;
prof_deriv3_wall3=prof_deriv3_wall;
prof_deriv1_cont3=prof_deriv1_cont;
prof_deriv2_cont3=prof_deriv2_cont;
prof_deriv1_top3=prof_deriv1_top;
prof_deriv2_top3=prof_deriv2_top;
prof_deriv3_top3=prof_deriv3_top;
'LE defined'
%%   EPS FIGURES in order to compare grids (just uncomment, if needed)
figure
axes1 = axes('Parent',figure);
hold(axes1,'on');
box(axes1,'on');
set(axes1,'YTick',...
[0 0.002 0.004 0.006 0.008 0.01 0.012 0.014 0.016 0.018 0.02 0.022]);
plot(s_JJ,deriv(s_JJ',1))
ylim(axes1,[0 0.022]);
xlim(axes1,[0 7.5]);

figure
axes1 = axes('Parent',figure);
hold(axes1,'on');
box(axes1,'on');
set(axes1,'YTick',...
[0 0.002 0.004 0.006 0.008 0.01 0.012 0.014 0.016 0.018 0.02 0.022]);
plot(deriv(s_JJ',1))
ylim(axes1,[0 0.022]);
xlim(axes1,[1 689]);
%%   Generate Distribution function (no input needed)
prof_deriv1_wall(1:size(Ic2,1)) = Poly6_s(size(Ic2,1),0,xTElt-0.5,prof_deriv1_wall1,prof_deriv1_wall2,0,0,0,0,0,'f');
prof_deriv1_wall(size(Ic2,1):size(Ic2,1)+size(Ic4,1)-1) = Poly6_s(size(Ic4,1),0,0.5-xq,prof_deriv1_wall2,prof_deriv1_wall3,0,0,0,0,0,'f');
prof_deriv1_wall(size(Ic2,1)+size(Ic4,1)-1:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2) = Poly6_s(size(Ic3,1),0,0.5-xq,prof_deriv1_wall3,prof_deriv1_wall2,0,0,0,0,0,'f');
prof_deriv1_wall(size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)+size(Ic1,1)-3) = Poly6_s(size(Ic1,1),0,xTEt-0.5,prof_deriv1_wall2,prof_deriv1_wall1,0,0,0,0,0,'f');
figure
plot(prof_deriv1_wall)
%
prof_deriv2_wall(1:size(Ic2,1)) = Poly6_s(size(Ic2,1),0,xTElt-0.5,prof_deriv2_wall1,prof_deriv2_wall2,0,0,0,0,0,'f');
prof_deriv2_wall(size(Ic2,1):size(Ic2,1)+size(Ic4,1)-1) = Poly6_s(size(Ic4,1),0,0.5-xq,prof_deriv2_wall2,prof_deriv2_wall3,0,0,0,0,0,'f');
prof_deriv2_wall(size(Ic2,1)+size(Ic4,1)-1:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2) = Poly6_s(size(Ic3,1),0,0.5-xq,prof_deriv2_wall3,prof_deriv2_wall2,0,0,0,0,0,'f');
prof_deriv2_wall(size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)+size(Ic1,1)-3) = Poly6_s(size(Ic1,1),0,xTEt-0.5,prof_deriv2_wall2,prof_deriv2_wall1,0,0,0,0,0,'f');
figure
plot(prof_deriv2_wall)
%
prof_deriv3_wall(1:size(Ic2,1)) = Poly6_s(size(Ic2,1),0,xTElt-0.5,prof_deriv3_wall1,prof_deriv3_wall2,0,0,0,0,0,'f');
prof_deriv3_wall(size(Ic2,1):size(Ic2,1)+size(Ic4,1)-1) = Poly6_s(size(Ic4,1),0,0.5-xq,prof_deriv3_wall2,prof_deriv3_wall3,0,0,0,0,0,'f');
prof_deriv3_wall(size(Ic2,1)+size(Ic4,1)-1:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2) = Poly6_s(size(Ic3,1),0,0.5-xq,prof_deriv3_wall3,prof_deriv3_wall2,0,0,0,0,0,'f');
prof_deriv3_wall(size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)+size(Ic1,1)-3) = Poly6_s(size(Ic1,1),0,xTEt-0.5,prof_deriv3_wall2,prof_deriv3_wall1,0,0,0,0,0,'f');
figure
plot(prof_deriv3_wall)
%
prof_deriv1_top(1:size(Ic2,1)) = Poly6_s(size(Ic2,1),0,xTElt-0.5,prof_deriv1_top1,prof_deriv1_top2,0,0,0,0,0,'f');
prof_deriv1_top(size(Ic2,1):size(Ic2,1)+size(Ic4,1)-1) = Poly6_s(size(Ic4,1),0,0.5-xq,prof_deriv1_top2,prof_deriv1_top3,0,0,0,0,0,'f');
prof_deriv1_top(size(Ic2,1)+size(Ic4,1)-1:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2) = Poly6_s(size(Ic3,1),0,0.5-xq,prof_deriv1_top3,prof_deriv1_top2,0,0,0,0,0,'f');
prof_deriv1_top(size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)+size(Ic1,1)-3) = Poly6_s(size(Ic1,1),0,xTEt-0.5,prof_deriv1_top2,prof_deriv1_top1,0,0,0,0,0,'f');
figure
plot(prof_deriv1_top)
%
prof_deriv2_top(1:size(Ic2,1)) = Poly6_s(size(Ic2,1),0,xTElt-0.5,prof_deriv2_top1,prof_deriv2_top2,0,0,0,0,0,'f');
prof_deriv2_top(size(Ic2,1):size(Ic2,1)+size(Ic4,1)-1) = Poly6_s(size(Ic4,1),0,0.5-xq,prof_deriv2_top2,prof_deriv2_top3,0,0,0,0,0,'f');
prof_deriv2_top(size(Ic2,1)+size(Ic4,1)-1:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2) = Poly6_s(size(Ic3,1),0,0.5-xq,prof_deriv2_top3,prof_deriv2_top2,0,0,0,0,0,'f');
prof_deriv2_top(size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)+size(Ic1,1)-3) = Poly6_s(size(Ic1,1),0,xTEt-0.5,prof_deriv2_top2,prof_deriv2_top1,0,0,0,0,0,'f');
figure
plot(prof_deriv2_top)
%
prof_deriv3_top(1:size(Ic2,1)) = Poly6_s(size(Ic2,1),0,xTElt-0.5,prof_deriv3_top1,prof_deriv3_top2,0,0,0,0,0,'f');
prof_deriv3_top(size(Ic2,1):size(Ic2,1)+size(Ic4,1)-1) = Poly6_s(size(Ic4,1),0,0.5-xq,prof_deriv3_top2,prof_deriv3_top3,0,0,0,0,0,'f');
prof_deriv3_top(size(Ic2,1)+size(Ic4,1)-1:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2) = Poly6_s(size(Ic3,1),0,0.5-xq,prof_deriv3_top3,prof_deriv3_top2,0,0,0,0,0,'f');
prof_deriv3_top(size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)+size(Ic1,1)-3) = Poly6_s(size(Ic1,1),0,xTEt-0.5,prof_deriv3_top2,prof_deriv3_top1,0,0,0,0,0,'f');
figure
plot(prof_deriv3_top)
%
prof_deriv1_cont(1:size(Ic2,1)) = Poly6_s(size(Ic2,1),0,xTElt-0.5,prof_deriv1_cont1,prof_deriv1_cont2,0,0,0,0,0,'f');
prof_deriv1_cont(size(Ic2,1):size(Ic2,1)+size(Ic4,1)-1) = Poly6_s(size(Ic4,1),0,0.5-xq,prof_deriv1_cont2,prof_deriv1_cont3,0,0,0,0,0,'f');
prof_deriv1_cont(size(Ic2,1)+size(Ic4,1)-1:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2) = Poly6_s(size(Ic3,1),0,0.5-xq,prof_deriv1_cont3,prof_deriv1_cont2,0,0,0,0,0,'f');
prof_deriv1_cont(size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)+size(Ic1,1)-3) = Poly6_s(size(Ic1,1),0,xTEt-0.5,prof_deriv1_cont2,prof_deriv1_cont1,0,0,0,0,0,'f');
figure
plot(prof_deriv1_cont)
%
prof_deriv2_cont(1:size(Ic2,1)) = Poly6_s(size(Ic2,1),0,xTElt-0.5,prof_deriv2_cont1,prof_deriv2_cont2,0,0,0,0,0,'f');
prof_deriv2_cont(size(Ic2,1):size(Ic2,1)+size(Ic4,1)-1) = Poly6_s(size(Ic4,1),0,0.5-xq,prof_deriv2_cont2,prof_deriv2_cont3,0,0,0,0,0,'f');
prof_deriv2_cont(size(Ic2,1)+size(Ic4,1)-1:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2) = Poly6_s(size(Ic3,1),0,0.5-xq,prof_deriv2_cont3,prof_deriv2_cont2,0,0,0,0,0,'f');
prof_deriv2_cont(size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)+size(Ic1,1)-3) = Poly6_s(size(Ic1,1),0,xTEt-0.5,prof_deriv2_cont2,prof_deriv2_cont1,0,0,0,0,0,'f');
figure
plot(prof_deriv2_cont)
%
control_dist(1:size(Ic2,1)) = Poly6_s(size(Ic2,1),0,xTElt-0.5,control_dist1,control_dist2,0,0,0,0,0,'f');
control_dist(size(Ic2,1):size(Ic2,1)+size(Ic4,1)-1) = Poly6_s(size(Ic4,1),0,0.5-xq,control_dist2,control_dist3,0,0,0,0,0,'f');
control_dist(size(Ic2,1)+size(Ic4,1)-1:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2) = Poly6_s(size(Ic3,1),0,0.5-xq,control_dist3,control_dist2,0,0,0,0,0,'f');
control_dist(size(Ic2,1)+size(Ic4,1)+size(Ic3,1)-2:size(Ic2,1)+size(Ic4,1)+size(Ic3,1)+size(Ic1,1)-3) = Poly6_s(size(Ic1,1),0,xTEt-0.5,control_dist2,control_dist1,0,0,0,0,0,'f');
figure
plot(control_dist)

prof_deriv1_wall_l1_BU = prof_deriv1_wall(1);
prof_deriv1_wall_u1_BU = prof_deriv1_wall(end);
'distribution function generated'
%% Final preparation for creating Blocks 
close all

%const_line=1;                       % Force number of points being on wall-optimised line
%const_line2=round(N_discr*0.1);     % Force number of points being on top-optimised line
%blend_C=round(Nc2*0.25)             % Number of points left & right of LE which are uesed for a blending

debug=false;
if debug
    figure
    plot(C_tot(:,1),C_tot(:,2))
    hold on
    plot(Ct_tot(:,1),Ct_tot(:,2))
    daspect([1 1 1])
end

% Create directory and safe Input file
safe_input='on'
formatOut = 'HH_MM_mm_dd_yy';
DateString = datestr(datetime,formatOut)
name_dir= ['Grid@',DateString]
mkdir(name_dir)
name_output= [name_dir,'/.']
name_input= [name_dir,'/Input.txt']
name_input_usr= [name_dir,'/Input_readable.txt']
dummy=1e99
status = system(['cp PolyGridWizZ_V*.m ', name_output])
if safe_input=='on'
    vars=[
        dummy %TE
        NTEw
        dummy %C1
        Nc1
        dscu
        dsTEu
        ddsTEu
        ddscu
        dddsTEu
        dddscu
        phi_center1
        dummy %C2
        Nc2
        dscl
        dsTEl
        ddsTEl
        ddscl
        dddsTEl
        dddscl
        phi_center2
        dummy %C3
        Nc3
        dsLE
        ddsLE
        dddsLE
        dummy %C4
        Nc4
        dsLE
        ddsLE
        dddsLE
        dummy %Aifoil extension
        Nbuffer
        Nw
        xTEc_save
        StretchW_
        Nk
        dsTEc
        ddsTEc
        dddsTEc
        dummy %Outlet
        dsOut
        ddsOut
        dddsOut
        dummy %TC1
        dsTEut
        ddsTEut
        dddsTEut
        dscut
        ddscut
        dddscut
        dummy %TC2
        dsTElt
        ddsTElt
        dddsTElt
        dsclt
        ddsclt
        dddsclt
        dummy %TC3
        N3tphic
        ddphi3t_LEc
        dLE_tan
        ddLE_tan
        dddLE_tan
        dummy %TC4
        N4tphic
        ddphi4t_LEc
        dummy %Eta gridlines
        dummy %TE
        Ny1
        Ny2
        NyTOT
        control_dist1
        control_dist2
        control_dist3
        prof_deriv1_top1
        prof_deriv2_top1
        prof_deriv3_top1
        prof_deriv1_wall1
        prof_deriv2_wall1
        prof_deriv3_wall1
        prof_deriv1_cont1
        prof_deriv2_cont1
        dummy %center
        prof_deriv1_top2
        prof_deriv2_top2
        prof_deriv3_top2
        prof_deriv1_wall2
        prof_deriv2_wall2
        prof_deriv3_wall2
        prof_deriv1_cont2
        prof_deriv2_cont2
        dummy %LE
        prof_deriv1_top3
        prof_deriv2_top3
        prof_deriv3_top3
        prof_deriv1_wall3
        prof_deriv2_wall3
        prof_deriv3_wall3
        prof_deriv1_cont3
        prof_deriv2_cont3
        dummy %Final preparation and additiona parameters
        const_line
        const_line2
        blend_C
        ]
    fid2 =fopen('Input_parameters.txt')
    fid = fopen(name_input_usr,'w');
    for i=1:size(vars,1)
        notes=fgetl(fid2);
        fprintf( fid, '%e \t %s\n', vars(i,1), notes);
    end
    fclose(fid2);
    fclose(fid);
    fid = fopen(name_input,'w');
    for i=1:size(vars,1)
        notes=fgetl(fid2);
        fprintf( fid, '%e \n', vars(i,1));
    end
    fclose(fid);
    fid = fopen('./Input.txt','w');
    for i=1:size(vars,1)
        notes=fgetl(fid2);
        fprintf( fid, '%e \n', vars(i,1));
    end
    fclose(fid);
end

'--- Ready for generating Block 2 ---'
stop
%% Loop trough every point in xi-direction (This can take a while)
%%      Test Loop for Block 2 

%Just copy the section Generate Block 2 and plot whatever you want
%You may want to comment out the section of interpolatin the spacing on the
%  gridlines (the interpolation is quite lame)
%Also you may switch off the writing out routine that is performed every loop
%--> I know that is very inefficient, but it was necessary to keep the RAM low.
%    Feel free to save all grid-points in an array and write the whole block at once. 

'Start testing Block 2'
step=10
figure
daspect([1 1 1])
for i=2:step:size(C_tot,1)-1
    i
    
    sig_phi=(Ct_tot(i,2)+1e-20)/(abs(Ct_tot(i,2))+1e-20);
    if(Ct_tot(i,2)==0)
        phi_tot(i)=pi/2;
    end
    smax=((C_tot(i,1)-Ct_tot(i,1))^2+(C_tot(i,2)-Ct_tot(i,2))^2)^0.5;

    s=linspace(0,smax,N_discr);
    
    f_pol=Poly6_sn(N_discr,0,N_discr,0,1,0,0,0,0,0,'f');
    f_blend=f_pol;
    f_blend2=1-fliplr(f_pol);
    clearvars f_pol
    
    xg1=C_tot(i,1)-sig_phi*s.*sin(phi_tot(i));
    yg1=C_tot(i,2)+sig_phi*s.*cos(phi_tot(i));
    xc1=linspace(C_tot(i,1),Ct_tot(i,1),N_discr);
    yc1=linspace(C_tot(i,2),Ct_tot(i,2),N_discr);
    
    xg=xg1+f_blend2.*(xc1-xg1);
    yg=yg1+f_blend2.*(yc1-yg1);
    
    xh2=Ct_tot(i,1)-sig_phi*(-s(end)+s).*sin(phit_tot(i));
    yh2=Ct_tot(i,2)+sig_phi*(-s(end)+s).*cos(phit_tot(i));
    xc2=linspace(C_tot(i,1),Ct_tot(i,1),N_discr);
    yc2=linspace(C_tot(i,2),Ct_tot(i,2),N_discr);
    
    xh=xc2+f_blend.*(xh2-xc2);
    yh=yc2+f_blend.*(yh2-yc2);
    
    x_fin=xg+(xh-xg).*f_blend2;
    y_fin=yg+(yh-yg).*f_blend2;
    
    GL(:,1)=x_fin;
    GL(:,3)=y_fin;
    
    s_GL=calc_s(GL(:,1),GL(:,3));
    
     %Version 1: Simplest version
    if Style==1
        s_J=Poly6_J(Ny1+Ny2-1,0,s_GL(end),prof_deriv1_wall(i),prof_deriv1_top(i),prof_deriv2_wall(i),prof_deriv2_top(i),prof_deriv3_wall(i),'f');
        %s_JJ=Poly6_J(Ny1+Ny2-1,0,s_GL(end),prof_deriv1_wall(i),dy_TEb,prof_deriv2_wall,prof_deriv2_top,prof_deriv3_top,'t');
    end
    %Version 2: Simple simple & no stretching at both ends (Cheaper)
    if Style==2
        s_J=Poly6_Jend(NyTOT,0,s_GL(end),prof_deriv1_wall(i),prof_deriv1_top(i),prof_deriv2_wall(i),prof_deriv2_top(i),prof_deriv2_top(i),'f');
    end
    % First draft for Version 3
    if Style==3
        s_J=Poly_J(Ny1,Ny2,0,control_dist(i),s_GL(end),prof_deriv1_wall(i),prof_deriv1_cont(i),prof_deriv1_top(i),prof_deriv2_wall(i),prof_deriv2_cont(i),prof_deriv2_top(i),'f');
    end
    
    GL_inter = Interp_onto(size(s_J,2),GL,s_GL,s_J,'f');

    plot(GL_inter(:,1),GL_inter(:,2))
    hold on
    daspect([1 1 1])
end

'Block 2 tested'
%%      Final loops
%%          Generate Block 2

'Start generating Block 2'
name_output= [name_dir,'/Bl2.dat']
fileID = fopen(name_output,'w');
%fprintf(fileID,'%d %d\n',size(C_tot,1)-2,Ny1+Ny2-1);
fprintf(fileID,'%d %d\n',size(C_tot,1),Ny1+Ny2-1);
%for i=2:1:size(C_tot,1)-1
for i=1:1:size(C_tot,1)
    sig_phi=(Ct_tot(i,2)+1e-20)/(abs(Ct_tot(i,2))+1e-20);
    if(Ct_tot(i,2)==0)
        phi_tot(i)=pi/2;
    end
    smax=((C_tot(i,1)-Ct_tot(i,1))^2+(C_tot(i,2)-Ct_tot(i,2))^2)^0.5;

    s=linspace(0,smax,N_discr);
    
    f_pol=Poly6_sn(N_discr,0,N_discr,0,1,0,0,0,0,0,'f');
    f_blend=f_pol;
    f_blend2=1-fliplr(f_pol);
    clearvars f_pol
    
    xg1=C_tot(i,1)-sig_phi*s.*sin(phi_tot(i));
    yg1=C_tot(i,2)+sig_phi*s.*cos(phi_tot(i));
    xc1=linspace(C_tot(i,1),Ct_tot(i,1),N_discr);
    yc1=linspace(C_tot(i,2),Ct_tot(i,2),N_discr);
    
    xg=xg1+f_blend2.*(xc1-xg1);
    yg=yg1+f_blend2.*(yc1-yg1);
    
    xh2=Ct_tot(i,1)-sig_phi*(-s(end)+s).*sin(phit_tot(i));
    yh2=Ct_tot(i,2)+sig_phi*(-s(end)+s).*cos(phit_tot(i));
    xc2=linspace(C_tot(i,1),Ct_tot(i,1),N_discr);
    yc2=linspace(C_tot(i,2),Ct_tot(i,2),N_discr);
    
    xh=xc2+f_blend.*(xh2-xc2);
    yh=yc2+f_blend.*(yh2-yc2);
 
    x_fin=xg+(xh-xg).*f_blend2;
    y_fin=yg+(yh-yg).*f_blend2;
    
    GL(:,1)=x_fin;
    GL(:,3)=y_fin;
        
    s_GL=calc_s(GL(:,1),GL(:,3));    

    %Version 1: Simplest version
    if Style==1
        s_J=Poly6_J(Ny1+Ny2-1,0,s_GL(end),prof_deriv1_wall(i),prof_deriv1_top(i),prof_deriv2_wall(i),prof_deriv2_top(i),prof_deriv3_wall(i),'f');
        %s_JJ=Poly6_J(Ny1+Ny2-1,0,s_GL(end),prof_deriv1_wall(i),dy_TEb,prof_deriv2_wall,prof_deriv2_top,prof_deriv3_top,'t');
    end
    %Version 2: Simple simple & no stretching at both ends (Cheaper)
    if Style==2
        s_J=Poly6_Jend(NyTOT,0,s_GL(end),prof_deriv1_wall(i),prof_deriv1_top(i),prof_deriv2_wall(i),prof_deriv2_top(i),prof_deriv2_top(i),'f');
    end
    % First draft for Version 3
    if Style==3
        s_J=Poly_J(Ny1,Ny2,0,control_dist(i),s_GL(end),prof_deriv1_wall(i),prof_deriv1_cont(i),prof_deriv1_top(i),prof_deriv2_wall(i),prof_deriv2_cont(i),prof_deriv2_top(i),'f');
    end
    %%%% end copy ---------------------------------------------------
    
    GL_inter = Interp_onto(size(s_J,2),GL,s_GL,s_J,'f');
    for j=1:size(GL_inter,1)
        fprintf(fileID,'%16.15f %16.15f %16.15f\n',GL_inter(j,1),0.0,GL_inter(j,2));
    end
    
end

fclose(fileID);
'Block 2 generated and written out'
%%          Generate Block 1 and 3 (blunt)
% Create Wake Profile at top
s_Ctwu=Poly6(Nbuffer,xTEu,Icwu(end,1),dsTEu,dsTEc,ddsTEu,ddsTEc,dddsTEc,'t');
s_Ctwl=Poly6(Nbuffer,xTEl,Icwu(end,1),dsTEl,dsTEc,ddsTEl,ddsTEc,dddsTEc,'f');

%Create Gridlines for first part
i_lim(1)=-1;

name_output1= [name_dir,'/Bl1.dat']
fileID1 = fopen(name_output1,'w');
fprintf(fileID1,'%d %d\n',Nw-1,Ny1+Ny2-1+NTEw/2-1);
fclose(fileID1);

name_output2= [name_dir,'/Bl3.dat']
fileID2 = fopen(name_output2,'w');
fprintf(fileID2,'%d %d\n',Nw-1,Ny1+Ny2-1+NTEw/2-1);
fclose(fileID2);

by2=linspace(yTE1,yTE3,N_discr);
dby2=by2(2)-by2(1)
N_discr_s=round((yTE1-yTE2)/dby2)
by3=linspace(yTE2,yTE1,N_discr_s);
dby3=by3(2)-by3(1)

clearvars f_blend
f_blend(1:10)=0;
f_pol=Poly6_s(Nbuffer-10-10+2,1,Nbuffer-10-10+2,1,0,0,0,0,0,0,'f');
f_blend(int32(10):Nbuffer-int32(10)+1)=fliplr(f_pol);
f_blend(Nbuffer-int32(10)+1:Nbuffer)=1;
bover=f_blend;
clearvars f_pol f_blend

prof_deriv1_wall(1)   = prof_deriv1_wall_l1_BU;
prof_deriv1_wall(end) = prof_deriv1_wall_u1_BU;

dsTE1=sqrt((Icwu(1,1)-Icwl(1,1))^2+(Icwu(1,2)-Icwl(1,2))^2)
dsTE2=sqrt((Icwu(end,1)-Icwl(end,1))^2+(Icwu(end,2)-Icwl(end,2))^2)

func_dsTE=Poly6(Nbuffer,1,0,0,0,0,0,0,'t');
func_dsTE=func_dsTE.^50;                     %<---- You could change here exponent to speed up transition from wall normal spacing at the TE to uniform spacing within the wake
figure
plot(func_dsTE)
for i=1:1:Nbuffer
    mfac=func_dsTE(i); % Blending factor
    prof_deriv1_wall_l1   = prof_deriv1_wall(1);
    prof_deriv1_wall_u1   = prof_deriv1_wall(end);
    s_TE                  = sqrt((Icwu(i,1)-Icwl(i,1))^2+(Icwu(i,2)-Icwl(i,2))^2);
    prof_deriv1_wall_l2   = s_TE/(NTEw-1);
    prof_deriv1_wall(1)   = prof_deriv1_wall_l1*mfac+prof_deriv1_wall_l2*(1-mfac);
    prof_deriv1_wall(end) = prof_deriv1_wall_u1*mfac+prof_deriv1_wall_l2*(1-mfac);

    %----------------------------------------------------------------------
    % GL1 - shape of gridline (upper part)
    sig_phi=1;

    smax=((Icwu(i,1)-s_Ctwu(1,i))^2+(Icwu(i,2)-yTE3)^2)^0.5;

    s=linspace(0,smax,N_discr);
    
    f_pol=Poly6_sn(N_discr,0,N_discr,0,1,0,0,0,0,0,'f');
    f_blend=f_pol;
    f_blend2=1-fliplr(f_pol);
    clearvars f_pol
    
    xg1=Icwu(i,1)-sig_phi*s.*sin(phi_wu(i));
    yg1=Icwu(i,2)+sig_phi*s.*cos(phi_wu(i));    
    xc1=linspace(Icwu(i,1),s_Ctwu(1,i),N_discr);
    yc1=linspace(Icwu(i,2),yTE3,N_discr);
    xg=xg1+f_blend2.*(xc1-xg1);
    yg=yg1+f_blend2.*(yc1-yg1);  
    
    xh2=s_Ctwu(1,i)-sig_phi*(-s(end)+s).*sin(phit_tot(end));
    yh2=yTE3+sig_phi*(-s(end)+s).*cos(phit_tot(end));
    xh=xc1+f_blend.*(xh2-xc1);
    yh=yc1+f_blend.*(yh2-yc1);    

    x_fin=xg+(xh-xg).*f_blend2;
    y_fin=yg+(yh-yg).*f_blend2;
   
    GL1(:,1)=x_fin;
    GL1(:,2)=y_fin;

    %----------------------------------------------------------------------
    % GL2 - shape of gridline (lower part)
    sig_phi=-1;
    
    smax=((Icwl(i,1)-s_Ctwl(1,i))^2+(Icwl(i,2)-yTE4)^2)^0.5;
 
    s=linspace(0,smax,N_discr);
    
    xg1=Icwl(i,1)-sig_phi*s.*sin(phi_wl(i));
    yg1=Icwl(i,2)+sig_phi*s.*cos(phi_wl(i));
    xc1=linspace(Icwl(i,1),s_Ctwl(1,i),N_discr);
    yc1=linspace(Icwl(i,2),yTE4,N_discr);
    xg=xg1+f_blend2.*(xc1-xg1);
    yg=yg1+f_blend2.*(yc1-yg1);
    
    xh2=s_Ctwl(1,i)-sig_phi*(-s(end)+s).*sin(phit_tot(1));
    yh2=yTE4+sig_phi*(-s(end)+s).*cos(phit_tot(1));
    xh=xc1+f_blend.*(xh2-xc1);
    yh=yc1+f_blend.*(yh2-yc1);
    
    x_fin=xg+(xh-xg).*f_blend2;
    y_fin=yg+(yh-yg).*f_blend2;

    GL2(:,1)=x_fin;
    GL2(:,2)=y_fin;

    %----------------------------------------------------------------------
    %Compute Spacing of upper part
    s_GL1=calc_s(GL1(:,1),GL1(:,2));
            
    %Version 1: Simplest version
    if Style==1
        s_J1=Poly6_J(Ny1+Ny2-1,0,s_GL1(end),prof_deriv1_wall(end),prof_deriv1_top(end),prof_deriv2_wall(end),prof_deriv2_top(end),prof_deriv3_wall(end),'f');
    end
    %Version 2: Simple simple & no stretching at both ends (Cheaper)
    if Style==2
        s_J1=Poly6_Jend(NyTOT,0,s_GL1(end),prof_deriv1_wall(end),prof_deriv1_top(end),prof_deriv2_wall(end),prof_deriv2_top(end),prof_deriv2_top(end),'f');
    end
    % First draft for Version 3
    if Style==3
        s_J1=Poly_J(Ny1,Ny2,0,control_dist(end),s_GL1(end),prof_deriv1_wall(end),prof_deriv1_cont(end),prof_deriv1_top(end),prof_deriv2_wall(end),prof_deriv2_cont(end),prof_deriv2_top(end),'f');
    end
    
    %Allocate Points along Gridlines
    GL1(:,3)=GL1(:,2);
    GL_inter1 = Interp_onto(size(s_J1,2),GL1,s_GL1,s_J1,'f');

    %----------------------------------------------------------------------
    %Compute Spacing of lower
    s_GL2=calc_s(GL2(:,1),GL2(:,2));

    %Version 1: Simplest version
    if Style==1
        s_J2=Poly6_J(Ny1+Ny2-1,0,s_GL2(end),prof_deriv1_wall(1),prof_deriv1_top(1),prof_deriv2_wall(1),prof_deriv2_top(1),prof_deriv3_wall(1),'f');
    end
    %Version 2: Simple simple & no stretching at both ends (Cheaper)
    if Style==2
        s_J2=Poly6_Jend(NyTOT,0,s_GL2(end),prof_deriv1_wall(1),prof_deriv1_top(1),prof_deriv2_wall(1),prof_deriv2_top(1),prof_deriv2_top(1),'f');
    end
    % First draft for Version 3
    if Style==3
        s_J2=Poly_J(Ny1,Ny2,0,control_dist(1),s_GL2(end),prof_deriv1_wall(1),prof_deriv1_cont(1),prof_deriv1_top(1),prof_deriv2_wall(1),prof_deriv2_cont(1),prof_deriv2_top(1),'f');
    end
    
    %Allocate Points along Gridlines
    GL2(:,3)=GL2(:,2);
    GL_inter2 = Interp_onto(size(s_J2,2),GL2,s_GL2,s_J2,'f');
    
    %Calculate first spacing and blend from constant to polynomial distribution
    dx3_u=(GL_inter1(2,1)-GL_inter1(1,1));
    dy3_u=(GL_inter1(2,2)-GL_inter1(1,2));
    ds3_u=sqrt(dx3_u^2+dy3_u^2);
    x31=linspace(Icwu(i,1),Icwu(i,1)+(Ny1+Ny2-1-1)*dx3_u,Ny1+Ny2-1);
    y31=linspace(Icwu(i,2),Icwu(i,2)+(Ny1+Ny2-1-1)*dy3_u,Ny1+Ny2-1);
    
    GL_inter1_(:,1)=x31;
    GL_inter1_(:,2)=y31;
    GL_inter1_(:,3)=y31;
    
    dx3_l=(GL_inter2(1,1)-GL_inter2(2,1));
    dy3_l=(GL_inter2(1,2)-GL_inter2(2,2));
    ds3_l=sqrt(dx3_l^2+dy3_l^2);
    x32=linspace(Icwl(i,1),Icwl(i,1)-(Ny1+Ny2-1-1)*dx3_l,Ny1+Ny2-1);
    y32=linspace(Icwl(i,2),Icwl(i,2)-(Ny1+Ny2-1-1)*dy3_l,Ny1+Ny2-1);
    
    GL_inter2_(:,1)=x32;
    GL_inter2_(:,2)=y32;
    GL_inter2_(:,3)=y32;
    
    %----------------------------------------------------------------------
    %Create connetction
    NTEwD=NTEw*2;
    x3_u=linspace(Icwu(i,1)-(NTEw-1)*dx3_u,Icwu(i,1),NTEwD);
    y3_u=linspace(Icwu(i,2)-(NTEw-1)*dy3_u,Icwu(i,2),NTEwD);
    x3_l=linspace(Icwl(i,1),Icwl(i,1)+(NTEw-1)*dx3_l,NTEwD);
    y3_l=linspace(Icwl(i,2),Icwl(i,2)+(NTEw-1)*dy3_l,NTEwD);
    x3_m=linspace(Icwl(i,1),Icwu(i,1),NTEwD);
    y3_m=linspace(Icwl(i,2),Icwu(i,2),NTEwD);
    s3_m=calc_s(x3_m',y3_m');
    
    clearvars f_blend
    nc=3;
    f_blend(1:nc)=0.0;
    f_pol=Poly6_s(NTEwD-nc,0.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,'f');
    f_blend(nc+1:size(f_pol,2)+nc-1)=f_pol(2:end);
    f_blend(size(f_pol,2)+nc:NTEwD)=1.0;
    clearvars f_pol nc    
   
    x3_i1=x3_l+f_blend.*(x3_m-x3_l);
    y3_i1=y3_l+f_blend.*(y3_m-y3_l);
    
    x3_i2=fliplr(x3_u)+f_blend.*(fliplr(x3_m)-fliplr(x3_u));
    y3_i2=fliplr(y3_u)+f_blend.*(fliplr(y3_m)-fliplr(y3_u));   
    
    x3_i2=fliplr(x3_i2);
    y3_i2=fliplr(y3_i2);

    clearvars f_blend
    nc=1;
    f_blend(1:nc)=0.0;
    f_pol=Poly6_s(NTEwD-nc,0.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,'f');
    f_blend(nc+1:size(f_pol,2)+nc-1)=f_pol(2:end);
    f_blend(size(f_pol,2)+nc:NTEwD)=1.0;
    clearvars f_pol nc
    x3_if=x3_i1+f_blend.*(x3_i2-x3_i1);
    y3_if=y3_i1+f_blend.*(y3_i2-y3_i1);
    
    GL3(:,1)=x3_if;
    GL3(:,2)=y3_if;
    GL3(:,3)=y3_if;
    s_GL3=calc_s(GL3(:,1),GL3(:,2));

    s_3if=Poly6_s(NTEw,0,NTEw,0,s_GL3(end),prof_deriv1_wall(1),prof_deriv1_wall(end),prof_deriv2_wall(1),prof_deriv2_wall(end),0,'f');
    plot(deriv(s_3if',1))
    
    GL_inter3 = Interp_onto(size(s_3if,2),GL3,s_GL3,s_3if,'f');

    fileID2 = fopen(name_output2,'a');
    for j=NTEw/2+1:size(GL_inter3,1)
        fprintf(fileID2,'%16.15f %16.15f %16.15f\n',GL_inter3(j,1),0.0,GL_inter3(j,2));
    end
    for j=2:size(GL_inter1,1)
        fprintf(fileID2,'%16.15f %16.15f %16.15f\n',GL_inter1(j,1),0.0,GL_inter1(j,2));
    end
    fclose(fileID2);
    
    fileID1 = fopen(name_output1,'a'); 
    for j=NTEw/2:-1:1
        fprintf(fileID1,'%16.15f %16.15f %16.15f\n',GL_inter3(j,1),0.0,GL_inter3(j,2));
    end
    for j=2:size(GL_inter2,1)
        fprintf(fileID1,'%16.15f %16.15f %16.15f\n',GL_inter2(j,1),0.0,GL_inter2(j,2));
    end
    fclose(fileID1);

end

GL_inter_Bl3(1:NTEw/2,1:2)=GL_inter3(NTEw/2+1:size(GL_inter3,1),1:2);
GL_inter_Bl3(NTEw/2:NTEw/2+Ny1+Ny2-1-1,1:2)=GL_inter1(1:end,1:2);
GL_inter_Bl1_(1:NTEw/2,1:2)=GL_inter3(1:NTEw/2,1:2);
GL_inter_Bl1=flipud(GL_inter_Bl1_);
GL_inter_Bl1(NTEw/2+1:NTEw/2+Ny1+Ny2-1-1,1:2)=GL_inter2(2:end,1:2);

GL_inter_Bl13=[flipud(GL_inter_Bl1);GL_inter_Bl3];

if s_Cwu2_end=='t'
    s_Cwu2=Poly6_end(Nw-Nbuffer,Icwu(end,1),xOut,dsTEc,dsOut,ddsTEc,ddsOut,0,'t');
else
    s_Cwu2=Poly6(Nw-Nbuffer,Icwu(end,1),xOut,dsTEc,dsOut,ddsTEc,ddsOut,0,'t');
end

for i=2:size(s_Cwu2,2)
    fileID2 = fopen(name_output2,'a');
    for j=1:size(GL_inter_Bl3,1)
        fprintf(fileID1,'%16.15f %16.15f %16.15f\n',s_Cwu2(1,i),0.0,GL_inter_Bl3(j,2));
    end
    fclose(fileID2);
    
    fileID1 = fopen(name_output1,'a');
    for j=1:size(GL_inter_Bl1,1)
        fprintf(fileID1,'%16.15f %16.15f %16.15f\n',s_Cwu2(1,i),0.0,GL_inter_Bl1(j,2));
    end
    fclose(fileID1);
end

daspect([1 1 1])
'Block 1&3 generated and written out'
stop
%%          Generate Block 1 and 3 (sharp)
% Create Wake Profile at top
StretchW=StretchW_;
s_Ctwu=Poly6(Nbuffer,xTEu,Icwu(end,1),dsTEu,dsTEc,ddsTEu,ddsTEc,dddsTEc,'f');
plot(s_Ctwu,deriv(s_Ctwu',1))
s_Ctwl=s_Ctwu;

prof_deriv1_wall(1)   = prof_deriv1_wall_l1_BU;
prof_deriv1_wall(end) = prof_deriv1_wall_u1_BU;
dy_Cwu=Poly6(size(s_Ctwu,2),prof_deriv1_wall(1),prof_deriv1_wall(1)*StretchW,0,0,0,0,0,'f');

%Create Gridlines for first part
i_lim(1)=-1;

name_output1= [name_dir,'/Bl1.dat']
fileID1 = fopen(name_output1,'w');

fprintf(fileID1,'%d %d\n',Nw-1,Ny1+Ny2-1);
fclose(fileID1);

name_output2= [name_dir,'/Bl3.dat']
fileID2 = fopen(name_output2,'w');
fprintf(fileID2,'%d %d\n',Nw-1,Ny1+Ny2-1);
fclose(fileID2);

by2=linspace(yTE1,yTE3,N_discr);
dby2=by2(2)-by2(1)
N_discr_s=round((yTE1-yTE2)/dby2)

clearvars f_blend
f_blend(1:10)=0;
f_pol=Poly6_s(Nbuffer-10-10+2,1,Nbuffer-10-10+2,1,0,0,0,0,0,0,'f');
f_blend(int32(10):Nbuffer-int32(10)+1)=fliplr(f_pol);
f_blend(Nbuffer-int32(10)+1:Nbuffer)=1;
bover=f_blend;
clearvars f_pol f_blend

prof_deriv1_wall(1)   = prof_deriv1_wall_l1_BU;
prof_deriv1_wall(end) = prof_deriv1_wall_u1_BU;
dsTE(1:2)=dsTEu;
for i=1:1:Nbuffer
    mfac=Poly6(size(s_Ctwu,2),1,0,0,0,0,0,0,'f');

    prof_deriv1_wall_l1   = prof_deriv1_wall(1);
    prof_deriv1_wall_u1   = prof_deriv1_wall(end);
%
    prof_deriv1_wall_l2   = dy_Cwu(i);
    prof_deriv1_wall_u2   = dy_Cwu(i);
%
    prof_deriv1_wall(1)   = prof_deriv1_wall_l1*mfac(i)+prof_deriv1_wall_l2*(1-mfac(i));
    prof_deriv1_wall(end) = prof_deriv1_wall_u1*mfac(i)+prof_deriv1_wall_l2*(1-mfac(i));
%
    % GL1 - upper part
    sig_phi=1;

    smax=((Icwu(i,1)-s_Ctwu(1,i))^2+(Icwu(i,2)-yTE3)^2)^0.5;

    s=linspace(0,smax,N_discr);
    
    f_pol=Poly6_sn(N_discr,0,N_discr,0,1,0,0,0,0,0,'f');
    f_blend=f_pol;
    f_blend2=1-fliplr(f_pol);
    clearvars f_pol
    
    xg1=Icwu(i,1)-sig_phi*s.*sin(phi_wu(i));
    yg1=Icwu(i,2)+sig_phi*s.*cos(phi_wu(i));    
    xc1=linspace(Icwu(i,1),s_Ctwu(1,i),N_discr);
    yc1=linspace(Icwu(i,2),yTE3,N_discr);
    xg=xg1+f_blend2.*(xc1-xg1);
    yg=yg1+f_blend2.*(yc1-yg1);  
    
    xh2=s_Ctwu(1,i)-sig_phi*(-s(end)+s).*sin(phit_tot(end));
    yh2=yTE3+sig_phi*(-s(end)+s).*cos(phit_tot(end));
    xh=xc1+f_blend.*(xh2-xc1);
    yh=yc1+f_blend.*(yh2-yc1);    

    x_fin=xg+(xh-xg).*f_blend2;
    y_fin=yg+(yh-yg).*f_blend2;
   
    GL1(:,1)=x_fin;
    GL1(:,2)=y_fin;

    % GL2 - lower part
    sig_phi=-1;
    
    smax=((Icwl(i,1)-s_Ctwl(1,i))^2+(Icwl(i,2)-yTE4)^2)^0.5;
 
    s=linspace(0,smax,N_discr);
    
    xg1=Icwl(i,1)-sig_phi*s.*sin(phi_wl(i));
    yg1=Icwl(i,2)+sig_phi*s.*cos(phi_wl(i));
    xc1=linspace(Icwl(i,1),s_Ctwl(1,i),N_discr);
    yc1=linspace(Icwl(i,2),yTE4,N_discr);
    xg=xg1+f_blend2.*(xc1-xg1);
    yg=yg1+f_blend2.*(yc1-yg1);
    
    xh2=s_Ctwl(1,i)-sig_phi*(-s(end)+s).*sin(phit_tot(1));
    yh2=yTE4+sig_phi*(-s(end)+s).*cos(phit_tot(1));
    xh=xc1+f_blend.*(xh2-xc1);
    yh=yc1+f_blend.*(yh2-yc1);
    
    x_fin=xg+(xh-xg).*f_blend2;
    y_fin=yg+(yh-yg).*f_blend2;

    GL2(:,1)=x_fin;
    GL2(:,2)=y_fin;

    %Compute Spacing upper
    s_GL1=calc_s(GL1(:,1),GL1(:,2));

    %Version 1: Simplest version
    if Style==1
        s_J1=Poly6_J(Ny1+Ny2-1,0,s_GL1(end),prof_deriv1_wall(end),prof_deriv1_top(end),prof_deriv2_wall(end),prof_deriv2_top(end),prof_deriv3_wall(end),'f');
    end
    %Version 2: Simple simple & no stretching at both ends (Cheaper)
    if Style==2
        s_J1=Poly6_Jend(NyTOT,0,s_GL1(end),prof_deriv1_wall(end),prof_deriv1_top(end),prof_deriv2_wall(end),prof_deriv2_top(end),prof_deriv2_top(end),'f');
    end
    % First draft for Version 3
    if Style==3
        s_J1=Poly_J(Ny1,Ny2,0,control_dist(end),s_GL1(end),prof_deriv1_wall(end),prof_deriv1_cont(end),prof_deriv1_top(end),prof_deriv2_wall(end),prof_deriv2_cont(end),prof_deriv2_top(end),'f');
    end
    
    
    %Allocate Points along Gridlines
    GL1(:,3)=GL1(:,2);
    GL_inter1 = Interp_onto(size(s_J1,2),GL1,s_GL1,s_J1,'f');

    %Compute Spacing lower
    s_GL2=calc_s(GL2(:,1),GL2(:,2));

    %Version 1: Simplest version
    if Style==1
        s_J2=Poly6_J(Ny1+Ny2-1,0,s_GL2(end),prof_deriv1_wall(1),prof_deriv1_top(1),prof_deriv2_wall(1),prof_deriv2_top(1),prof_deriv3_wall(1),'f');
    end
    %Version 2: Simple simple & no stretching at both ends (Cheaper)
    if Style==2
        s_J2=Poly6_Jend(NyTOT,0,s_GL2(end),prof_deriv1_wall(1),prof_deriv1_top(1),prof_deriv2_wall(1),prof_deriv2_top(1),prof_deriv2_top(1),'f');
    end
    % First draft for Version 3
    if Style==3
        s_J2=Poly_J(Ny1,Ny2,0,control_dist(1),s_GL2(end),prof_deriv1_wall(1),prof_deriv1_cont(1),prof_deriv1_top(1),prof_deriv2_wall(1),prof_deriv2_cont(1),prof_deriv2_top(1),'f');
    end
    
    %Allocate Points along Gridlines
    GL2(:,3)=GL2(:,2);
    GL_inter2 = Interp_onto(size(s_J2,2),GL2,s_GL2,s_J2,'f');
    
    %Calculate first spacing and blend from constant to polynomial distribution
    dx3_u=(GL_inter1(2,1)-GL_inter1(1,1));
    dy3_u=(GL_inter1(2,2)-GL_inter1(1,2));
    ds3_u=sqrt(dx3_u^2+dy3_u^2);
    x31=linspace(Icwu(i,1),Icwu(i,1)+(Ny1+Ny2-1-1)*dx3_u,Ny1+Ny2-1);
    y31=linspace(Icwu(i,2),Icwu(i,2)+(Ny1+Ny2-1-1)*dy3_u,Ny1+Ny2-1);
    
    GL_inter1_(:,1)=x31;
    GL_inter1_(:,2)=y31;
    GL_inter1_(:,3)=y31;
    
    dx3_l=(GL_inter2(1,1)-GL_inter2(2,1));
    dy3_l=(GL_inter2(1,2)-GL_inter2(2,2));
    ds3_l=sqrt(dx3_l^2+dy3_l^2);
    x32=linspace(Icwl(i,1),Icwl(i,1)-(Ny1+Ny2-1-1)*dx3_l,Ny1+Ny2-1);
    y32=linspace(Icwl(i,2),Icwl(i,2)-(Ny1+Ny2-1-1)*dy3_l,Ny1+Ny2-1);
    
    GL_inter2_(:,1)=x32;
    GL_inter2_(:,2)=y32;
    GL_inter2_(:,3)=y32;

    fileID2 = fopen(name_output2,'a');

    for j=1:size(GL_inter1,1)
        fprintf(fileID2,'%16.15f %16.15f %16.15f\n',GL_inter1(j,1),0.0,GL_inter1(j,2));
    end
    fclose(fileID2);
    
    fileID1 = fopen(name_output1,'a'); 

    for j=1:size(GL_inter2,1)
        fprintf(fileID1,'%16.15f %16.15f %16.15f\n',GL_inter2(j,1),0.0,GL_inter2(j,2));
    end
    fclose(fileID1);

end

GL_inter_Bl3(1:Ny1+Ny2-1,1:2)=GL_inter1(1:end,1:2);
GL_inter_Bl1(1:Ny1+Ny2-1,1:2)=GL_inter2(1:end,1:2);

GL_inter_Bl3=GL_inter1(1:end,1:2);
GL_inter_Bl1=GL_inter2(1:end,1:2);

GL_inter_Bl13=[flipud(GL_inter_Bl1);GL_inter_Bl3];

if s_Cwu2_end=='t'
    s_Cwu2=Poly6_end(Nw-Nbuffer,Icwu(end,1),xOut,dsTEc,dsOut,ddsTEc,ddsOut,0,'f');
else
    s_Cwu2=Poly6(Nw-Nbuffer,Icwu(end,1),xOut,dsTEc,dsOut,ddsTEc,ddsOut,0,'f');
end

for i=2:size(s_Cwu2,2)
    fileID2 = fopen(name_output2,'a');
    for j=1:size(GL_inter_Bl3,1)
        fprintf(fileID1,'%16.15f %16.15f %16.15f\n',s_Cwu2(1,i),0.0,GL_inter_Bl3(j,2));
    end
    fclose(fileID2);
    
    fileID1 = fopen(name_output1,'a');
    for j=1:size(GL_inter_Bl1,1)
        fprintf(fileID1,'%16.15f %16.15f %16.15f\n',s_Cwu2(1,i),0.0,GL_inter_Bl1(j,2));
    end
    fclose(fileID1);
end

daspect([1 1 1])
'Block 1&3 generated and written out'
stop
%% Finished

%% Testing (Delete first lines of Bl*.dat files)
%% Sharp TE
%% Load Grid
B1=load('Bl1.dat');
B2=load('Bl2.dat');
B3=load('Bl3.dat');
clearvars Bl13x Bl13y Blcx Blcy B2x_ B2y_ Bl13x Bl13y B1x B1y B3x B3y B1x_ B1y_ B3x_ B3y_
%global TE LE n2i n2j n1i n1j
TE=0 %Resolution of blunt Trailing edge
LE=size(Ic2,1)+size(Ic4,1)-2   %Position of Leading in Block 2 
n2i=size(C_tot,1) %-2
n2j=Ny1+Ny2-1
n1i=Nw-1
n1j=n2j%%%%%+TE/2-1
%Block 2
pos=1;
for i=1:n2i
    for j=1:n2j
        B2x(i,j)=B2(pos,1);
        B2y(i,j)=B2(pos,3);
        B2x_(i+2,j+2)=B2(pos,1);
        B2y_(i+2,j+2)=B2(pos,3);
        pos=pos+1;
    end
end

%Block 1
pos=1;
for i=1:n1i
    for j=1:n1j
        B1x(i,j)=B1(pos,1);
        B1y(i,j)=B1(pos,3);
        B1x_(i+2,j+2)=B1(pos,1);
        B1y_(i+2,j+2)=B1(pos,3);
        pos=pos+1;
    end
end
B2x_(n2i+2+2,n2j+2+2)=0;
B2y_(n2i+2+2,n2j+2+2)=0;
% % % B2x_(1,3:end-2)=B1x(end-2,TE/2:end);
% % % B2x_(2,3:end-2)=B1x(end-1,TE/2:end);
% % % B2y_(1,3:end-2)=B1y(end-2,TE/2:end);
% % % B2y_(2,3:end-2)=B1y(end-1,TE/2:end);
B2x_(1,3:end-2)=B1x(end-2,1:end);
B2x_(2,3:end-2)=B1x(end-1,1:end);
B2y_(1,3:end-2)=B1y(end-2,1:end);
B2y_(2,3:end-2)=B1y(end-1,1:end);
%Block 3
clearvars Bl13x Bl13y Blcx Blcy
Bl13x=fliplr(B1x);
Bl13y=fliplr(B1y);
pos=1
for i=1:n1i
    for j=1:n1j
        B3x(i,j)=B3(pos,1);
        B3y(i,j)=B3(pos,3);
        B3x_(i+2,j+2)=B3(pos,1);
        B3y_(i+2,j+2)=B3(pos,3);
        pos=pos+1;
    end
end
Bl13x(:,size(B1x,2)+1:size(B1x,2)+size(B3x,2))=B3x(:,:);
Bl13y(:,size(B1x,2)+1:size(B1x,2)+size(B3x,2))=B3y(:,:);
% % % Blcx(:,1:size(B2y,2))=flipud(B1x(1:end,TE/2:size(B1x,2)));
% % % Blcy(:,1:size(B2y,2))=flipud(B1y(1:end,TE/2:size(B1x,2)));
Blcx(:,1:size(B2y,2))=flipud(B1x(1:end,1:size(B1x,2)));
Blcy(:,1:size(B2y,2))=flipud(B1y(1:end,1:size(B1x,2)));
Blcx(size(B1x,1)+1:size(B1x,1)+size(B2x,1),1:end)=B2x;
Blcy(size(B1x,1)+1:size(B1x,1)+size(B2x,1),1:end)=B2y;
% Blcx(size(B1x,1)+size(B2x,1)+1:size(B1x,1)+size(B2x,1)+size(B3x,1),1:end)=B3x(1:end,TE/2:size(B3x,2));
% Blcy(size(B1x,1)+size(B2x,1)+1:size(B1x,1)+size(B2x,1)+size(B3x,1),1:end)=B3y(1:end,TE/2:size(B3x,2));
Blcx(size(B1x,1)+size(B2x,1)+1:size(B1x,1)+size(B2x,1)+size(B3x,1),1:end)=B3x(1:end,1:size(B3x,2));
Blcy(size(B1x,1)+size(B2x,1)+1:size(B1x,1)+size(B2x,1)+size(B3x,1),1:end)=B3y(1:end,1:size(B3x,2));

stop
%% Plot full Grid
figure
plot(B1x(1:end,1:end),B1y(1:end,1:end),'-r')
hold on
for n=1:n1i
plot(B1x(n,1:end),B1y(n,1:end),'-r')
end
plot(B2x(1:end,1:end),B2y(1:end,1:end),'k')
hold on
for n=1:n2i
plot(B2x(n,1:end),B2y(n,1:end),'k')
end
plot(B3x(1:end,1:end),B3y(1:end,1:end),'-b')
hold on
for n=1:n1i
plot(B3x(n,1:end),B3y(n,1:end),'-b')
end
daspect([1 1 1])

%% Blunt TE
%% Load Grid
B1=load('Bl1.dat');
B2=load('Bl2.dat');
B3=load('Bl3.dat');
%global TE LE n2i n2j n1i n1j
TE=NTEw %Resolution of blunt Trailing edge
LE=size(Ic2,1)+size(Ic4,1)-2   %Position of Leading in Block 2 
n2i=size(C_tot,1) %-2
n2j=Ny1+Ny2-1
n1i=Nw-1
n1j=n2j+TE/2-1
%Block 2
pos=1;
for i=1:n2i
    for j=1:n2j
        B2x(i,j)=B2(pos,1);
        B2y(i,j)=B2(pos,3);
        B2x_(i+2,j+2)=B2(pos,1);
        B2y_(i+2,j+2)=B2(pos,3);
        pos=pos+1;
    end
end


% fileID2 = fopen('Geo.dat','w');
%     fprintf(fileID2,'%d %d\n',1,B1y(1,1));
% for i=size(B2x,1):-5:1
%     fprintf(fileID2,'%d %d\n',B2x(i,1),B2y(i,1));
% end
%     fprintf(fileID2,'%d %d\n',1,B3y(1,1));
%     stop

%Block 1
pos=1;
for i=1:n1i
    for j=1:n1j
        B1x(i,j)=B1(pos,1);
        B1y(i,j)=B1(pos,3);
        B1x_(i+2,j+2)=B1(pos,1);
        B1y_(i+2,j+2)=B1(pos,3);
        pos=pos+1;
    end
end
B2x_(n2i+2+2,n2j+2+2)=0;
B2y_(n2i+2+2,n2j+2+2)=0;
B2x_(1,3:end-2)=B1x(end-2,TE/2:end);
B2x_(2,3:end-2)=B1x(end-1,TE/2:end);
B2y_(1,3:end-2)=B1y(end-2,TE/2:end);
B2y_(2,3:end-2)=B1y(end-1,TE/2:end);
%Block 3
clearvars Bl13x Bl13y Blcx Blcy
Bl13x=fliplr(B1x);
Bl13y=fliplr(B1y);
pos=1
for i=1:n1i
    for j=1:n1j
        B3x(i,j)=B3(pos,1);
        B3y(i,j)=B3(pos,3);
        B3x_(i+2,j+2)=B3(pos,1);
        B3y_(i+2,j+2)=B3(pos,3);
        pos=pos+1;
    end
end
Bl13x(:,size(B1x,2)+1:size(B1x,2)+size(B3x,2))=B3x(:,:);
Bl13y(:,size(B1x,2)+1:size(B1x,2)+size(B3x,2))=B3y(:,:);
Blcx(:,1:size(B2y,2))=flipud(B1x(1:end,TE/2:size(B1x,2)));
Blcy(:,1:size(B2y,2))=flipud(B1y(1:end,TE/2:size(B1x,2)));
Blcx(size(B1x,1)+1:size(B1x,1)+size(B2x,1),1:end)=B2x;
Blcy(size(B1x,1)+1:size(B1x,1)+size(B2x,1),1:end)=B2y;
Blcx(size(B1x,1)+size(B2x,1)+1:size(B1x,1)+size(B2x,1)+size(B3x,1),1:end)=B3x(1:end,TE/2:size(B3x,2));
Blcy(size(B1x,1)+size(B2x,1)+1:size(B1x,1)+size(B2x,1)+size(B3x,1),1:end)=B3y(1:end,TE/2:size(B3x,2));
stop
%% Plot full Grid
figure
plot(B1x(1:end,1:end),B1y(1:end,1:end),'-r')
hold on
for n=1:n1i
plot(B1x(n,1:end),B1y(n,1:end),'-r')
end
plot(B2x(1:end,1:end),B2y(1:end,1:end),'k')
hold on
for n=1:n2i
plot(B2x(n,1:end),B2y(n,1:end),'k')
end
plot(B3x(1:end,1:end),B3y(1:end,1:end),'-b')
hold on
for n=1:n1i
plot(B3x(n,1:end),B3y(n,1:end),'-b')
end
daspect([1 1 1])

%% Processor calculator
% This script suggests numbers of processors per block
%--------------------------------------------------------------------------
Node_in     =   60 % Nodes preferred 
%                    --> results will be filtered accordingly in output_new
npp         =   16 % Processors per node
%--------------------------------------------------------------------------
Node_range  =   round(Node_in*0.25); 
close all
% Block 1 
i1=(Nw-1);                  % Number of Points in xi
j1=(Ny1+Ny2-1+NTEw/2-1);    % Number of Points in eta
% Block 2
i2=(size(C_tot,1)-2);       % Number of Points in xi
j2=(Ny1+Ny2-1);             % Number of Points in eta
% Block 3
i3=(Nw-1);                  % Number of Points in eta
j3=(Ny1+Ny2-1+NTEw/2-1);    % Number of Points in eta

sum=2*j1*i1+j2*i2

% Start calculation
c3=i3*j3*k;
c1=i1*j1*k;
c2=i2*j2*k;

ctot=c1+c2+c3;
proc_min=(Node_in-Node_range)*npp;
proc_max=(Node_in+Node_range)*npp;
for i_p=proc_min:proc_max
    i_p;
    p1=i_p/ctot*c1;
    nj1=sqrt(p1/(i1/j1));
    ni1=(i1/j1)*nj1;
    
    p2=i_p/ctot*c2;
    nj2=sqrt(p1/(i2/j2));
    ni2=(i2/j2)*nj2;
    
    p3=i_p/ctot*c3;
    nj3=sqrt(p3/(i3/j3));
    ni3=(i3/j3)*nj3;
    
    ni1_r=round(ni1);
    nj1_r=round(nj1);
    ni2_r=round(ni2);
    nj2_r=round(nj2);
    entry=1;
    
    for ni1_a=ni1_r-10:ni1_r+10
        for ni2_a=ni1_r-10:ni2_r+10
            for nj2_a=nj1_r-10:nj2_r+10
                for nj1_a=nj1_r-10:nj1_r+10
                    if ni1_a*ni2_a*nj2_a*nj1_a>0 && (ni2_a*nj2_a+2*ni1_a*nj1_a)/npp>Node_in-Node_range
                        if mod((ni2_a*nj2_a+2*ni1_a*nj1_a),npp)==0
                            output(entry,1)=(ni2_a*nj2_a+2*ni1_a*nj1_a)/npp;    %1  Number of nodes
                            output(entry,2)=(ni2_a*nj2_a+2*ni1_a*nj1_a);        %2  Number of processors
                            output(entry,3)=(i1/ni1_a)/(j1/nj1_a);              %3  Aspect ratio B1 
                            output(entry,4)=(i2/ni2_a)/(j2/nj2_a);              %4  Aspect ratio B2
                            output(entry,5)=c1/(ni1_a*nj1_a);                   %5  Cells per processor B1
                            output(entry,6)=c2/(ni2_a*nj2_a);                   %6  Cells per processor B2
                            output(entry,7)=ni1_a;                              %7  Processors in I of Block 1
                            output(entry,8)=nj1_a;                              %8  Processors in J of Block 1
                            output(entry,9)=ni2_a;                              %9  Processors in I of Block 2
                            output(entry,10)=nj2_a;                             %10 Processors in J of Block 2
                            output(entry,11)=abs( ((c2/(ni2_a*nj2_a))-(c1/(ni1_a*nj1_a))));% / (c2/(ni2_a*nj2_a)) );           %11 devition of cells per processor
                            output(entry,12)=(abs((i1/ni1_a)/(j1/nj1_a)-1)*c1+abs((i2/ni2_a)/(j2/nj2_a)-1)*c2)/(c1+c2);     %12 averaged deviation from ideal aspect ratio aspect ratio
                            output(entry,13)=(i1/ni1_a);
                            output(entry,14)=(j1/nj1_a);
                            output(entry,15)=(i2/ni2_a);
                            output(entry,16)=(j2/nj2_a);
                            entry=entry+1;
                        end
                    end
                end
            end
        end
    end
end
output_new = output; %output((output(:,1)==Node_in),:);
output_new = sortrows(output_new,12);
output_new = sortrows(output_new,11);
'done'
%output = sort(output,1,'descend');
