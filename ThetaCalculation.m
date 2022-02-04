clc
clear all
%% Constant Data

a=6378.1370;
b=6356.7523;
e=sqrt((a^2-b^2)/(a^2));

%% Lat and Long Informations

%-------- Reference Point (Kashiwa) -------------------------------------
phi_R    = 35.8676*pi/180;
lambda_R = 139.9758*pi/180;

%----------- Rotating System -----------------------------------------

Point=3     % Select point: 1-> Tokyo, 2-> Itabashi, 3-> Adachi

if Point==1  % Tokyo
    phi_A    = 35.6762*pi/180;
    lambda_A = 139.6503*pi/180;
    Theta_R=10;
elseif Point==2 % Itabashi
    phi_A    = 35.7512*pi/180;
    lambda_A = 139.7093*pi/180;
    Theta_R=35.5;
elseif Point==3 % Adachi
    phi_A    = 35.7750*pi/180;
    lambda_A = 139.8044*pi/180;
    Theta_R=55.9;
end
%% Initialization

R_OA=[(a/sqrt(1-e^2*sin(phi_A)^2))*cos(phi_A)*cos(lambda_A);
    (a/sqrt(1-e^2*sin(phi_A)^2))*cos(phi_A)*sin(lambda_A);
    (a*(1-e^2)/sqrt(1-e^2*sin(phi_A)^2))*sin(phi_A) ];


R_OR=[(a/sqrt(1-e^2*sin(phi_R)^2))*cos(phi_R)*cos(lambda_R);
    (a/sqrt(1-e^2*sin(phi_R)^2))*cos(phi_R)*sin(lambda_R);
    (a*(1-e^2)/sqrt(1-e^2*sin(phi_R)^2))*sin(phi_R) ];

R_ON=[0;0;b];


R_AR=R_OR-R_OA;
R_AN=R_ON-R_OA;
%% Calculation of Theta

Method=3  % Select method: 1-> Spherical Earth Model, 2-> Ellipsoidal Earth Model , 3-> Ellipsoidal Earth Model (Approximation)


if Method==1  % Spherical Earth Model
    Theta=acos((cos(phi_A)*sin(phi_R)-cos(phi_R)*sin(phi_A)*cos(lambda_A-lambda_R)) /(sqrt(1- (cos(phi_A)*cos(phi_R)*cos(lambda_A-lambda_R)+sin(phi_A)*sin(phi_R))^2   ) ) )*180/pi
    Theta_N=Theta_R+Theta
    
elseif Method==2 % Ellipsoidal Earth Model (Approximation)
    n=[2*R_OA(1)/a^2; 2*R_OA(2)/a^2; 2*R_OA(3)/b^2];
    n=n/norm(n);
    Theta=acos(  (dot(cross(n,cross(R_AR,n)), cross(n,cross(R_AN,n)) ))/(norm(cross(n,cross(R_AR,n))) * norm(cross(n,cross(R_AN,n)) )) )*180/pi
    Theta_N=Theta_R+Theta
elseif Method==3 % Ellipsoidal Earth Model
    
    N=(cos(phi_A)^2+(1-e^2)^2*sin(phi_A)^2)*(cos(phi_R)^2+(1-e^2)^2*sin(phi_R)^2);
    Theta=acos(  (1-e^2)*( cos(phi_A)*sin(phi_R)- cos(phi_R)*sin(phi_A)*cos(lambda_A-lambda_R) ) / ( sqrt(N -(cos(phi_A)*cos(phi_R)*cos(lambda_A-lambda_R)+sin(phi_A)*sin(phi_R)*(1-e^2)^2)^2)  )   )*180/pi
    Theta_N=Theta_R+Theta
end






