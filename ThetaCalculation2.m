a=6378.1370;
b=6356.7523;
e=sqrt((a^2-b^2)/(a^2));


% phi_A=pi/8;
% phi_R=pi/7;
% 
% lambda_A=pi/4;
% lambda_R=pi/3+0.1;


phi_A=35.6762*pi/180;
phi_R=35.8676*pi/180;

lambda_A=139.6503*pi/180;
lambda_R=139.9758*pi/180;



%% 

R_OA=[(a/sqrt(1-e^2*sin(phi_A)^2))*cos(phi_A)*cos(lambda_A);
      (a/sqrt(1-e^2*sin(phi_A)^2))*cos(phi_A)*sin(lambda_A);
      (a*(1-e^2)/sqrt(1-e^2*sin(phi_A)^2))*sin(phi_A) ];

  
R_OR=[(a/sqrt(1-e^2*sin(phi_R)^2))*cos(phi_R)*cos(lambda_R);
      (a/sqrt(1-e^2*sin(phi_R)^2))*cos(phi_R)*sin(lambda_R);
      (a*(1-e^2)/sqrt(1-e^2*sin(phi_R)^2))*sin(phi_R) ];

R_ON=[0;0;b];  
  
  
R_AR=R_OR-R_OA;
R_AN=R_ON-R_OA;

n= R_OA/norm(R_OA) ;

acos(  (dot(cross(n,cross(R_AR,n)), cross(n,cross(R_AN,n)) ))/(norm(cross(n,cross(R_AR,n))) * norm(cross(n,cross(R_AN,n)) )) )*180/pi

%%

n=[2*R_OA(1)/a^2; 2*R_OA(2)/a^2; 2*R_OA(3)/b^2];

n=n/norm(n);
  
acos(  (dot(cross(n,cross(R_AR,n)), cross(n,cross(R_AN,n)) ))/(norm(cross(n,cross(R_AR,n))) * norm(cross(n,cross(R_AN,n)) )) )*180/pi


acos(  (dot(cross(R_OA,R_OR), cross(R_OA,R_ON) ))/(norm(cross(R_OA,R_OR)) * norm( cross(R_OA,R_ON))) )*180/pi





