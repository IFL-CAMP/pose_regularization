function plotPoseMatrix(U,scale)

N = size(U,3);
U2 = U;
r='r'; g='g'; b='b';

% figure;
for i=0:(N-1),
    x=U(1,4,i+1);
    y=U(2,4,i+1);
    z=U(3,4,i+1);
    p1=[x,y,z]+scale*U(1:3,1,i+1)'; 
    plot3([x,p1(1)],[y,p1(2)],[z,p1(3)],r,'LineWidth',1.5);
    hold on
    p2=[x,y,z]+scale*U(1:3,2,i+1)';
    plot3([x,p2(1)],[y,p2(2)],[z,p2(3)],g,'LineWidth',1.5);
    p3=[x,y,z]+scale*U(1:3,3,i+1)';
    plot3([x,p3(1)],[y,p3(2)],[z,p3(3)],b,'LineWidth',1.5);
end 
axis equal;
set(gca,'YTickLabel',{});set(gca,'ZTickLabel',{});
set(gca,'YTick',[0 eps]);set(gca,'ZTick',[0 eps]);
view(55,11);
