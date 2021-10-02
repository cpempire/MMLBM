% % % temp
% clc
% clear 
close all


% imagesc(LatticeC.Region);axis equal off
% hold on 

t =0:0.01:2*pi;
plot(8*cos(t),8*sin(t))
axis equal
grid on 




% x = 0:10; 
% y = 0:10;
% [gridy,gridx] = meshgrid(y,x);
% region = zeros(11);
% boundary =zeros(10,11,11);
% region(find((gridx-5).^2+(gridy-5).^2<13)) = 1;
% 
% LS2D.cx  = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];% lattice velocity -x
% LS2D.cy  = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1];% lattice velocity -y
% 
% for i = 1:11
%     for j = 1:11
% 
%         if region(i,j)==1 
%             for k = 1:9
%                 m = LS2D.cx(k);
%                 n = LS2D.cy(k);
%                 if i+m<=11&&i+m>=1&&j+n<=11&&j+n>=1
%                     if region(i+m,j+n)==0
%                         boundary(k,i,j)=1;
%                     end
%                     boundary(10,i,j) = boundary(10,i,j)+region(i+m,j+n);
%                 end
%             end
%             if boundary(10,i,j)<9
%                 boundary(10,i,j)=1;
%             else 
%                 boundary(10,i,j)=0;
%             end
%         end
%         
%     end
% end
