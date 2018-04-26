function dd=deriv(dnd,dim,per)
global bc4;

ds = 1.0;%np-1;
fact = ds/12;
order = [1 2 3];
sz = size(dnd);

if ~exist('per','var')
  per='false';
end
    
% Central differences
    for i=3:sz(1)-2
        dd(i,:)=(dnd(i-2,:)-dnd(i+2,:)+8.0*(dnd(i+1,:)-dnd(i-1,:)))*fact;
    end
    
    if(strcmp(per,'periodic'))
        
%         i = 1;
%         dd(i,:)=(dnd(end-1,:)-dnd(i+2,:)+8.0*(dnd(i+1,:)-dnd(end,:)))*fact;
%         i = 2;
%         dd(i,:)=(dnd(end,:)-dnd(i+2,:)+8.0*(dnd(i+1,:)-dnd(i-1,:)))*fact;
%         i = sz(1);
%         dd(i,:)=(dnd(i-2,:)-dnd(2,:)+8.0*(dnd(1,:)-dnd(i-1,:)))*fact;
%         i = sz(1)-1;
%         dd(i,:)=(dnd(i-2,:)-dnd(1,:)+8.0*(dnd(i+1,:)-dnd(i-1,:)))*fact;
%         
    else
% Carpenter at boundaries
        if(dim==2)
%             ji = 801;
%             je = 2600;
%             
%              ibM = 1;
%             for i=1:4
%                 ic = i+ibM-1;
%                 dd(ic,ji:je)=ds*(bc4(i,1)*dnd(ibM,ji:je)+bc4(i,2)*dnd(ibM+1,ji:je)+bc4(i,3)*dnd(ibM+2,ji:je)...
%                     +bc4(i,4)*dnd(ibM+3,ji:je)+bc4(i,5)*dnd(ibM+4,ji:je)+bc4(i,6)*dnd(ibM+5,ji:je));
%             end
%             ibP = sz(1);
%             for i=4:-1:1
%                 ic=ibP-i+1;
%                 dd(ic,:)=-ds*(bc4(i,1)*dnd(ibP,:)+bc4(i,2)*dnd(ibP-1,:)+bc4(i,3)*dnd(ibP-2,:)...
%                     +bc4(i,4)*dnd(ibP-3,:)+bc4(i,5)*dnd(ibP-4,:)+bc4(i,6)*dnd(ibP-5,:));
%             end
        else
            ibM = 1;
            for i=1:4
                ic = i+ibM-1;
                dd(ic,:)=ds*(bc4(i,1)*dnd(ibM,:)+bc4(i,2)*dnd(ibM+1,:)+bc4(i,3)*dnd(ibM+2,:)...
                    +bc4(i,4)*dnd(ibM+3,:)+bc4(i,5)*dnd(ibM+4,:)+bc4(i,6)*dnd(ibM+5,:));
            end
            ibP = sz(1);
            for i=4:-1:1
                ic=ibP-i+1;
                dd(ic,:)=-ds*(bc4(i,1)*dnd(ibP,:)+bc4(i,2)*dnd(ibP-1,:)+bc4(i,3)*dnd(ibP-2,:)...
                    +bc4(i,4)*dnd(ibP-3,:)+bc4(i,5)*dnd(ibP-4,:)+bc4(i,6)*dnd(ibP-5,:));
            end
        end
    
    end
% Differentiate accross the wake
%     if(dim == 2)
%         ji = 1;
%         je = 800;
%         for j=ji:je
%         dd(1,j)=(dnd(3,end-j+1)-dnd(1+2,j)+8.0*(dnd(1+1,j)-dnd(2,end-j+1)))*fact;
%         dd(2,j)=(dnd(2,end-j+1)-dnd(2+2,j)+8.0*(dnd(2+1,j)-dnd(1,j)))*fact;
%         dd(1,end-j+1)=(dnd(3,j)-dnd(1+2,end-j+1)+8.0*(dnd(1+1,end-j+1)-dnd(2,j)))*fact;
%         dd(2,end-j+1)=(dnd(2,j)-dnd(2+2,end-j+1)+8.0*(dnd(2+1,end-j+1)-dnd(1,end-j+1)))*fact;
%         end
%     end
    
    dd = permute(dd,order);
end