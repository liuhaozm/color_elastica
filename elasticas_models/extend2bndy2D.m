function fo=extend2bndy2D(fi)

fo=fi;

fo(1,:)=fo(2,:);
fo(end,:)=fo(end-1,:);
fo(:,1)=fo(:,2);
fo(:,end)=fo(:,end-1);



return