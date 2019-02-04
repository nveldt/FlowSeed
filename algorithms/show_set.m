function show_set(X,Nmap,Rfull,Rset,offset)

R = union(Rfull,Rset);
cmins = min(Nmap(:,R),[],2);
cmaxs = max(Nmap(:,R),[],2);

Sx = cmins(1)-offset:cmaxs(1)+offset;
Sy = cmins(2)-offset:cmaxs(2)+offset;
Sz = cmins(3)-offset:cmaxs(3)+offset;

RT = false(size(X));
RT(Rfull) = true;

RS = false(size(X));
RS(Rset) = true;


Z = double(X(Sx,Sy,Sz));
LT = RT(Sx,Sy,Sz);
LS = RS(Sx,Sy,Sz);

clf;
vol3d('CData',(Z/4000).*(3*LS+1.)/4,'Alpha',(50*LT+50*LS)/255)
colormap(1-gray);
alphamap(0:0.1:1)
