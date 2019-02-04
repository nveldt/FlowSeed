function show_set_iso(X,Nmap,Rfull,Rset,offset)

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

%Z = double(X(Sx,Sy,Sz));
LT = RT(Sx,Sy,Sz);
LS = RS(Sx,Sy,Sz);


clf;
data = smooth3(LT);
pdata = isosurface(data,0.5);
p = patch('Faces', pdata.faces, 'Vertices', pdata.vertices, 'FaceVertexCData', 0.5, ...
          'FaceColor', 'flat', 'FaceLighting', 'gouraud', ...
          'EdgeColor', 'none');
isonormals(data,p);
p.FaceColor = [128,128,255]/255;
%colormap(bone);
view(3); 

%assert(numel(Rset) == 0);

if numel(Rset) > 0
    p.FaceAlpha = 0.33;
    data = smooth3(LS);
    pdata = isosurface(data,0.5);
    p = patch('Faces', pdata.faces, 'Vertices', pdata.vertices, 'FaceVertexCData', 0.5, ...
          'FaceColor', 'flat', 'FaceLighting', 'gouraud', ...
          'EdgeColor', 'none');
    isonormals(data,p);
    p.FaceColor = [255,200,128]/255;
    p.FaceAlpha = 0.33;
end

camlight; lighting phong
axis off;
daspect([1 1 1]); axis tight; 
camroll(0)

