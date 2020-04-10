addpath /software/bioformats/5.3.4

[file, path] = uigetfile('*.*', 'Choose a file to open');
[dx,dy,dz,dt,data]=make5d([path,file]);
[mx,my,mz,mc,mt]=size(data);
cmap=[1,2];

%% Plot the data
data4d=uint16(squeeze(mean(data(:,:,:,:,:),1)));
bc=my_stretchlim(data4d)
playstack(data4d,bc,cmap)

%%
clf
threshold=2000;
stack=data(:,:,:,2,1)>threshold;
[xxx,yyy,zzz]=meshgrid([0:mx-1]*dx,[0:my-1]*dy,[0:mz-1]*dz);
subplot(1,3,1)
isosurface(xxx,yyy,zzz,stack,.5); title('Thresholded at 2000')
axis equal
axis ij
axis tight
%% And we can use bwareaopen to eliminate connnected regions smaller than 10 pixels
minvolume=10;
stack=bwareaopen(stack,minvolume);
subplot(1,3,2)
isosurface(xxx,yyy,zzz,bwareaopen(stack,minvolume)); title('Areas larger than 10 pixels')
axis equal
axis ij
axis tight
%% And then use bwconncomp to get the connected components
cc=bwconncomp(stack);
subplot(1,3,3)
imshow(label2rgb(max(labelmatrix(cc),[],3),'lines','k')); title('connected components');
axis ij

%% And then calculate the region properties and the centroids
props=regionprops(cc);
centroids=reshape([props(:).Centroid],3,[]);

%% And plot the results
figure(2)
clf
isosurface(xxx,yyy,zzz,stack,.5)
f=gcf()
p=f.Children(1).Children(2)
p.FaceAlpha=.2
hold on;
xcentroids=diag([dx,dy,dz])*(centroids-1);
plot3(xcentroids(1,:),xcentroids(2,:), xcentroids(3,:),'+k')
axis equal
axis ij
axis tight

%% And now we can import the data
[filecsv, pathcsv] = uigetfile('*.*', 'Choose trajectory data file to open');
trajectory_data=importfile([pathcsv, filecsv]);

%% Get unique timestamps
times=unique(trajectory_data.RelTimes(:));

%% Look at centroid positions at t=0
points=trajectory_data(trajectory_data.RelTimes==times(1),{'TrackID','ID','CentroidXm' 'CentroidYm' 'CentroidZm'});

%% And overlay them with centroids found by matlab
plot3(points.CentroidXm, points.CentroidYm, points.CentroidZm, '+')

%% Now fit surface
clf
background=mean(data(:,:,:,1,:),5);
hist(reshape(background,1,[]),1000)
%% Threshold at 6000
bwbackground=background>400;

%%
%% We can then take the first z value as the cell wall
for i=1:size(bwbackground,1)
    for j=1:size(background,2)
        k=find(bwbackground(i,j,:)>0,1);
        bwbackground(i,j,k+1:end)=0;
    end
end
clf
isosurface(xxx,yyy,zzz,bwbackground,[.9])
axis ij
axis equal

%% And then fit the data
[yy,xx,zz]=ind2sub(size(bwbackground),find(bwbackground));
xx=(xx-1)*dx;yy=(yy-1)*dy;zz=(zz-1)*dz;
[xData, yData, zData] = prepareSurfaceData( xx, yy, zz );
%% Calculate fit and plot result
ft = fittype( 'poly44' );
[fitresult, gof] = fit( [xData, yData], zData, ft );
clf
h = plot( fitresult);
hold on;
isosurface(xxx,yyy,zzz,bwbackground,[.9]) %permute(bwbackground,[2,1,3]),[.9])
axis ij
axis equal

%% Now find distance of each point from surface
[xx,yy]=meshgrid([0:mx-1]*dx,[0:my-1]*dy);
zz=fitresult(xx,yy);
for i=1:size(points,1)
    pos=[points.CentroidXm(i), points.CentroidYm(i), points.CentroidZm(i)];
    [d2,indx]=min((xx(:)-pos(1)).^2+(yy(:)-pos(2)).^2+(zz(:)-pos(3)).^2);
    intersect(i,:)=[xx(indx),yy(indx),zz(indx)];
    d(i)=sqrt(d2);
end
   
%% And plot results
hold on;
for i=1:size(points,1)
   plot3([intersect(i,1),points.CentroidXm(i)],[intersect(i,2), points.CentroidYm(i)],[intersect(i,3), points.CentroidZm(i)])
end

%% Now go through entire table and repeat analysis
for i=1:size(trajectory_data,1)
    pos=[trajectory_data.CentroidXm(i), trajectory_data.CentroidYm(i), trajectory_data.CentroidZm(i)];
    [d2,indx]=min((xx(:)-pos(1)).^2+(yy(:)-pos(2)).^2+(zz(:)-pos(3)).^2);
    intersect(i,:)=[xx(indx),yy(indx),zz(indx)];    
    d(i)=sqrt(d2);
end

%% Augment table with new data and write file
trajectory_data=[trajectory_data,array2table(intersect(:,:),'VariableNames',{'intersect_x','intersect_y','intersect_z'}),table(d(:),'VariableNames',{'Distance_to_Surface'})];
[outpath,outname,outext]=fileparts([pathcsv,filecsv]);
writetable(trajectory_data,[outpath,'/', outname,'_ext',outext]);