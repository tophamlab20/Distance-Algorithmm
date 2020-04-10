function [dx,dy,dz,dt,data5d]=make5d(file)
%addpath /Users/jcarrol5/Downloads/bfmatlab
data=bfopen(file);

data{1,4}
data{1,4}.getPixelsSizeX(0)
data{1,4}.getPixelsSizeX(0).getValue()

%% We can find out the number of x, y, z, c, and t pixels using
mx=data{1,4}.getPixelsSizeX(0).getValue();
my=data{1,4}.getPixelsSizeY(0).getValue();
mz=data{1,4}.getPixelsSizeZ(0).getValue();
mc=data{1,4}.getPixelsSizeC(0).getValue();
mt=data{1,4}.getPixelsSizeT(0).getValue();

%% And store the phyiscal size of each
dx=data{1,4}.getPixelsPhysicalSizeX(0).value.double()
dy=data{1,4}.getPixelsPhysicalSizeY(0).value.double()
dz=data{1,4}.getPixelsPhysicalSizeZ(0).value.double()
dt=data{1,4}.getPixelsTimeIncrement(0).value.double()

%% And get the pixel dimension order
order=data{1,4}.getPixelsDimensionOrder(0).getValue.toCharArray;

%% Note that images in Matlab are displayed as if they were matrices with M rows and N columns.
% So an MxN array appears as an image that is N pixels across and M pixels
% tall. So an array that is stored in YX order will appear correctly using
% Matlab's imaging routines.

% For this reason, each image plane read from the data file and stored in data, will be an array that
% is my x mx so that when imaged in Matlab, will appear as an mx x my picture.

% So even though the pixel dimension order might be 'XYCZT', when we
% combine all of the image planes together, it will technically be YXCZT.
% So we will swap the X and Y characters in the order character array.
order([find(order=='X'), find(order=='Y')]) = 'YX';


%% We then have to figure out how to permute the data to match our desired ordering
myorder=struct('Y',1,'X',2,'Z',3,'C',4,'T',5);      % Using a struct as a dictionary mapping letters to index order
mysizes=struct('X',mx,'Y',my,'Z',mz,'C',mc,'T',mt); % Another struct mapping letters to extent of dimension
perm=arrayfun(@(x)(getfield(myorder,x)),order)';    % construct permutation array based on the order of the data and the desired order
shape=arrayfun(@(x)(getfield(mysizes,x)),order)';   % construct a shape array

%% Then we can transform the data into our 5D array
data5d=permute(reshape([data{1,1}{:,1}],shape),perm); %transform data from cell array into 5D uint16 array
end
