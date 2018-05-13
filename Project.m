%--------------main program-------------

%encoder('Images/lena.512');
%I=imread('Images/lena.512');
%disp(I)

image_file = 'Images/barb.512';

encoder( image_file);

%decoder 'Outputs/DCTStream.txt',(image_file);

%{
k= 0;
im = zeros(rows,cols,'uint8');

fileID = fopen('Images/lena.512','r');

A = fread(fileID,'*ubit8');
A = A.';
for i =1:rows
    for j =1:cols
       if (i*cols)+j > size(A)
           break
       end
       im(i,j) = uint8(A((i*cols)+j));
    end
end

fclose(fileID);
%}
%imshow(I);
%imwrite(im,'test.png')
%imshow('test.png')