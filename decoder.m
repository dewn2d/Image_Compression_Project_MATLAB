function ret = decoder( file, Oimage )
    
	rows = 512;
	cols = 512;
    M = 8;
    N = 8;

	imagein = zeros(rows,cols,'uint8');
	image2 = zeros(rows,cols,'uint8');
    array =  zeros(rows*cols,2,'uint8');
    
    %disp(file);
	datafp = fopen(file, "r"); % open file of arithmatic decoded data
	%I = imshow("Images/output.512", "w");% open file for final image
    sc = fscanf(datafp, "%f "); % creat matrix of decoded data
    %sc = sc.';
    %disp(sc);
    temp = 1;
	for i = 1:size(sc)
        
            if i > 2*rows*cols
                break
            end
             
            if mod(i,2) == 0
                array(temp,2) = uint8(sc(i));
                temp=temp+1;
            else
                array(temp,1) = uint8(sc(i));
            end
            %disp(uint8(sc(i)));
    end
    arr = gen_arr(array,rows,cols);
    
    imagein = double(imagein);
	fclose(datafp);
    
	%**************decoder***************/
	for a=1:D:261442

        % Save true image of block
            curr_block = izigzag(arr(a:a+63));  
            %disp(curr_block);
            
            intemp = Iquantize(curr_block, 9,M,N); % inverse quantize matrix
        
        %disp(intemp);
            idct = IDCT(intemp,M,N); % inverse DCT trasnform matrix
            idct = idct + 128;
        
        %disp(image_mats(i,:,:));
%		
%		if (shift_rt==0)
%			printf("%d\n",i);
%			print_mat(image_mats[i]);
%		end
%		
            image2(a:a+M-1,k:k+M-1) = idct(:,:);
        
        end
    end   
	
	imwrite(image2, 'test.png' ); %write data to file
    %imshow('test.png');
    
    ret = 0;
end

function izig = izigzag( array)

    izig = zeros((8,8,'uint8')
k = 1;
    j = 1;
    i = 1;

	while(k<D)
	
		if k == 1
			izig(j,i)=array(k);
			k=1+k;
        end
    end
		
        
end

function im = gen_arr(arr,rows,cols)

    im = zeros(1,rows*cols,'uint8');
    [r,c] = size(arr);
    j=1;
    
    
    for i = 1:r
        
       if arr(i,1) == 0
           im(j) = arr(i,2);
           j=j+1;
       else
           j=j+arr(i,1);
           im(j) = arr(i,2);
       end
    end
        

end

function IQ = Iquantize( input, step_size,M,N)
	IQ = zeros(M,N,'uint8');
    for i=1:M
		for j=1:N
			IQ(i,j) = input(i,j)*step_size;
	 	end
    end
end

function output = IDCT( input,M,N)

    output = zeros(size(input));
    
	for k=0:N-1
		for  l=0:N-1
            
			temp = 0;
			for m=0:M-1
				for n =0:N-1
					temp = temp + (input(m+1,n+1)* cos(((2*k + 1)*m*pi)/(2*M))*cos(((2*l + 1)*n*pi)/(2*N))*alpha(m)*alpha(n));
                end
			end
		end
			
		output(k+1,l+1) = sqrt(2/N)*sqrt(2/M)*temp ; 
		%printf("%f, %f, %f\n",c(k)*c(l),temp, output[k][l]);

	end

end

function a = alpha( input)

	if input == 0
        a = 1/sqrt(2.0);
	else 
		a = 1;
	end
    
end

function PS = PSNR( Oimage, Nimage )
    PS = 0;
    O = imread(Oimage);
    N = imread(Nimage);
    for r = 1:rows
        for c = 1:cols
                PS = PS + O(r,c) - N(r,c)^2;
        end
    end
    
    PS = 10 * log10*((255^2)/(1/((rows+cols)/2)^2*PS));
    
end
