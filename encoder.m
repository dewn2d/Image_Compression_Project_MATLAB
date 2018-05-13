function en = encoder( image_file )
	
	% size of image
	rows = 512;
	cols = 512;
    M = 8;
    N = 8;
    D =(M*N);
    
	DC.size = 0;
    DC.coeff = zeros(1,D*D,'uint8');
    DC.prev_co = 0;
    
    %t = [139,144,149,153,155,155,155,155;144,151,153,156,159,156,156,156;150,155,160,163,158,156,156,156;159,161,162,160,160,159,159,159;159,160,161,162,162,155,155,155;161,161,161,161,160,157,157,157;162,162,161,163,162,157,157,157;162,162,161,161,163,158,158,158];
    %disp(t);
    %disp(DCT(t,M,N))
    
    imagefp = fopen(image_file,'r'); %open input image
	DCTStream = fopen( 'Outputs/DCTStream.txt', 'w'); % open file to store DCT coefficents DC and AC
	%fp = fopen("Outputs/mat.txt", "w");
	
    imagein = zeros(rows,cols,'uint8'); %must be type uint8 to create image
    im = fread(imagefp,'*ubit8');
    im = im.';
    
	for i =1:rows
        for j = 1:cols
            if (i*cols)+j > size(im)
                break
            end
            imagein(i,j) = uint8(im((i*cols)+j)); % form a matrix of the input values
            
        end
    end
    imagein = double(imagein) - 128;
    
    fclose(imagefp);
	%imshow(imagein);	
    image2 = zeros(size(imagein),'uint8');

	%*************encoding***************/
	%ampsizemat= zeros(D*D,M*N,3); 
   
	for a=1:M:rows
        for k=1:M:cols

        % Save true image of block
            curr_block = imagein(a:a+M-1,k:k+M-1);  
        
        
            dct = DCT(curr_block, M, N); %compute dct transform of matrix
            image2(a:a+M-1,k:k+M-1) = dct(:,:);
            quant= uni_quantizer(dct,9,M); %quantize matrix
            zscanmat = zzscan(quant,N,D); % zigzag scan matrix
            
            [arr DC] = ampsize( zscanmat,D, DC ); %get Size + amplitude representation of coefficients
            %disp(quant);
		
            for j=1:D
                
                if (arr(j,1) ~= 0) || (arr(j,2) ~= 0)
                    fprintf(DCTStream, "%.0f %.0f\n",arr(j,1),arr(j,2)); %save stream of DCT coefficents into file for the arithmatic coder
                end
                
                    %{
            if(i < 2) %printing for debugging
                disp(image_mats(i,:,:));
                disp(dct);
            
                for k=1:N*M
                    disp(zscanmat(k));
                end
                
            end
            %}
                end
			
            %	for j=0:D
            %		printf("(%.0f,%.0f)(%.0f), ",ampsizemat(i][j][0),ampsizemat[i][j][1],ampsizemat[i][j][2]);
            %	printf("\n");	
			
		
            
        end
    end
    
    imshow(image2)
	%printf("hellow\n");
	%for( i=0; i<DC.size; i++)
	%{	
		fprintf(DCTStream, "%.0f ", DC.coeff[i]);
		//printf("%d\n",i);
	%}

	%printf("%d\n",DC.size);

	fclose(DCTStream);
    en =quant;
    
end
    
function output = uni_quantizer( input, stp_size, M)

    output = zeros(M,M, 'uint8');
    
	for i=1:M
		for j=1:M
			output(i,j) = round(input(i,j)/stp_size);
			if(output(i,j) == -0)
				output(i,j) = 0;
			end
		end

	end
end

function scanArray = zzscan( input, N, D)

 	scanArray = zeros(1,D,'uint8'); 
	k = 1;
    j = 1;
    i = 1;

	while(k<D)
	
		if k == 1
			scanArray(:,k) = input(j,i);
			k=1+k;
        
		
        elseif i == N
		
			j=j+1;
			while (j<=N)
			
				scanArray(:,k) = input(j,i);
				%printf("a\t%d\t%d\t%d\t%.0f\n",k,i,j,scanArray[k]);
				k=k+1;
				if j~=N-1
					i=i-1;
					j=j+1;
				
                else
                    break;
				end
			end
			%printf("k%d\ti%d\tj%d\n\n",k,i,j);
         
		
         elseif j == N 
			
			i=i+1;
			while (i<=N)
				
				scanArray(:,k) = input(j,i);	
				%printf("b\t%d\t%d\t%d\t%.0f\n",k,i,j,scanArray[k]);
				k=k+1;
				if i~=N
					i=i+1;
					j=j-1;
				
                else
                    break;
				end
			end
			%printf("k%d\ti%d\tj%d\n\n",k,i,j);
            

        elseif j == 1
			
			i=i+1;

			while (i>=1)
				
				scanArray(:,k) = input(j,i);	
				%printf("c\t%d\t%d\t%d\t%.0f\n",k,i,j,scanArray[k]);
				k=k+1;
				if i~=1
					i=i-1;
					j=j+1;
				
                else
                    break;
				end
			end
            
			%printf("k%d\ti%d\tj%d\n\n",k,i,j);
		
		
        elseif(i == 1)
		
			j=j+1;
			while (j>=1)
			
				scanArray(:,k) = input(j,i);
				%printf("d\t%d\t%d\t%d\t%.0f\n",k,i,j,scanArray[k]);
				k=k+1;
				if j~=1
					i=j+1;
					j=j-1;
				
                else
                    break;
				end
			end
            
			%printf("k%d\ti%d\tj%d\n\n",k,i,j);
		end
	end
end

function scanArray = rasterscan( input, M, D)

	scanArray = zeros(1,D); 
    k =0;
    j =0;
	while(k<D)
		
		if mod(j,2) ==0
		
			for i =0:M
			
				scanArray(k) = input(i,j);
				k=k+1;
			end
		
        elseif mod(j,2) ==1
		
			for i =M:0
			
				
				scanArray(k) = input(i,j);
				k=k+1;
			end
		end
		j=j+1;
	end
end

function [ array, DC] = ampsize(input,D,DC)

	coeff = [2,2,3; 3,4,7; 4,8,15; 5,16,31; 6,32,63; 7,64,127; 8,128,255; 9,256,511; 10,512,1023;];
	lead0=0;
    array =  zeros(D,2,'uint8');

	for i=1:D
	
			if(input(i) == 0)
				lead0=lead0+1;

            elseif(input(i) == 1 || input(i) == -1)
			
				array(i,1) = lead0;
				array(i,2) = 1;
				lead0 = 0;
			
			else		
	
				for j=1:9
						
					if((input(i) >= coeff(j,2) && input(i) <= coeff(j,3)) ||(input(i) >= -coeff(j,3) && input(i) <= -coeff(j,2)))
					
						array(i,1) = lead0;
						array(i,2) = coeff(j,1);
						lead0 = 0;	
					end
                end
			
            end
    end
end

 function output = DCT(input,M,N)

    output = zeros(M,N,'uint8');

	for  k=0:7 
		for  l=0:7
		
			temp = 0;
			for  m=0:7			
				for n =0:7
					temp = temp + (((input(m+1,n+1))*cos(((2*m + 1)*k*pi)/16) * cos(((2*n + 1)*l*pi)/16)));
                    %disp(temp);     
                end
			end

			output(k+1,l+1) = lambda(k)*lambda(l)*(1/4)*temp;
		end
    end

 
end

function a = lambda( i)

	if( i == 0 )
		a = 1/sqrt(2);
	else
		a = 1;
	end
end
%{
function get_bitrate( void )

	bits = malloc(1000 * sizeof(char));
	bitfp = fopen( "map.art", "r" );
	int size;
	
	while(bits ~= EOF);
	
		bits =bits+ (char)fgetc(bitfp);
		size=size+1;
    end
end
	

void print_mat(float input[][N])
{
	int i, j;
	for( i=0; i<8;i++){
		for(j =0; j<8; j++){
			printf("%.1f\t", input[i][j]);
		}
		printf("\n");
	}

	printf("\n");
}
%}
