function diff_im = combinedAnisodiff(im, num_iter, delta_t, tao,sigmU,sigmL)

% Example usage:
% diffIm = combinedAnisodiff(image,322,1/4,5,2,0.5);
% imshow(diffIm,[]);

% I combined dynamic and ramp anisotrophic diffusion filter.


% c = 1/(1+(|grad(u)|/K)^alfa(x))
% c : coefficient function
% u : image
% K : calculated as in rampAnisodiff()
% alfa() : calculated as in dynamicAnisodiff();

% Matlab Implementation is based on
% Daniel Simoes Lopes's anisodiff2D().
% link: goo.gl/bDajuP

im = double(im);

% PDE (partial differential equation) initial condition.
diff_im = im;

% Center pixel distances.
dx = 1;
dy = 1;

% 2D convolution masks - finite differences.
hN = [0 1 0; 0 -1 0; 0 0 0];
hS = [0 0 0; 0 -1 0; 0 1 0];
hE = [0 0 0; 0 -1 1; 0 0 0];
hW = [0 0 0; 1 -1 0; 0 0 0];

cons1 = (2-sqrt(2))/4;
cons2 = (sqrt(2)-1)/2;
Udx = [cons1, cons2, cons1; 0, 0, 0; -cons1, -cons2, -cons1];
Udy = [cons1, 0, -cons1; cons2, 0, -cons2; cons1, 0, -cons1];


[row, colm] = size(diff_im);
uZeta = zeros(row,colm);
uNormal = zeros(row,colm);

% sigmU = 2;
% sigmL = 0.5;
sigm0 = sigmU;

% Anisotropic diffusion.
for t = 1:num_iter
        
    ux = imfilter(diff_im,Udx,'conv');
    uy = imfilter(diff_im,Udy,'conv');

    multXY = ux.*uy;
    
    for i = 1:row
        for j=1:colm
            
            %%%%%%%%%%%
            % nw n ne %
            % w  o e  %
            % sw s se %
            %%%%%%%%%%%
            uo = diff_im(i,j);
            if (j-1>=1)
                uw = diff_im(i,j-1);
            else
                uw = uo;
            end
            if (j+1<=colm)
                ue = diff_im(i,j+1);
            else
                ue = uo;
            end
            if (i-1>=1)
                un = diff_im(i-1,j);
            else
                un = uo;
            end
            if (i+1<=row)
                us = diff_im(i+1,j);
            else
                us = uo;
            end
            if (i-1>=1 && j-1>=1)
                unw = diff_im(i-1,j-1);
            elseif (i-1>=1 && j-1<1)
                unw = diff_im(i-1,j);
            elseif (i-1<1 && j-1>=1)
                unw = diff_im(i,j-1);
            else
                unw = uo;
            end
            if (i+1<=row && j-1>=1)
                usw = diff_im(i+1,j-1);
            elseif (i+1<=row && j-1<1)
                usw = diff_im(i+1,j);
            elseif (i+1>row && j-1>=1)
                usw = diff_im(i,j-1);
            else
                usw = uo;
            end
            if (i-1>=1 && j+1<=colm)
                une = diff_im(i-1,j+1);
            elseif (i-1>=1 && j+1>colm)
                une = diff_im(i-1,j);
            elseif (i-1<1 && j+1<=colm)
                une = diff_im(i,j+1);
            else
                une = uo;
            end
            if (i+1<=row && j+1<=colm)
                use = diff_im(i+1,j+1);
            elseif (i+1<=row && j+1>colm)
                use = diff_im(i+1,j);
            elseif (i+1>row && j+1<=colm)
                use = diff_im(i,j+1);
            else
                use = uo;
            end
            
            x = ux(i,j);
            y = uy(i,j);
            if(multXY(i,j)>0)                
                lamda = [-2*x^2-2*y^2+2*x*y, x^2-x*y, y^2-x*y, x*y, 0];
            else
                lamda = [-2*x^2-2*y^2-2*x*y, x^2+x*y, y^2+x*y, 0, -x*y];
            end
            uZeta(i,j) = -4*lamda(1)*uo+ lamda(2)*(ue+uw)+ ...
                lamda(3)*(un+us) + lamda(4)*(une+usw) + lamda(5)*(unw+use);
            
            x = -1*uy(i,j);
            y = ux(i,j);
            if(multXY(i,j)<0)                
                lamda = [-2*x^2-2*y^2+2*x*y, x^2-x*y, y^2-x*y, x*y, 0];
            else
                lamda = [-2*x^2-2*y^2-2*x*y, x^2+x*y, y^2+x*y, 0, -x*y];
            end
            uNormal(i,j) = -4*lamda(1)*uo+ lamda(2)*(ue+uw)+ ...
                lamda(3)*(un+us) + lamda(4)*(une+usw) + lamda(5)*(unw+use);
        end
    end
    
    E = abs(abs(uNormal) - abs(uZeta));
    maxE = max(max(E));
    E = E/maxE;
    kappa = tao*exp(-E);
        
    
    sigm = sigm0 + t*(sigmL-sigmU)/num_iter;
    Ng = ceil(6*sigm)+1;
    G = fspecial('gaussian',[Ng,Ng],sigm);
    nIm=diff_im/255;
    gi = imfilter(nIm,G,'conv');

    giN = imfilter(gi,hN,'conv');
    giS = imfilter(gi,hS,'conv');   
    giW = imfilter(gi,hW,'conv');
    giE = imfilter(gi,hE,'conv');

    k =0.1;
    alfaN = 2 - (2/(1+k*(norm(giN,1))^2));
    alfaS = 2 - (2/(1+k*(norm(giS,1))^2));
    alfaW = 2 - (2/(1+k*(norm(giW,1))^2));
    alfaE = 2 - (2/(1+k*(norm(giE,1))^2));

    nablaN = imfilter(diff_im,hN,'conv');
    nablaS = imfilter(diff_im,hS,'conv');   
    nablaW = imfilter(diff_im,hW,'conv');
    nablaE = imfilter(diff_im,hE,'conv');

    cN = 1./(1 + (abs(nablaN)./kappa).^alfaN);
    cS = 1./(1 + (abs(nablaS)./kappa).^alfaS);
    cW = 1./(1 + (abs(nablaW)./kappa).^alfaW);
    cE = 1./(1 + (abs(nablaE)./kappa).^alfaE);

    diff_im = diff_im + ...
              delta_t*(...
              (1/(dy^2))*cN.*nablaN + (1/(dy^2))*cS.*nablaS + ...
              (1/(dx^2))*cW.*nablaW + (1/(dx^2))*cE.*nablaE );


end
diff_im = uint8(diff_im);