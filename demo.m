% This package contains the code for Alternative Matting Laplacian as
% described in
%
% [Pitie16] An Alternative Matting Laplacian. F. Pitie (2016). In 
%           International Conference on Image Processing (ICCV'16)
%
% send an email to fpitie@mee.tcd.ie if you want more information

fprintf('\nAlternative Matting Laplacian demo based on:\n');
fprintf('  [Pitie16] An Alternative Matting Laplacian. F. Pitie (ICIP 2016).\n\n');

fprintf('... loading source image\n');
I = double(imread('GT04.png'))/255;

fprintf('... loading trimap\n');
trimap = double(imread('trimap-GT04.png'))/255;

fprintf('... loading initial alpha from sampling\n');
alpha0 = double(imread('alpha0-GT04.png'))/255;

fprintf('... loading initial confidence map for sampled alpha\n');
conf0 = double(imread('conf-GT04.png'))/255;

walpha0 = exp(log(conf0)*3);
walpha0(trimap < .1) = 2;
walpha0(trimap > .9) = 2;

fprintf('... solving for alpha and a for the entire picture\n');

[alpha, beta] = alternative_matting_laplacian_solver(...
    I, 'alpha0', alpha0, 'walpha0', walpha0, 'sigma_r', 1, 'T', .001);

fprintf('[ok]\n');

% display results
screensize = get(0,'ScreenSize');
sz = [576, 1024];
figure('Position', [ ceil((screensize(3)-sz(2))/2), ceil((screensize(4)-sz(1))/2), sz(2), sz(1)]);
subplot('Position',[0.01  0.4850 0.3200 .47]);
imshow(I); title('Original Image');

subplot('Position',[0.3400  0.4850 0.3200 .47]);
imshow(alpha0); title('initial alpha from sampling');

subplot('Position',[0.01 0.01 0.3200 .47]);
imshow(alpha); title('estimated alpha');

subplot('Position',[0.3400 0.01 0.3200 .47]);
imshow(beta(:,:,1:3)/5+.5); title('estimated a');

subplot('Position',[0.6700 0.01 0.3200 .47]);
imshow(beta(:,:,4)/5+.5); title('estimated b');

