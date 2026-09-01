% c0800_get_pfb_coeffs.m
% J.Ray    10/04/22
%
% Script to pull out the filter coefficients from the
% pfb_fir block in all of the c0800xNNNN designs.
%
% The following is a list of model files for the full 
% design and the pfb_fir black box:
%	for c0800x0032:
%	model = c0800x0032_x14_7_24t_095binw.slx
%	core  = codd_pfb_fir_0032ch_24t_06_core
%
%	for c0800x0064:
%	model = c0800x0064_x14_7_24t_095binw.slx
%	core  = codd_pfb_fir_0064ch_24t_095binw_core2
%
%	for c0800x0128:
%	model = c0800x0128_x14_7_24t_095binw.slx
%	core  = codd_pfb_fir_0128ch_24t_095binw_core2
%
%	for c0800x0256:
%	model = c0800x0256_x14_7_24t_095binw.slx
%	core  = codd_pfb_fir_0256ch_24t_095binw_core2
%
%	for c0800x0512:
%	model = c0800x0512_x14_7_24t_095binw.slx
%	core  = codd_pfb_fir_0512ch_24t_095binw_core
%
%	for c0800x1024:
%	model = c0800x1024_x14_7_24t_095binw.slx
%	core  = codd_pfb_fir_1024ch_24t_095binw_core2
%
%	for c0800x2048:
%	model = c0800x2048_x14_7_24t_095binw.slx
%	core  = codd_pfb_fir_2048ch_24t_095binw_core
%
%	for c0800x4096:
%	model = c0800x4096_x14_7_24t_095binw.slx
%	core  = codd_pfb_fir_4096ch_24t_095binw_core1
%
%
% Looking at the pfb_fir block in codd_pfb_fir_0032ch_24t_06_core.slx
% from pol1_coeffs_in1
% pfb_coeff_gen_calc(6, 24,'hamming',4, 0,0.94999999999999995559107901499374,1,false)
% pfb_coeff_gen_calc(6, 24,'hamming',4, 0,0.94999999999999995559107901499374,2,false)
%
% Looking through the different blocks, you can see the 5th and 7th arguments
% change depending on where in the design the block is.  So, I previously attempted
% to get these coeffs by using a nested for loop to cycle through the two arguments.
% This didn't work because I apparently don't understand the ordering.  No matter
% which argument I put in the inner loop vs outer loop, the results never looked
% correct.
%
% Consulted Dave MacMahon about this and he informed me that the pfb_coeff_gen_calc
% script has an option to output all coeffs in order, by setting the 7th argument to
% -1. So, this m-file exercises that option once for each design from 32 to 4096 channels
% and saves all the coeffs to separate files.
% 

% clear all variables
clear

% ====== c0800x0032 ======
pfb_coeffs = pfb_coeff_gen_calc(6, 24,'hamming',4, 0,0.94999999999999995559107901499374,-1,false);
% orig command:  pfb_coeff_gen_calc(6, 24,'hamming',4, 0,0.94999999999999995559107901499374,1,false)

% write coefficients to file
fileID = fopen('c0800x0032_x14_7_24t_095binw_get_pfb_coeffs.txt','w');
n_coeffs = 1536;  % 64pts * 24 taps
for z = 1:n_coeffs
    fprintf(fileID,'%15.14f\n',pfb_coeffs(z));
end
fclose(fileID);

% ====== c0800x0064 ======
pfb_coeffs = pfb_coeff_gen_calc(7, 24,'hamming',4, 0,0.94999999999999995559107901499374,-1,false);
% orig command:  pfb_coeff_gen_calc(7, 24,'hamming',4, 0,0.94999999999999995559107901499374,1,false)

% write coefficients to file
fileID = fopen('c0800x0064_x14_7_24t_095binw_get_pfb_coeffs.txt','w');
n_coeffs = 3072;  % 128pts * 24 taps
for z = 1:n_coeffs
    fprintf(fileID,'%15.14f\n',pfb_coeffs(z));
end
fclose(fileID);

% ====== c0800x0128 ======
pfb_coeffs = pfb_coeff_gen_calc(8, 24,'hamming',4, 0,0.94999999999999995559107901499374,-1,false);
% orig command:  pfb_coeff_gen_calc(8, 24,'hamming',4, 0,0.94999999999999995559107901499374,1,false)

% write coefficients to file
fileID = fopen('c0800x0128_x14_7_24t_095binw_get_pfb_coeffs.txt','w');
n_coeffs = 6144;  % 256pts * 24 taps
for z = 1:n_coeffs
    fprintf(fileID,'%15.14f\n',pfb_coeffs(z));
end
fclose(fileID);

% ====== c0800x0256 ======
pfb_coeffs = pfb_coeff_gen_calc(9, 24,'hamming',4, 0,0.94999999999999995559107901499374,-1,false);
% orig command:  pfb_coeff_gen_calc(9, 24,'hamming',4, 0,0.94999999999999995559107901499374,1,false)

% write coefficients to file
fileID = fopen('c0800x0256_x14_7_24t_095binw_get_pfb_coeffs.txt','w');
n_coeffs = 12288;  % 512pts * 24 taps
for z = 1:n_coeffs
    fprintf(fileID,'%15.14f\n',pfb_coeffs(z));
end
fclose(fileID);

% ====== c0800x0512 ======
pfb_coeffs = pfb_coeff_gen_calc(10, 24,'hamming',4, 0,0.94999999999999995559107901499374,-1,false);
% orig command:  pfb_coeff_gen_calc(10, 24,'hamming',4, 0,0.94999999999999995559107901499374,1,false)

% write coefficients to file
fileID = fopen('c0800x0512_x14_7_24t_095binw_get_pfb_coeffs.txt','w');
n_coeffs = 24576;  % 1024pts * 24 taps
for z = 1:n_coeffs
    fprintf(fileID,'%15.14f\n',pfb_coeffs(z));
end
fclose(fileID);

% ====== c0800x1024 ======
% **** NOTE THIS DESIGN ONLY HAS 12 TAPS!!!! ****
pfb_coeffs = pfb_coeff_gen_calc(11, 12,'hamming',4, 0,0.94999999999999995559107901499374,-1,false);
% orig command:  pfb_coeff_gen_calc(11, 12,'hamming',4, 0,0.94999999999999995559107901499374,1,false)

% write coefficients to file
fileID = fopen('c0800x1024_x14_7_24t_095binw_get_pfb_coeffs.txt','w');
n_coeffs = 24576;  % 2048pts * 12 taps
for z = 1:n_coeffs
    fprintf(fileID,'%15.14f\n',pfb_coeffs(z));
end
fclose(fileID);

% ====== c0800x2048 ======
pfb_coeffs = pfb_coeff_gen_calc(12, 24,'hamming',4, 0,0.94999999999999995559107901499374,-1,false);
% orig command:  pfb_coeff_gen_calc(12, 24,'hamming',4, 0,0.94999999999999995559107901499374,1,false)

% write coefficients to file
fileID = fopen('c0800x2048_x14_7_24t_095binw_get_pfb_coeffs.txt','w');
n_coeffs = 98304;  % 4096pts * 24 taps
for z = 1:n_coeffs
    fprintf(fileID,'%15.14f\n',pfb_coeffs(z));
end
fclose(fileID);

% ====== c0800x4096 ======
pfb_coeffs = pfb_coeff_gen_calc(13, 24,'hamming',4, 0,0.94999999999999995559107901499374,-1,false);
% orig command:  pfb_coeff_gen_calc(13, 24,'hamming',4, 0,0.94999999999999995559107901499374,1,false)

% write coefficients to file
fileID = fopen('c0800x4096_x14_7_24t_095binw_get_pfb_coeffs.txt','w');
n_coeffs = 196608;  % 8192pts * 24 taps
for z = 1:n_coeffs
    fprintf(fileID,'%15.14f\n',pfb_coeffs(z));
end
fclose(fileID);