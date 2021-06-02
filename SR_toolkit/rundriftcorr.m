function rundriftcorr(directory, segpara, image_size, pixel_size, bin_size)

addpath(directory);
fr = readtable('FitResults.txt');
A = table2array(fr(:, [17,16, 6]));

[~, finaldrift] = RCC(A, segpara, image_size, pixel_size, bin_size, 0.2);

t = 1:size(finaldrift, 1);
to_save = array2table([t(:) finaldrift]);
to_save.Properties.VariableNames(1:3) = {'Time', 'X', 'Y'};
writetable(to_save, strcat(directory, '\', 'RCC_Drift.txt'), 'Delimiter','\t');
exit