function [xhat, SNRr, SNR, RMSE] = Bandpass(~, ~, inp, ~, ~)
input_x = inp.time_series;
xhat = bandpass(input_x,[0.6,35],250);

%KF Quality
% SNR wrt to raw
SNRr = snr(input_x,abs(input_x-xhat));
% SNR wrt to filtered signal
SNR = snr(xhat,abs(input_x-xhat));
% MSE -  Mean Square Error
RMSE = sqrt(mean(input_x-xhat).^2);
end