function [ make_b0 ] = make_b0(S_hat)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
make_b0=S_hat(:,:,2,2).*S_hat(:,:,3,3)-S_hat(:,:,2,3).*S_hat(:,:,2,3);
end

