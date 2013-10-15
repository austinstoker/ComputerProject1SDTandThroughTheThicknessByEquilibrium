function [ make_b1 ] = make_b1(S_hat)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
make_b1=S_hat(:,:,2,3).*S_hat(:,:,1,3)-S_hat(:,:,1,2).*S_hat(:,:,3,3);
end

