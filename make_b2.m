function [ make_b2 ] = make_b2(S_hat)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
make_b2=S_hat(:,:,1,2).*S_hat(:,:,2,3)-S_hat(:,:,2,2).*S_hat(:,:,1,3);
end

