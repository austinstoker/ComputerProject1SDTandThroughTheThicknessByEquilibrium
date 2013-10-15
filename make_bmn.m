function [ make_bmn ] = make_bmn(Shat,b0,b1,b2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
make_bmn=Shat(:,:,1,1).*b0+Shat(:,:,1,2).*b1+Shat(:,:,1,3).*b2;
end

