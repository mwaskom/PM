function [ output_args ] = reconAllBatch( par )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
reconCommand = ['recon-all -i ' par.hiresimg ' -s ' par.substr ' -all'];
unix (reconCommand)
end

