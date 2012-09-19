function out = saveobj(obj)
%SAVEOBJ (overloaded)
warning('YALMIP objects cannot be saved in binary format. You will run into troubles if you try to load this file later. You should clear all YALMIP objects first, or avoid having any YALMIP objects in the list of variables which you save.');
out = [];
