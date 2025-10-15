function [ output_args ] = abs3D( input_args, dim )
%% abs3D  first dim is X Y Z default
% abs3D=sqrt((x.^2+y.^2+z.^2))
if nargin ==1
    dim=1;
end

output_args=sqrt(sum(input_args.*conj(input_args),dim));

end

