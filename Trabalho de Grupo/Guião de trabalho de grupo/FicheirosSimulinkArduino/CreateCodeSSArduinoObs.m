function [file_name] = CreateCodeSSArduino(fileName,SSDisc,Kr,n,m)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


fileID = fopen('fileName','w');
%fprintf(fileID,'%6s %12s\n','x','exp(x)');

fprintf(fileID,'%s\n\n', '//Definição de nº de estados e entradas');

fprintf(fileID,'#define STATES %d  // Number of states\n',n);
fprintf(fileID,'#define INPUTS %d  // Number of inputs\n',m);

fprintf(fileID,'%s\n\n', '//Definição da Matriz A');
 fprintf(fileID,'float A[STATES][STATES] = {\n');

for i=1:n
    fprintf(fileID,'{');
   for j=1:n-1
      fprintf(fileID,'%0.5g,', SSDisc.A(i,j));
   end
   
   if i<n 
    fprintf(fileID,'%0.5g},\n', SSDisc.A(i,j+1));
   else
    fprintf(fileID,'%0.5g}\n', SSDisc.A(i,j+1));
   end
end
 fprintf(fileID,'};\n');


  fprintf(fileID,'\n%s\n\n', '// MatrizB -');
  fprintf(fileID,'float B[STATES][INPUTS] = {\n');
for i=1:n
    fprintf(fileID,'{');
   for j=1:m-1
      fprintf(fileID,'%0.5g,', SSDisc.B(i,j));
   end

   if i<n 
    fprintf(fileID,'%0.5g},\n', SSDisc.B(i,j+1));
   else
    fprintf(fileID,'%0.5g}\n', SSDisc.B(i,j+1));
   end
 
end
 fprintf(fileID,'};\n');

 %ganhos Kr do LQR
 fprintf(fileID,'\n%s\n\n', '// Kr');
  fprintf(fileID,'float Kr[STATES] =  {');

 %float Kr[STATES] =  {0.0, 0.0, 0.2, 0.0, 0.0,0.0};
 for i=1:n
    
  if i<n 
    fprintf(fileID,'%0.5g,', Kr(i));
   else
    fprintf(fileID,'%0.5g};\n', Kr(i));
   end
 
end
%%
% for n=1:6
%       disp(sprintf('{%0.5g, %0.5g},', SSDisc.B(n,1)),SSDisc.B(n,2)));
% end

% 
% for n=1:n
%       disp(sprintf('{%0.5g, %0.5g, %0.5g, %0.5g, %0.5g, %0.5g},', SSDisc.A(n,1),SSDisc.A(n,2),SSDisc.A(n,3),SSDisc.A(n,4),SSDisc.A(n,5),SSDisc.A(n,6)));
% end
% %%
% for n=1:6
%       disp(sprintf('{%0.5g, %0.5g},', SSDisc.B(n,1)),SSDisc.B(n,2)));
%  end
fclose(fileID);

file_name = fileName;

end