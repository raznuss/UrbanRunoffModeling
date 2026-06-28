function CreateASCII(mat,Path_out,nanvalue,xll,yll,cellsize)
A = mat;
A(isnan(A))= nanvalue;
fid = fopen([Path_out '.asc'],'wt');
fprintf(fid, 'ncols   %.f\n', size(A,2));
fprintf(fid, 'nrows   %.f\n', size(A,1));
fprintf(fid, 'xllcorner   %.f\n', xll);
fprintf(fid, 'yllcorner   %.f\n', yll);
fprintf(fid, 'cellsize   %.3f\n', cellsize);
fprintf(fid, 'NODATA_value  %.f\n', nanvalue);


for ii = 1:size(A,1)
    fprintf(fid,'%g\t' ,A(ii,:));
%     fprintf(fid,'%.1f\t', A(ii, :));
    fprintf(fid,'\n');
end
fclose(fid);
end


%% DEM.refmat(3,2)-(0.05*size(DEM.Z,1))

