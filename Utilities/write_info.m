function write_info(file_name, str)
    fid = fopen(file_name, 'a+');
	fprintf(fid, '%s \n', str);
    fclose(fid);