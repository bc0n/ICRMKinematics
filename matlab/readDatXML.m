function d = readDatXML( datname, xmlname )
% setComputer;
% dnames = rdir('testSquareXYZ_i*.dat');
% datname = dnames(1).name;
% xnames = rdir('testSquareXYZ_i*.xml');
% xmlname = xnames(1).name;
    d.name = datname;

    fid = fopen(datname,'r','l');
    fseek(fid,0,'eof');
    fsize = ftell(fid);
    fseek(fid,0,'bof');

    data = fread(fid,fsize,'double'); %(1+16)*na + (1+5+3+3+1+1+1)*na + (1+5+16)*nb +3;
    fclose(fid);
    d.na = data(1); %1based
    d.nb = data(2);

    %commanded square
    a = 3; b = a + d.na*16 - 1;
    d.Hsq = permute( reshape(data(a:b),4,4,d.na),[3,2,1]);

    % converged data
    a = b+1; b = a + d.na;
    d.map = data(a:b);
    d.a2b{1} = 0 + 1:d.map(1);
    for i = 2:d.na+1;
        d.a2b{i} = d.a2b{i-1}(end) + (1:d.map(i));
    end
    a = b+1; b = a + d.na*5 - 1;
    d.qcv = reshape( data(a:b), 5, d.na)';
    a = b+1; b = a + d.na*3 - 1;
    d.xcm = reshape( data(a:b), 3, d.na)'; %IK xyz goal
    a = b+1; b = a + d.na*3 - 1;
    d.xcv = reshape( data(a:b), 3, d.na)'; %IK xyz converged
    a = b+1; b = a + d.na - 1;
    d.nrm = data(a:b);
    a = b+1; b = a + d.na - 1;
    d.ret = data(a:b);

    %measured
    a = b+1; b = a + d.nb*5 -1;
    d.qms = reshape( data(a:b), 5, d.nb)';
    a = b+1; b = a + d.nb*16 -1;
    d.Hms = permute( reshape(data(a:b),4,4,d.nb),[3,2,1]);
    
    %xml
    x = read_ni_xml_object( xmlread( xmlname ));
    d = mergestruct(d,x);
    
    d.col = rand(1,3); %random color
end
