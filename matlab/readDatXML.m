function d = readDatXML( datname, xmlname )
% setComputer;
% dnames = rdir('testSquareXYZ_i*.dat');
% datname = dnames(1).name;
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
    xml = xml2struct(xmlname );
    %nlparams
    d.nlParams.maxIts = str2double(xml.LVData.Cluster.Cluster{1}.DBL{1}.Val.Text);
    d.nlParams.maxTime = str2double(xml.LVData.Cluster.Cluster{1}.DBL{2}.Val.Text);
    d.nlParams.method = str2double( xml.LVData.Cluster.Cluster{1}.EW.Val.Text );
    d.nlParams.errTol = str2double(xml.LVData.Cluster.Cluster{1}.DBL{3}.Val.Text);
    d.nlParams.funTol = str2double(xml.LVData.Cluster.Cluster{1}.DBL{4}.Val.Text);
    d.nlParams.stepTol = str2double(xml.LVData.Cluster.Cluster{1}.DBL{5}.Val.Text);
    %kn0 used in square IK
    d.kn0.tx01 = str2double(xml.LVData.Cluster.Cluster{2}.DBL{1}.Val.Text);
    d.kn0.ty01 = str2double(xml.LVData.Cluster.Cluster{2}.DBL{2}.Val.Text);
    d.kn0.tz01 = str2double(xml.LVData.Cluster.Cluster{2}.DBL{3}.Val.Text);
    d.kn0.ry01 = str2double(xml.LVData.Cluster.Cluster{2}.DBL{4}.Val.Text);
    d.kn0.rz01 = str2double(xml.LVData.Cluster.Cluster{2}.DBL{5}.Val.Text);
    d.kn0.tx23 = str2double(xml.LVData.Cluster.Cluster{2}.DBL{6}.Val.Text);
    d.kn0.ry34 = str2double(xml.LVData.Cluster.Cluster{2}.DBL{7}.Val.Text);
    d.kn0.rz34 = str2double(xml.LVData.Cluster.Cluster{2}.DBL{8}.Val.Text);
    d.kn0.lCath = str2double(xml.LVData.Cluster.Cluster{2}.DBL{9}.Val.Text);
    d.kn0.ry45 = str2double(xml.LVData.Cluster.Cluster{2}.DBL{10}.Val.Text);
    d.kn0.rz45 = str2double(xml.LVData.Cluster.Cluster{2}.DBL{11}.Val.Text);
    %knPM used for invP around kn0
    d.knPM.tx01 = str2double(xml.LVData.Cluster.Cluster{3}.DBL{1}.Val.Text);
    d.knPM.ty01 = str2double(xml.LVData.Cluster.Cluster{3}.DBL{2}.Val.Text);
    d.knPM.tz01 = str2double(xml.LVData.Cluster.Cluster{3}.DBL{3}.Val.Text);
    d.knPM.ry01 = str2double(xml.LVData.Cluster.Cluster{3}.DBL{4}.Val.Text);
    d.knPM.rz01 = str2double(xml.LVData.Cluster.Cluster{3}.DBL{5}.Val.Text);
    d.knPM.tx23 = str2double(xml.LVData.Cluster.Cluster{3}.DBL{6}.Val.Text);
    d.knPM.ry34 = str2double(xml.LVData.Cluster.Cluster{3}.DBL{7}.Val.Text);
    d.knPM.rz34 = str2double(xml.LVData.Cluster.Cluster{3}.DBL{8}.Val.Text);
    d.knPM.lCath = str2double(xml.LVData.Cluster.Cluster{3}.DBL{9}.Val.Text);
    d.knPM.ry45 = str2double(xml.LVData.Cluster.Cluster{3}.DBL{10}.Val.Text);
    d.knPM.rz45 = str2double(xml.LVData.Cluster.Cluster{3}.DBL{11}.Val.Text);
    
    for i =1:5; d.qps0(i,1) = str2double(xml.LVData.Cluster.Array{1}.DBL{i}.Val.Text); end;
    for i =1:5; d.qps0dn(i,1) = str2double(xml.LVData.Cluster.Array{2}.DBL{i}.Val.Text); end;
    for i =1:5; d.qps0up(i,1) = str2double(xml.LVData.Cluster.Array{3}.DBL{i}.Val.Text); end;
    for i =1:5; d.vLim(i,1) = str2double(xml.LVData.Cluster.Array{4}.DBL{i}.Val.Text); end;
    
    d.jLim.qpUp = [str2double(xml.LVData.Cluster.Cluster{4}.Array{1}.DBL{1}.Val.Text);
              str2double(xml.LVData.Cluster.Cluster{4}.Array{1}.DBL{2}.Val.Text);
              str2double(xml.LVData.Cluster.Cluster{4}.Array{1}.DBL{3}.Val.Text); 
              str2double(xml.LVData.Cluster.Cluster{4}.Array{1}.DBL{4}.Val.Text);
              str2double(xml.LVData.Cluster.Cluster{4}.Array{1}.DBL{5}.Val.Text)];
    d.jLim.qpDn = [str2double(xml.LVData.Cluster.Cluster{4}.Array{2}.DBL{1}.Val.Text);
              str2double(xml.LVData.Cluster.Cluster{4}.Array{2}.DBL{2}.Val.Text);
              str2double(xml.LVData.Cluster.Cluster{4}.Array{2}.DBL{3}.Val.Text); 
              str2double(xml.LVData.Cluster.Cluster{4}.Array{2}.DBL{4}.Val.Text);
              str2double(xml.LVData.Cluster.Cluster{4}.Array{2}.DBL{5}.Val.Text)];
    
    %square
    nm = xml.LVData.Cluster.Cluster{5}.Name.Text;
    for i = 1:numel(xml.LVData.Cluster.Cluster{5}.DBL)
        field = xml.LVData.Cluster.Cluster{5}.DBL{i}.Name.Text;
        field = strrep(field, '[','');
        field = strrep(field, ']','');
        field = strrep(field, '/','per');
        d.(nm).(field) = str2double(xml.LVData.Cluster.Cluster{5}.DBL{i}.Val.Text);
    end
    for i = 1:numel(xml.LVData.Cluster.Cluster{5}.Array.DBL)
        d.(nm).qp0(i,1) = str2double(xml.LVData.Cluster.Cluster{5}.Array.DBL{i}.Val.Text);
    end
    
    %invP result
    for i = 1:5; d.invp.qps0(i,1) = str2double(xml.LVData.Cluster.Cluster{6}.Array.DBL{i}.Val.Text); end;
    for i = 1:11;
        field = xml.LVData.Cluster.Cluster{6}.Cluster.DBL{i}.Name.Text;
        d.invp.kn0.(field) = str2double(xml.LVData.Cluster.Cluster{6}.Cluster.DBL{i}.Val.Text);
    end
    d.invp.ret = str2double(xml.LVData.Cluster.Cluster{6}.I32.Val.Text);
    d.invp.fmin = str2double(xml.LVData.Cluster.Cluster{6}.DBL.Val.Text);
    
    d.tsms = str2double(xml.LVData.Cluster.DBL.Val.Text);
    d.col = [0,0,0];
    
end