% empty test
clear all
close all
clc
for k=10:45
    fid=fopen('cfg2321.dat','r+');
    while(1)
        cLine=fgetl(fid);
        if strcmpi(cLine,'algorithm free parameters:'), break, end
    end
    cLine=fgetl(fid);
    cLine(32:33)=num2str(k);
    fseek(fid,-(length(cLine)+2),'cof');
    fprintf(fid,'%s',cLine);
    fclose(fid);
    pause(1)
    !./a.exe
    beep
    pause(1)
end
beep;pause(0.3)
beep;pause(0.3)
beep;pause(0.3)