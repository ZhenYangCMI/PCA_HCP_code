close all
clear all
clc

b1=textread('subjectList_ascending.txt','%d');
[number,txt,raw]= xlsread('behaviorial.csv');
% [number,txt,raw]= xlsread('demographic.csv');
[r,c]=size(number);
[rt,xt]=size(txt);
[r1,c1]=size(b1);

headlines=rt-r;

for i=1:headlines
    str1='A';
    linenumber=sprintf('%s%d',str1,i);
    xlswrite('results2.xls', raw(i,:), 'sheet1',linenumber); 
end

bcounter=headlines+1;
for i=1:r
    for j=1:r1
        if number(i,1)==b1(j,1)
            str1='A';
            linenumber=sprintf('%s%d',str1,bcounter);
            xlswrite('results2.xls', raw(i+headlines,:), 'sheet1',linenumber); 
            bcounter=bcounter+1;
        end
    end
end
