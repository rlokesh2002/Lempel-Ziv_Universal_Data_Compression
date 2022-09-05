%   REKHA LOKESH
%   19EC10052
%   Digicomm Assignment 1

clear all;
close all;
clc;

src_path = "MnM_source_file.txt";   %   Path for the Source Text File
src_temp = readlines(src_path); %    Importing the text as string array

src_str = "";
for i=1:size(src_temp, 1)
    src_str = src_str + src_temp(i) + newline;   %   Writing the whole string array as single string where is for newline
end

src = char(src_str);   %Converting string to char vec

w = 2^14; %Length of the sliding window

%   --------    Preparing Dictionary   ---------
uniq = unique(src); %Finding the unique occurences of characters in the string

%Let the dictionary consists of the unique ASCII value characters in the
%string

M = size(uniq,2);   %Number of unique characters
key_vec = repmat(cellstr(char('a')), M, 1); %Creating the keyvec to have only the number of unique occurences
key_vec(1, :) = cellstr(string(uniq(1)));
% bin_size = ceil(log2(size(key_vec,1))); %length of the binary values in the fixed length codes of the dictionary
bin_size = log2(w);
value_vec = repmat(convertCharsToStrings(dec2bin(0, bin_size)), M, 1);
for i = 1:(M-1)
   key_vec(i+1, :) = cellstr(string(uniq(i+1)));
   value_vec(i+1, :) = (convertCharsToStrings(dec2bin(i, bin_size)));
end

dict = containers.Map(transpose(key_vec), transpose(value_vec));   %Preparing the dictionary is done


%   --------    ENCODING STARTS   ---------


%Encoding the first w symbols of src using fixed length encoding
encd = "";   %String that stores the encoded values
for i = 1:min([w, size(src, 2)], [], 'all')
    encd = encd + "1" + string(dict(string(src(i))));   %Fixed length encoding is done with n = 1
end

P = w;  %Pointer is set to w

disp("Size of src is: "+int2str(size(src, 2)));
src_size = size(src, 2);
while (P<src_size)
    [n, u] = findlargestmatch(src, P, w);
    if (n==1)
        encd = encd + "1" + string(dict(string(src(P+1))));   %Fixed length encoding is done
    else
        temp1 = dec2bin(n, 2*floor(log2(n))+1); %n codeword Unary Binary Encoding
        temp2 = dec2bin(u, log2(w)); %u codeword Fixed length Encoding
        encd = encd + string(temp1) + string(temp2);    %Encoded n, u
    end
    P = P + n;
    disp("Rem: "+int2str(src_size-P));    %To see the progress of the simulation
end

%Writing encoded data in a txt file
fileID = fopen("Encoded_Data_w="+int2str(w)+".txt", 'w');
fprintf(fileID, encd);
fclose(fileID);
comp_ratio = size(src, 2)*bin_size/strlength(encd);  %since Compression Ratio = Uncompressed Size/Compressed Size
comp_percent = (comp_ratio-1)/comp_ratio*100;
disp("w = " + int2str(w) + ", Uncompressed size = " + int2str(src_size*bin_size) + ", Compressed size = " + strlength(encd) + ", Compression Ratio is: "+sprintf("%.4f", comp_ratio)+", Compressed % = "+sprintf("%.4f", comp_percent));
%Encoding Ends

%   --------    DECODING STARTS   ---------
% Calling the Decoding Function
char_dec = lnz_decode(encd, dict);   %Calling the decoding functions whose arguments are encoded string and dictionary which returns the decoded string 

if(char_dec==src)   
    disp("Decoding completed successfully and correctly");
else
    disp("Decoding completed successfully and incorrectly");
end

%Printing the decoding symbols in text file
fileID = fopen("Decoded_text_file_w = "+int2str(w)+".txt", 'w');
m = 1;
sz = size(char_dec, 2);
while (m<=sz)
    q = 0;
    while ((m+q)<=sz && char_dec(m+q)~=char(newline))
       q = q + 1;   %Traversing the whole sentence until a newline character is detected
    end
    if(m+q+1 <= sz)
        fprintf(fileID, string(char_dec(m:(m+q-1)))+'\n');  %print newline if this is not the last string to be entered
    else
        fprintf(fileID, string(char_dec(m:(m+q-1))));
    end
    m = m+q+1;
end
fclose(fileID);



function [res1, res2] = findlargestmatch(src, P, w)
    sz = size(src, 2);
    if ((sz-P)<2)
        res1 = 1;
        res2 = 1;
    else
        res1 = 1;
        res2 = 1;
        for n = 2:(sz-P)
            substr = src(1, (P+1):(P+n));
            k = strfind((src(1, (P+2-w):(P+n-1))), substr);   %Since Range of u is 1 to w
            if(all(size(k)~= [0, 0], 'all') && (k(1, 1)+P-w)<=P)
               res1 = n;
               res2 = P-(k(1, 1)+P-w+1)+1;
            end
        end   
    end 
end

%   --------    DECODING Function   ---------
function dest_dec = lnz_decode(enc_src, dic)  
    enc = char(enc_src); %Converting encoded string to char to access easily
    res = '';  %Initializing decoded char array
    val = values(dic);
    val_sz = strlength(string(val(1))); %finding the values size from the dictionary
    dec_dict = containers.Map(values(dic), keys(dic));  %Inverse Mapping the dictionary values
    i = 1;
    while (i<=size(enc, 2))
        if(all(enc(i)=='1', "all"))   %Fixed Length encoding is detected
            %Then detect the next w/val_sz bits from where the symbol can
            %be detected by reverse mapping
            i = i + 1;
            tmp = enc(i:(i+val_sz-1));
            symb = dec_dict(tmp);   %decoding the symbol as string
            res = char(string(res)+symb); %Concatenating the symbol to the decoded string
            i = i + val_sz;
        else
           %It is a unary-binary code. So go on detecting zeroes until 1
           %appears. If k zeros detected then binary length of n is 2k+1
           k = 1; %number of zeros detected
           while (enc(i+k)=='0')
              k = k + 1; 
           end
           n = bin2dec(enc(i:(i+2*k)));
           i = i + 2*k+1;
           u = bin2dec(enc(i:(i+val_sz-1)));
           p = size(res, 2);
           st = p-u+1;
           en = st+n-1;
           if(en<=p)    %no overlapping between string and symbols
                res = char(string(res)+string(res(st:en)));%Appending the string 
           else         %overlapping between string and symbols
                symb = string(res(st:p));

                len = p-st+1;
                rem = n;
                stri = "";
                while (rem>=len)
                    stri = stri + symb;
                    rem = rem - len;
                end
                if(rem>0)
                   stri = stri + string(res(st:(st+rem-1))) ;
                end
                res = char(string(res)+stri);%Appending the string 
           end
           i = i + val_sz;
        end 
    end
    dest_dec = res;
end