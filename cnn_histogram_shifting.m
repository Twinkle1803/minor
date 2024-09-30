
function imgW = cnn_histogram_shifting(img, predictedImage, watermark, oddOrEvenNum)

%----------Histogram shift reversible watermarking method-----------

img = double(img);
[M,N]=size(img);
mn=nextpow2(M*N);

imgW = img;
watermarkLength = length(watermark);

oddOrEvenPlace = zeros(M,N);
locationMap = zeros(M,N);
predictedError = zeros(M,N);%Calculating prediction error
complexity = zeros(M,N);

for i=2:M-1
    for j=2:N-1
        if mod((i+j),2)==oddOrEvenNum %oddOrEvenNum is "1" for odd number embedding and "0" for even number embedding
            oddOrEvenPlace(i,j)=1;
            complexity(i, j)=calculate_complexity(imgW, i, j);%Calculating background complexity
        end
    end
end

inserablePlace = find((oddOrEvenPlace==1));
[~, indexComplexity] = sort(complexity(inserablePlace));%Embeddable location
minimumLsbLength = watermarkLength+6*(mn);%Preset the length of information to be replaced by the LSB algorithm
currentLength = minimumLsbLength; %To facilitate subsequent loops, define the current length, start position and end position
startPlace=1;
endPlace=currentLength; 
while(1)
    for k=startPlace:endPlace
        i = mod(inserablePlace(indexComplexity(k)),M);
        j = ceil(inserablePlace(indexComplexity(k))/M);
        x = imgW(i,j);
        xPredict = predictedImage(i,j);
        value = 2*x - xPredict;
        predictedError(i, j) = x - xPredict; %Calculating prediction error
        if (value>254)||(value<0)
            %Position marker 1 that cannot be embedded
            locationMap(i,j)=1; 
        end            
    end
    %Calculate the location map, number1 and number0 are used for arithmetic coding
    number1 = sum(sum(locationMap));
    number0 = M*N - number1;
    if number1==0
        compressedLocationMapBitlist = 0;
        number0 = 0;
    else
        compressedLocationMapBitlist = arithenco(locationMap(:)+1,[number0,number1]);
    end
    if currentLength - minimumLsbLength < length(compressedLocationMapBitlist)+6*mn
        currentLength = currentLength + 1000; %The current length increases by 1000 each time
        startPlace = endPlace + 1;
        endPlace=currentLength;
    else
       break;
    end
end
%Perform LSB replacement at locations with high background complexity, extracting the least significant bit first
extractedLsbBitlist = [];
for k=length(indexComplexity):-1:length(indexComplexity)-(length(compressedLocationMapBitlist)+6*mn)+1
    extractedLsbBitlist(length(indexComplexity)-k+1) = bitget(imgW(inserablePlace(indexComplexity(k))),1); %lsb_replace_bitlist is the lsb bit of the corresponding length
end
%The LocationMap to be embedded is composed of length + bit stream
lengthCompressedLocationMapBitlist = num2bitlist(length(compressedLocationMapBitlist), mn);
wholeCompressedLocationMapBitlist = [lengthCompressedLocationMapBitlist, compressedLocationMapBitlist'];

%Information to be embedded using histogram shift
messageToEmbed = [extractedLsbBitlist,watermark];
messageToEmbedLength = num2bitlist(length(messageToEmbed),mn);

%Calculate Tp and Tn for easy embedding
[Tp, Tn] = calculate_tp_tn(predictedError(inserablePlace(indexComplexity(1:endPlace))), length(messageToEmbed));

TpBitlist = num2bitlist(Tp, mn);
TnBitlist = num2bitlist(abs(Tn), mn);

number1Bitlist = num2bitlist(number1, mn);
number0Bitlist = num2bitlist(number0, mn);

%Bitstream information to be replaced by LSB£¬6*mn
lsbToReplaceBitlist = [wholeCompressedLocationMapBitlist,...
messageToEmbedLength, TpBitlist, TnBitlist, number1Bitlist, number0Bitlist];

%Perform LSB replacement
for k=length(indexComplexity):-1:length(indexComplexity)-length(lsbToReplaceBitlist)+1
    imgW(inserablePlace(indexComplexity(k))) = imgW(inserablePlace(indexComplexity(k))) ...
        - bitget(imgW(inserablePlace(indexComplexity(k))),1)...
        + lsbToReplaceBitlist(length(indexComplexity)-k+1); 
end

%Replace LSB and embed watermark information
indexMessage = 1;
for k=1:length(indexComplexity)
    i = mod(inserablePlace(indexComplexity(k)),M);
    j = ceil(inserablePlace(indexComplexity(k))/M);
    if locationMap(i,j)==0
        x = imgW(i,j);
        xPredict = predictedImage(i,j);
        dij= x - xPredict;
        if dij > Tp
            Dij = dij + Tp + 1;
            imgW(inserablePlace(indexComplexity(k)))=Dij + xPredict;
        elseif dij < Tn
            Dij = dij + Tn;
            imgW(inserablePlace(indexComplexity(k)))=Dij + xPredict;
        else
            if messageToEmbed(indexMessage)==1
                Dij = 2*dij + 1;
            else
                Dij = 2*dij + 0;
            end
            imgW(inserablePlace(indexComplexity(k)))=Dij + xPredict;
            indexMessage = indexMessage + 1;
            if indexMessage > length(messageToEmbed)
                break;
            end
        end
    end
end


end