function table(superheader,header,data,formats)
%TABLE Internal function to display tables

% Author Johan Löfberg
% $Id: table.m,v 1.2 2004-07-02 08:17:32 johanl Exp $


  [nheadersy,nheadersx] = size(header);
  [ndatay,ndatax] = size(data);
  
  datasizes = zeros(ndatay,ndatax);
  for i = 1:ndatay
    for j = 1:ndatax
      if isa(data{i,j},'double')
	data{i,j} = num2str(data{i,j});
      end
        datasizes(i,j) = length(data{i,j}); 
    end
  end
  
  headersizes = zeros(1,nheadersx);
  for j = 1:nheadersx
    if isa(header{j},'double')
      header{j} = num2str(header{j});
    end
    headersizes(1,j) = length(header{j}); 
  end
  
  if nargin<4
    for i = 1:ndatax
      formats{i}.header.just = 'right';
      formats{i}.data.just = 'right';
    end
  end
  
  datawidth = sum(datasizes,2);
  
  MaxWidth = max([headersizes;datasizes]);
  HeaderLine = ['|'];
  for i = 1:nheadersx
    HeaderLine = [HeaderLine ' ' strjust(fillstringRight(header{i},MaxWidth(i)+2),formats{i}.header.just) '|'];
  end
  HeaderLine = [HeaderLine ''];
  
  for j = 1:ndatay
    DataLine{j} = ['|'];
    for i = 1:ndatax
      DataLine{j} = [DataLine{j} ' '  strjust(fillstringRight(data{j,i},MaxWidth(i)+2),formats{i}.data.just) '|'];
    end
  end 
  if ~isempty(superheader)
    disp(char(repmat(double('+'),1,length(HeaderLine))))
    disp(['|' strjust(fillstringLeft(superheader{1},length(HeaderLine)-2),'center') '|'])
  end
  disp(char(repmat(double('+'),1,length(HeaderLine))))
  disp(HeaderLine)
  disp(char(repmat(double('+'),1,length(HeaderLine))))
  for i = 1:length(DataLine)
    disp(DataLine{i});
  end
  disp(char(repmat(double('+'),1,length(HeaderLine))))
  


function x= truncstring(x,n)
if length(x) > n
	x = [x(1:n-3) '...'];
end

function x = fillstringLeft(x,n)
x = [x blanks(n-length(x))];

function x = fillstringRight(x,n)
x = [blanks(n-length(x)) x];



