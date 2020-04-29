%1 Load the image
filem = mfilename('fullpath')

map = imread('D:\Data\ProgramMatlab\CurrentSetPath\spm8\toolbox\nirsHSJ_JT\GUI\Template\BebeROIDETECTOR.png')
figure
image(map(:,:,:))
datacursormode on
%save the cursor to the workspace


%Save the workspace in a mat file