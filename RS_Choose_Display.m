%add or subtract from presented number
function RS_Choose_Display()


%turn off screen test for testing
Screen('Preference', 'SkipSyncTests', 1)
clear
commandwindow

%save data?
save_data = true;
file_name = 'test';

require_response = false; %require response? make false to show trials quickly

%open psychtoolbox window
window.bgColor = [80 80 80];
window.dispScreen = max(Screen('Screens'));
[onScreen, rct] = Screen('OpenWindow', window.dispScreen, window.bgColor);
[xo, yo] = RectCenter(rct);
Screen('TextSize', onScreen, 40);


%experiment parameters
dot_num = 7; %number of dots in Grid
dot_size = 5;
display_type=[4,4]; %which tracking task is displayed on each side (0 =nothing, 1=spinner, 2=grid, 3=noisy spinner, 4=bouncy box)
fr = 60; %refresh rate Hz
duration = 6; %number of seconds to show displays for
fade_time = 1.5; %how long to fade dots off for (seconds)
%num_range = [10,90]; %range for number selection
nTrials = 48; %number of trials

%Spinner parameters
spin_speed = 350; %deg/s
spin_change = 1/75; %odds of changing spin direction each frame
spin_min_b4change = 50; %minimum number of spins before direction change
sx1 = xo/2 - 100;%spinner x locations
sx2 = 3*xo/2 + 100;
sy = yo; %spinner y locations
r = 200; %radius of spinners
%noisy spinner
if display_type(1)==3 || display_type(2)==3
    r = r*2/3;
end

%grid Parameters
gx2=[3*xo/2-100, 3*xo/2+100, 3*xo/2+300, 3*xo/2-100, 3*xo/2+100, 3*xo/2+300, 3*xo/2-100, 3*xo/2+100, 3*xo/2+300]; %x coordinates for right hand grid
gx1=[xo/2-300, xo/2-100, xo/2+100, xo/2-300, xo/2-100, xo/2+100, xo/2-300, xo/2-100, xo/2+100];     %x coordinates for left hand grid
gy=[yo+200, yo+200, yo+200, yo, yo, yo, yo-200, yo-200, yo-200]; %y coordinates for grids
perm_num= 25; %number of permutations grid goes through
tranmat=[2,4;1,3;2,6;1,7;0,0;3,9;4,8;7,9;6,8]; %every row shows allowed next transition spaces

%Bounce Box Parameters
box_locl=[xo/2-300,yo-200,xo/2+100,yo+200];
box_locr=[3*xo/2-100,yo-200,3*xo/2+300,yo+200];
closeness=dot_size*1.5;
bounce_speed=4; %number of steps moved each frame 

for cm = 1:2
   color_matrix{cm}= zeros(3,dot_num); %dot colors
end
trial_count = 0;


for trial_num = 1:nTrials
    
    trial_count = trial_count + 1;
    acc = [NaN, NaN];
    
    %initialize display variable so Draw_Display can be called without
    %error. These variables will be set if they are actually needed.
    movmatxl=[];
    movmatyl=[];
    movmatxr=[];
    movmatyr=[];
    z=[];
    bouncematxl=[];
    bouncematxr=[];
    bouncematyl=[];
    bouncematyr=[];
    
    %generate random frequency and amplitude varying sin function for noisy
    %spinner
    if display_type(1)==3 || display_type(2)==3
        x=linspace(0,1000,10000);
        y=rand(7,3);
        z1=4*y(1,1)*sin(10*y(1,2)*x+3*y(1,3));
        z2=4*y(2,1)*sin(8*y(2,2)*x+3*y(2,3));
        z3=4*y(3,1)*sin(7*y(3,2)*x+3*y(3,3));
        z4=2*y(4,1)*sin(3*y(4,2)*x+3*y(4,3));
        z5=2*y(5,1)*sin(y(5,2)*x+3*y(5,3));
        z6=3*y(6,1)*cos(y(6,2)*x+3*y(5,3));
        z7=y(7,1)*cos(7*y(3,2)*x+3*y(7,3));
        z=z1+z2+z3+z4+z5+z6+z7;
        z=z*r/(3*max(z))*1.3;
    end
    
    %generate matrix of legal permutations through the grid display
    if display_type(1)==2 || display_type(2)==2
        gidmatl=Generate_gidmat(dot_num, perm_num, tranmat); %left gidmat
        gidmatr=Generate_gidmat(dot_num,perm_num,tranmat); %right gidmat
        [movmatxl, movmatyl]=Generate_Movmat(gidmatl,gx1,gy,fr*duration); %left movmat
        [movmatxr, movmatyr]=Generate_Movmat(gidmatr, gx2, gy, fr*duration); %right movmat
    end
    
    %Generate matrix defining bounce box dot movement
    if display_type(1)==4
        [bouncematxl, bouncematyl]=Generate_Bouncemat(box_locl, dot_num, dot_size, fr*duration, bounce_speed, closeness);  
    end
    
    if display_type(2)==4
        [bouncematxr, bouncematyr]=Generate_Bouncemat(box_locr, dot_num, dot_size, fr*duration, bounce_speed, closeness);
    end
    
   %Display text info on number of trials completed periodically
   if rem(trial_count, 16) == 1
      Draw_Text(onScreen, [num2str(trial_count-1),'/',num2str(nTrials), ' Completed'], xo, yo, [0,0,0]);
      Screen('Flip', onScreen);
      qcnt = false;
      while ~qcnt
         [~,~,keyCode] = KbCheck();
         if keyCode(KbName('space'))
            qcnt = true;
         end
      end
   end
   
    ori = [randi(360), randi(360)]; %randomize inital orientation
    d = [randi(2), randi(2)]; %randomize initial spin directection
    d(d==2) = -1;
   
   
   targs = 1; %target dots
   fsldc = [0,0]; %frames since last direction change
   
   
   %allow blinking of target dots
   for cm = 1:2
      color_matrix_td{cm} = color_matrix{cm}(:,:);
      color_matrix_td{cm}(:,targs) = window.bgColor(:,1)';
   end
   
   %blink target dots before trial
   HideCursor;
   for i = 1:2
      fn=1;
      Draw_Display(display_type, color_matrix, dot_num, dot_size, onScreen,...
    fn, sx1, sx2, sy, ori, r, movmatxl, movmatyl, movmatxr, movmatyr, z, bouncematxl,...
    bouncematyl, bouncematxr, bouncematyr, box_locl, box_locr, closeness, bounce_speed)
      
      fixation(onScreen, xo, yo);
      Screen('Flip', onScreen);
      WaitSecs(.5);
      
      Draw_Display(display_type, color_matrix, dot_num, dot_size, onScreen,...
    fn, sx1, sx2, sy, ori, r, movmatxl, movmatyl, movmatxr, movmatyr, z, bouncematxl,...
    bouncematyl, bouncematxr, bouncematyr, box_locl, box_locr, closeness, bounce_speed)
      
     Draw_Display(display_type, color_matrix, dot_num, dot_size, onScreen,...
    fn, sx1, sx2, sy, ori, r, movmatxl, movmatyl, movmatxr, movmatyr, z, bouncematxl,...
    bouncematyl, bouncematxr, bouncematyr, box_locl, box_locr, closeness, bounce_speed)
      
      fixation(onScreen, xo, yo);
      Screen('Flip', onScreen);
      WaitSecs(.5);
   end
   
   %allow fading of target dots from bg_color to black
   fade_matrix = linspace(window.bgColor(1), 0, fr*fade_time);
   %number direction decision acc
   
   
   %spin pinwheels for desired time
   
   for fn = 1:duration*fr
      
      ori = ori + spin_speed/fr*d;
      
      %1/spin_change chance of changing spin direction
      cd = [randi(1/spin_change), randi(1/spin_change)];
      cd(fsldc<spin_min_b4change) = 1;
      d(cd==1/spin_change) = d(cd==1/spin_change) - sign(d(cd==1/spin_change))*2;
      fsldc = fsldc + 1;
      fsldc(cd==1/spin_change) = 0;
        
      
     Draw_Display(display_type, color_matrix, dot_num, dot_size, onScreen,...
    fn, sx1, sx2, sy, ori, r, movmatxl, movmatyl, movmatxr, movmatyr, z, bouncematxl,...
    bouncematyl, bouncematxr, bouncematyr, box_locl, box_locr, closeness, bounce_speed)
      
      %fade target dots from bg_color to black
      if fn < fr * fade_time
         for cm = 1:2
            color_matrix_td{cm}(:,targs) = repmat(fade_matrix(fn),3,1);
         end
         Draw_Display(display_type, color_matrix, dot_num, dot_size, onScreen,...
    fn, sx1, sx2, sy, ori, r, movmatxl, movmatyl, movmatxr, movmatyr, z, bouncematxl,...
    bouncematyl, bouncematxr, bouncematyr, box_locl, box_locr, closeness, bounce_speed)
      
      end
      
      fixation(onScreen, xo, yo);
      Screen('Flip', onScreen);
      
   end
   
   if require_response
      
      %response phase
      ShowCursor;
      resp1 = false;
      resp2 = false;
      if display_type(1) == 0
         resp1 = true;
      end
      if display_type(2) == 0
         resp2 = true;
      end
      
      for ds = 1:2
         color_matrix_resp{ds} = color_matrix{ds};
         color_matrix_final{ds} = color_matrix{ds};
      end
      
      for s = 1:2
         pcount = 0;
         for z1 = 0:dot_num-1
            pcount = pcount + 1;
            if s == 1
               xr{s}(pcount) = sx1 + r*cosd(ori(s)+(360/dot_num)*z1);
            elseif s == 2
               xr{s}(pcount) = sx2 + r*cosd(ori(s)+(360/dot_num)*z1);
            end
            yr{s}(pcount) = sy - r*sind(ori(s)+(360/dot_num)*z1);
         end
      end
      
      
      while ~(resp1 && resp2)
         [mx,my,buttons] = GetMouse(onScreen);
         if buttons(1) %if a click
            if ~resp1
               if display_type(1) == 1
                  dis = sqrt((mx-xr{1}).^2 + (my-yr{1}).^2); %assess distance of click
               elseif display_type(1) == 2
                  dis = sqrt((mx-movmatxl(:, fn)).^2 + (my-movmatyl(:, fn)).^2);
               end
               if any(dis < dot_size/2) %if on a dot
                  if any(targs==find(dis<dot_size))  %see if it's a target dot
                     
                     clicked(1) = find(dis<dot_size);
                     resp1 = true;
                     acc(1) = 1;
                     
                     color_matrix_resp{1}(:, clicked(1)) = window.bgColor(:,1)';
                     color_matrix_final{1}(:, clicked(1)) = [0,255,0]';
                     if ~resp2
                        %draw resp1 and blank resp2
                        Draw_Display(display_type, color_matrix, dot_num, dot_size, onScreen,...
    fn, sx1, sx2, sy, ori, r, movmatxl, movmatyl, movmatxr, movmatyr, z, bouncematxl,...
    bouncematyl, bouncematxr, bouncematyr, box_locl, box_locr, closeness, bounce_speed)
      
                        Draw_Display(display_type, color_matrix, dot_num, dot_size, onScreen,...
    fn, sx1, sx2, sy, ori, r, movmatxl, movmatyl, movmatxr, movmatyr, z, bouncematxl,...
    bouncematyl, bouncematxr, bouncematyr, box_locl, box_locr, closeness, bounce_speed)
      
                        fixation(onScreen, xo, yo);
                        Screen('Flip', onScreen);
                     end
                     
                  else %if a miss
                     clicked(1) = find(dis<dot_size);
                     resp1 = true;
                     acc(1) = 0;
                     color_matrix_resp{1}(:, clicked(1)) = window.bgColor(:,1)';
                     color_matrix_final{1}(:, clicked(1)) = [255,0,0]';
                     if ~resp2
                        %draw resp1 and blank resp2
                        Draw_Display(display_type, color_matrix, dot_num, dot_size, onScreen,...
    fn, sx1, sx2, sy, ori, r, movmatxl, movmatyl, movmatxr, movmatyr, z, bouncematxl,...
    bouncematyl, bouncematxr, bouncematyr, box_locl, box_locr, closeness, bounce_speed)
      
                        Draw_Display(display_type, color_matrix, dot_num, dot_size, onScreen,...
    fn, sx1, sx2, sy, ori, r, movmatxl, movmatyl, movmatxr, movmatyr, z, bouncematxl,...
    bouncematyl, bouncematxr, bouncematyr, box_locl, box_locr, closeness, bounce_speed)
      
                        fixation(onScreen, xo, yo);
                        Screen('Flip', onScreen);
                     end
                     
                  end
               end
            end
            if ~resp2
               if display_type(2) == 1
                  dis = sqrt((mx-xr{2}).^2 + (my-yr{2}).^2); %assess distance of click
               elseif display_type(2) == 2
                  dis = sqrt((mx-movmatxr(:, fn)).^2 + (my-movmatyr(:, fn)).^2);
               end
               if any(dis < dot_size/2) %if on a dot
                  if any(targs==find(dis<dot_size))  %see if it's a target dot
                     
                     clicked(2) = find(dis<dot_size);
                     resp2 = true;
                     acc(2) = 1;
                     color_matrix_resp{2}(:, clicked(2)) = window.bgColor(:,1)';
                     color_matrix_final{2}(:, clicked(2)) = [0,255,0]';
                     if ~resp1
                        %draw resp2 and blank resp1
                        Draw_Display(display_type, color_matrix, dot_num, dot_size, onScreen,...
    fn, sx1, sx2, sy, ori, r, movmatxl, movmatyl, movmatxr, movmatyr, z, bouncematxl,...
    bouncematyl, bouncematxr, bouncematyr, box_locl, box_locr, closeness, bounce_speed)
      
                        Draw_Display(display_type, color_matrix, dot_num, dot_size, onScreen,...
    fn, sx1, sx2, sy, ori, r, movmatxl, movmatyl, movmatxr, movmatyr, z, bouncematxl,...
    bouncematyl, bouncematxr, bouncematyr, box_locl, box_locr, closeness, bounce_speed)
      
                        fixation(onScreen, xo, yo);
                        Screen('Flip', onScreen);
                     end
                     
                  else %if a miss
                     clicked(2) = find(dis<dot_size);
                     resp2 = true;
                     acc(2) = 0;
                     color_matrix_resp{2}(:, clicked(2)) = window.bgColor(:,1)';
                     color_matrix_final{2}(:, clicked(2)) = [255,0,0]';
                     if ~resp1
                        %draw resp2 and blank resp1
                        Draw_Display(display_type, color_matrix, dot_num, dot_size, onScreen,...
    fn, sx1, sx2, sy, ori, r, movmatxl, movmatyl, movmatxr, movmatyr, z, bouncematxl,...
    bouncematyl, bouncematxr, bouncematyr, box_locl, box_locr, closeness, bounce_speed)
      
                        Draw_Display(display_type, color_matrix, dot_num, dot_size, onScreen,...
    fn, sx1, sx2, sy, ori, r, movmatxl, movmatyl, movmatxr, movmatyr, z, bouncematxl,...
    bouncematyl, bouncematxr, bouncematyr, box_locl, box_locr, closeness, bounce_speed)
      
                        fixation(onScreen, xo, yo);
                        Screen('Flip', onScreen);
                     end
                     
                  end
               end
            end
         end
      end
      
      
      Draw_Display(display_type, color_matrix, dot_num, dot_size, onScreen,...
    fn, sx1, sx2, sy, ori, r, movmatxl, movmatyl, movmatxr, movmatyr, z, bouncematxl,...
    bouncematyl, bouncematxr, bouncematyr, box_locl, box_locr, closeness, bounce_speed)
      
      Draw_Display(display_type, color_matrix, dot_num, dot_size, onScreen,...
    fn, sx1, sx2, sy, ori, r, movmatxl, movmatyl, movmatxr, movmatyr, z, bouncematxl,...
    bouncematyl, bouncematxr, bouncematyr, box_locl, box_locr, closeness, bounce_speed)
      
      fixation(onScreen, xo, yo);
      Screen('Flip', onScreen);
      WaitSecs(1);
      
   end
   
   %record data
   D.acc1(trial_num) = acc(1);
   D.acc2(trial_num) = acc(2);
   D.display_type = display_type;
   
   %save data
   if save_data
      save(file_name, 'D');
   end
   
end

WaitSecs(.5);
end


%Functions to Draw Text
function Draw_Text(onScreen, txt, x, y, text_color)
[xt,yt] = getTextCenter(txt,onScreen,x,y);
DrawFormattedText(onScreen, txt,xt, yt, text_color);
end

function [ sx, sy ] = getTextCenter( textStr, screenPtr, x,y )
%GETTEXTCENTER Summary of this function goes here
%   Detailed explanation goes here
% textStr is the text you want to center
% x,y is the point you want to center it on
% sx sy is where the pen should start drawing to center this
% text on this point
[normBoundsRect]=Screen('TextBounds', screenPtr, textStr);
boxcor = normBoundsRect;
startX = x - (boxcor(3) - boxcor(1))/2;
%endX = x + (boxcor(3) - boxcor(1))/2;
startY = y - (boxcor(4) - boxcor(2))/2;
%endY = y + (boxcor(4) - boxcor(2))/2;
sx = startX;
sy = startY;
end

%fixation
function fixation(onScreen, xo, yo)
% Screen('DrawLine', onScreen,  [0,0,0], xo-15, yo, xo+15, yo, 2);
% Screen('DrawLine', onScreen,  [0,0,0], xo, yo-15,xo,yo+15, 2);
Screen('DrawDots', onScreen, [xo;yo], 5, [0;0;0], [], 1);
end

%draw a grid of grid identifier vector
function Draw_Gid(onScreen, gx, gy, dot_size, color_matrix,gid)

for i=1:numel(gid)
   
   dot_locs(1,i)=gx(gid(i));
   dot_locs(2,i)=gy(gid(i));
end

Screen('DrawDots', onScreen, dot_locs, dot_size, color_matrix, [], 1)
end


%Generate a matrix of random legal grid permutations
function [gidmat]=Generate_gidmat(dot_num, perm_num, tranmat)
gidmat=zeros(perm_num-1,dot_num);
a=[1:9];
gid=[]; %grid identifier vector
for i=1:dot_num %Generate random starting position 'gid'
   x=a(randi(numel(a)));
   gid=[gid x];
   
   for j=1:numel(a)
      if a(j)==x
         a(j)=[];
         break
      end
   end
end
gidmat=[gid;gidmat]; %matrix representation of permutations
for m=2:perm_num %Generate gitmat
   check=1;
   while check==1
      check=2;
      gidmat(m,:)=zeros(1,dot_num);
      l=randi(dot_num); %move random initial dot
      while gidmat(m-1,l)==5
         l=randi(dot_num);
      end
      c=randi(2);
      gidmat(m,l)=tranmat(gidmat(m-1,l),c);
      n=[1:dot_num];
      n(l)=[];
      for i=1:dot_num-1     %move rest of dots in sequence
         if gidmat(m-1,n(i))==5 %move dot in 5 position last
            temp=n(i);
            n(i)=n(dot_num-1);
            n(dot_num-1)=temp;
         end
         if gidmat(m-1,n(i))==5 %move dot from 5 position
            c=1:9;
            c(5)=[];
            d=[];
            for j=1:length(c)
               for jj=1:dot_num
                  if c(j)==gidmat(m,jj)
                     d=[d 1];
                     break
                  end
                  if gidmat(m-1,jj)==c(j) && gidmat(m,jj)==5
                     d=[d 1];
                     break
                     
                  end
                  if jj==dot_num
                     d=[d 0];
                     
                  end
               end
            end
            c(logical(d))=[];
            gidmat(m,n(i))=c(randi(length(c)));
         else     %move dot from non 5 position
            c=randi(2);
            gidmat(m,n(i))=tranmat(gidmat(m-1,n(i)),c);
            if m>2
               if gidmat(m,n(i))==gidmat(m-2,n(i)) %preference given to new space
                  if c==1
                     gidmat(m,n(i))=tranmat(gidmat(m-1,n(i)),2);
                  else
                     gidmat(m,n(i))=tranmat(gidmat(m-1,n(i)),1);
                  end
               end
            end
            
            k=1:dot_num;
            k(n(i))=[];
            for ii=1:dot_num-1 % tests for legal move
               if gidmat(m,n(i))==gidmat(m,k(ii)) || (gidmat(m,n(i))==gidmat(m-1,k(ii)) && gidmat(m-1,n(i))==gidmat(m,k(ii)))
                  if c==1
                     gidmat(m,n(i))=tranmat(gidmat(m-1,n(i)),2);
                  else
                     gidmat(m,n(i))=tranmat(gidmat(m-1,n(i)),1);
                  end
                  for iii=1:dot_num-1 %test second option for legal move
                     if gidmat(m,n(i))==gidmat(m,k(iii)) || (gidmat(m,n(i))==gidmat(m-1,k(iii)) && gidmat(m-1,n(i))==gidmat(m,k(iii)))
                        gidmat(m,n(i))=5;
                        break
                     end
                  end
               end
            end
         end
         if i==dot_num-1 && gidmat(m-1,n(i))~=5 %  Give dot a chance to move to 5
            for ii=1:dot_num-2
               if gidmat(m,n(ii))==5
                  break
               end
               if ii==dot_num-2
                  c=randi(2);
                  if c==1
                     gidmat(m,n(i))=5;
                  end
               end
            end
         end
      end
      c=0;
      for i=1:dot_num
         if gidmat(m,i)==5
            c=c+1;
         end
         if c>=2
            check=1;
            break
         end
      end
   end
end
end

function [movmatx, movmaty]=Generate_Movmat(gidmat,gx,gy,totframenum)
movmatx=[];
movmaty=[];
for i=1:size(gidmat,2)
   xmov=[];
   ymov=[];
   for ii=1:(size(gidmat,1)-1)
      xmov=[xmov linspace(gx(gidmat(ii,i)),gx(gidmat(ii+1,i)),totframenum/(size(gidmat,1)-2)+1)];
      ymov=[ymov linspace(gy(gidmat(ii,i)),gy(gidmat(ii+1,i)),totframenum/(size(gidmat,1)-2)+1)];
      if ii~=(size(gidmat,1)-1)
         xmov(end)=[];
         ymov(end)=[];
      end
   end
   movmatx=[movmatx; xmov];
   movmaty=[movmaty; ymov];
   
end
end

function [bouncematx, bouncematy]=Generate_Bouncemat(box_loc, dot_num, dot_size, totframenum, bounce_speed, closeness)
try_again=true;
while try_again
    try_again=false;
    bouncematx=double(randi([box_loc(1)+dot_size box_loc(3)-dot_size]));
    bouncematy=double(randi([box_loc(2)+dot_size box_loc(4)-dot_size]));
    check_angle=true;
    while check_angle
        theta=randi([1 359],dot_num,1);      %Starting dot movement angles
        for i=1:dot_num
            if theta(i,1)==90 || theta(i,1)==180 || theta(i,1)==270
                break
            end
            if i==dot_num
                check_angle=false;
            end
        end
    end
    
    for i=2:dot_num       %Generate properly spaced starting x coordinate dot positions
        check_spacing=true;
        while check_spacing
            bouncematx(i,1)=randi([box_loc(1)+dot_size box_loc(3)-dot_size]);
            for ii=1:(size(bouncematx,1)-1)
                if abs(bouncematx(i,1)-bouncematx(ii,1))<closeness
                    break
                end
                if ii==size(bouncematx,1)-1
                    check_spacing=false;
                end
            end
        end
    end
    for i=2:dot_num       %Generate properly spaced starting y coordinate dot positions
        check_spacing=true;
        while check_spacing
            bouncematy(i,1)=randi([box_loc(2)+dot_size box_loc(4)-dot_size]);
            for ii=1:size(bouncematy,1)-1
                if abs(bouncematy(i,1)-bouncematy(ii,1))<closeness
                    break
                end
                if ii==size(bouncematy,1)-1
                    check_spacing=false;
                end
            end
        end
    end
    
    %generate next dot positions
    for i=2:totframenum
        bouncematx(:,i)=bouncematx(:,i-1)+bounce_speed*cosd(theta(:,1));
        bouncematy(:,i)=bouncematy(:,i-1)+bounce_speed*sind(theta(:,1));
        %check to see if dot touching wall, if so change angle
        for ii=1:dot_num
            if abs(bouncematx(ii,i)-box_loc(1))<bounce_speed || abs(bouncematx(ii,i)-box_loc(3))<bounce_speed
                theta(ii,1)=mod(-theta(ii,1)+180,360);
            end
            if abs(bouncematy(ii,i)-box_loc(2))<bounce_speed || abs(bouncematy(ii,i)-box_loc(4))<bounce_speed
                theta(ii,1)=mod(-theta(ii,1),360);
            end
        end
        %check to see if any dots are too close, and start over if they are
        for ii=1:dot_num-1
            for iii=ii+1:dot_num
                if abs(bouncematx(ii,i)-bouncematx(iii,i))<closeness && abs(bouncematy(ii,i)-bouncematy(iii,i))<closeness
                    try_again=true;
                    break
                end
            end
            if try_again
                break
            end
        end
        if try_again
            break
        end
    end
end

end
            


                
            

%function for dot rotation
function Spinner(onScreen, x, y, ori, r, dot_size, color_matrix, dot_num)
for dl = 1:dot_num
   dot_locs(1,dl) = x + r*cosd(ori + (360/dot_num)*(dl-1));
   dot_locs(2,dl) = y - r*sind(ori + (360/dot_num)*(dl-1));
end
% dot_locs = [x - r*cosd(ori), x + r*cosd(ori), x - r*cosd(ori+90), x + r*cosd(ori+90);...
%     y+r*sind(ori), y-r*sind(ori), y+r*sind(ori+90), y-r*sind(ori+90)];
Screen('DrawDots', onScreen, dot_locs, dot_size, color_matrix, [], 1);
end


%function for dot rotation with noisy sin wave
function Nspinner(onScreen, x, y, ori, r, dot_size, color_matrix, dot_num, z, fn)
for dl = 1:dot_num
    dot_locs(1,dl) = x +(r+z(417*(dl-1)+fn))*cosd(ori + (360/dot_num)*(dl-1));
    dot_locs(2, dl) = y - (r+z(417*(dl-1)+fn))*sind(ori + (360/dot_num)*(dl-1));
end
Screen('DrawDots', onScreen, dot_locs, dot_size, color_matrix, [], 1);
end

function Draw_Grid(onScreen, movmatx, movmaty, fn, dot_size, color_matrix)
for i=1:size(movmatx,1) 
   dot_locs(1,i)=movmatx(i,fn);
   dot_locs(2,i)=movmaty(i,fn);
end
Screen('DrawDots',onScreen, dot_locs, dot_size, color_matrix, [], 1)


end



function Draw_Bounce_Box(onScreen, bouncematx, bouncematy, fn, dot_size, color_matrix, box_loc)
for i=1:size(bouncematx,1)
    dot_locs(1,i)=bouncematx(i,fn);
    dot_locs(2,i)=bouncematy(i,fn);
end
dot_locs=double(dot_locs);
Screen('DrawDots', onScreen, dot_locs, dot_size, color_matrix, [], 1)
Screen('FrameRect', onScreen, [0,0,0], box_loc)
end



function Draw_Display(display_type, color_matrix, dot_num, dot_size, onScreen,...
    fn, sx1, sx2, sy, ori, r, movmatxl, movmatyl, movmatxr, movmatyr, z, bouncematxl,...
    bouncematyl, bouncematxr, bouncematyr, box_locl, box_locr, closeness, bounce_speed)
% Draw two displays simultaneously, kind depends on display_type
if display_type(1)==1 
    Spinner(onScreen, sx1, sy, ori(1), r, dot_size, color_matrix{1}, dot_num);
end
if display_type(2)==1
    Spinner(onScreen, sx2, sy, ori(2), r, dot_size, color_matrix{2}, dot_num);
end
if display_type(1)==2
    Draw_Grid(onScreen, movmatxl, movmatyl, fn, dot_size, color_matrix{1});
end
if display_type(2)==2
    Draw_Grid(onScreen, movmatxr, movmatyr, fn, dot_size, color_matrix{2});
end
if display_type(1)==3
    Nspinner(onScreen, sx1, sy, ori(1), r, dot_size, color_matrix{1}, dot_num, z, fn);
end
if display_type(2)==3
    Nspinner(onScreen, sx2, sy, ori(2), r, dot_size, color_matrix{2}, dot_num, z, fn);
end
if display_type(1)==4
    Draw_Bounce_Box(onScreen, bouncematxl, bouncematyl, fn, dot_size, color_matrix{1}, box_locl)
end
if display_type(2)==4
    Draw_Bounce_Box(onScreen, bouncematxr, bouncematyr, fn, dot_size, color_matrix{2}, box_locr)
end
    
end
