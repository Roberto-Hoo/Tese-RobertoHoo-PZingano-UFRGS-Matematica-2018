
fileID = fopen('u24h_serial_INPUT.dat','r');

A = fscanf(fileID,'%f',[1 inf]);   % <--- this command reads 
                                   %     the entire data file 
fclose(fileID);

M1 = A(1); M2 = A(2); 
N1 = A(3); N2 = A(4);

n = A(5);

p = A(6);
h = A(7);

cfl = A(8);

x_min = A( 9); x_max = A(10);
y_min = A(11); y_max = A(12);

previous_tF = A(13); tF = A(14);  dt_dump = A(15);
t0 = 0;

x1 = A(16); x2 = A(17); x3 = A(18);
y1 = A(19); y2 = A(20); y3 = A(21);

i1 = A(22); i2 = A(23); i3 = A(24);
j1 = A(25); j2 = A(26); j3 = A(27);

b1 = A(28); b2 = A(29); b3 = A(30);
b  = A(31);
expected_b = (b1+b2+b3)/3;

refv1 = A(32); refv2 = A(33); refv3 = A(34);

count = A(35);
no_runs = A(36);

i_ff1a = A(37); i_ff1b = A(38); 
i_ff2a = A(39); i_ff2b = A(40); 
i_ff3a = A(41); i_ff3b = A(42); 

j_ff1a = A(43); j_ff1b = A(44); 
j_ff2a = A(45); j_ff2b = A(46); 
j_ff3a = A(47); j_ff3b = A(48); 

length_i1 = A(49); length_i2 = A(50); length_i3 = A(51);
length_j1 = A(52); length_j2 = A(53); length_j3 = A(54);

ff1_value_u0 = A(55);
ff2_value_u0 = A(56);
ff3_value_u0 = A(57);

min_u0 = A(58); max_u0 = A(59);
mass_u0 = A(60);

t1_min = A(61); t2_min = A(62); t3_min = A(63);
t1_max = A(64); t2_max = A(65); t3_max = A(66);

pos_init = 67;

% Reading LOG summary variables:

DATE_start_LOG_previous  = zeros(no_runs+1,8);
DATE_finish_LOG_previous = zeros(no_runs+1,8);
tF_LOG_previous = zeros(no_runs+1,1);
ff1_value_LOG_previous = zeros(no_runs+1,1);
ff2_value_LOG_previous = zeros(no_runs+1,1);
ff3_value_LOG_previous = zeros(no_runs+1,1);
variation_mass_LOG_previous = zeros(no_runs+1,1);
variation_sup_LOG_previous  = zeros(no_runs+1,1);
min_u_LOG_previous = zeros(no_runs+1,1);
max_u_LOG_previous = zeros(no_runs+1,1);
elapsed_time_LOG_previous = zeros(no_runs+1,1);

for i = 1: no_runs + 1
    DATE_start_LOG_previous(i,:) = A(pos_init:pos_init+7);
    pos_init = pos_init + 8;
end
for i = 1: no_runs + 1
    DATE_finish_LOG_previous(i,:) = A(pos_init:pos_init+7);
    pos_init = pos_init + 8;
end

tF_LOG_previous = A(pos_init:pos_init+no_runs);
pos_init = pos_init + no_runs + 1;

ff1_value_LOG_previous = A(pos_init:pos_init+no_runs);
pos_init = pos_init + no_runs + 1;
ff2_value_LOG_previous = A(pos_init:pos_init+no_runs);
pos_init = pos_init + no_runs + 1;
ff3_value_LOG_previous = A(pos_init:pos_init+no_runs);
pos_init = pos_init + no_runs + 1;

variation_mass_LOG_previous = A(pos_init:pos_init+no_runs);
pos_init = pos_init + no_runs + 1;
variation_sup_LOG_previous  = A(pos_init:pos_init+no_runs);
pos_init = pos_init + no_runs + 1;

min_u_LOG_previous = A(pos_init:pos_init+no_runs);
pos_init = pos_init + no_runs + 1;
max_u_LOG_previous = A(pos_init:pos_init+no_runs);
pos_init = pos_init + no_runs + 1;

elapsed_time_LOG_previous = A(pos_init:pos_init+no_runs);
pos_init = pos_init + no_runs + 1;

% Reading meshgrid´points (x,y)
% and current solution values:

x = A(pos_init:pos_init+M2-M1);
pos_init = pos_init + M2-M1 + 1;

y = A(pos_init:pos_init+N2-N1);
pos_init = pos_init + N2-N1 + 1;

u = zeros(M2-M1+1, N2-N1+1);
for i = 1:M2-M1+1
    u(i,:) = A(pos_init:pos_init+N2-N1);
    pos_init = pos_init + N2-N1 + 1;
end

% Reading SOLUTION STATISTICS data:

new_length = A(pos_init);
previous_length = new_length;

ff1_value_u = zeros(new_length+1,1);
ff2_value_u = zeros(new_length+1,1);
ff3_value_u = zeros(new_length+1,1);

variation_mass = zeros(new_length+1,1);
variation_sup  = zeros(new_length+1,1);

min_u = zeros(new_length+1,1);
max_u = zeros(new_length+1,1);

time_per_cycle = zeros(new_length+1,1);

for i = 1: new_length + 1
    ff1_value_u(i) = A(pos_init+1);
    ff2_value_u(i) = A(pos_init+2);
    ff3_value_u(i) = A(pos_init+3);
    pos_init = pos_init + 3;
end

for i = 1: new_length + 1
    variation_mass(i) = A(pos_init+1);
    variation_sup(i)  = A(pos_init+2);
    pos_init = pos_init + 2;
end

for i = 1: new_length + 1
    min_u(i) = A(pos_init+1);
    max_u(i) = A(pos_init+2);
    pos_init = pos_init + 2;
end

for i = 1: new_length + 1
    time_per_cycle(i) = A(pos_init+1);
    pos_init = pos_init + 1;
end

disp('total number of data elements read:')
size_A = size(A)
pos_init 

disp('plotting...')

[X,Y] = meshgrid(x,y);

figure(10)
mesh(X',Y',u)
colormap(jet)
shading interp
xlabel('x')
ylabel('y')

if tF > t0 
    
   t = linspace(t0,tF,new_length+1); 
   % t = t0:dt_dump:tF;

   initial_ffv = ff1_value_u(1);
   figure(20)
   hold off
   plot(t,ff1_value_u,'-r',t,ff2_value_u,'-b',t,ff3_value_u,'-k')
   legend('ff1','ff2','ff3')
   hold on
   plot([t0 tF],[initial_ffv initial_ffv],':k')
   plot([t0 tF],[expected_b expected_b],':b')
   title('Far-field values: all 3 zones')
   xlabel('t')
   ylabel('ffv')

   min_var_mass = min(variation_mass);
   max_var_mass = max(variation_mass);
   if max_var_mass < 0 
      max_var_mass = 0;
   else
      max_var_mass = 1.10*max_var_mass; 
   end
   if min_var_mass > 0
      min_var_mass = 0;
   else
      min_var_mass = 1.10*min_var_mass; 
   end    
   figure(30)
   hold off
   plot(t,variation_mass,'-r')
   hold on
   plot([t0 tF],[0 0],':k')
   title('MASS of  u - uprevious')
   axis([t0 tF min_var_mass max_var_mass])
   xlabel('t')
   ylabel('variation mass')

   t0_picture = t0;
   max_var_sup = max(variation_sup);
   if ( tF > 7.5 ) 
      t0_picture = min([5,tF]);
      ii = max( find( t < t0_picture ) );
      max_var_sup = max(variation_sup(ii:end));
   end   
   figure(35)
   hold off
   plot(t,variation_sup,'-r')
   hold on
   plot([t0 tF],[0 0],':k')
   axis([t0_picture tF 0 1.10*max_var_sup])
   title('SUP of  u - uprevious')
   xlabel('t')
   ylabel('variation sup')

   figure(40)
   plot(t,min_u,'-r',t,max_u,'-b')
   axis([t0 tF min_u(1)-1 max_u(1)+1])
   title('MIN and MAX values of u over the grid')
   xlabel('t')
   ylabel('min & max')
   legend('min','max')

   figure(50)
   plot(t(2:end),time_per_cycle(2:end),'-r')
   title('Time per cycle (seconds)')
   xlabel('t')
   ylabel('elapsed time')

end 











