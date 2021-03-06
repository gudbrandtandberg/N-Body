%% Part 1 - Solar System (18 Bodies)
% Use positions, masses and velocities of all the planets in the solar
% system to write intial files. Has been written to write n = 3, 6, 13 18
% bodies to file. All the bodies are:
% 
% 1       Sun
% 2       Mercury
% 3       Venus
% 4       Earth
% 5       Mars
% 6       Jupiter
% 7       Saturn
% 8       Neptune
% 9       Uranus
% 10      Pluto
% 11      Moon
% 12      Phobos
% 13      Deimos
% 14      Io
% 15      Europa
% 16      Ganymede
% 17      Callisto
% 18      Halley's comet

% Unit scaling. We use astronomical units, weeks and set G = 1 in the
% siulations. This determines the characteristic mass; m_s. 
AU = 1.495978707E11;   % astronomical units
w = 60*60*24*7;        % weeks
G = 6.67E-11;          % Gravitational constant
m_s = AU^3/(G*w^2);    % characteristic mass

% File-writing options
N = 18;                % # of bodies in total
n = 11;                % # of bodies to write to file

% Masses of the 18 bodies in units of m_s 
m = 1/m_s*[1.9891E30 ... %sun
    3.302E23 4.8685E24 5.97219E24 6.4185E23 1.8986E27 5.6846E26...
    8.681E25 1.0243E26 1.3E22...%pluto
    7.3E22...%moon
    1.08E20...%phobos
    1.8E20...%deimos
    8.9E22...%io
    4.8E22...
    1.4819E23...
    1.0759E23...%callisto
    2.2E14 %Halley
    ];

% Initial positions in astronomical units
r_0 = 1000/AU*[0, 0, 0 ...                                              %sun
4.325741590616594E+07, -4.447598205227128E+07, -7.602793726150981E+06...%mercury
-1.065184862973610E+08,  1.426409267960166E+07,  6.342590226882618E+06...
1.484997959420039E+08,  1.966516392815287E+07, -1.414904880365577E+03...
7.561157938277872E+07, -1.979342346107423E+08, -6.003106128253528E+06 ...
-4.781071082144656E+08,  6.306239474601380E+08,  8.079395665674826E+06...
-8.681948970465746E+08, -1.205449819868287E+09,  5.552046389307594E+07 ...
2.902003974098086E+09,  7.358295383909013E+08, -3.484989767043513E+07...
4.101038351114639E+09, -1.811465478741618E+09, -5.719412717146279E+07...
1.063767728198170E+09, -4.775732770106762E+09,  2.032610886473450E+08... %pluto
1.484808788426344E+08,  1.928964235278643E+07,  3.084322958738403E+04... %moon
7.561951745378238E+07, -1.979307087677598E+08, -6.006860788618978E+06... %phobos  
7.559220083042181E+07, -1.979445171572549E+08, -5.994788581249264E+06... %deimos
-4.780994593819707E+08,  6.302018642958680E+08,  8.064344034494819E+06...%io
 -4.786656648485875E+08,  6.309900741640583E+08,  8.088050878001666E+06...
-4.791763379696921E+08,  6.307044425781475E+08,  8.068481493980590E+06... 
-4.762581559059567E+08,  6.303505565910144E+08,  8.095657797811498E+06...%callisto
-3.061126454498667E+09,  3.760634554053051E+09, -1.461631754870368E+09...%Halley
];

% Initial velocities in AU/week
v_0 = 1000*w/AU*[0, 0, 0 ...                                            %sun
2.525822913860079E+01,  3.627931518613809E+01,  6.469242537532837E-01...%mercury
-4.838066659924075E+00, -3.486533967772683E+01, -1.986418837879511E-01...
-4.406220350009290E+00,  2.941898165409254E+01, -1.333034560677749E-03...
2.355273541766074E+01,  1.072888188164977E+01, -3.532953957527511E-01...
-1.057863541952222E+01, -7.282677815690089E+00,  2.669653552828916E-01...
7.305583007781340E+00, -5.678067175058013E+00, -1.921159722403775E-01...
-1.730861423436294E+00,  6.275826934901269E+00,  4.597420467011978E-02...
2.152788365068528E+00,  4.996563864711828E+00, -1.527584614750138E-01...
5.389758769831626E+00,  5.017824639329278E-02, -1.578243168379842E+00...
-3.368658078218818E+00,  2.941691744785754E+01,  3.022524496536850E-02...%moon
2.293361231296165E+01,  1.269753491616835E+01,  1.228783448035231E-01...%phobos
2.404367682320214E+01,  9.520073139020946E+00, -7.039594982150323E-01...%deimos
6.725238795474569E+00, -6.911890224983383E+00,  5.430875771496093E-01...%io
-1.805503630424283E+01, -1.887336981626284E+01, -3.448311375533070E-01...
-1.137957128164578E+01, -1.810663574007972E+01, -1.515591072689991E-01...
-9.383472561750258E+00,  8.878262368908618E-01,  5.433379097068503E-01...%callisto
-1.236183883191173E-01,  1.621369229744805E+00, -3.015245553639735E-01...%Halley
];

bodies = zeros(N, 6);

for i = 1:N

bodies(i, 1:3) = r_0(3*i-2:3*i);
bodies(i, 4:6) = v_0(3*i-2:3*i); 

end

% bodies now contains ALL the information above. 
% the index vector I chooses which to write to file.

bodies = [m(1:N)' bodies];

if n == 3  % sun - earth - moon
    I = [1 4 11];
end
if n == 6  % sun - merc - venus - earth - mars - jupiter
    I = [1 2 3 4 5 6];
end
if n == 11
    I = [1 2 3 4 5 6 11 14 15 16 17];
end
if n == 13 % sun -> jupiter - moon - phobos - Galilleian moons
   I = [1 2 3 4 5 6 11 12 13 14 15 16 17];
end
if n == 18 % all bodies
    I = 1:18;
end

dlmwrite(sprintf('../initial_conditions/solarsystem%d.dat', n),...
    bodies(I,:), 'delimiter', ' ');
    
%% Part 2 - Open cluster
% Generate N random intial conditions for simulating cold collapse of a
% galaxy.

N = 1000;    % number of 'bodies' in simulation
R0 = 20;    % pre-collapse radius in AU

bodies = zeros(N, 7);

% Fill bodies with random masses and positions and zero velocities

for i = 1:N

    % Normal distributed masses (use solar masses as characteristic mass)
    m = normrnd(10, 1);
    
    % Uniformly distributed positions in sphere
    x = -R0 + 2*R0*rand();    
    ylim = sqrt(R0^2 - x^2);
    y = -ylim + 2*ylim*rand();    
    zlim = sqrt(R0^2 - y^2 - x^2);
    z = -zlim + 2*zlim*rand();
   
    bodies(i,:) = [m x y z 0 0 0];
    
end

dlmwrite(sprintf('../initial_conditions/cluster%d.dat', N), bodies, ...
    'delimiter', ' ');

                            