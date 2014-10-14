
AU = 1.495978707E11;
w = 60*60*24*7;
G = 6.67E-11;
m_s = AU^3/(G*w^2);
N = 11;

m = 1/m_s*[1.9891E30 ... %sun
    3.302E23 4.8685E24 5.97219E24 6.4185E23 1.8986E27 5.6846E26...
    8.681E25 1.0243E26 1.3E22...%pluto
    2.2E14 %Halley
    ];

r_0 = 1000/AU*[0 0 0 ...
    2.417458149944681E+07  3.961908823472586E+07  1.019152972912366E+06...  %sun
   -8.713874815025689E+07  6.267133929263873E+07  5.887624863078220E+06...
    -9.832893329273786E+07  1.098127067722844E+08 -3.257632516971247E+03...
    -2.437981627448902E+08  4.819784790635903E+07  6.994082893304829E+06...
    -2.331512428061480E+08  7.427122742511928E+08  2.132579853305100E+06...
    -1.012749850225338E+09 -1.076436562107142E+09  5.903530065123624E+07...
    2.935069311018662E+09  6.038275159130698E+08 -3.577124383069587E+07...
    4.054749746608424E+09 -1.915367119297712E+09 -5.398515950879561E+07...
    9.511822513264761E+08, -4.776354464417923E+09,  2.359614907267640E+08... %pluto
    -3.057669847992851E+09,  3.725674239897915E+09, -1.455042686200472E+09... %Halley
    ];

v_0 = 1000*w/AU*[0 0 0 ...
    -5.129420747164487E+01  2.732245601122905E+01  6.938641077916809E+00...  %sun
    -2.058556491270471E+01 -2.860628332870155E+01  7.959826095832389E-01...
    -2.267685408359372E+01 -1.998317369161212E+01  6.841972831700973E-04...
    -3.793312925984044E+00 -2.169904870654408E+01 -3.615448792248455E-01...
    -1.263517979410647E+01 -3.297975401056563E+00  2.964314434433838E-01...
    6.500277929864231E+00 -6.647002225408530E+00 -1.433612204947691E-01...
    -1.431512921791133E+00  6.348399029861570E+00  4.225470418730855E-02...
    2.275369579533704E+00  4.942024008970557E+00 -1.543036821871557E-01...
    5.436226188111175E+00 -3.630481576334465E-02 -1.568594107858565E+00... %pluto
    -1.813952903438697E-01  1.709843170334322E+00 -3.330679308312253E-01... %halley
    ];


bodies = zeros(N, 6);

for i = 1:N

bodies(i, 1:3) = r_0(3*i-2:3*i);
bodies(i, 4:6) = v_0(3*i-2:3*i); 

end

bodies = [m' bodies];

dlmwrite('./initial_conditions/solarsystem11.csv', bodies);
    
    
    
                           
                            