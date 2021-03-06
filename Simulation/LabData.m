% lab data
% Yuan Sun
% 2015-11-24

clc
clear

%% Real lab data 
% Step response with no control

% Sample   Time     Commanded Pos   Encoder 1 Pos   Encoder 2 Pos 

data = [    0    0.000            1000               0               0;
     1    0.009            1000              25               0;
     2    0.018            1000             101               0;
     3    0.027            1000             230               0;
     4    0.035            1000             399               0;
     5    0.044            1000             564              27;
     6    0.053            1000             642              41;
     7    0.062            1000             606              71;
     8    0.071            1000             534             105;
     9    0.080            1000             479             141;
    10    0.089            1000             451             176;
    11    0.097            1000             452             211;
    12    0.106            1000             486             244;
    13    0.115            1000             548             276;
    14    0.124            1000             633             307;
    15    0.133            1000             736             338;
    16    0.142            1000             848             367;
    17    0.151            1000             956             395;
    18    0.159            1000            1045             426;
    19    0.168            1000            1099             458;
    20    0.177            1000            1116             491;
    21    0.186            1000            1111             524;
    22    0.195            1000            1089             556;
    23    0.204            1000            1056             587;
    24    0.213            1000            1017             617;
    25    0.221            1000             973             644;
    26    0.230            1000             929             670;
    27    0.239            1000             886             695;
    28    0.248            1000             850             717;
    29    0.257            1000             823             733;
    30    0.266            1000             811             745;
    31    0.274            1000             816             762;
    32    0.283            1000             844             776;
    33    0.292            1000             887             787;
    34    0.301            1000             939             797;
    35    0.310            1000             995             805;
    36    0.319            1000            1050             811;
    37    0.328            1000            1101             816;
    38    0.336            1000            1143             819;
    39    0.345            1000            1174             821;
    40    0.354            1000            1192             822;
    41    0.363            1000            1196             822;
    42    0.372            1000            1185             822;
    43    0.381            1000            1159             821;
    44    0.390            1000            1122             819;
    45    0.398            1000            1076             816;
    46    0.407            1000            1026             812;
    47    0.416            1000             976             808;
    48    0.425            1000             929             803;
    49    0.434            1000             889             797;
    50    0.443            1000             862             789;
    51    0.452            1000             850             781;
    52    0.460            1000             853             771;
    53    0.469            1000             871             761;
    54    0.478            1000             897             750;
    55    0.487            1000             929             739;
    56    0.496            1000             964             727;
    57    0.505            1000             998             716;
    58    0.514            1000            1030             705;
    59    0.522            1000            1056             694;
    60    0.531            1000            1076             684;
    61    0.540            1000            1087             674;
    62    0.549            1000            1090             664;
    63    0.558            1000            1084             655;
    64    0.567            1000            1068             646;
    65    0.576            1000            1045             637;
    66    0.584            1000            1016             628;
    67    0.593            1000             983             620;
    68    0.602            1000             949             612;
    69    0.611            1000             914             605;
    70    0.620            1000             882             598;
    71    0.629            1000             853             590;
    72    0.638            1000             831             584;
    73    0.646            1000             818             577;
    74    0.655            1000             813             570;
    75    0.664            1000             817             563;
    76    0.673            1000             832             557;
    77    0.682            1000             853             551;
    78    0.691            1000             880             545;
    79    0.699            1000             910             540;
    80    0.708            1000             941             535;
    81    0.717            1000             971             530;
    82    0.726            1000             996             526;
    83    0.735            1000            1014             522;
    84    0.744            1000            1024             519;
    85    0.753            1000            1026             517;
    86    0.761            1000            1019             515;
    87    0.770            1000            1003             515;
    88    0.779            1000             981             514;
    89    0.788            1000             955             514;
    90    0.797            1000             927             514;
    91    0.806            1000             899             514;
    92    0.815            1000             873             514;
    93    0.823            1000             851             514;
    94    0.832            1000             835             514;
    95    0.841            1000             826             514;
    96    0.850            1000             825             514;
    97    0.859            1000             832             514;
    98    0.868            1000             847             514;
    99    0.877            1000             868             514;
   100    0.885            1000             893             514;
   101    0.894            1000             920             514;
   102    0.903            1000             947             514;
   103    0.912            1000             972             514;
   104    0.921            1000             992             514;
   105    0.930            1000            1004             514;
   106    0.939            1000            1010             514;
   107    0.947            1000            1009             514;
   108    0.956            1000             999             514;
   109    0.965            1000             983             514;
   110    0.974            1000             962             514;
   111    0.983            1000             939             514;
   112    0.992            1000             916             514;
   113    1.001            1000             893             514;
   114    1.009               0             852             514;
   115    1.018               0             766             514;
   116    1.027               0             634             514;
   117    1.036               0             466             514;
   118    1.045               0             285             510;
   119    1.054               0             170             497;
   120    1.063               0             184             470;
   121    1.071               0             266             435;
   122    1.080               0             344             397;
   123    1.089               0             403             360;
   124    1.098               0             439             323;
   125    1.107               0             448             287;
   126    1.116               0             428             251;
   127    1.124               0             380             216;
   128    1.133               0             307             181;
   129    1.142               0             216             147;
   130    1.151               0             112             114;
   131    1.160               0               3              81;
   132    1.169               0            -104              48;
   133    1.178               0            -202              16;
   134    1.186               0            -284             -16;
   135    1.195               0            -345             -49;
   136    1.204               0            -380             -81;
   137    1.213               0            -389            -112;
   138    1.222               0            -371            -144;
   139    1.231               0            -329            -175;
   140    1.240               0            -268            -205;
   141    1.248               0            -192            -234;
   142    1.257               0            -108            -262;
   143    1.266               0             -27            -285;
   144    1.275               0              25            -292;
   145    1.284               0              16            -308;
   146    1.293               0             -55            -316;
   147    1.302               0            -145            -323;
   148    1.310               0            -232            -326;
   149    1.319               0            -308            -329;
   150    1.328               0            -367            -331;
   151    1.337               0            -405            -333;
   152    1.346               0            -419            -334;
   153    1.355               0            -409            -335;
   154    1.364               0            -373            -336;
   155    1.372               0            -316            -336;
   156    1.381               0            -243            -336;
   157    1.390               0            -162            -336;
   158    1.399               0             -86            -333;
   159    1.408               0             -36            -326;
   160    1.417               0             -31            -315;
   161    1.426               0             -60            -300;
   162    1.434               0             -97            -284;
   163    1.443               0            -134            -269;
   164    1.452               0            -166            -253;
   165    1.461               0            -189            -237;
   166    1.470               0            -203            -222;
   167    1.479               0            -207            -208;
   168    1.488               0            -198            -194;
   169    1.496               0            -178            -181;
   170    1.505               0            -149            -168;
   171    1.514               0            -113            -156;
   172    1.523               0             -71            -144;
   173    1.532               0             -28            -133;
   174    1.541               0              14            -122;
   175    1.549               0              52            -111;
   176    1.558               0              83            -101;
   177    1.567               0             105             -92;
   178    1.576               0             116             -82;
   179    1.585               0             117             -73;
   180    1.594               0             106             -64;
   181    1.603               0              87             -56;
   182    1.611               0              61             -48;
   183    1.620               0              31             -42;
   184    1.629               0              -1             -35;
   185    1.638               0             -32             -30;
   186    1.647               0             -60             -26;
   187    1.656               0             -83             -22;
   188    1.665               0             -99             -19;
   189    1.673               0            -108             -17;
   190    1.682               0            -111             -15;
   191    1.691               0            -105             -15;
   192    1.700               0             -92             -15;
   193    1.709               0             -73             -15;
   194    1.718               0             -49             -15;
   195    1.727               0             -23             -15;
   196    1.735               0               3             -17;
   197    1.744               0              29             -18;
   198    1.753               0              52             -19;
   199    1.762               0              70             -20;
   200    1.771               0              83             -21;
   201    1.780               0              88             -21;
   202    1.789               0              88             -22;
   203    1.797               0              79             -22;
   204    1.806               0              64             -23;
   205    1.815               0              45             -23;
   206    1.824               0              22             -24;
   207    1.833               0              -1             -24;
   208    1.842               0             -23             -24;
   209    1.851               0             -44             -25;
   210    1.859               0             -60             -25;
   211    1.868               0             -71             -26;
   212    1.877               0             -77             -26;
   213    1.886               0             -79             -27;
   214    1.895               0             -72             -28;
   215    1.904               0             -60             -29;
   216    1.913               0             -44             -29;
   217    1.921               0             -25             -30;
   218    1.930               0              -6             -31;
   219    1.939               0              14             -31;
   220    1.948               0              31             -32;
   221    1.957               0              46             -32;
   222    1.966               0              56             -32;
   223    1.974               0              61             -33;
   224    1.983               0              61             -33;
   225    1.992               0              57             -33];

t = data(: ,2);
y1 = data(:, 4);
y2 = data(:, 5); % Encode2
u = data(:, 3); % Control Position as input
% System identification
zf = iddata(y2, u, 0.00884);
Tfd = tfest(zf,4,0);

figure
plot(t, y2, 'linewidth', 2);
grid;
%axis([0 0.2 0 400]);
