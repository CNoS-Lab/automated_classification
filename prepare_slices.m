function slices = prepare_slices()

% prepare all the slices we're intrested in template matching for

% Coronal slices: 1-15
% Axial slices: 16-37

% C 1-4, A 16 (1) 
L_cog_eval = {'Bilateral eyeball sitters','Bilateral space invaders shooters','Above the line (caudate)','Sad face antennae & flushed cheeks','X marks the spot'};
C_cog_eval = {{170,174,178,182,186,190},{136,140,144,148,152},{116,120,124,128,132},{76,80,84,88,92}};
A_cog_eval = {{22,26,30,34}};

% C 5-6, A 17-18 (2-3)
L_ling_proc = {'Rail shot coronal','Tears blown leftwards & eyebrows','Rail shot axial','Disappearing face (start right disappear left)'};
C_ling_proc = {{134,138,142,146,150,154},{46,52,60,68,76,84,92,100}};
A_ling_proc = {{64,68,72,76,80},{34,38,42,46,50,54}};

% A 19-22 (4-7)
L_int_attn = {'Left-lateralized upper triangle','Left-lateralized lower triangle','Right-handed crab claw','Found a peanut'};
A_int_attn = {{136,128,120,116},{100,104,108,112},{24,32,40,48},{72,78,82,86,90}};

% C 7-8, A 23-24 (8-9)
L_ext_attn = {'Jumping jack flash','Ape nostrils','Flexing hands','Wipe your mouth bear triple jam'};
C_ext_attn = {{136,144,152,160},{114,108,102,96}};
A_ext_attn = {{118,108,98,88,78},{38,42,46,50}};

% C 9, A 25-26 (10-11)
L_init = {'Raised eyebrows','When I''m 64','De Divina Proportione front guy'};
C_init = {{52,57,62,67,72}};
A_init = {{126,136,146},{132,127,122,117,112}};

% C 10-11, A 27-28 (12-13)
L_resp = {'Bat (one sided if one-handed response)','Thalamus kite surfer','Butterfly (one sided if one-handed response)','Compact crab claw'};
C_resp = {{116,126,136,146},{120,115,110,105,100}};
A_resp = {{144,134,124,114},{56,52,48,44}};

% C 12-14, A 29-31 (14-16)
L_dmn = {'Snowman nose vs mouth','Medial temporal dots - prominent vs muted','T bird vs stickman','Tripod vs baby dragon','Mandibles vs laughing clown','Angel wings - muted vs prominent'};
C_dmn = {{148,153,158,163},{112,116,120,124,128},{40,50,60,70,80,90}};
A_dmn = {{66,76,86,96,106},{25,30,35,40},{52,62,72,82,92}};

% A 32-33 (17-18)
L_aud = {'Thing 1','Thing 2'};
A_aud = {{50,60,70,80,90,100},{45,50,55,60,65}};

% A 34-36 (19-21)
L_aar = {'Happy 28th birthday long face/right angle','On fire','Small smile'};
A_aar = {{100,110,120,130,140},{65,75,85,95},{54,58,62,66,70}};

% C 15, A 37 (22)
L_fvf = {'Stay puft','Wishbone'};
C_fvf = {{37,42,47,52,57}};
A_fvf = {{57,62,67,72,77}};

slices(9) = struct();

slices(1).N = {'Cognitive Evaluation'};
slices(1).L = L_cog_eval;
slices(1).C = C_cog_eval;
slices(1).A = A_cog_eval;

slices(2).N = {'Lingustic Processing'};
slices(2).L = L_ling_proc;
slices(2).C = C_ling_proc;
slices(2).A = A_ling_proc;

slices(3).N = {'Internal Attention'};
slices(3).L = L_int_attn;
slices(3).C = [];
slices(3).A = A_int_attn;

slices(4).N = {'External Attention'};
slices(4).L = L_ext_attn;
slices(4).C = C_ext_attn;
slices(4).A = A_ext_attn;

slices(5).N = {'Initiation'};
slices(5).L = L_init;
slices(5).C = C_init;
slices(5).A = A_init;

slices(6).N = {'Response single vs double'};
slices(6).L = L_resp;
slices(6).C = C_resp;
slices(6).A = A_resp;

slices(7).N = {'Default Mode Networks Traditional vs Novel'};
slices(7).L = L_dmn;
slices(7).C = C_dmn;
slices(7).A = A_dmn;

slices(8).N = {'Primary Auditory'};
slices(8).L = L_aud;
slices(8).C = [];
slices(8).A = A_aud;

slices(9).N = {'Auditory Attention for Response'};
slices(9).L = L_aar;
slices(9).C = [];
slices(9).A = A_aar;

slices(10).N = {'Focus on Visual Features'};
slices(10).L = L_fvf;
slices(10).C = C_fvf;
slices(10).A = A_fvf;

end


