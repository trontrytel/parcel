% initial condition
outfile = 'test_dan.nc';
RH_init = .99;
T_init  = 300.;
p_init  = 101300.;

% recalculate RH to water vapor mixing ratio
r_init  = py.functions.rh_to_rv(RH_init, T_init, p_init);

% parcel model options
cmd = ['python parcel.py', ...
   ... % paramteres common to all microphysics
   ' --outfile ',       outfile, ...
   ' --p_0 ',           num2str(p_init), ...
   ' --T_0 ',           num2str(T_init), ...
   ' --r_0 ',           num2str(r_init), ...
   ' --dt ',            num2str(.1), ...
   ' --w ',             num2str(1.), ...
   ' --z_max ',         num2str(200), ...
   ' --sd_conc ',       num2str(100), ...
   ' --kappa ',         num2str(.5), ...
   ' --n_tot ',         num2str(160e6), ...
   ' --mean_r ',        num2str(.04e-6), ...
   ' --gstdev ',        num2str(1.4), ...
   ' --out_bin ',       '''{"aerosol": {"rght": 1e-3, "moms": [0,3], "drwt": "dry", "nbin": 50, "lnli": "log", "left": 1e-9}, "n0": {"rght": 1e-3, "moms": [0], "drwt": "wet", "nbin": 50, "lnli": "log", "left": 2e-6}  }''', ...
   ' --stop_at_RHmax ', num2str(1), ...
   ' 2>run_parcel.stderr.log' ...
];

disp(cmd);

% run parcel model
[status, result] = unix(cmd);
if (status ~= 0)
    unix(['echo "', cmd, '" >> fail.log']) % log commands which failed
    throw(MException('run_parcel:fatal', ['parcel.py failed with exit code ', num2str(status), ' any model output might be in ', outfile]));
end;
  
% read output:
rhod = ncreadatt(outfile, '/', 'rhod_RHmax'); %dry air density

%TODO - change parcel - now we save at RHmax + 1 timestep

% 1.  S_max [%]
RHmax = ncreadatt(outfile,'/', 'RHmax');
S_max = 100*(RHmax-1);
disp(S_max);

% 2.  N_0 [1/cm3]
tmp = ncread(outfile, 'n0_m0');  
N_0 = sum(tmp(:,2)) * rhod * 1e-6;
disp(N_0);

% 3.  N_1
disp('TODO');

% 4.  N_2 [1/cm3]
N_tmp = ncreadatt(outfile,'/', 'N_act_at_RH_max');
N_2 = N_tmp * rhod * 1e-6;
disp(N_2);

% 5.  N_3
disp('TODO');

% 6.  mass accommodation coeff
%     it's hardcoded in the library and is equal 1
disp(1);

% 7.  vertical velocity [m/s]
w = ncreadatt(outfile,'/','w');
disp(w);

% 8.  z_RHmax - z_LCL 
% TODO - add to attributes
% TODO - how is LCL defined?
z = ncread(outfile, 'z');
z_RHmax = z(2);
disp(z_RHmax);

% 9.  mass concentration of soluble aerosol  [ug of aerosol mass / m3 of air]
% TODO - add to attributes
% TODO - this will be the same as in the initial condition
tmp = ncread(outfile, 'aerosol_m3');
chem_rho = ncreadatt(outfile,'/', 'chem_rho');
mass_a = sum(tmp(:,2)) * rhod * chem_rho * 4/3. * pi * 1e9;
disp(mass_a);

% 10. aerosol number concentration [1/cm3]
%TODO - check converting o STP for initial condition for N_tot for aerosols
%TODO - add to attributes
%TODO - why is it not 60?
%TODO - this will be the same as initial condition?
tmp = ncread(outfile, 'aerosol_m0');  
N_a = sum(tmp(:,2)) * rhod * 1e-6;
disp(N_a);

%  B0 = kvp{2}(find(ismember(kvp{1}, 'S_max_B0')));
%  RH = kvp{2}(find(ismember(kvp{1}, 'S_max_RH')));
%  M0 = kvp{2}(find(ismember(kvp{1}, 'S_max_M0')));
%  rhod = kvp{2}(find(ismember(kvp{1}, 'S_max_rhod')));
