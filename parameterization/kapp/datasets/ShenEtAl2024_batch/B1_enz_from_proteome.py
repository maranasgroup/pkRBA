# finds max ENZSYN fluxes for all enzymes
# update model-specific settings in kapp_options.py
load_model = False
from kapp_options import *

path_gams = '../../../../GAMS/parameterization/enz_from_proteome/'
path_out = './enz_from_proteome/'

run_setting_file_from = './GAMS_setting_files/enz_from_proteome_GAMS_settings.txt'
run_setting_file_to = './enz_from_proteome/enz_from_proteome_GAMS_settings.txt'

#### Create directory and copy run settings
if os.path.isdir(path_out) == False:
    os.makedirs(path_out)
shutil.copy(run_setting_file_from, run_setting_file_to);

#### Process data
with open(os.path.join(path_gams, 'pro_and_enz.txt')) as f:
    pro_list = f.read().split('\n')
pro_list = pro_list[1:-1]
pro_list = [i[1:-1] for i in pro_list]
pro_list = [i for i in pro_list if i.split('-')[0] == 'PRO']
data_pro_name_matches = {}
for k in df_data['id'].to_list():
    new_id = vars(met_ID.from_name(f'{ptmpro_prefix}-{k}'))
    new_id.pop('PTMs', None)
    new_id = str(met_ID(**new_id)).removeprefix(f'{ptmpro_prefix}-')
    data_pro_name_matches.setdefault(new_id, []).append(k)

data = set(); pro_data = set(); pro_nodata = set()
for met in pro_list:
    _,sid = met.split('-', maxsplit=1)
    # match in df_data by sid or by uniprot id (if sid is not uniprot id)
    if sid in data_pro_name_matches:
        pro_data.add("'PROIN-" + sid + "'")
        data.add("'PROIN-" + sid + "' " + str(df_data.loc[df_data['id'].isin(data_pro_name_matches[sid]), 'vtrans (mmol/gDW/h)'].values[0]))
    else:
        pro_nodata.add("'PROIN-" + sid + "'")

data = sorted(data); pro_data = sorted(pro_data); pro_nodata = sorted(pro_nodata)
data = ['/'] + data + ['/']
pro_data = ['/'] + pro_data + ['/']
pro_nodata = ['/'] + pro_nodata + ['/']

# Write out run files
with open(os.path.join(path_out, 'proteome_data.txt'), 'w') as f:
    f.write('\n'.join(data))
with open(os.path.join(path_out, 'rxns_pro_data.txt'), 'w') as f:
    f.write('\n'.join(pro_data))
with open(os.path.join(path_out, 'rxns_pro_nodata.txt'), 'w') as f:
    f.write('\n'.join(pro_nodata))
    
#### Simulation
shutil.copy(os.path.join(path_gams, 'soplex.opt'),
            os.path.join(path_out, 'soplex.opt'));

modelstat_path = 'enz_from_proteome.modelStat.txt' 
fname = os.path.join(path_out, modelstat_path) 
# remove to prevent past runs from overriding status checks
if os.path.isfile(fname):
    os.remove(fname)

cmds = ['cd ' + path_out,
        'module load gams',
        'gams enz_from_proteome.gms' + output_redirect_str]
os.system('\n'.join(cmds))
stop_if_run_failed(fname)
