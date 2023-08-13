import sys
import yaml

print (sys.argv[0], len(sys.argv))

dataname = sys.argv[1]
config_file_path = sys.argv[2]
uname = sys.argv[3]
datapath = sys.argv[4]
codepath = sys.argv[5]

print('datapath', datapath)
print('codepath', codepath)

with open(config_file_path, 'r') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
    config['homedir'] = config['homedir'].replace('username', uname)
    config['dataname'] = dataname
    config['datapath'] = datapath
    config['codepath'] = codepath

# dump yaml
with open(f'{datapath}/build/config.yaml', 'w') as f:
    yaml.dump(config, f)

