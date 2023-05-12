import sys
import yaml

print (sys.argv[0], len(sys.argv))

dataname = sys.argv[1]
config_file_path = sys.argv[2]
uname = sys.argv[3]

with open(config_file_path, 'r') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
    config['homedir'] = config['homedir'].replace('username', uname)
    config['dataname'] = dataname

# dump yaml
with open('./build/config.yaml', 'w') as f:
    yaml.dump(config, f)

