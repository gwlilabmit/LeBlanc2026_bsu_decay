import sys

## path to directory with raw data to align
directory_file=sys.argv[1]

## name of config file to output
config_file=sys.argv[2]

def main():

	with open(directory_file, 'r') as f:
		info=f.read().strip().split()

	with open(config_file, 'w') as o:
		o.write('ArrayTaskID,FastqID\n')

		for i in enumerate(info):
			index=i[0]+1
			name=i[1]
			file_id=name.split('.')
			o.write(str(index)+','+file_id[0]+'\n')

if __name__=='__main__':
	main()