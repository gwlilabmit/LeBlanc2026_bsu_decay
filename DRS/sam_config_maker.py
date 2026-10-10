import sys

directory_file=sys.argv[1]
config_file=sys.argv[2]

def main():

	with open(directory_file, 'r') as f:
		info=f.read().strip().split()
		
		id_list=[]
		for i in info:
			if i[0]=='R':
				id_list.append(i)

	with open(config_file, 'w') as o:
		o.write('ArrayTaskID,SAMID\n')
		for i in enumerate(id_list):
			index=i[0]+1
			name=i[1]
			o.write(str(index)+','+name+'\n')

if __name__=='__main__':
	main()