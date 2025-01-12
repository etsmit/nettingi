import sys,os
import glob



infiles = glob.glob('*.cands')

for filenm in infiles:
	print(filenm)
	go=False
	out_temp = filenm+'.temp'
	if os.path.isfile(out_temp):
		os.system(f'rm {out_temp}')
	out = open(out_temp,'a')
	with open(filenm) as candfile:
		lines = candfile.readlines()
		line_to_read = 0
		while line_to_read < len(lines):
			line = lines[line_to_read]
			
			if line.startswith("#"):
				go = True
				print(f'found # on line {line_to_read}')
				
		
			if go:
				out.write(line)

			line_to_read += 1
				

	out.close()
	os.system(f'cp {out_temp} {filenm}')
	os.system(f'rm {out_temp}')

			
