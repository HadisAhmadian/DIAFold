########################################################
with open('./end2end_T20.msa', 'rb') as f:
#####################################################
    data = f.read()
    files = data.split(b'\0')
    print(len(files))
    for i, file in enumerate(files[:-1]):
        s=file[1:file.find(b'\n')].decode('utf-8')
#####################################################
        with open(f'./mmout1/'+s+'.a3m', 'wb') as new_file:
#####################################################
            new_file.write(file)