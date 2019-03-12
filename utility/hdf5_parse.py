import h5py
import numpy as np
import pandas as pd

#baseDir = '/Users/bhofmeister/Documents/Research/other/'
baseDir = '/lustre1/bth29393/share/'

f = h5py.File(baseDir+'imputed_snps_binary.hdf5','r')
positions = np.array(f['positions']) # ndarray (10709949,) dtype=int32
### Chromosome regions: 
### (chr1) [ 0, 2597825 ],
### (chr2) [ 2597825, 4466694 ],
### (chr3) [ 4466694, 6661059 ],
### (chr4) [ 6661059, 8428147 ],
### (chr5) [ 8428147, 10709949 ]
chr1 = np.full(2597825-0,'chr1',dtype='U4')
chr2 = np.full(4466694-2597825,'chr2',dtype='U4')
chr3 = np.full(6661059-4466694,'chr3',dtype='U4')
chr4 = np.full(8428147-6661059,'chr4',dtype='U4')
chr5 = np.full(10709949-8428147,'chr5',dtype='U4')
chrms = np.concatenate((chr1,chr2,chr3,chr4,chr5), axis=0)
accessions = np.array(f['accessions']) # shape (1135,) dtype=S5
accessions = accessions.astype('U5')
snps=np.array(f['snps']) # shape (10709949, 1135) dtype=int8

data = pd.DataFrame(snps)
data.columns = accessions
data['chrm'] = chrms
data['pos'] = positions
data.set_index(['chrm','pos'], inplace=True)

dataG = data.groupby(level='chrm')

for name,group in dataG:
	print(name)
	outFile = baseDir + name + '_snps.csv'
	group.to_csv(outFile)
