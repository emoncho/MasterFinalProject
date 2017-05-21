#PREDICTION OF GENE INTERACTIONS USING 'FASTING' DATASET (https://www.ncbi.nlm.nih.gov/sites/GDSbrowser / GDS5473: Short-term fasting effect on skeletal muscle: time course)

#GET ROW DATA AND CALCULATE CORRELATION MATRIX

from pandas import *
import numpy as np
from scipy.stats.stats import pearsonr
import itertools
from sklearn.preprocessing import Normalizer
import os

filefast100 = os.path.abspath('Fasting/fasting100gens.csv')   #ID from 100 first gens
fast100 = pandas.read_csv(filefast100, header=None, index_col=0)
    
filefast = os.path.abspath("Fasting/fasting.csv")             #genetic expression info all the genes
fast = pandas.read_csv(filefast, index_col = 0)

complet100 = fast[fast.index.isin(list(fast100.index))]       #genetic expression info 100 first genes 
complet100.head()
from sklearn.preprocessing import Normalizer
import numpy as np

df1 = complet100.T
df2 = np.log2(df1)
normalizeddf2 = Normalizer().fit_transform(df2)
matriucorrelacio = pandas.DataFrame(normalizeddf2).corr(method='pearson')
matriucorrelacio.index = [fast100.index]
matriucorrelacio.columns = [fast100.index]
matriucorrelacio.head()                                          #correlation matrix


#ACO MODIFIED (ACOGeneInteraction): It have been used increasing 'times' sequentially (from 1 to 20)
from ACOGeneInteraction import *
world1 = WorldGI(matriucorrelacio, times = 1)#Input is a correlation matrix. It also has a 'times' parameters to choose how many times each nodes is visited.
print(world1.nodes) #nodes are numbered
print(world1.node_names) #each gene in the correlation matrix has to be a unique name.

solver = SolverGI()
solution = solver.solve(world1)
print(solution.interactions) #List of interactions predicted (duplicated values have been removed)
print(solution.tour_interactions) #List of interactions predicted (duplicated values have NOT been removed)
print(solution.tour) #List of nodes visited in order
print(solution.distance) #Distance value of the complete optimal tour, we are not interested in this value.



#DATAFRAME WITH ALL THE PREDICTED INTERACTIONS ('times' = 1,..,20)
df = pandas.DataFrame() 
file20 = os.path.abspath('Fasting/times/20time.csv')
time20 = pandas.read_csv(file20, index_col=0, header=0)
df['time20'] = time20['Gen1'] + '_'+ time20['Gen2']
file19 = os.path.abspath('Fasting/times/19time.csv')
time19 = pandas.read_csv(file19, index_col=0, header=0)
df['time19'] = time19['Gen1'] + '_'+ time19['Gen2']
file18 = os.path.abspath('Fasting/times/18time.csv')
time18 = pandas.read_csv(file18, index_col=0, header=0)
df['time18'] = time18['Gen1'] + '_'+ time18['Gen2']
file17 = os.path.abspath('Fasting/times/17time.csv')
time17 = pandas.read_csv(file17, index_col=0, header=0)
df['time17'] = time17['Gen1'] + '_'+ time17['Gen2']
file16 = os.path.abspath('Fasting/times/16time.csv')
time16 = pandas.read_csv(file16, index_col=0, header=0)
df['time16'] = time16['Gen1'] + '_'+ time16['Gen2']
file15 = os.path.abspath('Fasting/times/15time.csv')
time15 = pandas.read_csv(file15, index_col=0, header=0)
df['time15'] = time15['Gen1'] + '_'+ time15['Gen2']
file14 = os.path.abspath('Fasting/times/14time.csv')
time14 = pandas.read_csv(file14, index_col=0, header=0)
df['time14'] = time14['Gen1'] + '_'+ time14['Gen2']
file13 = os.path.abspath('Fasting/times/13time.csv')
time13 = pandas.read_csv(file13, index_col=0, header=0)
df['time13'] = time13['Gen1'] + '_'+ time13['Gen2']
file12 = os.path.abspath('Fasting/times/12time.csv')
time12 = pandas.read_csv(file12, index_col=0, header=0)
df['time12'] = time12['Gen1'] + '_'+ time12['Gen2']
file11 = os.path.abspath('Fasting/times/11time.csv')
time11 = pandas.read_csv(file11, index_col=0, header=0)
df['time11'] = time11['Gen1'] + '_'+ time11['Gen2']
file10 = os.path.abspath('Fasting/times/10time.csv')
time10 = pandas.read_csv(file10, index_col=0, header=0)
df['time10'] = time10['Gen1'] + '_'+ time10['Gen2']
file9 = os.path.abspath('Fasting/times/9time.csv')
time9 = pandas.read_csv(file9, index_col=0, header=0)
df['time9'] = time9['Gen1'] + '_'+ time9['Gen2']
file8 = os.path.abspath('Fasting/times/8time.csv')
time8 = pandas.read_csv(file8, index_col=0, header=0)
df['time8'] = time8['Gen1'] + '_'+ time8['Gen2']
file7 = os.path.abspath('Fasting/times/7time.csv')
time7 = pandas.read_csv(file7, index_col=0, header=0)
df['time7'] = time7['Gen1'] + '_'+ time7['Gen2']
file6 = os.path.abspath('Fasting/times/6time.csv')
time6 = pandas.read_csv(file6, index_col=0, header=0)
df['time6'] = time6['Gen1'] + '_'+ time6['Gen2']
file5 = os.path.abspath('Fasting/times/5time.csv')
time5 = pandas.read_csv(file5, index_col=0, header=0)
df['time5'] = time5['Gen1'] + '_'+ time5['Gen2']
file4 = os.path.abspath('Fasting/times/4time.csv')
time4 = pandas.read_csv(file4, index_col=0, header=0)
df['time4'] = time4['Gen1'] + '_'+ time4['Gen2']
file3 = os.path.abspath('Fasting/times/3time.csv')
time3 = pandas.read_csv(file3, index_col=0, header=0)
df['time3'] = time3['Gen1'] + '_'+ time3['Gen2']
file2 = os.path.abspath('Fasting/times/2time.csv')
time2 = pandas.read_csv(file2, index_col=0, header=0)
df['time2'] = time2['Gen1'] + '_'+ time2['Gen2']
file1 = os.path.abspath('Fasting/times/1time.csv')
time1 = pandas.read_csv(file1, index_col=0, header=0)
df['time1'] = time1['Gen1'] + '_'+ time1['Gen2']
df


#GRAPH NUMBER OF INTERACTIONS PREDICTED USING DIFFERENT 'times'
import matplotlib.pyplot as plt
f = plt.figure()
plt.plot([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20], 
         [len(time1), len(time2),len(time3), len(time4),len(time5),len(time6),len(time7),
          len(time8),len(time9),len(time10),len(time11), len(time12), len(time13), len(time14),
          len(time15), len(time16), len(time17),len(time18), len(time19), len(time20)],'ro')
plt.title('Number of interactions predicted')
plt.xlabel('times = n')
plt.ylabel('Number of interactions')
plt.show()



#LIST NUMBER OF DIFFERENT INTERACTIONS PREDICTED SO FAR USING DIFFERENT 'times'
alltheinteractions1 = list(df['time1'])
alltheinteractions2 = list(df['time2']) + alltheinteractions1
alltheinteractions3 = list(df['time3']) + alltheinteractions2
alltheinteractions4 = list(df['time4']) + alltheinteractions3
alltheinteractions5 = list(df['time5']) + alltheinteractions4
alltheinteractions6 = list(df['time6']) + alltheinteractions5
alltheinteractions7 = list(df['time7']) + alltheinteractions6
alltheinteractions8 = list(df['time8']) + alltheinteractions7
alltheinteractions9 = list(df['time9']) + alltheinteractions8
alltheinteractions10 = list(df['time10']) + alltheinteractions9
alltheinteractions11 = list(df['time11']) + alltheinteractions10
alltheinteractions12 = list(df['time12']) + alltheinteractions11
alltheinteractions13 = list(df['time13']) + alltheinteractions12
alltheinteractions14 = list(df['time14']) + alltheinteractions13
alltheinteractions15 = list(df['time15']) + alltheinteractions14
alltheinteractions16 = list(df['time16']) + alltheinteractions15
alltheinteractions17 = list(df['time17']) + alltheinteractions16
alltheinteractions18 = list(df['time18']) + alltheinteractions17
alltheinteractions19 = list(df['time19']) + alltheinteractions18
alltheinteractions20 = list(df['time20']) + alltheinteractions19
print('Number of different interactions 1:',len(set(alltheinteractions1)))
print('Number of different interactions 1-2:',len(set(alltheinteractions2)))
print('Number of different interactions 1-3:',len(set(alltheinteractions3)))
print('Number of different interactions 1-4:',len(set(alltheinteractions4)))
print('Number of different interactions 1-5:',len(set(alltheinteractions5)))
print('Number of different interactions 1-5:',len(set(alltheinteractions6)))
print('Number of different interactions 1-5:',len(set(alltheinteractions7)))
print('Number of different interactions 1-8:',len(set(alltheinteractions8)))
print('Number of different interactions 1-9:',len(set(alltheinteractions9)))
print('Number of different interactions 1-10:',len(set(alltheinteractions10)))
print('Number of different interactions 1-11:',len(set(alltheinteractions11)))
print('Number of different interactions 1-12:',len(set(alltheinteractions12)))
print('Number of different interactions 1-13:',len(set(alltheinteractions13)))
print('Number of different interactions 1-14:',len(set(alltheinteractions14)))
print('Number of different interactions 1-15:',len(set(alltheinteractions15)))
print('Number of different interactions 1-16:',len(set(alltheinteractions16)))
print('Number of different interactions 1-17:',len(set(alltheinteractions17)))
print('Number of different interactions 1-18:',len(set(alltheinteractions18)))
print('Number of different interactions 1-19:',len(set(alltheinteractions19)))
print('Number of different interactions 1-20:',len(set(alltheinteractions20)))



#GRAPH NUMBER OF DIFFERENT INTERACTIONS PREDICTED USING 'times' SEQUENTIALLY (1,..,20)
import matplotlib.pyplot as plt
f = plt.figure()
plt.plot([20,19,18,17,16,15,14,13, 12,11,10,9,8,7,6,5,4,3,2,1], 
         [len(set(alltheinteractions20)),
          len(set(alltheinteractions19)),len(set(alltheinteractions18)),
          len(set(alltheinteractions17)),len(set(alltheinteractions16)),
          len(set(alltheinteractions15)),len(set(alltheinteractions14)),
          len(set(alltheinteractions13)),len(set(alltheinteractions12)),
          len(set(alltheinteractions11)),len(set(alltheinteractions10)),
          len(set(alltheinteractions9)),len(set(alltheinteractions8)),
          len(set(alltheinteractions7)),len(set(alltheinteractions6)),
          len(set(alltheinteractions5)),len(set(alltheinteractions4)), 
          len(set(alltheinteractions3)),len(set(alltheinteractions2)),
          len(set(alltheinteractions1))],'ro')
plt.title('Number of different interactions predicted')
plt.xlabel('times = n ')
plt.ylabel('Different interactions predicted')
plt.show()


#DATAFRAME INTERACTIONS PREDICTED SORTED BY HOW MANY TIMES THEY APPEAR WHEN RUNNING ACOGeneInteraction USING DIFFERENT 'times'
from collections import Counter
Counter(alltheinteractions17)
data = DataFrame(Counter(alltheinteractions20).most_common())
data2 = data.ix[1:]
data3 = data2[0].str.split('_')
data3

filenomgens = os.path.abspath('Fasting/gens200.csv')               #all the infor from 200 gens (among them there are the 100 gens we are working with)
nomgens = pandas.read_csv(filenomgens,index_col=7, header=None)
nomgens


col3 = []
col4 = []
for i, j in data3:
    col3.append(i)
    col4.append(j)
data4 = pandas.DataFrame()
data4 = pandas.DataFrame()
data4['N times'] = data2[1]
data4['Gene1'] = col3
data4['Gene2'] = col4
description1 =[]
for i in col3:
    x = nomgens.at[i,25]
    description1.append(x)
description2 = []
for i in col4:
    x = nomgens.at[i,25]
    description2.append(x)

data4['Description1'] = description1 
data4['Description2'] = description2
data4                                # gens interactions: nome of Gene1 and Gene2, plus a short escription of each gene


#NUMBER OF INTERACTIONS PREDICTED FOR EACH GENE (Genes with more interactions are more likely to have a more central role in the network)
data4
data5 = [data4['Gene1'],data4['Gene2']]
data5
data6 = DataFrame(Counter(data5[1]).most_common())
data6
description1 =[]
for i in data6[0]:
    x = nomgens.at[i,25]
    description1.append(x)

data6['Description'] = description1
data6



#COMPARING KEY GENES USING ACO-TSP AND ACOGeneInteraction
keygensTSP = time1[['Gen1','Gen2']]
keygensTSP
data4
x = data4.ix[:35]
x
keygenesGI = x[['Gene1','Gene2']]
keygenesGI
compKG = pandas.DataFrame()
compKG['keygenesTSP'] = keygensTSP['Gen1'] +'-'+ keygensTSP['Gen2']
compKG['keygenesGI'] = keygenesGI['Gene1'] +'-'+ keygenesGI['Gene2']
compKG
listtotal = list(compKG['keygenesTSP']) + list(compKG['keygenesGI'])
listtotal

listtotalsorted = DataFrame(Counter(listtotal).most_common())
listtotalsorted ##interactions that appear 2 times in the list are interactions share with keyinteractions ACOTSP and ACOGI, the others are inetarctions no shared




