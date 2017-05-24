#GET ROW DATA AND CALCULATE CORRELATION MATRIX

from pandas import *
import numpy as np
from scipy.stats.stats import pearsonr
import itertools
from sklearn.preprocessing import Normalizer
import os

file = os.path.abspath('Clioquinol/CQexpressiondata.csv')  #expression data
matr = pandas.read_csv(file, index_col=0, header=None)
df1 = matr.T
df2 = np.log2(df1)                              #log2 transformation
normalizeddf2 = Normalizer().fit_transform(df2) #normalization(-1,1)

matriucorrelacio = pandas.DataFrame(normalizeddf2).corr(method='pearson')
nomgensfile = os.path.abspath('Clioquinol/CQgenedata.csv')
nomgens = pandas.read_csv(nomgensfile, index_col=1, header=None)
gens100 = DataFrame(list(matr.index))
genename =[]
for i in gens100[0]:  
    #If the gene has an geneID this name will be used, 
    #otherwise the ID from the platform will be used
    x = nomgens.at[i,3]
    if x == '0' :
        genename.append(i)
    if x != '0':
        genename.append(x)
matriucorrelacio.columns = [genename]
matriucorrelacio.index = [genename]
matriucorrelacio



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
file20 = os.path.abspath('Clioquinol/times/20time.csv')
time20 = pandas.read_csv(file20, index_col=0, header=0)
df['time20'] = time20['0'] +'-'+ time20['1']
file19 = os.path.abspath('Clioquinol/times/19time.csv')
time19 = pandas.read_csv(file19, index_col=0, header=0)
df['time19'] = time19['0'] +'-'+ time19['1']
file18 = os.path.abspath('Clioquinol/times/18time.csv')
time18 = pandas.read_csv(file18, index_col=0, header=0)
df['time18'] = time18['0'] +'-'+ time18['1']
file17 = os.path.abspath('Clioquinol/times/17time.csv')
time17 = pandas.read_csv(file17, index_col=0, header=0)
df['time17'] = time17['0'] +'-'+ time17['1']
file16 = os.path.abspath('Clioquinol/times/16time.csv')
time16 = pandas.read_csv(file16, index_col=0, header=0)
df['time16'] = time16['0'] +'-'+ time16['1']
file15 = os.path.abspath('Clioquinol/times/15time.csv')
time15 = pandas.read_csv(file15, index_col=0, header=0)
df['time15'] = time15['0'] +'-'+ time15['1']
file14 = os.path.abspath('Clioquinol/times/14time.csv')
time14 = pandas.read_csv(file14, index_col=0, header=0)
df['time14'] = time14['0'] +'-'+ time14['1']
file13 = os.path.abspath('Clioquinol/times/13time.csv')
time13 = pandas.read_csv(file13, index_col=0, header=0)
df['time13'] = time13['0'] +'-'+ time13['1']
file12 = os.path.abspath('Clioquinol/times/12time.csv')
time12 = pandas.read_csv(file12, index_col=0, header=0)
df['time12'] = time12['0'] +'-'+ time12['1']
file11 = os.path.abspath('Clioquinol/times/11time.csv')
time11 = pandas.read_csv(file11, index_col=0, header=0)
df['time11'] = time11['0'] +'-'+ time11['1']
file10 = os.path.abspath('Clioquinol/times/10time.csv')
time10 = pandas.read_csv(file10, index_col=0, header=0)
df['time10'] = time10['0'] +'-'+ time10['1']
file9 = os.path.abspath('Clioquinol/times/9time.csv')
time9 = pandas.read_csv(file9, index_col=0, header=0)
df['time9'] = time9['0'] +'-'+ time9['1']
file8 = os.path.abspath('Clioquinol/times/8time.csv')
time8 = pandas.read_csv(file8, index_col=0, header=0)
df['time8'] = time8['0'] +'-'+ time8['1']
file7 = os.path.abspath('Clioquinol/times/7time.csv')
time7 = pandas.read_csv(file7, index_col=0, header=0)
df['time7'] = time7['0'] + '-'+ time7['1']
file6 = os.path.abspath('Clioquinol/times/6time.csv')
time6 = pandas.read_csv(file6, index_col=0, header=0)
df['time6'] = time6['0'] +'-'+ time6['1']
file5 = os.path.abspath('Clioquinol/times/5time.csv')
time5 = pandas.read_csv(file5, index_col=0, header=0)
df['time5'] = time5['0'] +'-'+ time5['1']
file4 = os.path.abspath('Clioquinol/times/4time.csv')
time4 = pandas.read_csv(file4, index_col=0, header=0)
df['time4'] = time4['0'] +'-'+ time4['1']
file3 = os.path.abspath('Clioquinol/times/3time.csv')
time3 = pandas.read_csv(file3, index_col=0, header=0)
df['time3'] = time3['0'] +'-'+ time3['1']
file2 = os.path.abspath('Clioquinol/times/2time.csv')
time2 = pandas.read_csv(file2, index_col=0, header=0)
df['time2'] = time2['0'] +'-'+ time2['1']
file1 = os.path.abspath('Clioquinol/times/1time.csv')
time1 = pandas.read_csv(file1, index_col=0, header=0)
df['time1'] = time1['0'] + '-'+ time1['1']
df.head()

#GRAPH NUMBER OF INTERACTIONS PREDICTED USING DIFFERENT 'times'
import matplotlib.pyplot as plt
f = plt.figure()
plt.plot(list(range(1,21), 
         [len(time1), len(time2),len(time3), len(time4),len(time5),len(time6),len(time7),
          len(time8),len(time9), len(time10),len(time11), len(time12), len(time13), len(time14),
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
import matplotlib.patheffects as path_effects

f = plt.figure()
plt.plot(list(reversed(range(1, 21))), 
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
          len(set(alltheinteractions1))],'g--', path_effects=[path_effects.SimpleLineShadow(),
                       path_effects.Normal()])
plt.title('NUMBER OF DIFFERENT INTERACTIONS PREDICTED', color= 'green')
plt.xlabel('times = n')
plt.ylabel('Different interactions predicted')
plt.annotate('max = 291', xy=(19, 291), xytext=(12, 210),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )
plt.show()


#DATAFRAME INTERACTIONS PREDICTED SORTED BY HOW MANY TIMES THEY APPEAR WHEN RUNNING ACOGeneInteraction USING DIFFERENT 'times'
from collections import Counter
Counter(alltheinteractions17)
data = DataFrame(Counter(alltheinteractions20).most_common())
data2 = data.ix[1:]
data3 = data2[0].str.split('-')
filenomgens = os.path.abspath('Clioquinol/CQgenedata.csv')
nomgens = pandas.read_csv(filenomgens,index_col=1, header=None)
nomgens
col1 =[]
for i in nomgens.index:
    x = nomgens.at[i,3]
    if x == '0' :
        col1.append(i)
    if x != '0':
        col1.append(x)

nomgens.index = col1
nomgens

col3 = []
col4 = []
for i, j in data3:
    col3.append(i)
    col4.append(j)
data4 = pandas.DataFrame()
data4['N times'] = data2[1]
data4['Gene1'] = col3
data4['Gene2'] = col4

description1 =[]
for i in col3:
    x = nomgens.at[i,4]
    description1.append(x)
description2 = []
for i in col4:
    x = nomgens.at[i,4]
    description2.append(x)

data4['Description1'] = description1 
data4['Description2'] = description2
data4


#NUMBER OF INTERACTIONS PREDICTED FOR EACH GENE (Genes with more interactions are more likely to have a more central role in the network)
from collections import Counter
data4
data5 = [data4['Gene1'],data4['Gene2']]
data5
data6 = DataFrame(Counter(data5[1]).most_common())
data6
description1 =[]
for i in data6[0]:
    x = nomgens.at[i,4]
    description1.append(x)

data6['Description'] = description1
data6


#COMPARING KEY GENES USING ACO-TSP AND ACOGeneInteraction
keygensTSP = time1[['0','1']]
keygensTSP
x = data4.ix[:32]
keygenesGI = x[['Gene1','Gene2']]
keygenesGI
compKG = pandas.DataFrame()
compKG['keygenesTSP'] = keygensTSP['0'] +'-'+ keygensTSP['1']
compKG['keygenesGI'] = keygenesGI['Gene1'] +'-'+ keygenesGI['Gene2']
compKG
listtotal = list(compKG['keygenesTSP']) + list(compKG['keygenesGI'])
listtotalsorted = DataFrame(Counter(listtotal).most_common())
listtotalsorted ##interactions that appear 2 times in the list are interactions share with keyinteractions ACOTSP and ACOGI, the others are inetarctions no shared




