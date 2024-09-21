import os, sys
import _analysis_scaling_relation
line= []
for config in _analysis_scaling_relation.config.keys():
    for i in range(len(_analysis_scaling_relation.config[config])):
        line+=[str(i) + ': ' + _analysis_scaling_relation.config[config][i]['name'].replace(config,'{}')]
    line+=['================================================', '']
    break
with open('inventory.txt', 'w') as f:
    for line_ in line:
        f.write(line_)
        f.write('\n')
