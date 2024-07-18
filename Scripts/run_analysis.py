import os
import functions

# steps of analysis to be done
# a - combine data, b - plot data, c - correlation of data, d - survival
steps = ['a','b','c','d','e']

# input files
data_folder     = '../Data/'
result_folder   = '../Results/'
data1           = '{}PetLit_Data_01.xlsx'.format(data_folder)
data2           = '{}PetLit_Data_02.xlsx'.format(result_folder)

# first step - combine data
if 'a' in steps:
    log = '{}a_combine.log'.format(result_folder)
    command = 'a_combine_data.py -in {} -out {} -log {}'.format(data1, data2, log)
    functions.run_python_tool(command)

# second step - plot
if 'b' in steps:
    log = '{}b_timeline.log'.format(result_folder)
    command = 'b_timeline.py -in {} -out {} -log {}'.format(data2, result_folder, log)
    functions.run_python_tool(command)

# third step - correlation
if 'c' in steps:
    log = '{}c_correlation.log'.format(result_folder)
    command = 'c_correlation.py -in {} -out {} -log {}'.format(data2, result_folder, log)
    functions.run_python_tool(command)

# fourth step - run survival analyses
if 'd' in steps:
    log = '{}d_survival.log'.format(result_folder)
    command = 'd_survival.py -in {} -out {} -log {}'.format(data2, result_folder, log)
    functions.run_python_tool(command)

# fifth step - collect data for venn diagram
if 'e' in steps:
    command = 'e_venn.py'
    functions.run_python_tool(command)

