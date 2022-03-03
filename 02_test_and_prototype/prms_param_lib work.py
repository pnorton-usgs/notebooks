# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.7
# ---

# %%
import prms_param_lib as prms

# %%
filename = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/PRMS_master/06267400/input/params/HaySep.params'
reload(prms)
pobj = prms.parameters(filename)

# %%
replist = [0.5, 0.9, 0.1, 1.0, 0.6, 0.7]
# .757
print pobj.get_var('emis_noppt')
pobj.replace_values('emis_noppt', 5.0)
print pobj.get_var('emis_noppt')

pobj.write_param_file('crap.param')

# %%
#pobj.check_var('gwflow_coef')
pobj.get_var('gwflow_coef')
#pobj.check_all_vars()
#pobj.check_var('dday_slope')

# %%

# %%
replist = [0.5, 0.9, 0.1, 1.0, 0.6, 0.7]
badlist = [0.5, 0.9, 0.1, 1.0, 0.6]

pobj.replace_values('gwflow_coef', badlist)
pobj.get_var('gwflow_coef')

# %%
rowdelim = '####'
catdelim = '**'     # Delimiter for categories of variables (e.g. Dimensions)

filename = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t1/06267400/input/params/HaySep.params'

infile = open(filename, 'r')
rawdata = infile.read().splitlines()
infile.close()

paramdict = {}

it = iter(rawdata)

for line in it:
    #print "Curr:", line
    
    if line == '$Id:$' or line[0:7] == 'Version':
        continue
    if line[0:2] == catdelim:
        if line == '** Dimensions **':
            # new category
            paramdict['Dimensions'] = []
            dimsection = True
        elif line == '** Parameters **':
            paramdict['Parameters'] = []
            dimsection = False

    elif line == rowdelim:
        nextVar = True

    else:
        # We're within a category section and dealing with variables
        if dimsection:
            # Deal with dimensions differently from parameters
            paramdict['Dimensions'].append({'name': line, 'value': next(it)})
        else:
            # We're in a parameter section            
            vardict = {}    # temporary to build variable info

            varname = line.split(' ')[0]
            print 'varname:', varname
            vardict['name'] = varname
            
            numdim = int(next(it))    # number of dimensions for this variable
            dimlist = []
            for dd in range(0,numdim):
                dimlist.append(next(it))
            vardict['dimnames'] = dimlist
            
            numval = int(next(it))
            valuetype = int(next(it))
            vardict['valuetype'] = valuetype
            
            vals = []
            
            for vv in range(0, numval):
                vals.append(next(it))
            vardict['values'] = vals
            
            paramdict['Parameters'].append(vardict)
print "Size of paramdict:", len(paramdict)
print 'Size of paramdict[Dimensions]:', len(paramdict['Dimensions'])
print 'Size of paramdict[Parameters]:', len(paramdict['Parameters'])

# %%
