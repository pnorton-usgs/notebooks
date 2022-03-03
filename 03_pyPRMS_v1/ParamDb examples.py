# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python [conda env:bandit_38]
#     language: python
#     name: conda-env-bandit_38-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
# import numpy as np
# from collections import OrderedDict
from pyPRMS.ParamDb import ParamDb

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11/paramdb_v11_default_CONUS'

# %% [markdown]
# ## Reading a Parameter Database

# %%
# #%prun 
# # %timeit -n1 -r1 
pdb = ParamDb(paramdb_dir=workdir, verbose=True, verify=True)

# using np.loadtxt: 1m25s +/- 4.36s
# using pandas.read_csv: 3.31s
# using iterator: 16s

# %% [markdown]
# ## Accessing parameters and dimensions
# The ParamDb object inherits from ParameterSet. Sub-classes based on ParameterSet have Parameters (parameters property) and Dimensions (dimensions property). The Parameters are composed of one or more Parameter objects. The Dimensions are composed of one or more dimension objects. Note the ParameterSet Dimensions object contains the dimensions used by all of the Parameters in the ParameterSet. Each Parameter will have the dimensions required for a that parameter.

# %% [markdown]
# ### List parameters
# The parameters property returns a Parameters object. The __getattr__ function defaults to returning the Parameters-object ordered dictionary of Parameters. The .keys() will return the name of each Parameter object.

# %%
list(pdb.parameters.keys())

# %% [markdown]
# ### List dimensions
# The parameters property returns a Dimensions object. The __getattr__ function defaults to returning the Dimensions-object ordered dictionary of Dimensions. Note that these are global dimensions; in other words, this is the aggregation of the dimensions used by the Parameters. 

# %%
list(pdb.dimensions.keys())

# %%
pdb.parameters.check()

# %%
pdb.parameters['poi_gage_id'].data

# %%
pdb.parameters['print_type'].data

# %%
pdb.write_paramdb('/Users/pnorton/tmp/check_paramdb_v11')

# %%
pdb.parameters['hru_area'].stats()

# %%
print(pdb.parameters.get('tmax_cbh_adj'))

# %%
pdb.parameters['hru_area'].units

# %%
bad_hrus = [21, 230, 250, 395, 521, 537, 539, 641, 743, 784, 805, 851, 858, 922, 953, 1042, 1083, 1213, 1217, 1274, 1278, 1326, 1389, 1433, 1528, 1837, 2041, 2107, 2116, 3052, 3291, 3966, 4348, 4385, 4408, 4455, 4614, 4679, 4863, 4945, 4953, 5359, 6291, 6542, 6549, 6562, 6653, 6741, 6753, 6791, 6945, 6968, 6983, 6988, 6992, 6993, 7479, 7733, 8337, 8608, 8840, 8841, 8888, 9157, 9447, 9786, 9801, 9802, 9806, 9823, 9868, 9894, 9912, 9932, 9933, 9936, 9962, 9988, 10014, 10018, 10020, 10032, 10105, 10129, 10183, 10248, 10331, 10343, 10407, 10422, 10471, 10475, 10477, 10526, 10600, 10637, 10682, 10729, 10757, 10762, 10824, 10872, 10879, 11268, 11298, 11423, 11696, 11725, 11826, 12034, 12035, 12036, 12038, 12047, 12048, 12053, 12054, 12055, 12057, 12059, 12067, 12092, 12110, 12113, 12146, 12149, 12172, 12203, 12298, 12370, 12529, 12532, 12625, 12663, 12676, 12753, 12782, 12823, 12830, 12839, 12842, 12857, 12910, 13237, 13365, 13377, 13402, 13439, 13451, 13460, 13534, 13666, 13739, 13808, 13900, 14067, 14428, 14467, 14537, 14587, 14696, 14726, 14727, 14730, 14736, 14742, 14743, 14762, 14831, 14882, 15035, 15270, 15405, 15406, 15551, 15617, 15842, 15875, 15905, 15993, 16003, 16090, 16220, 16259, 16303, 16513, 17261, 17274, 17298, 17480, 17585, 17990, 18780, 18838, 18844, 19003, 19120, 19334, 19403, 19461, 19689, 19799, 19822, 19900, 19934, 19944, 19965, 20021, 20036, 20072, 20266, 20723, 20851, 21161, 21256, 21283, 21312, 21359, 21444, 21461, 21530, 21567, 21589, 21651, 21755, 21877, 22037, 22149, 22302, 22363, 22411, 22421, 22556, 22569, 22594, 22595, 22664, 23012, 23124, 23239, 23489, 23928, 23929, 24027, 24144, 24152, 24158, 24652, 24761, 24932, 25032, 25276, 25830, 25884, 25924, 26414, 26794, 26843, 26850, 26852, 26859, 27840, 28238, 28809, 32044, 32735, 34022, 34114, 34137, 34344, 34424, 34469, 34497, 34604, 34944, 34957, 35011, 35655, 36509, 36517, 36556, 36672, 36687, 36728, 36837, 36940, 36951, 36957, 37017, 37034, 37082, 37171, 37175, 37186, 37201, 37222, 37350, 37382, 37419, 37560, 37656, 37683, 37821, 37893, 38448, 38458, 38612, 38763, 38891, 39102, 39104, 39116, 39205, 39224, 39229, 39404, 39495, 39590, 39653, 39686, 39706, 39712, 39791, 39838, 39871, 39877, 39942, 40002, 40014, 40025, 40049, 40076, 40140, 40157, 40159, 40161, 40162, 40169, 40176, 40177, 40189, 40216, 40217, 40231, 40251, 40257, 40265, 40304, 40335, 40421, 40422, 40889, 40935, 40952, 40974, 40988, 41041, 41086, 41092, 41103, 41113, 41126, 41154, 41182, 41187, 41190, 41228, 41233, 41267, 41287, 41292, 41341, 41347, 41424, 41497, 41510, 41581, 41587, 41625, 41717, 41738, 42020, 42049, 42150, 42152, 42177, 42186, 42190, 42270, 42293, 42372, 42803, 43009, 43066, 43243, 43517, 43649, 43918, 43994, 44071, 44141, 44341, 44417, 44550, 44588, 44613, 44664, 44685, 44738, 44747, 44791, 44806, 44824, 44835, 44841, 44844, 44922, 45527, 45624, 45635, 45662, 45672, 45726, 45729, 45732, 45753, 45758, 45761, 45781, 45794, 45836, 45850, 45857, 45863, 45869, 45884, 45894, 45912, 45923, 45933, 45935, 45936, 45937, 45940, 45941, 45943, 45944, 45945, 45948, 45956, 45957, 45958, 45970, 45991, 46025, 46035, 46044, 46088, 46121, 46168, 46197, 46218, 46240, 46243, 46249, 46339, 46341, 46352, 46353, 46376, 46430, 46462, 46471, 46508, 46514, 46532, 46570, 46612, 46645, 46666, 46697, 46700, 46702, 46704, 46710, 46714, 46720, 46742, 46765, 46768, 46783, 46811, 46821, 46823, 46851, 46859, 46861, 46883, 46889, 46890, 46894, 46900, 46902, 46905, 46909, 46930, 47144, 47177, 48335, 48433, 48451, 48473, 49081, 49098, 49114, 49217, 49227, 49272, 49378, 49432, 49449, 49478, 49483, 49516, 49534, 49537, 49561, 49562, 49629, 50344, 50406, 50504, 50632, 50644, 50658, 50663, 50666, 50702, 50712, 50819, 50913, 50927, 50956, 50982, 50984, 51042, 51053, 51116, 51148, 51164, 51192, 51193, 51194, 51197, 51198, 51200, 51270, 51282, 51736, 51808, 51919, 51960, 52011, 52052, 52082, 52127, 52249, 52254, 52262, 52302, 52303, 52326, 52337, 52348, 52361, 52392, 52402, 52432, 52434, 52441, 52451, 52485, 52490, 52508, 52509, 52513, 52554, 52571, 52588, 52651, 52669, 52671, 52679, 52684, 52702, 52718, 52736, 52737, 52748, 52749, 52750, 52773, 52822, 52836, 52842, 52844, 52861, 52884, 52913, 52932, 52962, 53002, 53004, 53015, 53017, 53024, 53098, 53221, 53682, 53790, 53797, 53800, 53844, 53857, 53879, 53883, 53913, 53925, 53962, 53966, 53976, 54000, 54002, 54011, 54015, 54032, 54033, 54052, 54058, 54084, 54109, 54129, 54131, 54138, 54148, 54152, 54157, 54168, 54181, 54186, 54189, 54203, 54215, 54220, 54248, 54255, 54264, 54293, 54318, 54342, 54361, 54374, 54419, 54502, 54536, 54552, 54567, 54659, 54736, 54898, 54903, 54974, 54991, 55145, 55444, 55549, 55768, 55966, 56088, 56171, 56518, 56798, 56803, 56819, 56923, 56944, 56949, 56951, 56964, 56990, 57002, 57099, 57258, 57377, 57382, 57393, 57442, 57462, 57465, 57492, 57520, 57606, 57637, 57651, 57691, 57693, 57794, 57795, 57799, 57803, 57815, 57825, 57836, 57853, 57858, 57859, 57873, 57885, 57892, 57947, 57948, 57949, 57961, 57970, 57973, 58399, 58696, 58730, 58812, 58818, 58923, 59084, 59091, 59409, 59450, 59542, 59735, 59775, 59888, 59934, 60079, 60110, 60194, 60221, 60256, 60324, 60372, 60380, 60389, 60511, 60562, 60569, 60605, 60619, 60636, 60656, 60677, 60700, 60719, 60737, 60747, 60749, 60753, 60768, 61164, 61254, 61259, 61282, 61327, 61345, 61358, 61386, 61596, 61676, 61784, 61905, 61914, 61931, 62021, 62085, 62094, 62134, 62194, 62218, 62233, 62257, 62270, 62271, 62341, 62353, 62380, 62395, 62406, 62418, 62434, 62549, 62569, 62641, 62643, 62669, 62769, 62784, 62797, 62837, 62879, 62978, 63119, 63121, 63125, 63130, 63141, 63824, 63835, 63840, 63911, 63951, 64058, 64091, 64118, 64142, 64168, 64274, 64288, 64325, 64543, 64596, 64633, 64668, 64835, 64852, 64912, 64956, 64967, 64987, 65020, 65024, 65150, 65218, 65256, 65264, 65277, 65317, 65340, 65388, 65408, 65411, 65424, 65589, 65833, 65834, 65866, 65965, 66099, 66119, 66223, 66347, 66375, 66409, 66486, 66607, 66619, 66631, 66664, 66770, 66905, 66966, 66979, 67033, 67061, 67065, 67072, 67106, 67151, 67154, 67224, 67231, 67249, 67259, 67269, 67270, 67274, 67291, 67292, 67294, 67296, 67308, 67315, 67316, 67322, 67324, 67331, 67342, 67355, 67359, 67373, 67379, 67390, 67392, 67404, 67405, 67419, 67425, 67441, 67443, 67463, 67488, 67496, 67497, 67521, 67537, 67548, 67569, 67572, 67589, 67594, 67596, 67613, 67617, 67623, 67667, 67706, 67716, 68816, 68829, 68858, 69080, 69158, 69232, 69296, 69524, 69554, 69577, 69588, 69611, 69759, 69992, 70693, 70731, 70852, 70946, 70968, 70971, 71054, 71071, 71100, 71146, 71188, 71275, 71326, 71341, 71381, 71388, 71426, 71474, 71487, 71517, 71535, 71543, 71579, 71638, 71642, 71643, 71705, 71733, 71742, 71784, 71791, 71810, 71812, 71838, 71887, 71900, 71922, 71960, 71980, 72000, 72004, 72024, 72031, 72086, 72097, 72099, 72105, 72125, 72135, 72136, 72150, 72318, 72344, 72355, 72408, 72511, 72541, 72555, 72557, 72569, 72580, 72609, 72618, 72620, 72660, 72692, 72694, 72723, 72729, 72730, 72734, 72747, 72749, 72774, 72781, 72783, 72792, 72819, 72829, 72873, 72877, 72893, 72909, 72930, 72962, 72985, 72992, 73005, 73018, 73034, 73056, 73086, 73096, 73155, 73156, 73197, 73198, 73200, 73201, 73221, 73246, 73273, 73296, 73449, 73660, 76575, 77239, 77248, 77266, 77269, 77271, 77274, 77287, 77288, 77293, 77296, 77299, 77753, 77801, 78225, 78438, 79002, 79435, 79446, 79454, 79471, 79476, 79488, 79582, 80495, 80524, 80527, 80538, 80578, 80650, 80662, 80690, 80737, 80823, 80879, 80898, 80899, 80903, 80904, 80905, 80906, 80907, 80914, 80958, 81068, 81072, 81074, 81082, 81084, 81088, 81089, 81091, 81095, 81116, 81117, 81119, 81140, 81143, 81157, 81159, 81160, 81168, 81171, 81173, 81175, 81177, 81185, 81188, 81190, 81192, 81193, 81194, 81197, 81203, 81205, 81206, 81210, 81212, 81217, 81218, 81221, 81226, 81227, 81246, 81249, 81254, 81257, 81261, 81270, 81276, 81277, 81279, 81281, 81283, 81286, 81290, 81291, 81292, 81304, 81306, 81341, 81373, 81489, 81506, 81518, 81525, 81541, 81584, 81608, 81620, 81624, 81626, 81630, 81631, 81634, 81635, 81636, 81638, 81644, 81653, 81654, 81657, 81663, 81670, 81678, 81694, 81705, 81706, 81795, 81843, 81844, 81845, 81846, 81848, 81849, 81850, 81852, 81853, 81856, 81863, 81872, 81877, 81878, 81886, 81904, 81916, 81944, 81960, 81976, 81987, 82000, 82025, 82042, 82060, 82072, 82076, 82078, 82083, 82107, 82116, 82132, 82140, 82143, 82147, 82158, 82170, 82177, 82179, 82185, 82187, 82189, 82192, 82213, 82215, 82217, 82219, 82229, 82248, 82278, 82314, 82344, 82365, 82435, 82447, 82478, 82493, 82512, 82553, 82590, 82612, 82628, 82646, 82649, 82664, 82700, 82702, 82705, 82719, 82804, 82884, 82886, 82890, 82891, 82895, 82905, 82913, 82920, 82922, 82925, 82936, 82939, 82940, 82941, 82945, 82947, 82950, 82971, 82974, 82976, 82977, 82978, 82980, 82981, 82984, 82985, 82986, 82991, 82992, 82993, 82994, 82995, 82997, 82998, 83001, 83002, 83004, 83005, 83006, 83011, 83013, 83014, 83032, 83041, 83042, 83045, 83070, 83074, 83075, 83125, 83133, 83163, 83273, 83284, 83305, 83309, 83374, 83404, 83427, 83442, 83455, 83460, 83481, 83542, 83544, 83578, 83582, 83615, 83651, 83655, 83662, 83710, 83712, 83724, 83786, 83791, 83828, 83859, 83865, 83901, 83916, 83928, 83986, 84013, 84015, 84019, 84025, 84064, 84115, 84162, 84297, 84320, 84335, 84675, 84842, 85196, 85204, 85218, 85260, 85286, 85294, 85352, 85355, 85384, 85620, 85945, 85965, 85993, 86006, 86074, 86075, 86078, 86086, 86116, 86135, 86152, 86153, 86188, 86194, 86224, 86248, 86275, 86280, 86425, 86570, 86907, 86915, 86921, 86944, 86961, 86962, 86966, 86974, 86977, 86993, 87004, 87010, 87011, 87012, 87017, 87027, 87063, 87066, 87068, 87092, 87109, 87119, 87122, 87130, 87143, 87149, 87162, 87197, 87200, 87228, 87236, 87270, 87277, 87278, 87284, 87314, 87404, 87425, 87427, 87433, 87434, 87461, 87468, 87472, 87497, 87521, 87553, 87568, 87578, 87610, 87621, 87639, 87647, 87671, 87674, 87793, 87800, 87819, 87822, 87867, 87889, 87895, 87916, 87917, 87919, 87936, 87948, 87950, 87952, 87953, 87974, 88005, 88030, 88065, 88080, 88091, 88098, 88109, 88112, 88165, 88228, 88242, 88266, 88302, 88305, 88313, 88314, 88333, 88437, 88486, 88546, 88597, 88633, 88643, 88794, 88867, 88873, 88887, 88889, 88942, 88943, 88957, 88973, 88978, 88998, 89005, 89012, 89024, 89034, 89054, 89057, 89075, 89089, 89099, 89100, 89111, 89122, 89139, 89154, 89166, 89180, 89232, 89233, 89238, 89259, 89272, 89276, 89300, 89321, 89339, 89346, 89350, 89351, 89383, 89446, 89463, 89514, 89545, 89586, 89587, 89593, 89602, 89647, 89649, 89651, 89665, 89710, 89718, 89784, 89837, 89838, 89866, 89876, 89960, 89966, 90008, 90010, 90022, 90031, 90033, 90227, 90280, 90350, 90353, 90359, 90361, 90366, 90375, 90377, 90378, 90383, 90385, 90391, 90393, 90409, 90410, 90421, 90431, 90434, 90439, 90443, 90450, 90451, 90460, 90477, 90490, 90496, 90527, 90538, 90553, 90556, 90569, 90589, 90599, 90607, 90613, 90623, 90643, 90652, 90655, 90661, 90662, 90664, 90673, 90716, 90719, 90720, 90727, 90743, 90767, 90777, 90802, 90806, 90810, 90812, 90816, 90818, 90830, 90868, 90885, 90890, 90910, 90911, 90913, 90915, 90918, 90927, 90929, 90932, 90939, 90944, 90950, 90954, 90955, 90959, 90960, 90969, 90973, 90978, 91000, 91001, 91016, 91025, 91037, 91052, 91079, 91091, 91092, 91098, 91142, 91162, 91188, 91193, 91237, 91251, 91263, 91293, 91311, 91343, 91353, 91366, 91375, 91376, 91441, 91516, 91535, 91555, 91636, 91640, 91641, 91650, 91700, 91722, 91723, 91734, 91736, 91738, 91742, 91780, 91785, 91809, 91813, 91830, 91848, 91879, 91883, 91892, 91902, 91904, 91908, 91909, 91915, 91929, 91945, 91966, 91975, 91978, 91981, 92008, 92016, 92077, 92078, 92086, 92096, 92106, 92115, 92116, 92124, 92128, 92135, 92155, 92182, 92195, 92202, 92216, 92281, 92286, 92300, 92313, 92320, 92351, 92361, 92373, 92422, 92426, 92430, 92484, 92486, 92543, 92554, 92593, 92622, 92664, 92710, 92725, 92746, 92747, 92771, 92893, 93028, 93043, 93116, 93131, 93198, 93267, 93271, 93278, 93336, 93409, 93422, 93450, 93455, 93456, 93502, 93513, 93586, 93590, 93639, 93672, 93737, 93800, 93831, 93874, 93880, 93886, 93891, 93901, 93927, 93993, 94014, 94029, 94053, 94068, 94073, 94131, 94156, 94158, 94159, 94186, 94218, 94219, 94237, 94272, 94283, 94304, 94325, 94375, 94399, 94406, 94438, 94522, 94641, 94673, 94874, 95339, 95384, 95416, 95420, 95503, 95517, 95643, 95654, 95791, 95796, 95868, 95891, 95895, 95898, 95917, 96021, 96087, 96104, 96257, 96279, 96325, 96468, 96472, 96516, 96568, 96751, 96857, 96883, 96890, 96941, 96946, 97146, 97171, 97244, 97271, 97273, 97314, 97553, 97796, 97821, 97829, 97831, 97834, 97835, 97839, 97862, 97890, 97893, 97895, 97897, 97898, 97912, 97927, 97939, 97954, 97958, 97966, 97973, 97986, 97994, 98009, 98038, 98054, 98148, 98188, 98189, 98245, 98283, 98285, 98342, 98367, 98376, 98389, 98443, 98508, 98641, 98795, 98830, 98848, 98988, 99224, 99696, 100038, 100142, 100228, 100275, 100337, 100351, 100365, 100388, 100459, 100531, 101022, 101080, 101093, 101094, 101241, 101331, 101498, 101522, 101653, 101722, 102983, 103097, 103101, 103347, 104115, 104116, 104117, 104132, 104136, 104143, 104160, 104185, 104204, 104234, 104236, 104244, 104258, 104335, 104361, 104401, 104418, 104421, 104442, 104446, 104523, 104525, 104550, 104622, 104668, 104696, 104754, 104755, 104786, 104793, 104800, 104812, 104817, 104867, 104895, 104897, 104907, 104951, 104975, 104993, 104995, 104997, 105001, 105029, 105038, 105044, 105056, 105072, 105073, 105079, 105088, 105109, 105119, 105124, 105130, 105155, 105160, 105162, 105183, 105213, 105226, 105255, 105266, 105281, 105304, 105339, 105348, 105349, 105354, 105405, 105410, 105451, 105504, 105546, 105565, 105618, 105993, 106042, 106115, 106142, 106195, 106227, 106263, 106264, 106278, 106319, 106374, 106375, 106490, 106503, 106586, 106602, 106626, 106699, 106948, 107437, 107525, 107617, 107673, 107796, 108040, 108043, 108116, 108234, 108241, 108277, 108305, 108335, 108414, 108430, 108435, 108441, 108477, 108701, 109022, 109293, 109436, 109525, 109532, 109562, 109565, 109689, 109753, 109755, 109756, 109771, 109807, 109840, 109960, 109981, 109982, 110005, 110012, 110018, 110034, 110050, 110052, 110089, 110163, 110320, 110469, 111499, 112026, 112107, 114487, 114957]

# %%
hru_area = pdb.parameters['hru_area'].data
bad_areas = []
for xx in bad_hrus:
    bad_areas.append(hru_area[xx-1])
    print(f'{xx}: {hru_area[xx-1]}')

# %%
print(f'min area: {min(bad_areas)}')
print(f'max area: {max(bad_areas)}')

# %%

# %%
len(bad_areas)

# %%
import pandas as pd

# %%
curr_file = f'{workdir}/seg_elev.csv'

tmp_data = pd.read_csv(curr_file, skiprows=0, usecols=[1], squeeze=True)

# %%
tmp_data

# %%
df2 = pd.read_csv(curr_file, skiprows=0, usecols=[1]).squeeze('columns')
df2

# %%

# %%

# %%
