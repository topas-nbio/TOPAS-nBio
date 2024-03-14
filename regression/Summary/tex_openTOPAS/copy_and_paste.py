import os

regression_tests = ['DBSCAN/','LET/','FrickeIRT/','GvalueStepByStep/','Gvalue_LET-SBS/','GvalueIRT/',
                    'Gvalue_LET-IRT/','GvalueIRT_H2O2/','GvalueIRT_H/','GvalueIRT-Temperature/','NanodosimetryI/',
                    'NanodosimetryII/','NanodosimetryIII/','DNASSBPulsed/']

tex_dir = './openTOPAS/'

directory_check = os.path.isdir(tex_dir)
if directory_check == True:
    os.system('rm -rf ' + tex_dir)
    print('Image folder found, deleting and recreating')
    os.system('mkdir ' + tex_dir)
else:
    os.system('mkdir ' + tex_dir)

for test in regression_tests:
    reg_dir = '../../{}/results/'.format(test)

    if test == 'DBSCAN/':
        folder_path = tex_dir + test
        os.system('mkdir ' + folder_path)
        os.system('cp ' + reg_dir + 'mainTopas/DBSCAN1.pdf ' + folder_path + 'DBSCAN1_TsEmDNAPhysics.pdf')
        os.system('cp ' + reg_dir + 'mainTopas/DBSCAN2.pdf ' + folder_path + 'DBSCAN2_TsEmDNAPhysics.pdf')
        os.system('cp ' + reg_dir + 'mainOpt2/DBSCAN1.pdf ' + folder_path + 'DBSCAN1_g4em-dna_opt2.pdf')
        os.system('cp ' + reg_dir + 'mainOpt2/DBSCAN2.pdf ' + folder_path + 'DBSCAN2_g4em-dna_opt2.pdf')
        os.system('cp ' + reg_dir + 'mainOpt4/DBSCAN1.pdf ' + folder_path + 'DBSCAN1_g4em-dna_opt4.pdf')
        os.system('cp ' + reg_dir + 'mainOpt4/DBSCAN2.pdf ' + folder_path + 'DBSCAN2_g4em-dna_opt4.pdf')
        os.system('cp ' + reg_dir + 'mainOpt6/DBSCAN1.pdf ' + folder_path + 'DBSCAN1_g4em-dna_opt6.pdf')
        os.system('cp ' + reg_dir + 'mainOpt6/DBSCAN2.pdf ' + folder_path + 'DBSCAN2_g4em-dna_opt6.pdf')

    elif test == 'LET/':
        folder_path = tex_dir + test
        os.system('mkdir ' + folder_path)
        os.system('cp ' + reg_dir + 'mainTopas/LET.pdf ' + folder_path + 'LET_TsEmDNAPhysics.pdf')
        os.system('cp ' + reg_dir + 'mainOpt2/LET.pdf ' + folder_path + 'LET_g4em-dna_opt2.pdf')
        os.system('cp ' + reg_dir + 'mainOpt4/LET.pdf ' + folder_path + 'LET_g4em-dna_opt4.pdf')
        os.system('cp ' + reg_dir + 'mainOpt6/LET.pdf ' + folder_path + 'LET_g4em-dna_opt6.pdf')

    elif test == 'FrickeIRT/':
        folder_path = tex_dir + test
        os.system('mkdir ' + folder_path)
        os.system('cp ' + reg_dir + 'mainTopas/Gvalue.pdf ' + folder_path + 'Gvalue.pdf')

    elif test == 'GvalueStepByStep/':
        folder_path = tex_dir + test
        os.system('mkdir ' + folder_path)
        os.system('cp ' + reg_dir + 'mainTopas/Gvalue.pdf ' + folder_path + 'Gvalue.pdf')

    elif test == 'Gvalue_LET-SBS/':
        folder_path = tex_dir + test
        os.system('mkdir ' + folder_path)
        os.system('cp ' + reg_dir + 'runMain/Gvalue_LET-SBS.pdf ' + folder_path + 'Gvalue_LET-SBS.pdf')

    elif test == 'GvalueIRT/':
        folder_path = tex_dir + test
        os.system('mkdir ' + folder_path)
        os.system('cp ' + reg_dir + 'mainTopas/GvalueIRT.pdf ' + folder_path + 'GvalueIRT.pdf')

    elif test == 'Gvalue_LET-IRT/':
        folder_path = tex_dir + test
        os.system('mkdir ' + folder_path)
        os.system('cp ' + reg_dir + 'runMain/Gvalue_LET-IRT.pdf ' + folder_path + 'Gvalue_LET-IRT.pdf')

    elif test == 'GvalueIRT_H2O2/':
        folder_path = tex_dir + test
        os.system('mkdir ' + folder_path)
        os.system('cp ' + reg_dir + 'mainPython/TimeEvolution.pdf ' + folder_path + 'TimeEvolution.pdf')

    elif test == 'GvalueIRT_H/':
        folder_path = tex_dir + test
        os.system('mkdir ' + folder_path)
        os.system('cp ' + reg_dir + 'mainPython/TimeEvolution.pdf ' + folder_path + 'TimeEvolution.pdf')

    elif test == 'GvalueIRT-Temperature/':
        folder_path = tex_dir + test
        os.system('mkdir ' + folder_path)
        os.system('cp ' + reg_dir + 'mainPython/TimeEvolution.pdf ' + folder_path + 'TimeEvolution.pdf')
        os.system('cp ' + reg_dir + 'mainPython/TemperatureEvolution.pdf ' + folder_path + 'TemperatureEvolution.pdf')

    elif test == 'NanodosimetryI/':
        folder_path = tex_dir + test
        os.system('mkdir ' + folder_path)
        os.system('cp ' + reg_dir + 'mainTopas/IDDistribution.pdf ' + folder_path + 'IDDistribution_TsEmDNAPhysics.pdf')
        os.system('cp ' + reg_dir + 'mainOpt2/IDDistribution.pdf ' + folder_path + 'IDDistribution_opt2.pdf')
        os.system('cp ' + reg_dir + 'mainOpt4/IDDistribution.pdf ' + folder_path + 'IDDistribution_opt4.pdf')
        os.system('cp ' + reg_dir + 'mainOpt6/IDDistribution.pdf ' + folder_path + 'IDDistribution_opt6.pdf')

    elif test == 'NanodosimetryII/':
        folder_path = tex_dir + test
        os.system('mkdir ' + folder_path)
        os.system('cp ' + reg_dir + 'mainTopas/NanodosimetryII.pdf ' + folder_path + 'nanoII_TsEmDNAPhysics.pdf')
        os.system('cp ' + reg_dir + 'mainOpt2/NanodosimetryII.pdf ' + folder_path + 'nanoII_opt2.pdf')
        os.system('cp ' + reg_dir + 'mainOpt4/NanodosimetryII.pdf ' + folder_path + 'nanoII_opt4.pdf')
        os.system('cp ' + reg_dir + 'mainOpt6/NanodosimetryII.pdf ' + folder_path + 'nanoII_opt6.pdf')

    elif test == 'NanodosimetryIII/':
        folder_path = tex_dir + test
        os.system('mkdir ' + folder_path)
        os.system('cp ' + reg_dir + 'mainTopas/NanodosimetryIII.pdf ' + folder_path + 'nanoIII_TsEmDNAPhysics.pdf')
        os.system('cp ' + reg_dir + 'mainOpt2/NanodosimetryIII.pdf ' + folder_path + 'nanoIII_opt2.pdf')
        os.system('cp ' + reg_dir + 'mainOpt4/NanodosimetryIII.pdf ' + folder_path + 'nanoIII_opt4.pdf')
        os.system('cp ' + reg_dir + 'mainOpt6/NanodosimetryIII.pdf ' + folder_path + 'nanoIII_opt6.pdf')

    elif test == 'DNASSBPulsed/':
        folder_path = tex_dir + test
        os.system('mkdir ' + folder_path)
        os.system('cp ' + reg_dir + 'mainPython/SSB_vs_DMSO.pdf ' + folder_path + 'SSB_vs_DMSO.pdf')
