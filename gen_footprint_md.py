'''
This is a script that will parse the configuration file and generate a
markdown file with different tables populated with the different entries in
the configuration file.
'''

try:
    import ConfigParser
except ImportError:
    import configparser as ConfigParser

def gen_markdown():

    md_file = 'config_test.md'
    f = open(md_file, 'w')

    fn = 'footprint.cfg'
    backgrounds, cmb, cmbpol, lss = process_config(fn)

    f.write('This file is automatically constructed from the configuration ')
    f.write('file and is meant as an easier way to see what surveys are in ')
    f.write('the configuration file')
    f.write('\n')

#   Write the background table
    f.write('# Backgrounds\n')
    f.write('\n')
    f.write('This table lists all the different backgrounds in the ')
    f.write('configuration file that can be added to the footprint figure ')
    f.write('by inputting a string instead of a Healpix map.\n')
    f.write('\n')
    f.write('| Experiment | Name | Description |\n')
    f.write('|------------|------|-------------|\n')

    for tmp in backgrounds:
        line = '| ' + tmp[0] + ' | ' + tmp[1] + ' | ' + tmp[2] + ' |\n'
        line = line.replace('"', '')
        f.write(line)
    
    f.write('\n')

#   Write the CMB intensity table
    f.write('# CMB Intensity Surveys\n')
    f.write('\n')
    f.write('This next table lists all the surveys associated with CMB ')
    f.write('intensity experiments\n')
    f.write('\n')
    f.write('| Experiment | Survey | Description |\n')
    f.write('|------------|--------|-------------|\n')
    
    for tmp in cmb:
        line = '| ' + tmp[0] + ' | ' + tmp[1] + ' | ' + tmp[2] + ' |\n'
        line = line.replace('"', '')
        f.write(line)

    f.write('\n')

#   Write the CMB polarization table
    f.write('# CMB Polarization Surveys\n')
    f.write('\n')
    f.write('This next table lists all the surveys associated with CMB ')
    f.write('polarization experiments\n')
    f.write('\n')
    f.write('| Experiment | Survey | Description |\n')
    f.write('|------------|--------|-------------|\n')
    
    for tmp in cmbpol:
        line = '| ' + tmp[0] + ' | ' + tmp[1] + ' | ' + tmp[2] + ' |\n'
        line = line.replace('"', '')
        f.write(line)

    f.write('\n')

#   Write the LSS table
    f.write('# Large-Scale Structury Surveys\n')
    f.write('\n')
    f.write('This next table lists all the surveys associated with LSS ')
    f.write('surveys\n')
    f.write('\n')
    f.write('| Experiment | Survey | Description |\n')
    f.write('|------------|--------|-------------|\n')
    
    for tmp in lss:
        line = '| ' + tmp[0] + ' | ' + tmp[1] + ' | ' + tmp[2] + ' |\n'
        line = line.replace('"', '')
        f.write(line)

    f.write('\n')


    f.close()

def process_config(fn):

    config = ConfigParser.ConfigParser()
    config.read(fn)

    output_file = 'readable_config.md'

#   Add the entries to the correct table
    background = []
    cmb = []
    cmbpol = []
    lss = []

    for section in config.sections():
        tableval = config.get(section,'survey_type')
        descr = config.get(section,'description')
        citation = config.get(section,'citation')
        expr_name = config.get(section,'experiment')
        if tableval == 'background':
            background.append([expr_name,section,descr,citation])
        elif tableval == 'cmb':
            cmb.append([expr_name,section,descr,citation])
        elif tableval == 'cmbpol':
            cmbpol.append([expr_name,section,descr,citation])
        elif tableval == 'lss':
            lss.append([expr_name,section,descr,citation])

    return background, cmb, cmbpol, lss

if __name__ == '__main__':
    gen_markdown()
