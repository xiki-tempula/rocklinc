from rocklinc.lipid import StepLipid

lipid = StepLipid([80, 80, 80], [49, 49, 49])
lipid.write_dx_input([6.22602388e-04, 18.5969598, 1.98288109,
                                 24.1998191, 1.91347013],
                      [13.80515, 2, 80])

lipid.write_apbs_input()
lipid.run_APBS(apbs_exe='/Users/grte2001/opt/miniconda3/bin/apbs')
charge, diel, lipid = lipid.read_APBS()