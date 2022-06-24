echo $'************************************************************************'
echo $'*                                                                      *'
echo $'*                         TESTING APE-GEN                              *'
echo $'*                                                                      *'
echo $'************************************************************************'
echo $'\n  *****           Starting with no MODELLER key                  *****'
echo $'\n\n  *    no PTM |  known allotype structure  | no openMM | verbose     *\n'
python New_APE-Gen.py ARSEDEVILS HLA-A*11:01 --verbose
echo $'\n\n  *    no PTM |  known allotype structure  |   openMM  | verbose     *\n'
python New_APE-Gen.py ARSEDEVILS HLA-A*11:01 --verbose --score_with_openmm
echo $'\n\n  *      PTMs |  known allotype structure  |   openMM  | verbose     *\n'
python New_APE-Gen.py ARpSEpTEVIpYS HLA-A*11:01 --verbose --score_with_openmm
echo $'\n\n  *    no PTM | unknown allotype structure |   openMM  | verbose     *'
echo $'  *    error expected (no MODELLER key)                              *\n'
python New_APE-Gen.py ARSEDEVILS HLA-A*11:02 --verbose --score_with_openmm
echo $'\n\n  *      PTMs | unknown allotype structure |   openMM  | verbose     *'
echo $'  *      error expected (no MODELLER key)                            *\n'
python New_APE-Gen.py ARpSEpTEVIpYS HLA-A*11:02 --verbose --score_with_openmm
echo $'\n\n\n\n'

echo $'\n\n  *    no PTM |  known allotype structure  | no openMM |  silent     *\n'
python New_APE-Gen.py ARSEDEVILS HLA-A*11:01 
echo $'\n\n  *    no PTM |  known allotype structure  |   openMM  |  silent     *\n'
python New_APE-Gen.py ARSEDEVILS HLA-A*11:01 --score_with_openmm
echo $'\n\n  *      PTMs |  known allotype structure  |   openMM  |  silent     *\n'
python New_APE-Gen.py ARpSEpTEVIpYS HLA-A*11:01 --score_with_openmm
echo $'\n\n  *    no PTM | unknown allotype structure |   openMM  |  silent     *'
echo $'  *    error expected (no MODELLER key)                              *\n'
python New_APE-Gen.py ARSEDEVILS HLA-A*11:02 --score_with_openmm
echo $'\n\n  *      PTMs | unknown allotype structure |   openMM  |  silent     *'
echo $'  *      error expected (no MODELLER key)                            *\n'
python New_APE-Gen.py ARpSEpTEVIpYS HLA-A*11:02 --score_with_openmm
echo $'\n\n\n\n'


echo $'  *****                  Getting MODELLER key                    *****'
rm /conda/envs/apegen/lib/modeller-10.2/modlib/modeller/config.py
cp testing/modeller_key.py /conda/envs/apegen/lib/modeller-10.2/modlib/modeller/config.py
echo $'\n\n  *    no PTM |  known allotype structure  | no openMM | verbose     *\n'
python New_APE-Gen.py ARSEDEVILS HLA-A*11:01 --verbose
echo $'\n\n  *    no PTM |  known allotype structure  |   openMM  | verbose     *\n'
python New_APE-Gen.py ARSEDEVILS HLA-A*11:01 --verbose --score_with_openmm
echo $'\n\n  *      PTMs |  known allotype structure  |   openMM  | verbose     *\n'
python New_APE-Gen.py ARpSEpTEVIpYS HLA-A*11:01 --verbose --score_with_openmm
echo $'\n\n  *    no PTM | unknown allotype structure |   openMM  | verbose     *\n'
python New_APE-Gen.py ARSEDEVILS HLA-A*11:02 --verbose --score_with_openmm
echo $'\n\n  *      PTMs | unknown allotype structure |   openMM  | verbose     *\n'
python New_APE-Gen.py ARpSEpTEVIpYS HLA-A*11:02 --verbose --score_with_openmm
echo $'\n\n\n\n'

echo $'\n\n  *    no PTM |  known allotype structure  | no openMM |  silent     *\n'
python New_APE-Gen.py ARSEDEVILS HLA-A*11:01 
echo $'\n\n  *    no PTM |  known allotype structure  |   openMM  |  silent     *\n'
python New_APE-Gen.py ARSEDEVILS HLA-A*11:01 --score_with_openmm
echo $'\n\n  *      PTMs |  known allotype structure  |   openMM  |  silent     *\n'
python New_APE-Gen.py ARpSEpTEVIpYS HLA-A*11:01 --score_with_openmm
echo $'\n\n  *    no PTM | unknown allotype structure |   openMM  |  silent     *\n'
python New_APE-Gen.py ARSEDEVILS HLA-A*11:02 --score_with_openmm
echo $'\n\n  *      PTMs | unknown allotype structure |   openMM  |  silent     *\n'
python New_APE-Gen.py ARpSEpTEVIpYS HLA-A*11:02 --score_with_openmm


echo $'\n\n\n\n***********************************************************************'
echo $'*                                                                     *'
echo $'*                         TESTING COMPLETE                            *'
echo $'*                                                                     *'
echo $'***********************************************************************'