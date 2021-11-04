waLBerla_link_files_to_builddir(*.prm)
waLBerla_link_files_to_builddir(*.obj)

waLBerla_generate_target_from_python(NAME        generate_cumulant_lbm
                                     FILE        generate_cumulant_lbm.py
                                     OUT_FILES   CumulantMRTSweep.cpp CumulantMRTSweep.h
			                         CumulantMRTPackInfo.cpp CumulantMRTPackInfo.h
                        			 InitialPDFsSetter.cpp InitialPDFsSetter.h
                                     #CumulantMRTSimpleUBB.cpp CumulantMRTSimpleUBB.h 
                                     CumulantMRTDynamicUBB.cpp CumulantMRTDynamicUBB.h
                                     CumulantMRTOutflow.cpp CumulantMRTOutflow.h 
                                     CumulantMRTNoSlip.cpp CumulantMRTNoSlip.h 
                                    )

waLBerla_add_executable(NAME    aneurysm_simulation
                        FILES   main.cpp
                        DEPENDS generate_cumulant_lbm geometry domain_decomposition blockforest boundary core field lbm mesh stencil timeloop vtk
                       )

