# Add executable called "helloDemo" that is built from the source files
# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.
add_executable (makeCorrelationRunTree makeCorrelationRunTree.cxx)
add_executable (makeSingleEventCorrelationTree makeSingleEventCorrelationTree.cxx)

# Link the executable to the Hello library. Since the Hello library has
# public include directories we will use those link directories when building
# helloDemo
target_link_libraries (makeCorrelationRunTree LINK_PUBLIC AnitaCorrelator AnitaEvent AntarcticaRoot)
  
install (TARGETS makeCorrelationRunTree DESTINATION ${UTIL_BIN_DIR})

target_link_libraries (makeSingleEventCorrelationTree LINK_PUBLIC AnitaCorrelator AnitaEvent AntarcticaRoot)
  
install (TARGETS makeSingleEventCorrelationTree DESTINATION ${UTIL_BIN_DIR})

option ( FORCE_OLD_GPP_ABI " Force old g++ ABI; this might be necessary if using new g++ with ROOT compiled with older g++ or other similar situations" OFF ) 
if (FORCE_OLD_GPP_ABI) 
  add_definitions( -D_GLIBCXX_USE_CXX11_ABI=0 ) 
endif() 
 
