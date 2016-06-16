# Add executable called "helloDemo" that is built from the source files
# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.
add_executable (makeCorrelationRunTree makeCorrelationRunTree.cxx)

# Link the executable to the Hello library. Since the Hello library has
# public include directories we will use those link directories when building
# helloDemo
target_link_libraries (makeCorrelationRunTree LINK_PUBLIC AnitaCorrelator ${ANITA_EVENT} )
  
install (TARGETS makeCorrelationRunTree DESTINATION ${UTIL_BIN_DIR})
