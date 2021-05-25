if(BUILD_GLO_OPT)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../ThirdParty/)

  add_executable(glopt "main.cpp" )
  target_link_libraries(glopt PUBLIC BranchNBound)
endif()
