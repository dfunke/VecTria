if (TESTS)
    add_executable(test
            main.cpp
            ../src/Predicates.cpp
            ../src/Logger.cpp)

    if (GTest_FOUND)
        target_link_libraries(test ${GTEST_BOTH_LIBRARIES} 'pthread')
    endif ()

    if (Vc_FOUND)
        target_link_libraries(test
                ${Vc_LIBRARIES})
    endif ()
endif ()