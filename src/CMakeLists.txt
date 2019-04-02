add_executable(main
        main.cpp
        Predicates.cpp
        )

if(VTUNE_FOUND)
    target_link_libraries(main
            ${VTUNE_LIBRARIES})
endif()