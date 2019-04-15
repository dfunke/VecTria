add_executable(main
        main.cpp
        Predicates.cpp
        )

if(VTUNE_FOUND)
    target_link_libraries(main
            ${VTUNE_LIBRARIES})
endif()

if(Advisor_FOUND)
    target_link_libraries(main
            ${Advisor_LIBRARIES})
endif()

if(Vc_FOUND)
    target_link_libraries(main
            ${Vc_LIBRARIES})
endif()