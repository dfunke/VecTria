add_executable(main
        main.cpp
        Predicates.cpp
        )

if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	target_link_libraries(main
	        stdc++fs
	)
endif()

if (VTUNE_FOUND)
    target_link_libraries(main
            ${VTUNE_LIBRARIES})
endif ()

if (Advisor_FOUND)
    target_link_libraries(main
            ${Advisor_LIBRARIES})
endif ()

if (Vc_FOUND)
    target_link_libraries(main
            ${Vc_LIBRARIES})
endif ()

option(STATS "Detailed statistics" OFF)
if (STATS)
    message(WARNING "STATS is enabled")
    add_definitions(-DENABLE_STATS)
endif()
