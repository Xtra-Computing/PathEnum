SET(CYCLE_ENUMERATOR_SRC_FILE_LIST main.cpp cycle_enumerator.cpp)
ADD_EXECUTABLE(CycleEnumerator.out ${CYCLE_ENUMERATOR_SRC_FILE_LIST})
TARGET_LINK_LIBRARIES(CycleEnumerator.out graph log io)