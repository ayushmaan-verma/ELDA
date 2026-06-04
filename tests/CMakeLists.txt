add_executable(test_main test_main.cpp)
target_link_libraries(test_main PRIVATE elda)
add_test(NAME ELDA_Unit_Tests COMMAND test_main)
