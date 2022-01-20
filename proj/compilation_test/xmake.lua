includes("src/common")

target("compilation_tests")
    set_kind("binary")
    add_files("*.h","*.cpp","*.cu")