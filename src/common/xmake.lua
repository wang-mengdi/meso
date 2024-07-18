set_languages("cxx17")

add_requires("eigen >=3.4.0")
add_requires("fmt =8.1.1")
add_requires("nlohmann_json >=3.10.5")
add_requires("boost =1.78.0", {all=true})
add_requires("cuda", {system=true, configs={utils={"cublas","cusparse","cusolver"}}})


target("common")
    set_kind("static")
    add_headerfiles("*.h")
    add_files("*.cpp","*.cu","*.cxx")
    add_includedirs(".",{public=true})
    add_packages("cuda",{public=true})
    add_packages("eigen",{public=false})
    add_packages("fmt",{public=true})
    add_packages("nlohmann_json",{public=true})
    add_packages("boost",{public=true})
    add_cugencodes("native","compute_75")
    add_cuflags("-extended-lambda --std=c++17")
    if is_plat("windows") then
        add_cuflags("-Xcompiler /bigobj -Xptxas=\"-v\" -rdc=true")
    end
