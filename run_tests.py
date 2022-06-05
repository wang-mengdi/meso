import os
import sys
import platform

def Run_Project(proj_name):
    bin_dir=os.path.join('..','..','bin',proj_name)
    build_dir=os.path.join('..','..','build',proj_name)
    lua_file=os.path.join('proj',proj_name,'xmake.lua')
    proj_dir=os.path.join('proj',proj_name)
    clean_cmd="xmake c -P {} -a".format(proj_dir)
    config_cmd="xmake f -o {} -P {} -y -m release".format(bin_dir, proj_dir)
    #project_cmd="xmake project -k vsxmake -v -a \"x64\" -P {} {}".format(proj_dir, build_dir)
    compile_cmd="xmake build -P {} {}".format(proj_dir,proj_name)
    run_cmd="xmake run -P {} {}".format(proj_dir,proj_name)
    #print(clean_cmd)
    #os.system(clean_cmd)
    print(config_cmd)
    ret=os.system(config_cmd)
    assert (ret==0), "Config failed"
    print(compile_cmd)
    ret=os.system(compile_cmd)
    assert (ret==0), "Compilation failed"
    print(run_cmd)
    ret=os.system(run_cmd)
    assert (ret==0), "Tests failed"

if __name__=='__main__':
    Run_Project("_tests_kernel")
    Run_Project("_tests_dec_system")
    Run_Project("_tests_grid_algorithm")