import os
import sys
import platform

if __name__=='__main__':
    assert (len(sys.argv)==2), 'Error: must pass 1 parameter(s). Example: python make_project.py _tests_reservoir'
    proj_name=sys.argv[1]
    bin_dir=os.path.join('..','..','bin',proj_name)
    build_dir=os.path.join('..','..','build',proj_name)
    lua_file=os.path.join('proj',proj_name,'xmake.lua')
    proj_dir=os.path.join('proj',proj_name)
    config_cmd="xmake f -o {}".format(bin_dir)
    project_cmd="xmake project -k vsxmake -v -a \"x64\" {}".format(build_dir)
    os.chdir(proj_dir)
    print("Current Working Directory" , os.getcwd())
    print(config_cmd)
    os.system(config_cmd)
    print("Current Working Directory" , os.getcwd())
    print(project_cmd)
    os.system(project_cmd)